#define ADAPT 0
#include "phase_change_code/centered-pc.h" // change the projection step
#include "phase_change_code/two-phase-pc.h"
#include "phase_change_code/phase-change.h"
#include "tension.h"
#include "navier-stokes/conserving.h" // This file implements momentum-conserving VOF advection of the velocity components for the two-phase Navier-Stokes solver.
#include "navier-stokes/perfs.h" // record velocity/pressure statistic 
#include "view.h"

#define REDUCED 1
#if REDUCED
# include "reduced.h"
#endif

#define level 6

double maxruntime = HUGE;

#define F_ERR 1e-10

#define T_sat 373.15// could not be zero, because in Lee model for denominator
#define T_sup 10 // difference between wall and saturated temperatures
#define T_wall (T_sat + T_sup)

#define cp_1 4216
#define cp_2 2080
#define lambda_1 0.68
#define lambda_2 0.025
#define D_L lambda_1/(rho1*cp_1)
#define D_V lambda_2/(rho2*cp_2)
#define L_h 2.26e6
#define L0 0.01

face vector u_l[];    // liquid face velocity
face vector uf_new[];
scalar div_1[];       // one-fluid velocity divergence
scalar div_2[];       // liquid velocity divergence
vector ul[];     // liquid velocity divergence

/*********************************** temperature tracers **************************************/
scalar T[], *tracers = {T}; // vapor and liquid temp and use tracer.h to advect them

/************************************ boundary conditions **************************************/
// outflow boundary conditions on heated wall
p[right] = dirichlet(0);
pf[right] = dirichlet(0.);
u.n[right] = neumann(0.);
p_new[right] = dirichlet(0.);
ul.n[right] = dirichlet(0.);

T[left] = dirichlet(T_wall);
pf[left] = dirichlet(0);
p_new[left] = dirichlet(0.);
u.n[left] = dirichlet(0.);
ul.n[left] = neumann(0.);

int main() {
  size (L0); // square domain 1m each side
  init_grid(1<<level); 
  rho1 = 958;
  rho2 = 0.6;
  mu1 = 2.82e-4;
  mu2 = 1.23e-5;
  f.sigma = 0.059;
  TOLERANCE = 1e-3;
  #if REDUCED
    G.x = 0;
  #endif
  run();
}

//#define H0 L0/(1>>level)
#define H0 3.225e-4
#define plane(x, y, H) (x-H0)
#define rhoo (1/rho2 - 1/rho1)
event init (t = 0) {
  fraction (f,  plane(x,y,H0));

  foreach()
    T[] = (1-f[])*(T_wall-T_sup/H0*x) + f[]*T_sat;
  foreach_face()
    uf.x[] = 0.;
  boundary({T, uf});
  }

/************************ define mass source *******************************/
scalar velocity[];     // velocity magnitude
scalar m_dot[];
scalar div_pc[];
face vector evp[];
event stability(i++)   // executed before vof and ns events
{
  scalar delta_s[];
  delta_magnini(f,delta_s);
  foreach()
    {
      T.tr_eq = T_sat;
      T.rho = rho1*f[]+rho2*(1-f[]);
      T.cp = cp_1*f[]+cp_2*(1-f[]);
      T.lambda = lambda_1*f[]+lambda_2*(1-f[]);
    }
  /* if you calculate mass flux, we can use several mass transfer models, just change the function names, I have
  already make the parameters in the function same. */
  lee_model(T,f,m_dot,L_h); 
  /* when calculating the divergenc of velocity, we have three methods:
  1. using mass diffusion method, this method is for the interface instable case, for example: large density ratios
  2. just multiply the delta_s function with the mass flux, the predition is delta_s function can be computed
  3. calculte the divergence after we get the new vof value.
  */
  foreach()
    div_pc[] = m_dot[]*delta_s[];
  boundary({div_pc});
  //mass_diffusion(div_pc,f,m_dot); //Hardt method for large density ratios
  // compute velocity magnitude
  foreach()
    velocity[] = sqrt(sq(u.x[]) + sq(u.y[]));
  boundary({velocity});
}



event vof(i++)
{
  //double dt_1 = 1/(D_L*4/sq(L0/(1<<level)));
/******************* step-2: shift interface accounting for phase-change *******************/  
  foreach()
  {
   /* Note that VOF function is clipped when the interface displacement extends beyond the cell boundary.
   * This could be solved by addig the clipped value to neighbouring cells, which is automatically solved by Basilisk vof solver*/
  if(interfacial(point,f))
    {
      // f_old = f[];
      coord n = interface_normal( point, f); // interface normal
      double alpha = plane_alpha(f[],n); // distance from original point to interface 
      alpha -= m_dot[]*dt/(Delta*rho1)*sqrt(sq(n.x)+sq(n.y)); // phase-change shifting distance 
      f[] = line_area(n.x,n.y,alpha); // cell volume fraction for n+1
      //div[] = (f_old - f[])/dt*(rho1/rho1 - 1); // this is the third method to calculate divergence of velocity
    }
  }
  boundary({f,T});
}

scalar dd[];
event tracer_diffusion(i++){
  double T_p;
  T.tr_eq = T_sat;
  foreach()
    { 
      dd[] = 1 - f[];// volume fraction of vapor
      coord n  = interface_normal( point, dd) , p, q; 
      q.x = -n.x, q.y = -n.y;
      double alpha  = plane_alpha (dd[], q);
      line_center (q, alpha,dd[],&p);
      if(interfacial(point,f))
        T_p = interpolate_2 (point, T, p);// we calculate the temperature in the barycenter of vapor portion in interfacial cells
      if(interfacial(neighborp(-1,-1),f)&&f[]==1)
        T[] = 2*T.tr_eq - T_p;
      T[] = clamp(T[],T.tr_eq,T_wall);
      T.D = D_V*f[] + D_L*(1-f[]);
      T.rho = rho1*f[]+rho2*(1-f[]);
      T.cp = cp_1*f[]+cp_2*(1-f[]);
      T.lambda = lambda_1*f[]+lambda_2*(1-f[]);
    }
    heat_source (T, f, div_pc, L_h);// add heat source term in diffusion equation
}

/* first approximation projection method after prediction of face velocity to compute advection velocity, 
this step also considers mass source due to projection step
*/
event advection_term (i++,last)
{
  if (!stokes) {
    mgpf = project_pc (uf, pf, alpha, dt/2., mgpf.nrelax,div_pc,rhoo);
  }
}


event projection(i++){
  
  project_pc(uf,p,alpha,dt,4,div_pc,rhoo); // projection step for entire domain
  centered_gradient (p, g);
  correction (dt);

  // copy velocity to u_l
  foreach_face()
    u_l.x[] = uf.x[];
  boundary((scalar *){u_l});

  // compute de divergence of u_l (should be exactly equal to div_pc)
  foreach(){
    div_1[] = 0;
    foreach_dimension()
      div_1[] += (u_l.x[1] - u_l.x[0])/Delta;
  }
  boundary({div_1}); 
  
  // solve pressure equation with phase change velocity divergence
  scalar divv[];
  foreach(){
    divv[] = div_pc[]*rhoo;
    foreach_dimension()
      divv[] += 0;
    divv[] /= dt;
  }
  poisson (p_new, divv, alpha, tolerance = TOLERANCE/sq(dt), nrelax = 4);
  
  // add phase change correction to u_l
  foreach_face()
    {
      uf_new.x[] = -dt*alpha.x[]*face_gradient_x (p_new, 0);
      u_l.x[] += uf_new.x[];
      ul.x[] = (uf_new.x[] + uf_new.x[1])/(fm.x[] + fm.x[1] + SEPS);
    }
  boundary ((scalar *){u_l, uf_new, ul});
  // compute again velocity divergence of u_l (should be zero)
  foreach(){
    div_2[] = 0;
    foreach_dimension()
      div_2[] += (u_l.x[1] - u_l.x[0])/Delta;
  }
  boundary({div_2});  
}

#if ADAPT
event adapt (i++) {
  adapt_wavelet ({f}, (double[]){1e-6}, 8,6);
}
#endif


void mg_print (mgstats mg)
{
  if (mg.i > 0 && mg.resa > 0.)
    printf ("%d %g %g %g %d ", mg.i, mg.resb, mg.resa,
	    mg.resb > 0 ? exp (log (mg.resb/mg.resa)/mg.i) : 0,
	    mg.nrelax);
}


event logfile (t = 0; t <= 0.1; t += 0.01) {

  double xb = 0., vx = 0.,vy = 0., sb = 0.,yb = 0., nu = 0.;
  foreach(reduction(+:xb) reduction(+:vx) reduction(+:sb) reduction(+:yb) reduction(+:vy) reduction(+:nu)) {
    double dv = (1. - f[])*dv();
    xb += x*dv;
    yb += y*dv;
    vx += u.x[]*dv;
    vy += u.y[]*dv;
    sb += dv; 
  }
  printf (" %g %g %g %g %g %g %g %g %g %g ", 
	  t, sb, -1.,  xb/sb, yb/sb,vx/sb, vy/sb,dt, perf.t, perf.speed);
  mg_print (mgp);
  mg_print (mgpf);
  mg_print (mgu);
  putchar ('\n');
  fflush (stdout);
}

// output
event snap (t=0;t += 0.01)
 {
   char name[80];
   sprintf (name, "snapshot-%g.gfs",t);
   output_gfs (file = name);
 }


event movies (t = 0; t <= 0.1; t += 0.01) {
  
  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f});

  char legend[1000];
  sprintf(legend, "t = %0.2g", t);
  clear();
// interface
  view (tx = 0,ty=-0.5);
  draw_vof ("f",lc = {0,0,0},lw = 2);
  squares("T");
   draw_string(legend, 1, size = 30., lw = 2.);
    box (notics=true,lw=2);
     mirror({-1}){
    draw_vof("f",lc = {0,0,0},lw = 2);
    draw_string(legend, 1, size = 30., lw = 2.);
     box (notics=true,lw=2);
}
save("interface.mp4");
  }