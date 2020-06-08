#define ADAPT 1
//#define FILTERED// for large density and viscosity ratios
#include "phase_change_code/centered-pc.h" // change the projection step
//#include "phase_change_code/double-projection-pc.h"
#include "phase_change_code/two-phase-pc.h"
#include "phase_change_code/mass-transfer.h" // mass transfer models
#include "phase_change_code/heat-transfer.h" // heat transfer models
#include "tension.h"
#include "phase_change_code/conserving-pc.h" // This file implements momentum-conserving VOF advection of the velocity components for the two-phase Navier-Stokes solver.
#include "navier-stokes/perfs.h" // record velocity/pressure statistic 
#include "view.h"

#define REDUCED 1
#if REDUCED
# include "reduced.h"
#endif
#define maxlevel 6
#define level 6
#define minlevel 5
#define F_ERR 1e-10

#define T_sat 373.15// could not be zero, because in Lee model for denominator
#define T_sup 10 // difference between wall and saturated temperatures
#define T_wall (T_sat + T_sup)

#define cp_1 4216
#define cp_2 2030
#define lambda_1 0.68
#define lambda_2 0.025
#define D_L lambda_1/(rho1*cp_1)//1.68e-7
#define D_V lambda_2/(rho2*cp_2)//2.05e-5
#define L_h 2.26e6
#define L0 0.01

face vector u_l[];    // liquid face velocity
face vector uf_new[];
scalar div_1[];       // one-fluid velocity divergence
scalar div_2[];       // liquid velocity divergence
scalar p_new[];
vector ul[];     // 

/*********************************** temperature tracers **************************************/
scalar T[], *tracers = {T}; // vapor and liquid temp and use tracer.h to advect them

/************************************ boundary conditions **************************************/
// outflow boundary conditions on liquid side wall
T[right] = dirichlet(T_sat);
p[right] = dirichlet(0);
pf[right] = dirichlet(0.);
u.n[right] = neumann(0.);
ul.n[right] = neumann(0.);
// stationary wall boundary condition on left superheated wall
T[left] = dirichlet(T_wall);
pf[left] = dirichlet(0.);
p_new[left] = dirichlet(0.);
u.n[left] = dirichlet(0.);
ul.n[left] = dirichlet(0.);

u_l.n[top] = 0;
u_l.n[bottom] = 0;
u_l.n[left] = 0;

int main() {
  size (L0); // square domain 1m each side
  init_grid(1<<level); 
  rho1 = 958;
  rho2 = 0.6;
  mu1 = 2.82e-4;
  mu2 = 1.23e-5;
  f.sigma = 0.059;
  TOLERANCE = 1e-4;
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
  //mask(y>4*L0/(1<<minlevel)?top:none);
  fraction (f,  plane(x,y,H0));

  foreach()
    T[] = (1-f[])*(T_wall-T_sup/H0*x) + f[]*T_sat;
  foreach_face()
    uf.x[] = 0.;
  boundary({T, uf});

  //unrefine(y>4*L0/(1<<minlevel));
  }

/****************************************** analytical solution **********************************************/
// Newton interpolation method
#define TOLER 1e-6
double coeff()
{
  double x = 0.05, x_old;  
  double f, df;
  do
  {
    x_old = x;
    f = x*exp(sq(x))*erf(x) - cp_2*T_sup/(L_h*sqrt(pi)); 
    df = 2*sq(x)*exp(sq(x))*erf(x) + exp(sq(x))*erf(x) + 2*x/sqrt(pi);
    x = x_old - f/df;
  }while (fabs(x - x_old) > TOLER);
  return x;
}

double exact(double delta_x, double t)
{
  double a = coeff();
  delta_x = 2*a*sqrt(D_V*t);

  return delta_x;
}



/************************ define mass source *******************************/
scalar velocity[];     // velocity magnitude
scalar m_dot[];
scalar div_pc[],ff[];
face vector evp[];scalar temp[];
event mass_flux(i++)   
{
  scalar delta_s[];
  delta_magnini(f,delta_s);
  foreach()
    {
      f[] = clamp(f[],0,1);
      T.tr_eq = T_sat;
      T.rho = rho1*f[]+rho2*(1-f[]);
      T.cp = cp_1*f[]+cp_2*(1-f[]);
      T.lambda = lambda_1*f[]+lambda_2*(1-f[]);
      T.rho_cp = f[]*rho1*cp_1 + (1-f[])*rho2*cp_2;
    }
    
  /* if you calculate mass flux, we can use several mass transfer models, just change the function names, I have
  already make the parameters in the function same. */
  /*************************************************************************************************/
  /*********************************** mass transfer models ***************************************/
  /*************************************************************************************************/

  sharp_simple_model(T,f,m_dot,L_h); // !!!!! change this to different phase change models
  //sun_model(T,f,div_pc,L_h); // !!!!! change this to different phase change models
  //zhang_model(T,f,m_dot,L_h); //note:!!! this should be combined with mass_diffusion(), not (m_dot[]*delta_s[])
  //zhang_model_1(T,f,m_dot,L_h); //note:!!! this should be combined with mass_diffusion(), not (m_dot[]*delta_s[])
  //zhang_model_2(T,f,m_dot,L_h); //note:!!! this should be combined with mass_diffusion(), not (m_dot[]*delta_s[])
  //malan_model(T,f,m_dot,L_h); //note:!!! this should be combined with mass_diffusion(), not (m_dot[]*delta_s[])
  //tanasawa_model(T,f,m_dot,L_h); //note:!!! this should be combined with mass_diffusion(), not (m_dot[]*delta_s[])
  //lee_model(T,f,m_dot,L_h);
  //rattner_model(T,f,m_dot,L_h,dt); // this should be checked again, the resource paper is 'NUMERICAL ANALYSIS ON HEAT TRANSFER CHARACTERISTICS OF LOOPED MINICHANNEL USING PHASE-CHANGE VOF METHOD'
  
  /* when calculating the divergenc of velocity, we have three methods:
  1. using mass diffusion method, this method is for the interface instable case, for example: large density ratios
  2. just multiply the delta_s function with the mass flux, the precondition is delta_s function can be computed
  3. calculte the divergence after we get the new vof value.
  */

 /* note: two different methods to calculate volumetric mass source terms:
 */

// method 1: 
foreach()
  temp[] = m_dot[]*delta_s[];
boundary({temp});
mass_diffusion(div_pc,temp);


//scalar temp_1[];
// method 2:
  //mass_diffusion(temp_1,temp);
  //foreach()
    //div_pc[] = temp_1[]*delta_s[];
  //boundary({div_pc});

  
foreach()
  ff[] = clamp(ff[],0,1);
boundary({ff});
volume_diffusion(ff,f); //Hardt method

// compute velocity magnitude
  foreach()
    velocity[] = sqrt(sq(u.x[]) + sq(u.y[]));
  boundary({velocity});
}

event vof(i++)
{

/******************* step-1: construct entire divergence-free domain *******************/  
  scalar divv[];
  foreach(){
    divv[] = div_pc[]*rhoo;
    //divv[] = div_pc[]; // for volumetric mass source
    divv[] /= dt;
  }
  poisson (p_new, divv, alpha, tolerance = TOLERANCE/sq(dt), nrelax = 4);
  
  // add phase change correction to u_l
  foreach_face()
    {  
      uf_new.x[] = -dt*alpha.x[]*face_gradient_x (p_new, 0);
      u_l.x[] = uf.x[] + uf_new.x[];
      ul.x[] = (u_l.x[] + u_l.x[1])/(fm.x[] + fm.x[1] + SEPS);
    }
  boundary ((scalar *){u_l, uf_new, ul});
  // compute again velocity divergence of u_l (should be zero)
  foreach(){
    div_2[] = 0;
    foreach_dimension()
      div_2[] += (u_l.x[1] - u_l.x[0])/Delta;
  }
  
  foreach_boundary(left)
    div_2[] = 0;
  boundary({div_2});  

/******************* step-2: shift interface accounting for phase-change *******************/  
  foreach()
  {
   /* Note that VOF function is clipped when the interface displacement extends beyond the cell boundary.
   * This could be solved by addig the clipped value to neighbouring cells, which is automatically solved by Basilisk vof solver*/
  if(f[]>1e-12&&f[]<1-1e-12)
    {
      //double f_old = f[];
      coord n = interface_normal( point, f); // interface normal
      double alpha = plane_alpha(f[],n); // distance from original point to interface 
      //alpha -= m_dot[]*dt/(Delta*rho1)*sqrt(sq(n.x)+sq(n.y)); // phase-change shifting distance 
      alpha -= div_pc[]*dt/rho1*sqrt(sq(n.x)+sq(n.y)); // or use volumetric mass source in vof reconstruction directly
      f[] = line_area(n.x,n.y,alpha); // cell volume fraction for n+1
      //div_pc[] = (f_old - f[])/dt*(rho1/rho2 - 1); // for volumetric mass source, a little large for velocity divergence, this will cause large new velocity in extended domain
        /*
      the difference between mass flux(\dot{m}) and volumetric source term (S_h):
      S_h = \dot{m} |\nabla f| = \dot{m} \delta_s
      the unit is (kg/(m^2 s)) and (kg/(m^3 s))

      Therefore, we need to change the Malan expression to
        alpha -= div_pc[]*dt/rho1;
        not alpha -= div_pc[]*dt/(rho1*Delta);

        note this!
      */
     //f[] -= div_pc[]*dt/rho1;
    }
  }
  boundary({f});
}
/*
event tracer_advection(i++){
  foreach_face()
    uf.x[] = 0;
  boundary((scalar *){uf});
}
*/
scalar tt[];
event tracer_diffusion(i++){
  T.tr_eq = T_sat;
  scalar dd[],t_av[]; // volume fraction of vapor
  double T_p;
  vertex scalar t_cor[];
  foreach_vertex()
    t_cor[] = (T[] + T[-1,0] + T[0,-1] + T[-1,-1])/4;
  boundary({t_cor});
  foreach()
    { 
      /* step-1: determine ghost cell in liquid.
      In Zhang's paper, the ghost cell is defined as:
      1. the mixed cell whose center is outside the interface (use the distance from cell center to interface to determine this)
      2. the pure liquid cell which has a corner neighbor belonging to the mixed cells but not belonging to the ghost cells
      -> but here we use a special treatment, we calculate the vapor portion barycenter temperature in interface cell. 
      Then using the coner pure liquid cell to balance the expression: T[] = 2*T.tr_eq - T_p to keep the pure liquid cell saturated state.
      */
      dd[] = 1 - f[];// volume fraction of vapor
      coord n  = interface_normal( point, f) , p, q; // interface normal, barycenter in vapor potion, negative interface normal
      q.x = -n.x, q.y = -n.y; 
      double alpha  = plane_alpha (dd[], q); // cell center to interface but for vapor
      double alpha_l  = plane_alpha (f[], n); // cell center to interface but for vapor
      line_center (q, alpha,dd[],&p); // barycenter in vapor portion 
      coord pp,qq;
      //line_length_center (n, alpha_l, &pp); // interface center coordinate
      double distance_1;
      pp.x = (alpha_l + 0.5*(n.x + n.y))*n.x + 0.5*(n.x + n.y)*n.x;
      pp.y = (alpha_l + 0.5*(n.x + n.y))*n.y + 0.5*(n.x + n.y)*n.y;
      distance_1 = fabs(alpha_l + 0.5*(n.x + n.y) + 0.5*(n.x + n.y));
      qq.x = pp.x + distance_1*n.x;
      qq.y = pp.y + distance_1*n.y;
      
      // if the ghost cell in mixed cell
      if(interfacial(point,f)&&alpha_l>0.5*(n.x + n.y))
        //T_p = interpolate_1 (point, T, p);// we calculate the barycenter of interfacial vapor portion temperature
        {
          t_av[] = (t_cor[] + t_cor[1,0] + t_cor[0,1] + t_cor[1,1])/4;
          T_p =  interpolate_2(point,t_av,qq);
          //T_p =  interpolate_2(point,T,qq);
          //T[] = T_p*(1-f[]) + f[]*T_sat;
          tt[] = T_p;
        }

      // ghost cell method for different situations of interface normal directions, combined with heat_source() function
      if(n.x==0&&n.y>0)
        {
          if((interfacial(neighborp(0,1),f)&&f[]==1)||(interfacial(point,f)&&alpha_l<0.5*(n.x + n.y)))
            T[] = 2*T.tr_eq - T_p;
        }
      else if(n.x==0&&n.y<0)
        {
          if((interfacial(neighborp(0,-1),f)&&f[]==1)||(interfacial(point,f)&&alpha_l<0.5*(n.x + n.y)))
            T[] = 2*T.tr_eq - T_p;
        }
      else if(n.x>0&&n.y==0)
        {
          if((interfacial(neighborp(1,0),f)&&f[]==1)||(interfacial(point,f)&&alpha_l<0.5*(n.x + n.y)))
            T[] = 2*T.tr_eq - T_p;
        }
      else if(n.x<0&&n.y==0)
        {
          if((interfacial(neighborp(-1,0),f)&&f[]==1)||(interfacial(point,f)&&alpha_l<0.5*(n.x + n.y)))
            T[] = 2*T.tr_eq - T_p;
        }

      else if(n.x<0&&n.y<0)
        {
          if((interfacial(neighborp(-1,-1),f)&&f[]==1)||(interfacial(point,f)&&alpha_l<0.5*(n.x + n.y)))
            T[] = T.tr_eq - T_p;
        }
      else if(n.x<0&&n.y>0)
        {
          if((interfacial(neighborp(-1,1),f)&&f[]==1)||(interfacial(point,f)&&alpha_l<0.5*(n.x + n.y)))
            T[] = 2*T.tr_eq - T_p;
        }
      else if(n.x>0&&n.y<0)
        {
          if((interfacial(neighborp(1,-1),f)&&f[]==1)||(interfacial(point,f)&&alpha_l<0.5*(n.x + n.y)))
            T[] = 2*T.tr_eq - T_p;
        }
      else if(n.x>0&&n.y>0)
        {
          if((interfacial(neighborp(1,1),f)&&f[]==1)||(interfacial(point,f)&&alpha_l<0.5*(n.x + n.y)))
            T[] = 2*T.tr_eq - T_p;
        }

      T[] = clamp(T[],T.tr_eq,T_wall);
      if(f[] == 1)
        T[] = T_sat;
        
      if(interfacial(point,f))
        T.D = D_L*f[] + D_V*(1-f[]);
      else 
        T.D = D_V;
      T.rho = rho1*f[]+rho2*(1-f[]);
      T.cp = cp_1*f[]+cp_2*(1-f[]);
      T.lambda = lambda_1*f[]+lambda_2*(1-f[]);
      T.rho_cp = f[]*rho1*cp_1 + (1-f[])*rho2*cp_2;
    }
  /*
  three different methods to treat with energy diffusion, at this moment, just try the first one, it should work with most mass transfer models
  */
  heat_source (T, f, div_pc, L_h);// add heat source term in diffusion equation
  //dirichlet_diffusion(T,f,level,dt,10);
  //rattner_diffusion(T,f,div_pc,dt,L_h);
}


/* first approximation projection method after prediction of face velocity to compute advection velocity, 
this step also considers mass source due to projection step
*/

event advection_term (i++,last)
{
  if (!stokes) {
    mgpf = project_pc (uf, pf, alpha, dt/2., mgpf.nrelax,div_pc,rhoo);
    //// mgpf = project_liquid(uf, pf, alpha, dt/2., mgpf.nrelax,div_pc); // for volumetric mass source
  }
}

event acceleration (i++) {
  face vector av = a;
  foreach_face()
    {
      if (interfacial(point,f)) 
        av.x[] += alpha.x[]/fm.x[]*((1/rho2-1/rho1)*sq(m_dot[])*(f[]-f[-1])/Delta); //recoil pressure due to mass transfer
    }
  boundary((scalar *){av});
}


event projection(i++){
/*
we add the volumetric mass source in poisson equation projection step after solving intermediate fluid velocity.
This is for solving the fluid velocity and pressure with phase change.It is similar to the methodology without phase change.
*/
  project_pc(uf,p,alpha,dt,4,div_pc,rhoo); // projection step for entire domain
  //project_liquid(uf,p,alpha,dt,4,div_pc); // for volumetric mass source
  centered_gradient (p, g);// calculate centered acceleration
  correction (dt); // correct face velocity on centered one

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

}

#if ADAPT
event adapt (i++) {
  adapt_wavelet ({f}, (double[]){1e-6}, maxlevel,minlevel);
}
#endif

void mg_print (mgstats mg)
{
  if (mg.i > 0 && mg.resa > 0.)
    printf ("%d %g %g %g %d ", mg.i, mg.resb, mg.resa,
	    mg.resb > 0 ? exp (log (mg.resb/mg.resa)/mg.i) : 0,
	    mg.nrelax);
}

event logfile (t = 0; t <= 10; t += 0.01) {
  double aaa;
  double delta_interface = exact(aaa,t);

  double xb = 0., vx = 0.,vy = 0., sb = 0.,yb = 0., nu = 0.;
  foreach(reduction(+:xb) reduction(+:vx) reduction(+:sb) reduction(+:yb) reduction(+:vy) reduction(+:nu)) {
    double dv = (1. - f[])*dv();
    xb += x*dv;//interface position
    yb += y*dv;
    vx += u.x[]*dv;
    vy += u.y[]*dv;
    sb += dv; 
  }
  printf (" %g %g %g %g %g %g %g %g %g %g", 
	  t, sb, delta_interface,  2*xb/sb, yb/sb,vx/sb, vy/sb,dt, perf.t, perf.speed);
  mg_print (mgp);
  mg_print (mgpf);
  mg_print (mgu);
  putchar ('\n');
  fflush (stdout);
}

/*************************************** output ************************************************/
event snap (t+=0.1)
 {
   char name[80];
   sprintf (name, "snapshot-%g.gfs",t);
   output_gfs (file = name);
 }

event interface_position (t += 0.01) {
  scalar pos[];
  position (f, pos, {1,0});
  double max = statsf(pos).max;
  char name[80];
  sprintf (name, "in-pos");
  static FILE * fp = fopen (name, "w");
  fprintf (fp, "%g %g\n", t, max);   
  fflush (fp);
}

event temperature_value (t+=0.1){

  char *name = NULL;
  name = (char *) malloc(sizeof(char) * 256);
  sprintf (name, "temperature/termperature-%.1f",t);
  FILE*fp1 = fopen (name, "w");

  for (double x = 0; x < L0; x += 0.001)
      fprintf (fp1, " %g %g\n", x,
          interpolate (T, x, 0.005)); // generate velocity vector field, at least 12 grids space in liquid thin film(\delta)----reference from "Co-current flow effects on a rising Taylor bubble"
      fclose(fp1);
}

event velocity_magnitude (t+=0.1){

  char *name = NULL;
  name = (char *) malloc(sizeof(char) * 256);
  sprintf (name, "velocity/velocity-%.1f",t);
  FILE*fp1 = fopen (name, "w");

  for (double x = 0; x < L0; x += 0.001)
      fprintf (fp1, " %g %g %g\n", x,
          interpolate (u.x, x, 0.005),interpolate (u.y, x, 0.005)); // generate velocity vector field, at least 12 grids space in liquid thin film(\delta)----reference from "Co-current flow effects on a rising Taylor bubble"
      fclose(fp1);
}

event movies (t = 0; t <= 10; t += 0.1) {
  
  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f});

  char legend[1000];
  sprintf(legend, "t = %0.2g", t);
  clear();
// interface
  view (tx = -0.5,ty=-0.5);
  draw_vof ("f",lc = {0,0,0},lw = 2);
  squares("T");
   draw_string(legend, 1, size = 30., lw = 2.);
    box (notics=true,lw=2);
save("interface.mp4");
  }

/*
*** gnuplot ******
reset
set size square
set title 'Stefan problem'
set xlabel 'Time t'
set ylabel 'delta'
set key left
p [0:10][0:0.002]'out' every 10 u ($1+0.3):4 w p t 'Basilisk simulation',\
'reference' u 1:2 w l lw 2 lc rgb 'red' t 'Analytical solution'
*/