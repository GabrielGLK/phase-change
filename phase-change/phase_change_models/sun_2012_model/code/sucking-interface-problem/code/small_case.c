#define ADAPT 0
#define REDUCED 1

/***********************************************************************************************************************/
//#define FILTERED// smeared fluid properties for large ratios, in 'two-phase-pc.h'
#include "phase_change_code/centered-pc.h" // solve momentum equation 
#include "phase_change_code/two-phase-pc.h" // two-phase flow model
#include "phase_change_code/tracer-pc.h" // passive tracer (advection equation)
#include "tension.h" // surface tension model
#include "phase_change_code/heat-transfer.h" // diffusion equations
#include "phase_change_code/mass_transfer/mass_diffusion.h" // use Hardt method to do volume diffusion or mass diffusion
/*********************************************** phase-chnage models *********************************************/
/********************************************* mass transfer models options **********************************************/
#include "phase_change_code/mass_transfer/sun_model.h"
/***********************************************************************************************************************/
/********************************************* heat transfer models options **********************************************/
#include "phase_change_code/heat_transfer/zhang_diffusion.h"
/***********************************************************************************************************************/
#include "phase_change_code/conserving-pc.h" // This file implements momentum-conserving VOF advection of the velocity components for the two-phase Navier-Stokes solver.

#if REDUCED
# include "reduced.h"  // reduced gravity model
#endif

/******************************************* output ***************************************************/
#include "navier-stokes/perfs.h" // record velocity/pressure statistic 
#include "view.h" // post-processing
/***********************************************************************************************************************/

// grid level
#define maxlevel 6
#define level 7
#define minlevel 5

/*************** physical properties ****************************/
#define T_sat 373.15// could not be zero, because in Lee model for denominator
#define T_sup 12 // difference between wall and saturated temperatures
#define T_wall (T_sat + T_sup)
#define cp_l 10 // liquid heat capacity
#define cp_v 10 

#define lambda_l 0.0015 // liquid conductivity coefficient

#define lambda_v 0.0035
#define D_l lambda_l/(rho1*cp_l)//1.68e-7 // liquid diffusion coeffcient, used in 'diffusion.h'
#define D_v lambda_v/(rho2*cp_v)//2.05e-5
#define L_h 100 // latent heat 
#define L0 1 // domain size
#define rhoo (1/rho2 - 1/rho1) // used in mass source term 


/************************* face vector declaration, cell face ***************************/
face vector u_l[];    // liquid face velocity after add phase-change
face vector uf_new[]; // subdomain face velocity

/************************* scalar declaration, cell center ******************************/
scalar m_dot[];  // entire mass transfer flux
scalar div_pc[]; // volumetric mass source without rhoo
scalar delta_s[]; // delta_s function to make the mass flux just occurs at interface
// using general temperature gradient to calcuate mass flux
scalar m_dot_l[]; // liquid side temperature gradient
scalar m_dot_v[]; // vapor side temperature gradient
scalar div_1[];  // one-fluid velocity divergence
scalar div_2[];  // liquid velocity divergence
scalar p_new[];  // subdomain pressure    
scalar T[];      // tempearture
scalar velocity[];     // velocity magnitude without phase-change
scalar velocity_pc[];  // velocity magnitude with phase-change


// vector declaration, cell center
vector ul[];     // subdomain cell center velocity
/***********************************************************************************************/
scalar *tracers = {T}; // scalar tracers list

/************************************ boundary conditions **************************************/
// outflow boundary conditions on liquid side wall
T[right] = dirichlet(T_wall);
p[right] = dirichlet(0);
pf[right] = dirichlet(0.);
u.n[right] = neumann(0.);

T[left] = dirichlet(T_sat);
//e[left] = dirichlet(T_wall*cp_v*rho2);
pf[left] = dirichlet(0.);
p_new[left] = dirichlet(0.);
u.n[left] = dirichlet(0.);

u_l.n[top] = 0;
u_l.n[bottom] = 0;
u_l.n[left] = neumann(0.);
/***********************************************************************************************/

// automatically create folders for output files
void createFolder(const char* folder)

{
  char *str = NULL;
  str = (char *) malloc(sizeof(char) *1024);
  sprintf(str, "mkdir %s", folder); 
  system(str);   
}


int main() {
  size (L0); 
  init_grid(1<<level); 
  // fluid densities
  rho1 = 2.5;
  rho2 = 0.25;
  // fluid viscosity
  mu1 = 0.098;
  mu2 = 0.007;
  // fluid surface tension
  f.sigma = 0.059;
  TOLERANCE = 1e-8; // scalar value residual tolerance, used in 'possion.h'

  #if REDUCED
    G.x = 0;
  #endif
  
  run();
}

/***************************************** Stefan problem geometry configuration ************************************/
#define H0 0.1 // initial vapor thickness
#define stefan(x, y, H) (x-H0) // initial configuration
event init (t = 0) {

  createFolder("temperature");
  createFolder("velocity");
/***************************************************************/
  fraction (f,  stefan(x,y,H0)); // interface reconstruction, vof method

/**************** define initial temperature using volume weighted method *************************/
  foreach()
    T[] = (1-f[])*T_sat + f[]*T_wall; // initial temperature
  boundary({T});
  }


/************************ mass flux calculation *******************************/
event mass_flux(i++)
{
  T.tr_eq = T_sat; // define saturated temperature which assumes at interface
  delta_magnini(f,delta_s); // delta_s function calculation

  sun_model(T,f,div_pc,L_h);// directly obtain volumetric mass source term

  scalar ff[];
  foreach()
    ff[] = clamp(ff[],0,1);
  boundary({ff});
  volume_diffusion(ff,f); // volume fraction diffusion whith constructs one uniform domain to be capatible with AMR

// compute velocity magnitude
  foreach()
    {
      velocity[] = sqrt(sq(u.x[]) + sq(u.y[]));
      velocity_pc[] = sqrt(sq(ul.x[]) + sq(ul.y[]));
    }
  boundary({velocity,velocity_pc});
}

/********************************** vof advection with mass source **********************/
event vof(i++)
{
/******************* step-1: construct entire divergence-free domain *******************/  
  scalar divv[];
  foreach()
  {
    divv[] = div_pc[]*rhoo;
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

/******************* step-2: shift interface accounting for phase-change *******************/  
  foreach()
  {
/* Note that VOF function is clipped when the interface displacement extends beyond the cell boundary.
  * This could be solved by addig the clipped value to neighbouring cells, which is automatically solved by Basilisk vof solver*/
  if(f[]>1e-12&&f[]<1-1e-12)
    {
// first method: from PLIC method
      coord n = interface_normal( point, f); // interface normal
      double alpha = plane_alpha(f[],n); // distance from original point to interface 
      //alpha -= m_dot[]*dt/(Delta*rho1)*sqrt(sq(n.x)+sq(n.y)); // phase-change shifting distance 
      alpha -= div_pc[]*dt/rho1*sqrt(sq(n.x)+sq(n.y)); // or use volumetric mass source in vof reconstruction directly
      f[] = line_area(n.x,n.y,alpha); // cell volume fraction for n+1

// second method: from interface advection
    //f[] -= div_pc[]*dt/rho1;
    }
  }
  boundary({f});
}
/******************************************** energy diffusion ********************************/
event tracer_diffusion(i++)
{
/*
  The basic idea is to give a dirichlet boundary condition at interface (saturated temperature condition). Several models can be chosen here.
*/
  foreach()
  {
    T.rho = rho1*f[]+rho2*(1-f[]);
    T.cp = cp_l*f[]+cp_v*(1-f[]);
    T.rho_cp = f[]*rho1*cp_l + (1-f[])*rho2*cp_v;
    T.lambda = lambda_l*f[]+lambda_v*(1-f[]);
    T.D = D_v;
  }


  face vector diffusion_coef[];
  foreach_face()
    diffusion_coef.x[] = fm.x[]*T.D;
  boundary((scalar *){diffusion_coef});

  scalar heat_s[];
  foreach()
    heat_s[] = div_pc[]*L_h/T.rho_cp;
  boundary({heat_s});
/*************************************************heat transfer models**********************************************************/
  zhang_diffusion_liquid(T,f);
/***********************************************************************************************************/
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
        av.x[] -= alpha.x[]/fm.x[]*((1/rho2-1/rho1)*sq(m_dot[])*(f[]-f[-1])/Delta); //recoil pressure due to mass transfer
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


/****************************************** analytical solution **********************************************/
// Newton interpolation method
#define TOLER 1e-6
double coeff()
{
  double x = 1, x_old; // guess value
  double f, df;
  do
  {
    x_old = x;
    f = exp(sq(x))*erf(x)*(x - (T_wall - T_sat)*cp_v*lambda_l*sqrt(D_v)*exp(-sq(x)*(sq(rho2)*D_v)/(sq(rho1)*D_l))/
    (L_h*lambda_v*(sqrt(pi*D_l))*erfc(x*(rho2*sqrt(D_v)/(rho1*sqrt(D_v)))))) - cp_v*(T_wall - T_sat)/(L_h*sqrt(pi)); 
    df = (5497*x*exp(239*sq(x)/5000000)/(1026250000*erfc(0.007*x)) - 161*exp(-239*sq(x)/5000000)/(205250*sqrt(pi)*sq(erfc(0.007*x))) + 1)*exp(sq(x))*erf(x)
    + 2*x*(x - 46*(exp(-239*sq(x)/5000000))/(821*erfc(0.007*x)))*exp(sq(x))*erf(x) + 2*(x -  46*(exp(-239*sq(x)/5000000))/(821*erfc(0.007*x)))/(sqrt(pi));
    x = x_old - f/df;
  }while (fabs(x - x_old) > TOLER);
  return x;
}

// interface position
double exact(double t)
{
  double a = coeff();
  double delta_x = 2*a*sqrt(D_l*t);
  return delta_x;
}

// analytical temperature

double temper(double t, double x)
{

  double a = coeff();
  double temp_x = T_wall - ((T_sat - T_wall))/erf(a)*erf(x/(2*sqrt(t*D_v)));

  return temp_x;
}
event temperature_analy(t = 0.29; t += 1)
{
  
  char *name = NULL;
  name = (char *) malloc(sizeof(char) * 256);
  sprintf (name, "temperature/temperature-analy-%.1f",t);
  FILE*fp1 = fopen (name, "w");

  double xx = exact(t);
  for (double x = 0; x <= xx; x += 0.00001)
      fprintf (fp1, " %g %g\n", x, temper(t,x)); 
  for (double x = xx; x <= L0; x += 0.001)
      fprintf (fp1, " %g %g\n", x, 373.15); 
  fclose(fp1);
}

// analytical velocity
double velo(double t)
{
  double a = coeff();
  double velo_x = a*sqrt(D_v/t)*(1-rho2/rho1);
  return velo_x;
}
/*************************************** output ************************************************/
void mg_print (mgstats mg)
{
  if (mg.i > 0 && mg.resa > 0.)
    printf ("%d %g %g %g %d ", mg.i, mg.resb, mg.resa,
	    mg.resb > 0 ? exp (log (mg.resb/mg.resa)/mg.i) : 0,
	    mg.nrelax);
}

event logfile (t = 0; t <= 100; t += 1) { 
  double delta_interface = exact(t);
  double delta_velocity = velo(t);
  double xb = 0., vx = 0.,vy = 0., sb = 0.,yb = 0., nu = 0., mass = 0, sbb = 0;
  foreach(reduction(+:xb) reduction(+:vx) reduction(+:sb) reduction(+:yb) reduction(+:vy) reduction(+:nu) reduction(+:mass) reduction(+:sbb)) {
    double dv = (1. - f[])*dv();
    double dvv = f[]*dv();
    xb += x*dv;//interface position
    yb += y*dv;
    vx += u.x[]*dv;
    vy += u.y[]*dv;
    sb += dv; 
    sbb += dvv;
    mass += ul.x[]*dvv;
  }
  printf (" %g %g %g %g %g %g %g %g %g %g %g %g", 
	  t, sb, delta_interface, 2*xb/sb,delta_velocity,mass/sbb,  yb/sb,vx/sb, vy/sb,dt, perf.t, perf.speed);
  mg_print (mgp);
  mg_print (mgpf);
  mg_print (mgu);
  putchar ('\n');
  fflush (stdout);
}


event snap (t+=10)
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

event temperature_value (t+=10){
  char *name = NULL;
  name = (char *) malloc(sizeof(char) * 256);
  sprintf (name, "temperature/temperature-%.1f",t);
  FILE*fp1 = fopen (name, "w");

  for (double x = 0; x < L0; x += 0.0001)
      fprintf (fp1, " %g %g\n", x,
          interpolate (T, x, 0.005)); // generate velocity vector field, at least 12 grids space in liquid thin film(\delta)----reference from "Co-current flow effects on a rising Taylor bubble"
      fclose(fp1);
}

event velocity_magnitude (t+=10){

  char *name = NULL;
  name = (char *) malloc(sizeof(char) * 256);
  sprintf (name, "velocity/velocity-%.1f",t);
  FILE*fp1 = fopen (name, "w");

  for (double x = 0; x < L0; x += 0.001)
      fprintf (fp1, " %g %g %g\n", x,
          interpolate (u.x, x, 0.005),interpolate (u.y, x, 0.005)); // generate velocity vector field, at least 12 grids space in liquid thin film(\delta)----reference from "Co-current flow effects on a rising Taylor bubble"
      fclose(fp1);
}

event movies (t = 0; t <= 100; t += 1) {
  
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
set key right bottom
p [0:10][0:0.002]'in-pos' every 20 u ($1+0.29):2 w p pt 'o'  t 'Basilisk-64',\
'out' u 1:3 w l lw 2 lc rgb 'red' t 'Analytical solution'
*/