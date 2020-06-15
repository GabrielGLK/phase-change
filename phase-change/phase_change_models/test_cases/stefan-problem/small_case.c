#define ADAPT 0
#define REDUCED 1
/*******************************mass-transfer model options **************************/
#define MALAN_MASS 0 // Leon MALAN's paper mass transfer model
#define LEE_MASS 0 // Lee model
#define TANASAWA_MASS 0 // tanasawa model
#define SUN_MASS 1 // Dongliang Sun's paper model
#define RATTNER_MASS 0 // Rattner's paper model
#define ZHANG_MASS 0 // Jie zhang's paper model and modification models
/***********************************************************************************************************************/
/*******************************heat-transfer model options **************************/
#define ZHANG_DIFFUSION 1// Jie zhang's paper heat transfer model
#define NO_FLUX_DIFFUSION 0 // Quentin's sandbox for no diffusion flux excess the interface from vapor
#define ORIGIN_DIFFUSION 1 // original diffusion equations
#define MALAN_DIFFUSION 0 // Leon MALAN's paper heat ransfer model
#define RATTNER_DIFFUSION 0
#define DIRICHLET_ARTIFICIAL_DIFFUSION 0 //give an artifical temperature boundary condition at interface
/***********************************************************************************************************************/
// Solving advection equaiton for energy diffusion, if the tracer is temperature, we use 'tracers'. Otherwise we use f.tracers for thermal energy advection
#define THERMAL_ENERGY 0
#define TEMPERATURE_TRACER 1
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
#if MALAN_MASS 
#include "phase_change_code/mass_transfer/malan_model.h"
#endif

#if LEE_MASS
#include "phase_change_code/mass_transfer/lee_model.h"
#endif

#if TANASAWA_MASS
#include "phase_change_code/mass_transfer/tanasawa_model.h"
#endif

#if SUN_MASS
#include "phase_change_code/mass_transfer/sun_model.h"
#endif

#if RATTNER_MASS
#include "phase_change_code/mass_transfer/rattner_model.h"
#endif

#if ZHANG_MASS
#include "phase_change_code/mass_transfer/zhang_model.h"
#endif
/***********************************************************************************************************************/
/********************************************* heat transfer models options **********************************************/
#if MALAN_DIFFUSION
#include "phase_change_code/heat_transfer/malan_diffusion.h"
#endif

#if NO_FLUX_DIFFUSION
#include "phase_change_code/heat_transfer/no_flux_diffusion.h"
#endif

#if RATTNER_DIFFUSION
#include "phase_change_code/heat_transfer/rattner_diffusion.h"
#endif

#if DIRICHLET_ARTIFICIAL_DIFFUSION
#include "phase_change_code/heat_transfer/dirichlet_artificial_diffusion.h"
#endif

#if ZHANG_DIFFUSION
#include "phase_change_code/heat_transfer/zhang_diffusion.h"
#endif
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
#define T_sup 10 // difference between wall and saturated temperatures
#define T_wall (T_sat + T_sup)
#define cp_l 4216 // liquid heat capacity
#define cp_v 2030 

#if SUN_MASS
#define lambda_l 0. // liquid conductivity coefficient
#else
#define lambda_l 0.68
#endif

#define lambda_v 0.025
#define D_l lambda_l/(rho1*cp_l)//1.68e-7 // liquid diffusion coeffcient, used in 'diffusion.h'
#define D_v lambda_v/(rho2*cp_v)//2.05e-5
#define L_h 2.26e6 // latent heat 
#define L0 0.01 // domain size
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
scalar e[];      // thermal energy 
scalar T[];      // tempearture
scalar velocity[];     // velocity magnitude without phase-change
scalar velocity_pc[];  // velocity magnitude with phase-change


// vector declaration, cell center
vector ul[];     // subdomain cell center velocity
/***********************************************************************************************/
#if TEMPERATURE_TRACER
scalar *tracers = {T}; // scalar tracers list
#endif 

#if THERMAL_ENERGY
scalar *tracers = NULL; //no give scalar tracers, due to we use f.tracers to instead of this
#endif

/************************************ boundary conditions **************************************/
// outflow boundary conditions on liquid side wall
T[right] = dirichlet(T_sat);
p[right] = dirichlet(0);
pf[right] = dirichlet(0.);
u.n[right] = neumann(0.);

T[left] = dirichlet(T_wall);
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
  rho1 = 958;
  rho2 = 0.6;
  // fluid viscosity
  mu1 = 2.82e-4;
  mu2 = 1.23e-5;
  // fluid surface tension
  f.sigma = 0.059;
  TOLERANCE = 1e-4; // scalar value residual tolerance, used in 'possion.h'

  #if REDUCED
    G.x = 0;
  #endif
  
  run();
}

/***************************************** Stefan problem geometry configuration ************************************/
#define H0 3.225e-4 // initial vapor thickness
#define stefan(x, y, H) (x-H0) // initial configuration
event init (t = 0) {

  createFolder("temperature");
  createFolder("velocity");
/***************************************************************/
  fraction (f,  stefan(x,y,H0)); // interface reconstruction, vof method

/**************** define initial temperature using volume weighted method *************************/
  scalar Cp[];
  e.inverse = true;
  foreach()
    {
      Cp[] = (rho1*f[]*cp_l + rho2*(1-f[])*cp_v); // volumetric heat capacity
      T[] = (1-f[])*(T_wall - T_sup/H0*x) + f[]*T_sat; // initial temperature
      e[] = Cp[]*T[]; // initial thermal energy
    }
  boundary({Cp,T,e});
  }


/************************ mass flux calculation *******************************/
event mass_flux(i++)
{
  T.tr_eq = T_sat; // define saturated temperature which assumes at interface
  delta_magnini(f,delta_s); // delta_s function calculation

  #if MALAN_MASS // using Malan's paper method
  //sharp_simple_model_vapor(T,f,m_dot,L_h);//simplified fluid normal temperature
  malan_model_liquid(T,f,m_dot_l,L_h);//liquid normal temperature gradient
  malan_model_vapor(T,f,m_dot_v,L_h); // vapor normal temperature gradient
// mass flux and volumetric mass flux without rhoo
  foreach()
  {
    m_dot[] = (m_dot_v[] + m_dot_l[]);
    div_pc[] = m_dot[]*delta_s[];
  }
  boundary({m_dot,div_pc});
  #endif

  #if LEE_MASS
  lee_model(T,f,m_dot,L_h);
  foreach()
    div_pc[] = m_dot[]*delta_s[];
  boundary({div_pc});
  #endif

  #if TANASAWA_MASS
  tanasawa_model(T,f,m_dot,L_h);
  foreach()
    div_pc[] = m_dot[]*delta_s[];
  boundary({div_pc});
  #endif

  #if SUN_MASS
  sun_model(T,f,div_pc,L_h);// directly obtain volumetric mass source term
  #endif

  #if RATTNER_MASS
  rattner_model(T,f,m_dot,L_h);
  foreach()
    div_pc[] = m_dot[]*delta_s[];
  boundary({div_pc});
  #endif

  #if ZHANG_MASS
  zhang_model_liquid(T,f,m_dot_l,L_h);
  zhang_model_vapor(T,f,m_dot_v,L_h);
  //zhang_model_1(T,f,m_dot,L_h);
  //zhang_model_2(T,f,m_dot,L_h);
  foreach()
  {
    m_dot[] = (m_dot_v[] - m_dot_l[]);
    div_pc[] = m_dot[]*delta_s[];
  }
  boundary({m_dot,div_pc});
  #endif
  
  /*
  scalar temp[];
  foreach()
    temp[] = div_pc[];
  boundary({temp});
  mass_diffusion(div_pc,temp);
  */

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
#if THERMAL_ENERGY
static scalar * interfaces_save = NULL;
#endif
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

//We set the list of interfaces to NULL so that the default vof() event does nothing (otherwise we would transport f twice).
  #if THERMAL_ENERGY
  e.gradient = minmod2;; 
  f.tracers = {e}; // using same approach of vof for thermal energy, multi-dimentional advection 
  boundary({f, e});
  
  vof_advection ({f}, i);
  interfaces_save = interfaces;
  interfaces = NULL;
  #endif
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

event tracer_advection (i++) {
  #if THERMAL_ENERGY
  interfaces = interfaces_save;
  #endif
}
scalar tt[];
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
    //T.D = D_V*(1-f[]) + D_L*f[];
    #if SUN_MASS
    T.D = D_v;
    #else
    T.D = D_v*(1-f[]) + D_l*f[];
    #endif
  }

   // we have the intermediate thermal energy in the cell (after the sweep) and can calculate the intermediate temperature
  #if THERMAL_ENERGY
  foreach()
    T[] = e[]/T.rho_cp;
  boundary({T});
  #endif

  face vector diffusion_coef[];
  foreach_face()
    diffusion_coef.x[] = fm.x[]*T.D;
  boundary((scalar *){diffusion_coef});

  scalar heat_s[];
  foreach()
    heat_s[] = div_pc[]*L_h/T.rho_cp;
  boundary({heat_s});
/*************************************************heat transfer models**********************************************************/
  #if ZHANG_DIFFUSION
  //foreach()
    //T[] = clamp(T[],T.tr_eq,T_wall);
  //boundary({T});
  zhang_diffusion_liquid(T,f, tt);
  //zhang_diffusion_vapor(T,f, tt);
  #endif

  
  #if MALAN_DIFFUSION
  diffusion_1(T,dt,diffusion_coef,heat_s);
//diffusion_midpoint(T,dt,diffusion_coef,heat_source); //2nd order accuracy
  #endif

  #if NO_FLUX_DIFFUSION
  //foreach()
    //T[] = clamp(T[],T.tr_eq,T_wall);
  //boundary({T});
  no_flux_diffusion(T,f,dt);
  #endif

  #if ORIGIN_DIFFUSION
  heat_source (T, f, div_pc);// add heat source term in diffusion equation
  #endif

  #if RATTNER_DIFFUSION
  foreach()
    T[] = clamp(T[],T.tr_eq,T_wall);
  boundary({T});
  rattner_diffusion (T, f, div_pc, dt, L_h);// add heat source term in diffusion equation
  #endif

  #if DIRICHLET_ARTIFICIAL_DIFFUSION
  foreach()
    T[] = clamp(T[],T.tr_eq,T_wall);
  boundary({T});
  dirichlet_diffusion (T, f, level, dt, 10);
  #endif

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
    f = x*exp(sq(x))*erf(x) - cp_v*T_sup/(L_h*sqrt(pi)); 
    df = 2*sq(x)*exp(sq(x))*erf(x) + exp(sq(x))*erf(x) + 2*x/sqrt(pi);
    x = x_old - f/df;
  }while (fabs(x - x_old) > TOLER);
  return x;
}

double exact(double t)
{
  double a = coeff();
  double delta_x = 2*a*sqrt(D_v*t);
  return delta_x;
}

/*************************************** output ************************************************/
// a small event to check height function value in each direction
#if MALAN_DIFFUSION
event timeseries (t += 0.01 ) 
{
  double maxy = - HUGE,maxx = - HUGE;;
  foreach()
  if ((hh.y[] != nodata) && (hh.x[] != nodata)) 
  {
    double yi = y + height(hh.y[])*Delta;
    double xi = x + height(hh.x[])*Delta;
    if (yi > maxy)
        maxy = yi;
    if (xi > maxx)
        maxx = xi;
  }
  char s[80];
  sprintf (s, "hmax");
  static FILE * fp0 = fopen (s, "w");
  fprintf (fp0, "%g %g %g\n", t, maxx, maxy);
  fflush (fp0);
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
  double delta_interface = exact(t);

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


event snap (t+=0.01)
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
set key right bottom
p [0:10][0:0.002]'in-pos' every 20 u ($1+0.29):2 w p pt 'o'  t 'Basilisk-64',\
'out' u 1:3 w l lw 2 lc rgb 'red' t 'Analytical solution'
*/