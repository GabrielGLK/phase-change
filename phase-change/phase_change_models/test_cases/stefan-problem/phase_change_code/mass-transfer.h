#include "tracer-pc.h"
#include "vof-pc.h"
#include "diffusion-pc.h"
#include "mydiffusion.h"
#include "curvature.h"
#include "my_function.h"

attribute {
  double D;
  double tr_eq;
  double lambda;
  double cp;
  double rho;
  double rho_cp;
}

/*
 interpolation method to set neighboring cell of unsaturated fluid of interface to be saturated temperature
First interpolation method
*/
double interpolate_1 (Point point, scalar s, coord p)
{
  int i = sign(p.x), j = sign(p.y);
  if (f[i] && f[0,j] && f[i,j])
    // bilinear interpolation when all neighbors are defined
    return ((s[]*(1. - fabs(p.x)) + s[i]*fabs(p.x))*(1. - fabs(p.y)) + 
	    (s[0,j]*(1. - fabs(p.x)) + s[i,j]*fabs(p.x))*fabs(p.y));
  else {
    // linear interpolation with gradients biased toward the
    // cells which are defined
    double val = s[];
    foreach_dimension() {
      int i = sign(p.x);
      if (f[i])
	val += fabs(p.x)*(s[i] - s[]);
      else if (f[-i])
	val += fabs(p.x)*(s[] - s[-i]);
    }
    return val;
  }
}

// second interpolation method
static inline double interpolate_2 (Point point, scalar col,coord p) {
  struct _interpolate _r = { col, x + p.x*Delta, y + p.y*Delta, z + p.z*Delta };
  return interpolate_linear (point, _r);
}

// constant mass source scalar function
void mass_source(scalar f, scalar m)
{
  foreach()
  {
    if(f[]>1e-12&&f[]<1 - 1e-12)
      m[] = 0.05*cm[];
    else
      m[] = 0;
  }
  boundary({f,m});
}

// delta_s function
/*
The delta_s function is inspired from Magnini's thesis:
'CFD Modeling of Two-Phase Boiling Flows in the Slug Flow Regime with an Interface Capturing Technique'
page-61
*/
void delta_magnini(scalar f, scalar del) {

  foreach()  {
    f[] = clamp(f[], 0., 1.);
    double df, df2;
    df2 = 0.;
    df = 0;
    foreach_dimension () {
      df = (f[-1,1] + 2*f[0,1] + f[1,1] - f[-1,-1] - 2*f[0,-1] - f[1,-1])/(8*Delta);
      df2 += sq(df);
    }
    del[] = sqrt(df2);
  }
  boundary({del});
}
// smear volume fraction with the vertex-average of *f*
void smoother(scalar f, scalar ff) {
  #if dimension <= 2
  foreach()
    ff[] = (4.*f[] + 
	    2.*(f[0,1] + f[0,-1] + f[1,0] + f[-1,0]) +
	    f[-1,-1] + f[1,-1] + f[1,1] + f[-1,1])/16.;
#else // dimension == 3
  foreach()
    ff[] = (8.*f[] +
	    4.*(f[-1] + f[1] + f[0,1] + f[0,-1] + f[0,0,1] + f[0,0,-1]) +
	    2.*(f[-1,1] + f[-1,0,1] + f[-1,0,-1] + f[-1,-1] + 
		f[0,1,1] + f[0,1,-1] + f[0,-1,1] + f[0,-1,-1] +
		f[1,1] + f[1,0,1] + f[1,-1] + f[1,0,-1]) +
	    f[1,-1,1] + f[-1,1,1] + f[-1,1,-1] + f[1,1,1] +
	    f[1,1,-1] + f[-1,-1,-1] + f[1,-1,-1] + f[-1,-1,1])/64.;
#endif
}


// diffusion the rude mass source to the neighbouring cells
/* using this method, we don't need delts_s function at interface, as Zhang's paper did:
'Direct numerical simulations of incompressible multiphase magnetohydrodynamics with phase change'
*/
struct Mass_Diffusion{
  scalar tr;
  scalar m;
};

mgstats mass_diffusion (struct Mass_Diffusion p){

  scalar m = p.m, tr = p.tr;
  scalar mass_source[], b[], c[];
  face vector mass_diffusion_coef[];

foreach() 
  {
    mass_source[] = m[];
    b[] = 0;
    c[] = 0;
  }
boundary({mass_source, b , c});

foreach_face()
  mass_diffusion_coef.x[] = fm.x[]*sq(1.5*Delta);// control the mass diffusion amplitude
boundary((scalar *){mass_diffusion_coef});

return diffusion (tr, dt, D = mass_diffusion_coef, r = mass_source, beta = b, theta = c);
}

mgstats mass_poisson (struct Mass_Diffusion p){

  scalar m = p.m, tr = p.tr;
  face vector mass_diffusion_coef[];
  scalar mm[],b[],c[];
 foreach()
    {
      mm[] = m[];
      c[] = 0;
    }
  boundary ({tr,mm,b, c});

foreach_face()
  mass_diffusion_coef.x[] = fm.x[]*sq(2*Delta);// control the mass diffusion amplitude
boundary((scalar *){mass_diffusion_coef});

//return poisson(tr,mm,mass_diffusion_coef,lambda = b, tolerance = 1e-3, nrelax = 2);
return poisson (a = tr, b = mm, alpha = mass_diffusion_coef, lambda = c);
//return diffusion (tr, dt, D = mass_diffusion_coef, r = mm, theta = c);
}

/* if we want to use adaptive mesh refinement, the better way is to construct one scalar covering volume fraction.
This method is similar with mass diffusion.

*/
struct Volume_Diffusion{
  scalar tr;
  scalar f;
};

mgstats volume_diffusion (struct Volume_Diffusion p){

  scalar tr = p.tr, f = p.f;
  scalar mass_source[], b[], c[];
  face vector mass_diffusion_coef[];

 foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f, tr});

foreach() 
  {
    mass_source[] = f[];
    b[] = -1;
    c[] = 0;
  }
boundary({mass_source, b , c});

foreach_face()
  mass_diffusion_coef.x[] = fm.x[]*sq(1*Delta);// control the mass diffusion amplitude
boundary((scalar *){mass_diffusion_coef});

return diffusion (tr, dt, D = mass_diffusion_coef, r = mass_source, beta = b, theta = c);
}

/******************************************* mass transfer models**********************************/
/*
Here we give several different phase-change models, including calculate mass flux or volumetric mass source.

*/

/******************************* The sharp interface model ***************************************/
/*
The introduction of this model in my Github:
https://github.com/GabrielGLK/phase-change/blob/master/phase-change/README.MD
*/
void sharp_simple_model (scalar tr, scalar f, scalar temp, double L_h) 
{
  face vector v_pc[];
  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f, tr});
  
  face vector gtr[];
  foreach_face()
    gtr.x[] = (tr[] - tr[-1])/Delta;
  boundary((scalar*){gtr});

  vector n[];
  compute_normal (f, n);
  
  foreach_face() {
    v_pc.x[] = 0.;

    if (interfacial(point,f) || interfacial(neighborp(-1), f)) {
      coord nf;
      foreach_dimension()
        nf.x = 0.;
      if (interfacial(point, f)) {
        foreach_dimension()
          nf.x += n.x[];
      }
      if (interfacial(neighborp(-1), f)) {
        nf.x += n.x[-1];
        nf.y += n.y[-1];
      }
   
      double norm = 0.;
      foreach_dimension()
        norm += fabs(nf.x);
      foreach_dimension()
        nf.x /= norm;
      
      if (nf.x >= 0.) {
        v_pc.x[] = (fabs(nf.x)*gtr.x[1, 0] + fabs(nf.y)*(nf.y >= 0. ? gtr.x[1, 1] : gtr.x[1, -1]));
      }
      else if (nf.x < 0.) {
        v_pc.x[] = (fabs(nf.x)*gtr.x[-1, 0]+ fabs(nf.y)*(nf.y >= 0. ? gtr.x[-1, 1] : gtr.x[-1, -1]));
      }
    }
  } 
  boundary((scalar *){v_pc});
  scalar dd[];
  face vector f_v[];
  foreach_face()
    {
      dd[] = 1 - f[];
      f_v.x[] = (dd[1] - dd[-1])/(2*Delta);
    }
  boundary((scalar *){f_v});

  foreach()
    {
      temp[] = 0;
      foreach_dimension()
        temp[] += cm[]*2*0.025*v_pc.x[]*n.x[]/(sqrt(sq(n.x[]) + sq(n.y[]))*L_h); // unsaturated thermo-conductivity, film-boiling is 1
    }
  boundary({temp});
}
/******************************* The sun simplified model ***************************************/
void sun_model (scalar tr, scalar f, scalar temp, double L_h) 
{
  face vector v_pc[];
  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f, tr});
  
  face vector gtr[];
  foreach_face()
    gtr.x[] = (tr[] - tr[-1])/Delta;
  boundary((scalar*){gtr});

  vector n[];
  compute_normal (f, n);
  
  foreach_face() {
    v_pc.x[] = 0.;

    if (interfacial(point,f) || interfacial(neighborp(-1), f)) {
      coord nf;
      foreach_dimension()
        nf.x = 0.;
      if (interfacial(point, f)) {
        foreach_dimension()
          nf.x += n.x[];
      }
      if (interfacial(neighborp(-1), f)) {
        nf.x += n.x[-1];
        nf.y += n.y[-1];
      }
   
      double norm = 0.;
      foreach_dimension()
        norm += fabs(nf.x);
      foreach_dimension()
        nf.x /= norm;
      
      if (nf.x >= 0.) {
        v_pc.x[] = (fabs(nf.x)*gtr.x[1, 0] + fabs(nf.y)*(nf.y >= 0. ? gtr.x[1, 1] : gtr.x[1, -1]));
      }
      else if (nf.x < 0.) {
        v_pc.x[] = (fabs(nf.x)*gtr.x[-1, 0]+ fabs(nf.y)*(nf.y >= 0. ? gtr.x[-1, 1] : gtr.x[-1, -1]));
      }
    }
  } 
  boundary((scalar *){v_pc});
  scalar dd[];
  face vector f_v[];
  foreach_face()
    {
      dd[] = 1 - f[];
      f_v.x[] = (dd[1] - dd[-1])/(2*Delta);
    }
  boundary((scalar *){f_v});

  foreach()
    {
      temp[] = 0;
      foreach_dimension()
        temp[] += cm[]*2*0.025*v_pc.x[]*f_v.x[]/L_h; // unsaturated thermo-conductivity, film-boiling is 1
    }
  boundary({temp});
}

/****************** Lee-model****************************/
/*
The introduction of this model in my Github:
https://github.com/GabrielGLK/phase-change/blob/master/phase-change/README.MD
*/
void lee_model(scalar tr, scalar f, scalar m, double L_h){

  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f, tr});
  scalar r_i[], dd[];
  foreach()
  {
    m[] = 0;
    dd[] = 1 - f[];// vapor volume fraction in one cell
    coord n  = interface_normal( point, dd) , p, q;
    q.x = -n.x, q.y = -n.y;
    double alpha  = plane_alpha (dd[], q);
    line_center (q, alpha,dd[],&p);
    double T_p = 0;
    if(interfacial(point,f))
      T_p = interpolate_2 (point, tr, p); 
    if(interfacial(neighborp(-1,-1),f)&&f[]==1)
      { 
      r_i[] = 0;
      foreach_dimension()
        r_i[] += 0.68*tr.tr_eq/(L_h*(0.5+0.5*f[])*Delta*n.x*f[]*rho1);//mass transfer intensity factor expression // saturated thermo-conductivity, film-boiling is 40,but this does not work. it's good for stefan problem
      m[] = 2*cm[]*r_i[]*f[]*rho1*(T_p - tr.tr_eq)/tr.tr_eq;
      }
  }
  boundary({m});
}

/****************** Tanasawa-model****************************/
void tanasawa_model (scalar tr, scalar f, scalar m_dot,  double L_h)
{
  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary({f});
  //double a = 0.2;// because we don't know the liquid molecular weight, so we give the guess of mass transfer intensity multipy liquid molecular weight
  double b = 0.05;// change this value for different phase-change problems
  foreach()
  {
    if(interfacial(point,f))
      {
        if(f[] != f[1])
          m_dot[] = cm[]*(2*b/(2-b))*sqrt(18/(2*pi*8.314))*(rho2*L_h*(tr[] - tr.tr_eq))/(pow(tr.tr_eq,3/2));// molecular weight is for water
        else
          m_dot[] = 0;
      }
    }
  boundary({m_dot});
}

/****************** rattner-model****************************/ 
// should use with ratter_diffusion() energy diffusion
void rattner_model (scalar tr, scalar f, scalar m_dot,  double L_h, scalar h_s)
{
  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary({f});
  scalar mass_min_1[], mass_min_2[];
  foreach()
  {
  if(interfacial(point,f))
    {
      //mass_min_1[] = min(tr.rho_cp*(tr[] - tr.tr_eq)/dt,f[]*958*L_h/dt);//the parameters for stefan problem
      //mass_min_2[] = min(mass_min_1[], L_h/(dt*(1/0.6 - 1/958)));
    
      m_dot[] = -h_s[]/L_h;  // combined with rattner model
    }
  }
  boundary({m_dot});
}

/****************** Zhang-model and extended models****************************/
// combined with mass_diffusion() function
void zhang_model (scalar tr, scalar f, scalar m, double L_h)
{
  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary({f});

  scalar f_v[];
  foreach()
    {
      m[] = 0;
      f_v[] = 1 - f[];
      coord n  = interface_normal( point, f_v) , s;
      double alpha  = plane_alpha (f_v[], n);
      plane_area_center(n, alpha, &s);
      double delta_d = 1.85;// according to the paper, this value is between 1 and 2
      double xc = s.x-n.x*delta_d;
      double yc = s.y-n.y*delta_d;
      coord q;
      q.x = xc, q.y = yc;
      double t_center = interpolate_2(point,tr,q);
      if(interfacial(point,f))
        m[] = tr.lambda*(t_center - tr.tr_eq)/(L_h*delta_d*Delta);
    }
  boundary({m});
}

// Simplified mass transfer inspired from Jie zhang's paper
// not use interpolation method, using averaging method
void zhang_model_1 (scalar tr, scalar f, scalar m, double L_h)
{
  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary({f});

  scalar f_v[];
  vertex scalar t_cor[];
  foreach_vertex()
    t_cor[] = (tr[] + tr[-1,0] + tr[0,-1] + tr[-1,-1])/4;
  boundary({t_cor});
  scalar t_av[];
  foreach()
    {
      m[] = 0;
      f_v[] = 1 - f[];
      coord n  = interface_normal( point, f) , s;
      double alpha  = plane_alpha (f[], n);
      plane_area_center(n, alpha, &s);
      double xc = 0, yc = 0;
      double delta_d = 1.5;// according to the paper, this value is between 1 and 2
      if(interfacial(point,f))
        {
          xc += s.x*Delta-n.x*delta_d*Delta;
          yc += s.y*Delta-n.y*delta_d*Delta;
          t_av[] = (t_cor[] + t_cor[1,0] + t_cor[0,1] + t_cor[1,1])/4;
          double d1 = fabs(n.x*xc + n.y*yc - alpha*Delta - Delta*0.5*(n.x+n.y))/sqrt(sq(n.x) + sq(n.y));// distance for temperature gradient
          m[] =tr.lambda*(t_av[] - tr.tr_eq)/(L_h*d1);
        }
    }
  boundary({m});
}


// Simplified mass transfer inspired from Jie zhang's paper
// using barycenter value in interfacial cell
void zhang_model_2 (scalar tr, scalar f, scalar m, double L_h)
{
  foreach()
    {
      if(f[] == 1)
        tr[] = tr.tr_eq;
      f[] = clamp(f[], 0., 1.);
    }
  boundary({tr,f});

  scalar f_v[];
  vertex scalar t_cor[];
  foreach_vertex()
    t_cor[] = (tr[] + tr[-1,0] + tr[0,-1] + tr[-1,-1])/4;
  boundary({t_cor});
  face vector s[];
  foreach()
    {
      m[] = 0;
      f_v[] = 1 - f[];
      coord n  = facet_normal( point, f, s) , p, q, k;
      q.x = -n.x, q.y = -n.y;
      double alpha_new  = plane_alpha (f_v[], q);
      double alpha  = plane_alpha (f[], n); 
      line_center (q, alpha_new,f_v[],&p);// barycenter of vapor portion in interfacial cell
      line_center (n, alpha,f[],&k);// barycenter of liquid portion in interfacial cell
      double xc = 0, yc = 0;
      double xk = 0, yk = 0;
      if(f[]<1-1e-12&&f[]>1e-12)
        {
          xc += p.x*Delta;
          yc += p.y*Delta;
          xk += k.x*Delta;
          yk += k.y*Delta;
          //coord pp;
          //pp.x = xc, pp.y = yc;// mirror point in local corner neighboring cell
          double d1 = fabs(q.x*xc + q.y*yc - alpha_new*Delta - Delta*0.5*(q.x+q.y))/sqrt(sq(q.x) + sq(q.y));// distance for temperature gradient
          //double d1 = sqrt(sq(xk - xc) + sq(yk - yc));
          //double tt = interpolate_1 (point, tr, pp); // using interpolation scheme for temperpoint p
          double tt = interpolate_2(point,tr,p);
          m[] = tr.lambda*(tt - tr.tr_eq)/(L_h*d1);// new method here
        }
    }
  boundary({m});
}

/****************** Leon Malan's model****************************/
void malan_model (scalar tr, scalar f, scalar m, double L_h)
{
  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary({f});

  double d[3][4];
  double eps[3][4];
  double gamma[3][4];
  double phi[3][4];
  double mm[3][4];
  //double m_sum[4];
  double tt[3][4]; 
  
  foreach()
  {
    coord p[3] = {{-1,1},{-1,0},{0,1}};
    coord q[3] = {{-1,-1},{-1,0},{0,-1}};
    coord r[3] = {{1,1},{1,0},{0,1}};
    coord s[3] = {{1,-1},{1,0},{0,-1}};
    coord n = interface_normal(point,f);
    normalize(&n);
    double alpha  = plane_alpha (f[], n) + 0.5*(n.x + n.y);  
    if(f[]>1e-12&&f[]<1-1e-12)
    {
      tt[0][0] = tr[-1,1];
      tt[1][0] = tr[-1,0];
      tt[2][0] = tr[0,1];
      tt[0][1] = tr[-1,-1];
      tt[1][1] = tr[-1,0];
      tt[2][1] = tr[0,-1];
      tt[0][2] = tr[1,1];
      tt[1][2] = tr[1,0];
      tt[2][2] = tr[0,1];
      tt[0][3] = tr[1,-1];
      tt[1][3] = tr[1,0];
      tt[2][3] = tr[0,-1];

      if(n.x<0&&n.y>0)
        {
          d[0][0] = fabs(-alpha + p[0].x*n.x + p[0].y*n.y)/sqrt(sq(n.x) + sq(n.y));
          eps[0][0] = fabs(p[0].x*n.x + p[0].y*n.y);
          gamma[0][0] = eps[0][0]/(sq(p[0].x) + sq(p[0].y));
          phi[0][0] = gamma[0][0]/(gamma[0][0] + gamma[1][0] + gamma[2][0]);
          mm[0][0] = phi[0][0]*(tt[0][0] - tr.tr_eq)/(Delta*d[0][0]);

          d[1][0] = fabs(-alpha + p[1].x*n.x + p[1].y*n.y)/sqrt(sq(n.x) + sq(n.y));
          eps[1][0] = fabs(p[1].x*n.x + p[1].y*n.y);
          gamma[1][0] = eps[1][0]/(sq(p[1].x) + sq(p[1].y));
          phi[1][0] = gamma[1][0]/(gamma[0][0] + gamma[1][0] + gamma[2][0]);
          mm[1][0] = phi[1][0]*(tt[1][0] - tr.tr_eq)/(Delta*d[1][0]);

          d[2][0] = fabs(-alpha + p[2].x*n.x + p[2].y*n.y)/sqrt(sq(n.x) + sq(n.y));
          eps[2][0] = fabs(p[2].x*n.x + p[2].y*n.y);
          gamma[2][0] = eps[2][0]/(sq(p[2].x) + sq(p[2].y));
          phi[2][0] = gamma[2][0]/(gamma[0][0] + gamma[1][0] + gamma[2][0]);
          mm[2][0] = phi[2][0]*(tt[2][0] - tr.tr_eq)/(Delta*d[2][0]);

          m[] = mm[0][0] + mm[1][0] + mm[2][0];//for stefan problem, if we want to check the direction of interface normal vector, we can make this term disappear
        }
        
      if(n.x<=0&&n.y<=0)
        {
          d[0][1] = fabs(-alpha + q[0].x*n.x + q[0].y*n.y)/sqrt(sq(n.x) + sq(n.y));
          eps[0][1] = fabs(q[0].x*n.x + q[0].y*n.y);
          gamma[0][1] = eps[0][1]/(sq(q[0].x) + sq(q[0].y));
          phi[0][1] = gamma[0][1]/(gamma[0][1] + gamma[1][1] + gamma[2][1]);
          mm[0][1] = phi[0][1]*(tt[0][1] - tr.tr_eq)/(Delta*d[0][1]);
        
          d[1][1] = fabs(-alpha + q[1].x*n.x + q[1].y*n.y)/sqrt(sq(n.x) + sq(n.y));
          eps[1][1] = fabs(q[1].x*n.x + q[1].y*n.y);
          gamma[1][1] = eps[1][1]/(sq(q[1].x) + sq(q[1].y));
          phi[1][1] = gamma[1][1]/(gamma[0][1] + gamma[1][1] + gamma[2][1]);
          mm[1][1] = phi[1][1]*(tt[1][1] - tr.tr_eq)/(Delta*d[1][1]);
        
          d[2][1] = fabs(-alpha + q[2].x*n.x + q[2].y*n.y)/sqrt(sq(n.x) + sq(n.y));
          eps[2][1] = fabs(q[2].x*n.x + q[2].y*n.y);
          gamma[2][1] = eps[2][1]/(sq(q[2].x) + sq(q[2].y));
          phi[2][1] = gamma[2][1]/(gamma[0][1] + gamma[1][1] + gamma[2][1]);
          mm[2][1] = phi[2][1]*(tt[2][1] - tr.tr_eq)/(Delta*d[2][1]);

          m[] = mm[0][1] + mm[1][1] + mm[2][1];
        }
        
      
      if(n.x>=0&&n.y>=0)
        {
          d[0][2] = fabs(-alpha + r[0].x*n.x + r[0].y*n.y)/sqrt(sq(n.x) + sq(n.y));
          eps[0][2] = fabs(r[0].x*n.x + r[0].y*n.y);
          gamma[0][2] = eps[0][2]/(sq(r[0].x) + sq(r[0].y));
          phi[0][2] = gamma[0][2]/(gamma[0][2] + gamma[1][2] + gamma[2][2]);
          mm[0][2] = phi[0][2]*(tt[0][2] - tr.tr_eq)/(Delta*d[0][2]);
        
          d[1][2] = fabs(-alpha + r[1].x*n.x + r[1].y*n.y)/sqrt(sq(n.x) + sq(n.y));
          eps[1][2] = fabs(r[1].x*n.x + r[1].y*n.y);
          gamma[1][2] = eps[1][2]/(sq(r[1].x) + sq(r[1].y));
          phi[1][2] = gamma[1][2]/(gamma[0][2] + gamma[1][2] + gamma[2][2]);
          mm[1][2] = phi[1][2]*(tt[1][2] - tr.tr_eq)/(Delta*d[1][2]);
        
          d[2][2] = fabs(-alpha + r[2].x*n.x + r[2].y*n.y)/sqrt(sq(n.x) + sq(n.y));
          eps[2][2] = fabs(r[2].x*n.x + r[2].y*n.y);
          gamma[2][2] = eps[2][2]/(sq(r[2].x) + sq(r[2].y));
          phi[2][2] = gamma[2][2]/(gamma[0][2] + gamma[1][2] + gamma[2][2]);
          mm[2][2] = phi[2][2]*(tt[2][2] - tr.tr_eq)/(Delta*d[2][2]);

          m[] = mm[0][2] + mm[1][2] + mm[2][2];
        }

      if(n.x>0&&n.y<0)
        {
          d[0][3] = fabs(-alpha + s[0].x*n.x + s[0].y*n.y)/sqrt(sq(n.x) + sq(n.y));
          eps[0][3] = fabs(s[0].x*n.x + s[0].y*n.y);
          gamma[0][3] = eps[0][3]/(sq(s[0].x) + sq(s[0].y));
          phi[0][3] = gamma[0][3]/(gamma[0][3] + gamma[1][3] + gamma[2][3]);
          mm[0][3] = phi[0][3]*(tt[0][3] - tr.tr_eq)/(Delta*d[0][3]);
       
          d[1][3] = fabs(-alpha + s[1].x*n.x + s[1].y*n.y)/sqrt(sq(n.x) + sq(n.y));
          eps[1][3] = fabs(s[1].x*n.x + s[1].y*n.y);
          gamma[1][3] = eps[1][3]/(sq(s[1].x) + sq(s[1].y));
          phi[1][3] = gamma[1][3]/(gamma[0][3] + gamma[1][3] + gamma[2][3]);
          mm[1][3] = phi[1][3]*(tt[1][3] - tr.tr_eq)/(Delta*d[1][3]);
        
          d[2][3] = fabs(-alpha + s[2].x*n.x + s[2].y*n.y)/sqrt(sq(n.x) + sq(n.y));
          eps[2][3] = fabs(s[2].x*n.x + s[2].y*n.y);
          gamma[2][3] = eps[2][3]/(sq(s[2].x) + sq(s[2].y));
          phi[2][3] = gamma[2][3]/(gamma[0][3] + gamma[1][3] + gamma[2][3]);
          mm[2][3] = phi[2][3]*(tt[2][3] - tr.tr_eq)/(Delta*d[2][3]);

          m[] = mm[0][3] + mm[1][3] + mm[2][3];
        }
        
    }
    m[] *= tr.lambda/L_h;
  }
  boundary({m});
}

/*************************  Hard mass-transfer model ***********************/
// maybe this not works because I didn't check the guess value here !!!!! not checked!!!!!
void mass_hard (scalar tr, scalar f, scalar m, scalar m_0, scalar m_1, double L_h)
{

  scalar delta_s[];
  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary({f});

   /* step-1: initial mass flux
  */
  scalar J_evp[]; // evaporation mass flux
  /*
  the delta_s function defined here is different from div_pc[], we obtain this from 
  'CFD Modeling of Two-Phase Boiling Flows in the Slug Flow Regime with an Interface Capturing Technique' p-60
  */
  double R_int = 6.38e-8;

  // delta_s function (|\nabla c|)

  //delta_magnini(f,delta_s);
   face vector d_u[];
  // delta_s function (|\nabla c|)
  /* gradient of volume fraction*/
  foreach_face()
    //d_u.x[] = (f[-1,1] + 2*f[0,1] + f[1,1] - f[-1,-1] - 2*f[0,-1] - f[1,-1])/(2*Delta);
    d_u.x[] = (f[] - f[-1])/Delta;
  boundary((scalar *){d_u});
/* calculation of (|\nabla c|) */
  foreach()
    delta_s[] = sqrt(sq(d_u.x[]*fm.x[]) + sq(d_u.y[]*fm.x[]));
  boundary({delta_s});


  double N1 = 0, N_all = 0, N = 0;
  foreach(reduction(+:N_all) reduction(+:N1))
    {
      N1 += delta_s[]*dv();
      N_all += delta_s[]*f[]*dv();

      if (N_all > 1e-99)
        N = N1/N_all;
    }
    
  foreach()
    {
      J_evp[] = (tr[] - tr.tr_eq)/(R_int*L_h); // should be (tr[] - tr.sat)/R_int -> this is what we need to consider, not use gradient of temperature, much easier for mass flux
      m_0[] = N*J_evp[]*f[]*delta_s[];
    }
  boundary({J_evp, m_0});

/*
  step-2: smearing the mass flux 
*/

mass_diffusion(m_1,m_0);

/*
  step-3: From m[] the complete source-term incorporating evap-
orative mass transfer between the phases is determined as
*/
  double N_e = 0, N_v_all = 0, N_l_all = 0, N_v = 0, N_l = 0;
  foreach(reduction(+:N_e) reduction(+:N_v_all) reduction(+:N_l_all))
  {
    N_e += m_0[]*dv();
    N_v_all += (1-f[])*m_0[]*dv();
    N_l_all += f[]*m_0[]*dv();
    if(N_v_all > 1e-99)
      N_v = N_e/N_v_all;
    if(N_l_all > 1e-99)
      N_l = N_e/N_l_all;
  }

  foreach()
  {
    if(f[] < 1e-12)
      m[] = N_v*(1-f[])*m_1[];
    else if(f[] > 1 - 1e-12)
      m[] = -N_l*f[]*m_1[];
    else
      m[] = 0;
  }
  boundary({m});

}