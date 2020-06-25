
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
    b[] = -1;
    c[] = 0;
  }
boundary({mass_source, b , c});

foreach_face()
  mass_diffusion_coef.x[] = fm.x[]*sq(1*Delta);// control the mass diffusion amplitude
boundary((scalar *){mass_diffusion_coef});

return diffusion (tr, dt, D = mass_diffusion_coef, r = mass_source, beta = b, theta = c);
}

mgstats mass_poisson (struct Mass_Diffusion p){

  scalar m = p.m, tr = p.tr;
  face vector mass_diffusion_coef[];
  scalar mm[],c[];
  foreach()
    mm[] = m[];
  boundary ({tr,mm,c});

foreach_face()
  mass_diffusion_coef.x[] = fm.x[]*sq(2*Delta);// control the mass diffusion amplitude
boundary((scalar *){mass_diffusion_coef});

//return poisson(tr,mm,mass_diffusion_coef,lambda = b, tolerance = 1e-3, nrelax = 2);
return poisson (a = tr, b = mm, alpha = mass_diffusion_coef);
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