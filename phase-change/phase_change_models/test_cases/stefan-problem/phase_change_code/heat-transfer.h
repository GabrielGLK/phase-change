/********************************** energy equation ***********************************/
#include "diffusion-pc.h"
#include "curvature.h"
#include "my_function.h"

#define F_ERR 1e-10
attribute {
  double D;
  double tr_eq;
  double lambda;
  double cp;
  double rho;
  double rho_cp;
}

/****** zhang's ghost cell method to keep interface saturated temperature *****************/
struct Heat_Source {
  scalar tr;
  scalar f;
  scalar m;
  double L_h;
};

mgstats heat_source (struct Heat_Source p)
{

  scalar tr = p.tr, f = p.f, m = p.m;
  double L_h = p.L_h;
  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f, tr});
  scalar dirichlet_source_term[], volumic_metric[], b[], convec[];
  face vector diffusion_coef[];
  foreach()
    {
      if(interfacial(point,f))
        convec[] = m[]*tr.cp*tr[];
      volumic_metric[] = cm[];
      //dirichlet_source_term[] = m[]*L_h/(tr.cp*tr.rho) - m[]*2000*(tr[] - tr.tr_eq)/(tr.cp*tr.rho);
      dirichlet_source_term[] = -m[]*L_h/(tr.rho*tr.cp);
      //dirichlet_source_term[] = m[]*L_h/(tr.rho*tr.cp);
      b[] = 0;
    }
  boundary({volumic_metric,dirichlet_source_term,b});
  
  foreach_face()
    diffusion_coef.x[] = fm.x[]*tr.D;
  boundary((scalar *){diffusion_coef});
  
  return diffusion (tr, dt, D = diffusion_coef, r = dirichlet_source_term, beta = b, theta = volumic_metric);
}