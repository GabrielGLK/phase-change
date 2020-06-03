/********************************** energy equation ***********************************/
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
  scalar dd[]; // volume fraction of vapor
  double T_p;
  foreach()
    { 
      /* step-1: determine ghost cell in liquid
      1. the mixed cell whose center is outside the interface (use the distance from cell center to interface to determine this)
      2. the pure liquid cell which has a corner neighbor belonging to the mixed cells but not belonging to the ghost cells
      */
      dd[] = 1 - f[];// volume fraction of vapor
      coord n  = interface_normal( point, f) , p, q; // interface normal, barycenter in vapor potion, negative interface normal
      q.x = -n.x, q.y = -n.y; 
      double alpha  = plane_alpha (dd[], q); // cell center to interface but for vapor
      double alpha_l  = plane_alpha (f[], n); // cell center to interface but for vapor
      line_center (q, alpha,dd[],&p); // barycenter in vapor portion 
      // if the ghost cell in mixed cell
      if(interfacial(point,f))
        //T_p = interpolate_1 (point, T, p);
        T_p =  interpolate_2(point,tr,p);
      if(interfacial(neighborp(-1,-1),f)&&f[]==1)
        tr[] = 2*tr.tr_eq - T_p; 
    }
    scalar dirichlet_source_term[], volumic_metric[], b[];
    face vector diffusion_coef[];
    foreach()
      {
        volumic_metric[] = cm[];
        //dirichlet_source_term[] = m[]*L_h/(tr.cp*tr.rho) - m[]*2000*(tr[] - tr.tr_eq)/(tr.cp*tr.rho);
        dirichlet_source_term[] = m[]*L_h/(tr.rho_cp);
        b[] = 0;
      }
    boundary({volumic_metric,dirichlet_source_term,b});
    foreach_face()
      diffusion_coef.x[] = fm.x[]*tr.D;
    boundary((scalar *){diffusion_coef});
  
  return diffusion (tr, dt, D = diffusion_coef, r = dirichlet_source_term, beta = b, theta = volumic_metric);
}


/****************** equilibrium saturated temperature on interface *************************/
/*************** Quentin's method to keep the interface saturated temmperature ****************/
struct Dirichlet_Diffusion {
  // mandatory
  scalar tr;
  scalar f;
  int level;
  double dt;
  double time_factor;
  // optional
  scalar tr_op; // default uniform 1.
};

mgstats dirichlet_diffusion (struct Dirichlet_Diffusion p) {
  
  scalar tr = p.tr, f = p.f, tr_op = automatic (p.tr_op);
  int level = p.level;
  double time_factor = p.time_factor;

  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f, tr});
  if (p.tr_op.i)
    boundary ({tr_op});

  scalar volumic_metric[], dirichlet_source_term[], dirichlet_feedback_term[];
  face vector diffusion_coef[];

  foreach() {
    volumic_metric[] = cm[];
    dirichlet_source_term[] = cm[]*tr.tr_eq*tr.D*time_factor*(1-f[])*sq((1<<level)/L0);
    dirichlet_feedback_term[] = - cm[]*tr.D*time_factor*(1-f[])*sq((1<<level)/L0);
  }
  foreach_face()
    diffusion_coef.x[] = fm.x[]*tr.D;

  boundary({volumic_metric, dirichlet_source_term, dirichlet_feedback_term, diffusion_coef});
  
  return diffusion (tr, dt, D = diffusion_coef, r = dirichlet_source_term,
			              beta = dirichlet_feedback_term, theta = volumic_metric);
}

/****************** equilibrium saturated temperature on interface *************************/
/*************** Rattner's method to keep the interface saturated temmperature ****************/
// one simple article used here: NUMERICAL ANALYSIS ON HEAT TRANSFER CHARACTERISTICS OF LOOPED MINICHANNEL USING PHASE-CHANGE VOF METHOD
struct Rattner_Diffusion {
  // mandatory
  scalar tr;
  scalar f;
  scalar m;
  double dt;
  double L_h;
};

mgstats rattner_diffusion (struct Rattner_Diffusion p) {
  
  scalar tr = p.tr, f = p.f, m = p.m;
  double dt = p.dt, L_h = p.L_h;
  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f, tr});

  scalar volumic_metric[], dirichlet_source_term[], dirichlet_feedback_term[];
  face vector diffusion_coef[];

  foreach() {
    volumic_metric[] = cm[];
    dirichlet_source_term[] = tr.tr_eq/dt;
    dirichlet_feedback_term[] = -1/dt;
  }
  foreach_face()
    diffusion_coef.x[] = fm.x[]*tr.D;

  boundary({volumic_metric, dirichlet_source_term, dirichlet_feedback_term, diffusion_coef});
  
  return diffusion (tr, dt, D = diffusion_coef, r = dirichlet_source_term,
			              beta = dirichlet_feedback_term, theta = volumic_metric);
}
/****************** equilibrium saturated temperature on interface *************************/
/*************** Malan's method to keep the interface saturated temmperature ****************/
// to be added