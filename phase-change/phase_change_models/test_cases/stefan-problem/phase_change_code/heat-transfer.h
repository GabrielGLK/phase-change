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
  scalar dirichlet_source_term[], volumic_metric[], b[], convec[];
  face vector diffusion_coef[];
  foreach()
    {
      if(interfacial(point,f))
        convec[] = m[]*tr.cp*tr[];
      volumic_metric[] = cm[];
      //dirichlet_source_term[] = m[]*L_h/(tr.cp*tr.rho) - m[]*2000*(tr[] - tr.tr_eq)/(tr.cp*tr.rho);
      dirichlet_source_term[] = m[]*L_h/(tr.rho*tr.cp) + convec[]/(tr.rho*tr.cp);
      //dirichlet_source_term[] = m[]*L_h/(tr.rho*tr.cp);
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
  //int level = p.level;
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
  scalar heat_source;
  double dt;
  double L_h;
};

mgstats rattner_diffusion (struct Rattner_Diffusion p) {
  
  scalar tr = p.tr, f = p.f, heat_source = p.heat_source;
  double dt = p.dt, L_h = p.L_h;
  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f, tr});

  scalar volumic_metric[], dirichlet_feedback_term[];
  face vector diffusion_coef[];
  foreach() {
    double mass_min_1 = min(tr.rho_cp*tr.tr_eq/(dt*tr.rho*tr.cp),f[]*958*L_h/(tr.rho*tr.cp*dt));//the parameters for stefan problem
    double mass_min_2 = min(mass_min_1, L_h/(tr.rho*tr.cp*dt)/(1/0.6 - 1/958));
    volumic_metric[] = cm[];
    heat_source[] = -mass_min_2;
    dirichlet_feedback_term[] = -1/(tr.rho*tr.cp*dt);
  }
  foreach_face()
    diffusion_coef.x[] = fm.x[]*tr.D;

  boundary({volumic_metric, heat_source, dirichlet_feedback_term, diffusion_coef});
  
  return diffusion (tr, dt, D = diffusion_coef, r = heat_source,
			              beta = dirichlet_feedback_term, theta = volumic_metric);
}
/****************** equilibrium saturated temperature on interface *************************/
/*************** Malan's method to keep the interface saturated temmperature ****************/
// to be added