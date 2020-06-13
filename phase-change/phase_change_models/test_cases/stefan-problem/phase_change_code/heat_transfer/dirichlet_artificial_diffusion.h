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
