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