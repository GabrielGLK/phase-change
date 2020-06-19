/******************************* The sun simplified model ***************************************/
/* reference paper:
Development of a vapor–liquid phase change model for volume-of-fluid method in FLUENT
*/
void sun_model (scalar tr, scalar f, scalar temp, double L_h) 
{
  face vector grad_t[]; // gradient of temperature
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
    grad_t.x[] = 0.;

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
        grad_t.x[] = (fabs(nf.x)*gtr.x[1,0] + fabs(nf.y)*(nf.y >= 0. ? gtr.x[1, 1] : gtr.x[1, -1]));
      }
      else if (nf.x < 0.) {
        grad_t.x[] = (fabs(nf.x)*gtr.x[-1,0]+ fabs(nf.y)*(nf.y >= 0. ? gtr.x[-1, 1] : gtr.x[-1, -1]));
      }
    }
  } 
  boundary((scalar *){grad_t});
  scalar dd[];
  face vector f_v[];
  foreach()
    dd[] = 1 - f[];
  boundary({dd});

  foreach_face()
    f_v.x[] = (dd[] - dd[-1])/Delta;
  boundary((scalar *){f_v});

  foreach()
    {
      temp[] = 0;
      foreach_dimension()
        temp[] += 0.025*grad_t.x[]*f_v.x[]/L_h; // unsaturated thermo-conductivity, film-boiling is 1
    }
  boundary({temp});
}

void sharp_simple_model_vapor (scalar tr, scalar f, scalar temp, double L_h) 
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

    if (interfacial(point,f) ) {
      coord nf;
      foreach_dimension()
        nf.x = 0.;
      if (interfacial(point, f)) {
        foreach_dimension()
          nf.x += n.x[];
      }
   
      double norm = 0.;
      foreach_dimension()
        norm += fabs(nf.x);
      foreach_dimension()
        nf.x /= norm;
      
      if (nf.x > 0.) {
        v_pc.x[] = (fabs(nf.x)*gtr.x[1,0] + fabs(nf.y)*(nf.y >= 0. ? gtr.x[1, 1] : gtr.x[1, -1]));
      }
      else if (nf.x < 0.) {
        v_pc.x[] = (fabs(nf.x)*gtr.x[]+ fabs(nf.y)*(nf.y >= 0. ? gtr.x[-1, 1] : gtr.x[-1, -1]));
      }
    }
  } 
  boundary((scalar *){v_pc});
  foreach()
    {
      temp[] = 0;
      foreach_dimension()
        temp[] += 0.025*v_pc.x[]*n.x[]/L_h; // unsaturated thermo-conductivity, film-boiling is 1
    }
  boundary({temp});
}

void sharp_simple_model_liquid (scalar tr, scalar f, scalar temp, double L_h) 
{
  face vector v_pc[];
  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f, tr});
  
  face vector gtr[];
  foreach_face()
    gtr.x[] = (tr[] - tr[-1])/Delta;
  boundary((scalar*){gtr});
  scalar dd[];
  foreach()
    dd[] = 1 - f[];
  boundary({dd});
  vector n[];
  compute_normal (dd, n);
  
  
  foreach_face() {
    v_pc.x[] = 0.;

    if (interfacial(point,f) ) {
      coord nf;
      foreach_dimension()
        nf.x = 0.;
      if (interfacial(point, f)) {
        foreach_dimension()
          nf.x += n.x[];
      }
   
      double norm = 0.;
      foreach_dimension()
        norm += fabs(nf.x);
      foreach_dimension()
        nf.x /= norm;
      
      if (nf.x >= 0.) {
        v_pc.x[] = (fabs(nf.x)*gtr.x[1,0] + fabs(nf.y)*(nf.y >= 0. ? gtr.x[1, 1] : gtr.x[1, -1]));
      }
      else if (nf.x < 0.) {
        v_pc.x[] = (fabs(nf.x)*gtr.x[-1,0]+ fabs(nf.y)*(nf.y >= 0. ? gtr.x[-1, 1] : gtr.x[-1, -1]));
      }
    }
  } 
  boundary((scalar *){v_pc});
  foreach()
    {
      temp[] = 0;
      foreach_dimension()
        temp[] += 0.68*v_pc.x[]*n.x[]/L_h; // unsaturated thermo-conductivity, film-boiling is 1
    }
  boundary({temp});
}