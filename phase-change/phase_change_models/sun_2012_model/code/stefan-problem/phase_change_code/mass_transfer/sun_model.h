/******************************* The sun simplified model ***************************************/
/* reference paper:
Development of a vapor–liquid phase change model for volume-of-fluid method in FLUENT
*/
void sun_model_center (scalar tr, scalar f, scalar temp, double L_h) 
{
  face vector grad_t[]; // gradient of temperature
  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f, tr});

  scalar dd[];
  foreach()
    dd[] = 1 - f[];
  boundary({dd});

  scalar f_v[],gd[];
  foreach()
    {
      f_v[] = 0;
      gd[] = 0;
      foreach_dimension()
        {
          f_v[] += sq((dd[1] - dd[-1])/(2*Delta));
          gd[] += sq((tr[1] - tr[-1])/(2*Delta));
        }
        f_v[] = sqrt(f_v[]);
        gd[] = sqrt(gd[]);
    }
  boundary({f_v,gd});
/*
  foreach()
    {
      if(interfacial(point,f)&&f[-1] == 0)
        temp[] = 2*0.025*(gd[] + gd[1])*f_v[]/L_h;
      else
        temp[] = 0;
    }
  boundary({temp});
  */
  foreach()
    temp[] = 2*0.025*(gd[] + gd[1])*f_v[]/L_h;// the mass flux should be some value, we can put it in the cenered cell 
  boundary({temp});
}

void sun_model_simple (scalar tr, scalar f, scalar temp, double L_h) 
{
  
  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f, tr});

  vector m[];
  compute_normal(f,m);
  scalar gtr[];
  scalar tt[];
  foreach()
  {
    gtr[] = 0;
    if(interfacial(point,f))
    {
      coord n = interface_normal(point,f), p, q;
      double alpha = plane_alpha (f[], n);
      //double alpha_all = alpha + 0.5*(n.x + n.y);
      line_length_center (n, alpha, &p); // interface fragment center coordinate

      //q.x = p.x  + n.x;
      //q.y = p.y + n.y;
      tt[] = interpolate_2(point,tr,p);//using interpolation temperature in interfacial cell
      //tt[-1] = interpolate_2(point,tr,q);

      foreach_dimension()
        gtr[] += sq((2*tr.tr_eq - tr[-1] - tt[])/(2*Delta)*(2/(1 + alpha + fabs(alpha))));// give a coefficient to make sure the temperature gradient consistent with physical problem, new here
        //gtr.x[] += (tr.tr_eq - tt[-1])/(2*Delta);
      gtr[] = sqrt(gtr[]);
    }
  }
  boundary({tt,gtr});

  scalar dd[];
  foreach()
    dd[] = 1 - f[];
  boundary({dd});

  scalar f_v[];
  foreach()
    {
      f_v[] = 0;
      foreach_dimension()
        f_v[] += sq((dd[1] - dd[-1])/(2*Delta));//((dd[1] - dd[])/Delta + (dd[] - dd[-1])/Delta)/2
      f_v[] = sqrt(f_v[]);
    }
  boundary({f_v});

  // just simple use face temperature to get temperature gradient
  foreach()
    temp[] = 2*0.025*gtr[]*f_v[]/L_h; // Eq.(20), 0.025 is the unsaturated thermo conductivity
  boundary({temp});
}


// reference paper: A DIRECT NUMERICAL SIMULATION FOR NUCLEATE BOILING BY THE VOSET METHOD
void sun_model_face (scalar tr, scalar f, scalar temp, double L_h) 
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
