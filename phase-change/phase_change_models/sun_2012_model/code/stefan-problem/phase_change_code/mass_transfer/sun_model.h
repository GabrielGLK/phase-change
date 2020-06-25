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
          f_v[] += (dd[1] - dd[-1])/(2*Delta);
          gd[] += (tr[1] - tr[-1])/(2*Delta);
        }
    }
  boundary({f_v,gd});
/* //the segregation method, to avoid interface smearing while declaring mass transfer
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
  

  vector gtr[];
  foreach()
    {
      gtr.x[] = 0;
      foreach_dimension()
        gtr.x[] += (face_value(tr,1) - face_value(tr,0))/Delta;//face_value(a,i) for calculate face value
    }
  boundary((scalar*){gtr});

  scalar dd[];
  foreach()
    dd[] = 1 - f[];
  boundary({dd});

  scalar f_v[];
  foreach()
    {
      f_v[] = 0;
      foreach_dimension()
        f_v[] += (face_value(dd,1) - face_value(dd,0))/Delta;//((dd[1] - dd[])/Delta + (dd[] - dd[-1])/Delta)/2
    }
  boundary({f_v});

  // just simple use face temperature to get temperature gradient
  foreach()
  {
    temp[] = 0;
    foreach_dimension()
      temp[] += 2*0.025*(gtr.x[]+gtr.x[1])*f_v[]/L_h; // Eq.(20), 0.025 is the unsaturated thermo conductivity
  }
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

    if (interfacial(point,f)||interfacial(neighborp(-1), f)) {
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
      
      if (nf.x > 0.) {
        grad_t.x[] = (fabs(nf.x)*(gtr.x[1,0]) + fabs(nf.y)*(nf.y >= 0. ? gtr.x[1, 1] : gtr.x[1, -1]));
      }
      else if (nf.x < 0.) {
        grad_t.x[] = (fabs(nf.x)*gtr.x[-1,0]+ fabs(nf.y)*(nf.y >= 0. ? gtr.x[-1, 1] : gtr.x[-1, -1]));
      }
    }
  } 
  boundary((scalar *){grad_t});

  scalar dd[];
  foreach()
    dd[] = 1 - f[];
  boundary({dd});

  face vector f_v[];
  foreach()
    f_v.x[] = (dd[] - dd[-1])/Delta;
  boundary((scalar *){f_v});

  double mm;
  foreach()
    {
      mm = 0;
      foreach_dimension()
        mm += 0.025*grad_t.x[]*f_v.x[]/L_h;// the mass flux should be some value, we can put it in the cenered cell 
      temp[] = mm;// transfer the neigboring cell tempearture gradient to interfacial cell center
    }
  boundary({temp});
}