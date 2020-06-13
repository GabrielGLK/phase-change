vector hh[];
scalar h_x[];
void flux_diffusion (scalar c, face vector kappa, face vector F, scalar dc) {

  boundary ({c});
  scalar dd[];
  for(scalar f in interfaces)
    {
      foreach()
        dd[] = 1 - f[];
      heights (dd, hh);
    }
  boundary({dd});
  
  foreach()
    {
      h_x[] = Delta*height(hh.x[]);
      //h_x[] = clamp(height(h_x[]),1e-3,1);
    }
  boundary({h_x});
  
  foreach_face()
    { 
      coord n  = interface_normal( point, f);
      if(n.x<0&&n.y==0)
      {
        if(f[]>=0.5)
          {        
            F.x[] = kappa.x[]*(c.tr_eq - c[])/(Delta - height(hh.x[])*Delta ); 
            F.x[1] = kappa.x[1]*(c[1] - c[])/Delta; 
          }
        else if(f[]<0.5)
          {
            F.x[1] = kappa.x[1]*(c[] - c.tr_eq)/(height(hh.x[])*Delta ); 
            F.x[] = kappa.x[]*(c[] - c[-1])/Delta; 
          }
      }

      
    }
  boundary_flux ({F});  //flux on levels
  foreach() {
    dc[] = 0;
    foreach_dimension() //Rotates over the dimensions
      {
        if(f[]>=0.5)
          dc[] +=  (F.x[1] - F.x[])/(0.5*((Delta - height(hh.x[])*Delta ) + Delta)); 
        else if(f[]<0.5)
          dc[] +=  (F.x[1] - F.x[])/(0.5*((height(hh.x[])*Delta) + Delta));
      } 
  }
  boundary({dc});
}

void advance (scalar c, scalar dc, double dt, scalar r) {
  foreach()
    c[] += dt*(dc[] + r[]);
}


void diffusion_1 (scalar c, double dt, face vector kappa, scalar r) {
  face vector F[];
  scalar dc[];
  flux_diffusion (c, kappa, F, dc);
  advance (c, dc, dt,r);
}


void diffusion_midpoint (scalar c, double dt, face vector kappa, scalar r) {
  face vector F[];
  scalar dc[], c_temp[];
  foreach()
    c_temp[] = c[];                 //create a scratch
  flux_diffusion(c_temp, kappa, F, dc); 
  advance (c_temp, dc, dt/2.,r);     //advance to the mid point.
  flux_diffusion(c_temp, kappa, F, dc); //Re-estimate the fluxes at the mid point,
  advance (c, dc, dt,r);             //update the solution
}