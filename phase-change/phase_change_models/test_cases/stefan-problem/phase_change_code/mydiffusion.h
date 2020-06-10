vector hh[];
void flux_diffusion (scalar c, (const) face vector kappa, face vector F, scalar dc) {

  boundary ({c});

  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f});
  heights (f, hh);
  
  foreach_face()
    { 
      clamp(height(hh.x[]),1e-3,1);
      if(f[]>=0.5)
        {        
          F.x[1] = kappa.x[1]*(c.tr_eq - c[])/(Delta*height(hh.x[])); //2nd order accurate gradient *estimation*
          F.x[] = kappa.x[]*(c[] - c[-1])/Delta; //2nd order accurate gradient *estimation*
        }
      else if(f[]<0.5)
      {
          F.x[1] = kappa.x[1]*(c[] - c.tr_eq)/(Delta*height(hh.x[])); //2nd order accurate gradient *estimation*
          F.x[] = kappa.x[]*(c[] - c[-1])/Delta; //2nd order accurate gradient *estimation*
        }
    }
  boundary_flux ({F});  //flux on levels
  foreach() {
    dc[] = 0;
    foreach_dimension() //Rotates over the dimensions
      dc[] +=  (F.x[1] - F.x[])/(0.5*(height(hh.x[]) + Delta)); 
  }
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