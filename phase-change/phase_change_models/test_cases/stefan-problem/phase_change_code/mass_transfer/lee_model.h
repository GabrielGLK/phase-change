/****************** Lee-model****************************/
/*
The introduction of this model in my Github:
https://github.com/GabrielGLK/phase-change/blob/master/phase-change/README.MD
*/
void lee_model(scalar tr, scalar f, scalar m, double L_h){

  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f, tr});
  scalar r_i[], dd[];
  foreach()
  {
    m[] = 0;
    dd[] = 1 - f[];// vapor volume fraction in one cell
    coord n  = interface_normal( point, dd) , p, q;
    q.x = -n.x, q.y = -n.y;
    double alpha  = plane_alpha (dd[], q);
    line_center (q, alpha,dd[],&p);
    double T_p = 0;
    if(interfacial(point,f))
      T_p = interpolate_2 (point, tr, p); 
    r_i[] = 0.025*tr.tr_eq/(L_h*(0.5+0.5*f[])*Delta*f[]*rho1);//mass transfer intensity factor expression // saturated thermo-conductivity, film-boiling is 40,but this does not work. it's good for stefan problem
    if(interfacial(point,f))
      m[] = cm[]*r_i[]*f[]*rho1*(tr[-1] - tr.tr_eq)/tr.tr_eq;
  }
  boundary({m});
}