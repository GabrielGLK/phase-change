/****************** Tanasawa-model****************************/
void tanasawa_model (scalar tr, scalar f, scalar m_dot,  double L_h)
{
  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary({f});
  //double a = 0.2;// because we don't know the liquid molecular weight, so we give the guess of mass transfer intensity multipy liquid molecular weight
  double b = 0.05;// change this value for different phase-change problems
  foreach()
  {
    if(interfacial(point,f))
      {
        if(f[] != f[1])
          m_dot[] = cm[]*(2*b/(2-b))*sqrt(18/(2*pi*8.314))*(rho2*L_h*(tr[] - tr.tr_eq))/(pow(tr.tr_eq,3/2));// molecular weight is for water
        else
          m_dot[] = 0;
      }
    }
  boundary({m_dot});
}