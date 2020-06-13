/****************** rattner-model****************************/ 
// should use with ratter_diffusion() energy diffusion
// reference paper: SIMPLE MECHANISTICALLY CONSISTENT FORMULATION FOR VOLUME-OF-FLUID BASED COMPUTATIONS OF CONDENSING FLOWS
void rattner_model (scalar tr, scalar f, scalar m_dot,  double L_h, scalar h_s)
{
  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary({f});
  scalar mass_min_1[], mass_min_2[];
  foreach()
  {
  if(interfacial(point,f))
    {
      //mass_min_1[] = min(tr.rho_cp*(tr[] - tr.tr_eq)/dt,f[]*958*L_h/dt);//the parameters for stefan problem
      //mass_min_2[] = min(mass_min_1[], L_h/(dt*(1/0.6 - 1/958)));
    
      m_dot[] = -h_s[]/L_h;  // combined with rattner model
    }
  }
  boundary({m_dot});
}