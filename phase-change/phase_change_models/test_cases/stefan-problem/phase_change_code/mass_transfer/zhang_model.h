/****************** Zhang-model and extended models****************************/
// combined with mass_diffusion() function
void zhang_model (scalar tr, scalar f, scalar m, double L_h)
{
  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary({f});

  scalar f_v[];
  foreach()
    {
      m[] = 0;
      coord n  = interface_normal( point, f) , s;
      double alpha  = plane_alpha (f[], n);
      plane_area_center(n, alpha, &s);
      double delta_d = 1.8;// according to the paper, this value is between 1 and 2
      double xc = s.x+n.x*delta_d;
      double yc = s.y+n.y*delta_d;
      coord q;
      q.x = xc, q.y = yc;

      double t_center = interpolate_1(point,tr,q);
      if(interfacial(point,f))
        m[] = 2*0.025*(t_center - tr.tr_eq)/(L_h*delta_d*Delta);
        
    }
  boundary({m});
}

// Simplified mass transfer inspired from Jie zhang's paper
// not use interpolation method, using averaging method
void zhang_model_1 (scalar tr, scalar f, scalar m, double L_h)
{
  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary({f});

  scalar f_v[];
  vertex scalar t_cor[];
  foreach_vertex()
    t_cor[] = (tr[] + tr[-1,0] + tr[0,-1] + tr[-1,-1])/4;
  boundary({t_cor});
  scalar t_av[];
  foreach()
    {
      m[] = 0;
      f_v[] = 1 - f[];
      coord n  = interface_normal( point, f) , s;
      double alpha  = plane_alpha (f[], n);
      plane_area_center(n, alpha, &s);
      double xc = 0, yc = 0;
      double delta_d = 1.5;// according to the paper, this value is between 1 and 2
      if(interfacial(point,f))
        {
          xc += s.x*Delta-n.x*delta_d*Delta;
          yc += s.y*Delta-n.y*delta_d*Delta;
          t_av[] = (t_cor[] + t_cor[1,0] + t_cor[0,1] + t_cor[1,1])/4;
          double d1 = fabs(n.x*xc + n.y*yc - alpha*Delta - Delta*0.5*(n.x+n.y))/sqrt(sq(n.x) + sq(n.y));// distance for temperature gradient
          m[] =tr.lambda*(t_av[] - tr.tr_eq)/(L_h*d1);
        }
    }
  boundary({m});
}


// Simplified mass transfer inspired from Jie zhang's paper
// using barycenter value in interfacial cell
void zhang_model_2 (scalar tr, scalar f, scalar m, double L_h)
{
  foreach()
    {
      if(f[] == 1)
        tr[] = tr.tr_eq;
      f[] = clamp(f[], 0., 1.);
    }
  boundary({tr,f});

  scalar f_v[];
  vertex scalar t_cor[];
  foreach_vertex()
    t_cor[] = (tr[] + tr[-1,0] + tr[0,-1] + tr[-1,-1])/4;
  boundary({t_cor});
  face vector s[];
  foreach()
    {
      m[] = 0;
      f_v[] = 1 - f[];
      coord n  = facet_normal( point, f, s) , p, q, k;
      q.x = -n.x, q.y = -n.y;
      double alpha_new  = plane_alpha (f_v[], q);
      double alpha  = plane_alpha (f[], n); 
      line_center (q, alpha_new,f_v[],&p);// barycenter of vapor portion in interfacial cell
      line_center (n, alpha,f[],&k);// barycenter of liquid portion in interfacial cell
      double xc = 0, yc = 0;
      double xk = 0, yk = 0;
      if(f[]<1-1e-12&&f[]>1e-12)
        {
          xc += p.x*Delta;
          yc += p.y*Delta;
          xk += k.x*Delta;
          yk += k.y*Delta;
          //coord pp;
          //pp.x = xc, pp.y = yc;// mirror point in local corner neighboring cell
          double d1 = fabs(q.x*xc + q.y*yc - alpha_new*Delta - Delta*0.5*(q.x+q.y))/sqrt(sq(q.x) + sq(q.y));// distance for temperature gradient
          //double d1 = sqrt(sq(xk - xc) + sq(yk - yc));
          //double tt = interpolate_1 (point, tr, pp); // using interpolation scheme for temperpoint p
          double tt = interpolate_2(point,tr,p);
          m[] = tr.lambda*(tt - tr.tr_eq)/(L_h*d1);// new method here
        }
    }
  boundary({m});
}
