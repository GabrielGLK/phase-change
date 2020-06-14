struct Zhang_diffusion {
  scalar tr;
  scalar f;
  scalar m;
  double L_h;
};

mgstats zhang_diffusion_vapor (struct Zhang_diffusion p)
{

  scalar tr = p.tr, f = p.f, m = p.m;
  double L_h = p.L_h;
  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f, tr});

  scalar dd[],t_av[]; // volume fraction of vapor
  double T_p;
  vertex scalar t_cor[];
  foreach_vertex()
    t_cor[] = (tr[] + tr[-1,0] + tr[0,-1] + tr[-1,-1])/4;
  boundary({t_cor});
  foreach()
    { 
      dd[] = 1 - f[];// volume fraction of vapor
      coord n  = interface_normal( point, f) , p, q; // interface normal, barycenter in vapor potion, negative interface normal
      q.x = -n.x, q.y = -n.y; 
      double alpha  = plane_alpha (dd[], q); // cell center to interface but for vapor
      double alpha_l  = plane_alpha (f[], n); // cell center to interface but for vapor
      line_center (q, alpha,dd[],&p); // barycenter in vapor portion 
      coord pp,qq;
      //line_length_center (n, alpha_l, &pp); // interface center coordinate
      double distance_1;
      pp.x = (alpha_l + 0.5*(n.x + n.y))*n.x + 0.5*(n.x + n.y)*n.x;
      pp.y = (alpha_l + 0.5*(n.x + n.y))*n.y + 0.5*(n.x + n.y)*n.y;
      distance_1 = fabs(alpha_l + 0.5*(n.x + n.y) + 0.5*(n.x + n.y));
      qq.x = pp.x + distance_1*n.x;
      qq.y = pp.y + distance_1*n.y;
      
      // if the ghost cell in mixed cell
      if(interfacial(point,f)&&alpha_l>0.5*(n.x + n.y))
        {
          t_av[] = (t_cor[] + t_cor[1,0] + t_cor[0,1] + t_cor[1,1])/4;
          T_p =  interpolate_2(point,t_av,qq);
          //T_p = interpolate_1 (point, T, p);// we calculate the barycenter of interfacial vapor portion temperature
          //tr[] = T_p*(1-f[]) + f[]*T_sat;
        }

      // ghost cell method for different situations of interface normal directions, combined with heat_source() function
      if(n.x==0&&n.y>0)
        {
          if((f[0,-1]==1)||(interfacial(point,f)&&alpha_l<0.5*(n.x + n.y)))
            tr[] = 2*tr.tr_eq - T_p;
        }
      else if(n.x==0&&n.y<0)
        {
          if((interfacial(neighborp(0,1),f)&&f[]==1)||(interfacial(point,f)&&alpha_l<0.5*(n.x + n.y)))
            tr[] = 2*tr.tr_eq - T_p;
        }
      else if(n.x>0&&n.y==0)
        {
          if((interfacial(neighborp(-1,0),f)&&f[]==1)||(interfacial(point,f)&&alpha_l<0.5*(n.x + n.y)))
            tr[] = 2*tr.tr_eq - T_p;
        }
      else if(n.x<0&&n.y==0)
        {
          if((interfacial(neighborp(1,0),f)&&f[]==1)||(interfacial(point,f)&&alpha_l<0.5*(n.x + n.y)))
            tr[1,0] = 2*tr.tr_eq - T_p;
        }

      else if(n.x<0&&n.y<0)
        {
          if((interfacial(neighborp(1,1),f)&&f[]==1)||(interfacial(point,f)&&alpha_l<0.5*(n.x + n.y)))
            tr[] = tr.tr_eq - T_p;
        }
      else if(n.x<0&&n.y>0)
        {
          if((interfacial(neighborp(1,-1),f)&&f[]==1)||(interfacial(point,f)&&alpha_l<0.5*(n.x + n.y)))
            tr[] = 2*tr.tr_eq - T_p;
        }
      else if(n.x>0&&n.y<0)
        {
          if((interfacial(neighborp(-1,1),f)&&f[]==1)||(interfacial(point,f)&&alpha_l<0.5*(n.x + n.y)))
            tr[] = 2*tr.tr_eq - T_p;
        }
      else if(n.x>0&&n.y>0)
        {
          if((interfacial(neighborp(-1,-1),f)&&f[]==1)||(interfacial(point,f)&&alpha_l<0.5*(n.x + n.y)))
            tr[] = 2*tr.tr_eq - T_p;
        }
    }


  face vector diffusion_coef[];
  foreach_face()
    diffusion_coef.x[] = fm.x[]*tr.D;
  boundary((scalar *){diffusion_coef});
  
  return diffusion (tr, dt, D = diffusion_coef);
}