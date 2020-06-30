

// reference paper: A DIRECT NUMERICAL SIMULATION FOR NUCLEATE BOILING BY THE VOSET METHOD
void ling_model (scalar tr, scalar f, scalar temp, double L_h) 
{
  face vector grad_t[]; // gradient of temperature
  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f, tr});
  
  
  scalar t_inter_l[]; // interpolation temperature in interfacial cell
  scalar t_inter_v[]; 
  scalar A[];
  scalar B[];
  scalar gtr_l[],gtr_v[];
  //scalar gtr_l[]; // scalar centered value of temperature gradient

  scalar gtr_s[];
  scalar t_e[];

  vector n[];
  compute_normal (f, n);

  foreach() { 
      A[] = 0;
      B[] = 0;
      temp[] = 0;
      t_inter_l[] = 0;
      t_inter_v[] = 0;
      gtr_l[] = 0;
      gtr_v[] = 0;
      coord nf = normal(point,f),p;
      double alpha = plane_alpha(f[],nf); //interface distance to centered point in one cell
      double alpha_1 = plane_alpha(f[],nf) + 0.5*(nf.x + nf.y); 
      double delta_l = line_length_center (nf, alpha, &p)*Delta;
      if (interfacial(point,f))
      { 
        if(alpha > 0)
          {
            A[] = fabs(alpha*Delta)*(fabs(nf.x)/Delta*(nf.x < 0?tr[1]:tr[-1]) + fabs(nf.y)/Delta*(nf.y >= 0? tr[0,-1]:tr[0,1])) + tr.tr_eq;
            B[] = fabs(alpha*Delta)*(fabs(nf.x)/Delta + fabs(nf.y)/Delta) + 1;
            t_inter_l[] = A[]/B[];
            gtr_l[] = (t_inter_l[] - tr.tr_eq)/fabs(alpha*Delta);
            //gtr_l[] = ((nf.x < 0?tr[1]:tr[-1]) - t_inter_l[])/Delta*fabs(nf.x) + ((nf.y >= 0? tr[0,-1]:tr[0,1]) - t_inter_l[])/Delta*fabs(nf.y);
            temp[] = (0.68*gtr_l[])/L_h*delta_l/dv(); // liquid side
            }
          else
          {
            //a better accuracy
             /* coord p = {-1,0}, q = {-1,1}, s = {0,1};// coordinates of cells around interface
              double alpha_a = fabs(nf.y*s.y - alpha_1)/(sqrt(sq(nf.x) + sq(nf.y)));
              double alpha_b = fabs(nf.x*q.x + nf.y*q.y - alpha_1)/(sqrt(sq(nf.x) + sq(nf.y)));
              double alpha_c = fabs(nf.x*p.x - alpha_1)/(sqrt(sq(nf.x) + sq(nf.y)));
              if(fabs(nf.x) <= 1/2)// interface vector for x direction 
                gtr_v[] = ((nf.x<0?tr[0,1]:tr[0,-1]) - tr.tr_eq)/fabs(alpha_a*Delta);
              else if(fabs(nf.x) <= sqrt(3)/2&&fabs(nf.x)>1/2)
                gtr_v[] = ((nf.x<0&&nf.y>0)?tr[-1,1]:(nf.x<0&&nf.y<0)?tr[-1,-1]:(nf.x>0&&nf.y>0)?tr[1,1]:tr[1,-1]
                - tr.tr_eq)/fabs(alpha_b*Delta);
              else if(fabs(nf.x) > sqrt(3)/2)
                gtr_v[] = ((nf.x<0?tr[-1]:tr[1]) - tr.tr_eq)/fabs(alpha_c*Delta);// stefan problem using this term
             */ 
            A[] = fabs(alpha*Delta)*(fabs(nf.x)/Delta * (nf.x < 0?tr[-1]:tr[1]) + fabs(nf.y)/Delta * (nf.y >= 0? tr[0,1]:tr[0,-1])) + tr.tr_eq;
            B[] = fabs(alpha*Delta)*(fabs(nf.x)/Delta + fabs(nf.y)/Delta) + 1;
            t_inter_v[] = A[]/B[];
            gtr_v[] = (t_inter_v[] - tr.tr_eq)/fabs(alpha*Delta);
            //gtr_v[] = ((nf.x < 0?tr[-1]:tr[1]) - t_inter_l[])/Delta*fabs(nf.x) + ((nf.y > 0? tr[0,-1]:tr[0,1]) - t_inter_l[])/Delta*fabs(nf.y);
            temp[] = (0.025*gtr_v[])/L_h*delta_l/dv(); // vapor side
            }
    }

  }
  boundary({A, B, t_inter_v,t_inter_l,gtr_v,gtr_l,gtr_s, temp});
}

/* In this modified version, we ignore the effect of liquid side, because the liquid side is saturated condition,
Therefore, we just need to calculate vapor side using two differernt types as Ling's paper.
*/

void ling_model_modify (scalar tr, scalar f, scalar temp, double L_h) 
{
  face vector grad_t[]; // gradient of temperature
  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f, tr});
  
  
  scalar t_inter_l[]; // interpolation temperature in interfacial cell
  scalar t_inter_v[]; 
  scalar A[];
  scalar B[];
  scalar gtr_l[],gtr_v[];
  //scalar gtr_l[]; // scalar centered value of temperature gradient

  scalar gtr_s[];
  scalar t_e[];

  vector n[];
  compute_normal (f, n);

  foreach() { 
      A[] = 0;
      B[] = 0;
      temp[] = 0;
      t_inter_l[] = 0;
      t_inter_v[] = 0;
      gtr_l[] = 0;
      gtr_v[] = 0;
      coord nf = normal(point,f),p;
      double alpha = plane_alpha(f[],nf); //interface distance to centered point in one cell
      double alpha_1 = plane_alpha(f[],nf) + 0.5*(nf.x + nf.y); 
      double delta_l = line_length_center (nf, alpha, &p)*Delta;
      if (interfacial(point,f))
      { 
        A[] = fabs(alpha*Delta)*(fabs(nf.x)/Delta * (nf.x < 0?tr[-1]:tr[1]) + fabs(nf.y)/Delta * (nf.y >= 0? tr[0,1]:tr[0,-1])) + tr.tr_eq;
        B[] = fabs(alpha*Delta)*(fabs(nf.x)/Delta + fabs(nf.y)/Delta) + 1;
        t_inter_v[] = A[]/B[];
        gtr_v[] = (t_inter_v[] - tr.tr_eq)/fabs(alpha*Delta);
        //gtr_v[] = ((nf.x < 0?tr[-1]:tr[1]) - t_inter_l[])/Delta*fabs(nf.x) + ((nf.y > 0? tr[0,-1]:tr[0,1]) - t_inter_l[])/Delta*fabs(nf.y);
      /*
      coord p = {-1,0}, q = {-1,1}, s = {0,1};// coordinates of cells around interface
      double alpha_a = fabs(nf.y*s.y - alpha_1)/(sqrt(sq(nf.x) + sq(nf.y)));
      double alpha_b = fabs(nf.x*q.x + nf.y*q.y - alpha_1)/(sqrt(sq(nf.x) + sq(nf.y)));
      double alpha_c = fabs(nf.x*p.x - alpha_1)/(sqrt(sq(nf.x) + sq(nf.y)));
      if(fabs(nf.x) <= 1/2)// interface vector for x direction 
        gtr_v[] = ((nf.x<0?tr[0,1]:tr[0,-1]) - tr.tr_eq)/fabs(alpha_a*Delta);
      else if(fabs(nf.x) <= sqrt(3)/2&&fabs(nf.x)>1/2)
        gtr_v[] = ((nf.x<0&&nf.y>0)?tr[-1,1]:(nf.x<0&&nf.y<0)?tr[-1,-1]:(nf.x>0&&nf.y>0)?tr[1,1]:tr[1,-1]
            - tr.tr_eq)/fabs(alpha_b*Delta);
      else if(fabs(nf.x) > sqrt(3)/2)
        gtr_v[] = ((nf.x<0?tr[-1]:tr[1]) - tr.tr_eq)/fabs(alpha_c*Delta);// stefan problem using this term      
        */
        temp[] = (0.025*gtr_v[])/L_h*delta_l/dv(); // vapor side
    }

  }
  boundary({A, B, t_inter_v,t_inter_l,gtr_v,gtr_l,gtr_s, temp});
}