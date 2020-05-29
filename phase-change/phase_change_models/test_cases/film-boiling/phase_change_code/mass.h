void Malan (scalar tr, scalar f, scalar m, double L_h)
{
  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary({f});

  double d[3][4];
  double eps[3][4];
  double gamma[3][4];
  double phi[3][4];
  double mm[3][4];
  //double m_sum[4];
  double tt[3][4]; 
  
  foreach()
  {
    coord p[3] = {{-1,1},{-1,0},{-0,1}};
    coord q[3] = {{-1,-1},{-1,0},{0,-1}};
    coord r[3] = {{1,1},{1,0},{0,1}};
    coord s[3] = {{1,-1},{1,0},{0,-1}};
    coord n = interface_normal(point,f);
    normalize(&n);
    double alpha  = plane_alpha (f[], n) + 0.5*(n.x + n.y);  
    if(f[]>1e-12&&f[]<1-1e-12)
    {
      tt[0][0] = tr[-1,1];
      tt[1][0] = tr[-1,0];
      tt[2][0] = tr[0,1];
      tt[0][1] = tr[-1,-1];
      tt[1][1] = tr[-1,0];
      tt[2][1] = tr[0,-1];
      tt[0][2] = tr[1,1];
      tt[1][2] = tr[1,0];
      tt[2][2] = tr[0,1];
      tt[0][3] = tr[1,-1];
      tt[1][3] = tr[1,0];
      tt[2][3] = tr[0,-1];

      if(n.x<0&&n.y>0)
        {
          d[0][0] = fabs(-alpha + p[0].x*n.x + p[0].y*n.y)/sqrt(sq(n.x) + sq(n.y));
          eps[0][0] = fabs(p[0].x*n.x + p[0].y*n.y);
          gamma[0][0] = eps[0][0]/(sq(p[0].x) + sq(p[0].y));
          phi[0][0] = gamma[0][0]/(gamma[0][0] + gamma[1][0] + gamma[2][0]);
          mm[0][0] = phi[0][0]*(tt[0][0] - tr.tr_eq)/(Delta*d[0][0]);

          d[1][0] = fabs(-alpha + p[1].x*n.x + p[1].y*n.y)/sqrt(sq(n.x) + sq(n.y));
          eps[1][0] = fabs(p[1].x*n.x + p[1].y*n.y);
          gamma[1][0] = eps[1][0]/(sq(p[1].x) + sq(p[1].y));
          phi[1][0] = gamma[1][0]/(gamma[0][0] + gamma[1][0] + gamma[2][0]);
          mm[1][0] = phi[1][0]*(tt[1][0] - tr.tr_eq)/(Delta*d[1][0]);

          d[2][0] = fabs(-alpha + p[2].x*n.x + p[2].y*n.y)/sqrt(sq(n.x) + sq(n.y));
          eps[2][0] = fabs(p[2].x*n.x + p[2].y*n.y);
          gamma[2][0] = eps[2][0]/(sq(p[2].x) + sq(p[2].y));
          phi[2][0] = gamma[2][0]/(gamma[0][0] + gamma[1][0] + gamma[2][0]);
          mm[2][0] = phi[2][0]*(tt[2][0] - tr.tr_eq)/(Delta*d[2][0]);

          m[] = mm[1][0] + mm[2][0];
        }
        
      if(n.x<=0&&n.y<=0)
        {
          d[0][1] = fabs(-alpha + q[0].x*n.x + q[0].y*n.y)/sqrt(sq(n.x) + sq(n.y));
          eps[0][1] = fabs(q[0].x*n.x + q[0].y*n.y);
          gamma[0][1] = eps[0][1]/(sq(q[0].x) + sq(q[0].y));
          phi[0][1] = gamma[0][1]/(gamma[0][1] + gamma[1][1] + gamma[2][1]);
          mm[0][1] = phi[0][1]*(tt[0][1] - tr.tr_eq)/(Delta*d[0][1]);
        
          d[1][1] = fabs(-alpha + q[1].x*n.x + q[1].y*n.y)/sqrt(sq(n.x) + sq(n.y));
          eps[1][1] = fabs(q[1].x*n.x + q[1].y*n.y);
          gamma[1][1] = eps[1][1]/(sq(q[1].x) + sq(q[1].y));
          phi[1][1] = gamma[1][1]/(gamma[0][1] + gamma[1][1] + gamma[2][1]);
          mm[1][1] = phi[1][1]*(tt[1][1] - tr.tr_eq)/(Delta*d[1][1]);
        
          d[2][1] = fabs(-alpha + q[2].x*n.x + q[2].y*n.y)/sqrt(sq(n.x) + sq(n.y));
          eps[2][1] = fabs(q[2].x*n.x + q[2].y*n.y);
          gamma[2][1] = eps[2][1]/(sq(q[2].x) + sq(q[2].y));
          phi[2][1] = gamma[2][1]/(gamma[0][1] + gamma[1][1] + gamma[2][1]);
          mm[2][1] = phi[2][1]*(tt[2][1] - tr.tr_eq)/(Delta*d[2][1]);

          m[] = mm[0][1] + mm[1][1] + mm[2][1];
        }
        

      if(n.x>=0&&n.y>=0)
        {
          d[0][2] = fabs(-alpha + r[0].x*n.x + r[0].y*n.y)/sqrt(sq(n.x) + sq(n.y));
          eps[0][2] = fabs(r[0].x*n.x + r[0].y*n.y);
          gamma[0][2] = eps[0][2]/(sq(r[0].x) + sq(r[0].y));
          phi[0][2] = gamma[0][2]/(gamma[0][2] + gamma[1][2] + gamma[2][2]);
          mm[0][2] = phi[0][2]*(tt[0][2] - tr.tr_eq)/(Delta*d[0][2]);
        
          d[1][2] = fabs(-alpha + r[1].x*n.x + r[1].y*n.y)/sqrt(sq(n.x) + sq(n.y));
          eps[1][2] = fabs(r[1].x*n.x + r[1].y*n.y);
          gamma[1][2] = eps[1][2]/(sq(r[1].x) + sq(r[1].y));
          phi[1][2] = gamma[1][2]/(gamma[0][2] + gamma[1][2] + gamma[2][2]);
          mm[1][2] = phi[1][2]*(tt[1][2] - tr.tr_eq)/(Delta*d[1][2]);
        
          d[2][2] = fabs(-alpha + r[2].x*n.x + r[2].y*n.y)/sqrt(sq(n.x) + sq(n.y));
          eps[2][2] = fabs(r[2].x*n.x + r[2].y*n.y);
          gamma[2][2] = eps[2][2]/(sq(r[2].x) + sq(r[2].y));
          phi[2][2] = gamma[2][2]/(gamma[0][2] + gamma[1][2] + gamma[2][2]);
          mm[2][2] = phi[2][2]*(tt[2][2] - tr.tr_eq)/(Delta*d[2][2]);

          m[] = mm[0][2] + mm[1][2] + mm[2][2];
        }

      if(n.x>0&&n.y<0)
        {
          d[0][3] = fabs(-alpha + s[0].x*n.x + s[0].y*n.y)/sqrt(sq(n.x) + sq(n.y));
          eps[0][3] = fabs(s[0].x*n.x + s[0].y*n.y);
          gamma[0][3] = eps[0][3]/(sq(s[0].x) + sq(s[0].y));
          phi[0][3] = gamma[0][3]/(gamma[0][3] + gamma[1][3] + gamma[2][3]);
          mm[0][3] = phi[0][3]*(tt[0][3] - tr.tr_eq)/(Delta*d[0][3]);
       
          d[1][3] = fabs(-alpha + s[1].x*n.x + s[1].y*n.y)/sqrt(sq(n.x) + sq(n.y));
          eps[1][3] = fabs(s[1].x*n.x + s[1].y*n.y);
          gamma[1][3] = eps[1][3]/(sq(s[1].x) + sq(s[1].y));
          phi[1][3] = gamma[1][3]/(gamma[0][3] + gamma[1][3] + gamma[2][3]);
          mm[1][3] = phi[1][3]*(tt[1][3] - tr.tr_eq)/(Delta*d[1][3]);
        
          d[2][3] = fabs(-alpha + s[2].x*n.x + s[2].y*n.y)/sqrt(sq(n.x) + sq(n.y));
          eps[2][3] = fabs(s[2].x*n.x + s[2].y*n.y);
          gamma[2][3] = eps[2][3]/(sq(s[2].x) + sq(s[2].y));
          phi[2][3] = gamma[2][3]/(gamma[0][3] + gamma[1][3] + gamma[2][3]);
          mm[2][3] = phi[2][3]*(tt[2][3] - tr.tr_eq)/(Delta*d[2][3]);

          m[] = mm[0][3] + mm[1][3] + mm[2][3];
        }
        
    }
    m[] *= tr.lambda/L_h*f[];
  }
  boundary({m});
}