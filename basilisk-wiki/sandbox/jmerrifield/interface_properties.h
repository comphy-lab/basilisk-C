/**
# Capillary rise function
First, we define the type and the function for data recuperation of the meniscus of the capillary rise.
*/

typedef struct
{
  double h_int;
  double theta_right;
  double theta_left;
  double V;
} CapData;

CapData cap_properties (scalar f,
                        scalar cs,
                        double beta,
                        double R)
{
  CapData d = {nodata, nodata, nodata, 0.};

  double ct = cos(beta*pi/180.);
  double st = sin(beta*pi/180.);
  double xc = L0/2.;
  
  d.h_int = -HUGE;
  
  foreach() 
  {

    d.V += f[]*cs[]*dv();

    if (!interfacial(point,f))
      continue;

/**
We first use the PLIC reconstruction to reconstruct the liquid-gaz interface.
*/

    coord nf = interface_normal(point, f);
    double alpha = plane_alpha(f[], nf);
    
    coord p[2];
    if (facets(nf, alpha, p) != 2)
      continue;

    double nmag = sqrt(sq(nf.x) + sq(nf.y));
    if (nmag == 0.)
      continue;
    nf.x /= nmag;
    nf.y /= nmag;

/**
The facet vertices are converted to global coordinates.
*/
    double x1 = x + Delta*p[0].x;
    double y1 = y + Delta*p[0].y;

    double x2 = x + Delta*p[1].x;
    double y2 = y + Delta*p[1].y;

/**
We pass inthe capillary's reference frame.
*/

    double xp1 = (x1 -xc)*ct - y1*st;
    double yp1 = (x1 - xc)*st + y1*ct;

    double xp2 = (x2 - xc)*ct - y2*st;
    double yp2 = (x2 - xc)*st + y2*ct;
    
/**
We search for the intersection of the facets with the axis of the capillary
*/
    if (xp1*xp2 <= 0. && fabs(xp2 - xp1) > 1e-14) 
    {

      double lambda = -xp1/(xp2 - xp1);
      double yp = yp1 + lambda*(yp2 - yp1);

      if (yp > d.h_int)
        d.h_int = yp;
    }

/**
The contact angle is computed with a scalar product between the fluid normal nf of the PLIC reconstruction and the solid normal with its inclinaison.
*/

    double s = (x - xc)*ct - y*st;
    if (fabs(s - R) < Delta/2.)
    {
      coord ns = {ct, -st};
	
      double dot = -(nf.x*ns.x + nf.y*ns.y);
      dot = clamp(dot, -1., 1.);
		
      d.theta_right = acos(dot)*180./pi;
    }
      
    if (fabs(s + R) < Delta/2.)
    {
      coord ns = {-ct, st};
      
      double dot = -(nf.x*ns.x + nf.y*ns.y);
      dot = clamp(dot, -1., 1.);
      
      d.theta_left = acos(dot)*180./pi;
    }
  }
  return d;
}


/**
# Sessil spreading function
Now, we define the type and the function for data recuperation of the sessile spreading.
*/

typedef struct
{
  double x_left;
  double x_right;

  double U_left;
  double U_right;
  
  double theta_left;
  double theta_right;

  double height;
} SessileData;

SessileData contact_properties (scalar f,
                                scalar cs,
                                face vector fs,
                                double ywall)
{
  SessileData c = {HUGE, -HUGE,
                   nodata, nodata,
                   nodata, nodata,
                   0.};

  vector n[];
  scalar alpha[];

  reconstruction (f, n, alpha);

  foreach()
  {
    if (!interfacial(point,f))
      continue;

    if (!interfacial(point,cs))
      continue;

    coord nf;

    foreach_dimension()
      nf.x = n.x[];

    double nmag = sqrt(sq(nf.x) + sq(nf.y));

    if (nmag == 0.)
      continue;

    nf.x /= nmag;
    nf.y /= nmag;

/**
We set the solid normal for a horizontal plane
*/
    coord ns = {0.,1.};

    double a = alpha[];
    
    coord p[2];

    if (facets(nf,a,p) != 2)
      continue;

    double x1 = x + Delta*p[0].x;
    double y1 = y + Delta*p[0].y;
    
    double x2 = x + Delta*p[1].x;
    double y2 = y + Delta*p[1].y;

    double yp = ywall;

    if ((y1-yp)*(y2-yp) <= 1e-10 &&
        fabs(y2-y1) > 1e-14)
    {

      double lambda = (yp-y1)/(y2-y1);
      
      double X = x1 + lambda*(x2-x1);

      double theta;

      // left side
      if (X < c.x_left)
      {
        c.x_left = X;


/**
For the left side, we reverse the normal for it to be the outward normal of the fluid
*/

        coord nfg = {-nf.x,-nf.y};
        double dot = -(nfg.x*ns.x + nfg.y*ns.y);
        dot = clamp(dot,-1.,1.);
        theta = acos(dot)*180./pi;      

        c.theta_left = theta;


        double eps = 1e-6;
        c.U_left = interpolate(u.x,
                               X,
                               ywall + eps);
      }


      // right side
      if (X > c.x_right)
      {
        c.x_right = X;
        double dot = nf.x*ns.x + nf.y*ns.y;
        dot = clamp(dot,-1.,1.);
        theta = acos(dot)*180./pi;
        
        c.theta_right = theta;
        
        double eps = 1e-6;
        
        c.U_right = interpolate(u.x,
                                X,
                                ywall + eps);
      }
    }
  }

/**
Height of the sessile
*/
  foreach()
  {
    if (!interfacial(point, f))
      continue;

    coord nf;
    double a;

    foreach_dimension()
      nf.x = n.x[];

    double nmag = sqrt(sq(nf.x) + sq(nf.y));
    if (nmag == 0.)
      continue;

    nf.x /= nmag;
    nf.y /= nmag;

    a = alpha[];
    coord p[2];

    if (facets(nf, a, p) != 2)
      continue;

    //double x1 = x + Delta*p[0].x;
    double y1 = y + Delta*p[0].y;

    //double x2 = x + Delta*p[1].x;
    double y2 = y + Delta*p[1].y;

    double h = max(y1, y2) - ywall;

    if (h > c.height)
      c.height = h;
  }
  return c;
}

//fixme : these two functions could be assembled in one if we remove the dependency with the geometrical parameters