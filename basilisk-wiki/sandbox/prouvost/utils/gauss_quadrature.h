
/**
Integration, Gauss quadrature and total error computation
============================================

At the origin, we where interested in computing accurately the total error by comparing a numerical solution with a known analytical solution.
Thus, we implemented some [Gauss quadrature](https://en.wikipedia.org/wiki/Gaussian_quadrature) rules to compute the integral in $\mathcal{L}^p$ norm of the difference between a scalar field and a known analytical function.
Incidentally, the integral of a scalar field (in any $\mathcal{L}^p$ norm) can also be computed.

The Gauss integrations are written in 2D, and some of them (but not all of them) in 3D.

Source for the coefficients of the Gauss quadratures : Krylov V. I., *Approximate calculation of integrals*, Dover publications, 2005
	
*/



/**
Usefull functions
-----------------
	
We had some odd results with basilisk interpolation function in particular cases.
Thus we propose an interpolation function.

*/

#if dimension==2	
typedef struct { int x, y;} pseudo_v_i;
#endif
#if dimension==3
typedef struct { int x, y, z;} pseudo_v_i;
#endif

double interpolate_in_cell (Point point, scalar p, coord psi) {

   double val;

   //shift if required
   pseudo_v_i s;
   foreach_dimension () 
      s.x = -(psi.x < 0);

   //redefinition of the relative coordinates
   foreach_dimension () 
      psi.x -= s.x;
 
#if dimension == 1 
   val = p[s.x]*(1. - psi.x) 
      + p[s.x+1]*psi.x;
#endif

#if dimension == 2
   val = p[s.x,s.y]     *(1. - psi.x)*(1. - psi.y)
      + p[s.x+1,s.y]   *psi.x*(1. - psi.y)
      + p[s.x,s.y+1]   *(1. - psi.x)*psi.y
      + p[s.x+1,s.y+1] *psi.x*psi.y;
#endif

#if dimension == 3
   val = p[s.x,s.y,s.z]       *(1. - psi.x)*(1. - psi.y)*(1. - psi.z)
      + p[s.x+1,s.y,s.z]     *psi.x*(1. - psi.y)*(1. - psi.z)
      + p[s.x,s.y+1,s.z]     *(1. - psi.x)*psi.y*(1. - psi.z)
      + p[s.x+1,s.y+1,s.z]   *psi.x*psi.y*(1. - psi.z)
      + p[s.x,s.y,s.z+1]     *(1. - psi.x)*(1. - psi.y)*psi.z
      + p[s.x+1,s.y,s.z+1]   *psi.x*(1. - psi.y)*psi.z
      + p[s.x,s.y+1,s.z+1]   *(1. - psi.x)*psi.y*psi.z
      + p[s.x+1,s.y+1,s.z+1] *psi.x*psi.y*psi.z;
#endif

   return val;
}


/**
A convenient structure and a convenient function to estimate the interpolated value of a scalar or the difference between the interpolated value of a scalar and a known function at prescribed locations.

*/

struct RealError {
  scalar f;
  int Lnorm;
  double (* func) (double x, double y, double z);
};

double point_error (Point point, struct RealError re, coord psi) {

   coord xl;

   xl.x = x + Delta*psi.x;
   xl.y = y + Delta*psi.y;
   xl.z = z + Delta*psi.z;

   if (re.func)
    return fabs(re.func(xl.x, xl.y, xl.z) - interpolate_in_cell(point, re.f, psi));
   else
    return fabs(interpolate_in_cell(point, re.f, psi));

}

/**
Some non Gauss integrations
---------------------------
	
"Child-cell" integrations (integration at the location of the center of child cells).

First, integration of a function.

*/

struct IntegrateFun {
  int Lnorm;
  double (* func) (double x, double y, double z);
};

double integrate_in_cell (Point point, struct IntegrateFun p, int sublevel) {

  int Lnorm = p.Lnorm;
  int np = pow(2,sublevel);
  double val = 0.;
  coord xc;
  xc.x = x; xc.y = y; xc.z = z;

#if dimension == 1
  for (int i =0; i < np; i++) {
        xc.x = x + Delta*(-0.5 + i/((double) np)); 
        val += pow(p.func(xc.x, xc.y, xc.z) ,Lnorm)*dv()/pow(np,dimension);
  }
#endif

#if dimension == 2
  for (int i =0; i < np; i++) 
    for (int j =0; j < np; j++) {
        xc.x = x + Delta*(-0.5 + i/((double) np)); 
        xc.y = y + Delta*(-0.5 + j/((double) np)); 
        val += pow(p.func(xc.x, xc.y, xc.z) ,Lnorm)*dv()/pow(np,dimension);
    }
#endif

#if dimension == 3
  for (int i =0; i < np; i++) 
    for (int j =0; j < np; j++) 
      for (int k =0; k < np; k++) {
        xc.x = x + Delta*(-0.5 + i/((double) np)); 
        xc.y = y + Delta*(-0.5 + j/((double) np)); 
        xc.z = z + Delta*(-0.5 + j/((double) np)); 
        val += pow(p.func(xc.x, xc.y, xc.z) ,Lnorm)*dv()/pow(np,dimension);
      }
#endif

  return val;
}

double normfromChilds (struct IntegrateFun p, int sublevel)
{

  int Lnorm = p.Lnorm;
   double val = 0., volume = 0.;

   foreach(reduction(+:val) reduction(+:volume))
      if (dv() > 0.) {
        volume += dv();
        val += pow(integrate_in_cell(point, p, sublevel),Lnorm)*dv();
      }

   return (volume ? pow(val,1./Lnorm) : 0.);
}


/**
"Child"-integration to compute the total error (for scalar (and function)).

*/


double norm_child (struct RealError re)
{

  int Lnorm = re.Lnorm;

   double val = 0., volume = 0.;
   int np = pow(2,dimension);

#if dimension == 1
   static coord int_points[2] = {
				 {-0.25,0.,0.},{0.25,0.,0.}
   };
#endif
#if dimension == 2
   static coord int_points[4] = {
				 {-0.25,0.25,0.},{0.25,0.25,0.},{-0.25,-0.25,0.},{0.25,-0.25,0.}
   };
#endif
#if dimension == 3
   static coord int_points[8] = {
				 {-0.25,0.25,0.25},{0.25,0.25,0.25},{-0.25,-0.25,0.25},{0.25,-0.25,0.25},
				 {-0.25,0.25,-0.25},{0.25,0.25,-0.25},{-0.25,-0.25,-0.25},{0.25,-0.25,-0.25}
   };
#endif

   foreach(reduction(+:val) reduction(+:volume)) {
      if (dv() > 0.) {
	 volume += dv();

	 for (int i =0; i < np; i++) {
	   val += pow(point_error(point, re, int_points[i]), Lnorm)*dv()/np;
	}

      }
   }

   return (volume ? pow(val/volume,1./Lnorm) : 0.);
}

/**

Gauss-quadratures
-----------------

### 3 points Gauss quadrature 

Only in 2D for now.

*/

#if dimension == 2
double norm_gauss_3p (struct RealError re)
{

  int Lnorm = re.Lnorm;

  double val = 0., volume = 0.;
  int np = 9;

   double p1 = -sqrt(3./5.)/2.;
   double p3 = sqrt(3./5.)/2.;
   coord int_points[9] = {
			  {p1,p1,0.},{p1,0.,0.},{p1,p3,0.},
			  {0.,p1,0.},{0.,0.,0.},{0.,p3,0.},
			  {p3,p1,0.},{p3,0.,0.},{p3,p3,0.}
   };

   double coef_points[9] = {
			    0.3086419753086420,
			    0.4938271604938272,    
			    0.3086419753086420,    
			    0.4938271604938272,    
			    0.7901234567901235,    
			    0.4938271604938272,    
			    0.3086419753086420,    
			    0.4938271604938272,    
			    0.3086419753086420  
   };


   foreach(reduction(+:val) reduction(+:volume))
      if (dv() > 0.) {
	 volume += dv();

	 for (int i =0; i < np; i++) 
	    val += pow(point_error(point, re, int_points[i]), Lnorm)*coef_points[i] *dv()/4.;

      }

   return (volume ? pow(val/volume,1./Lnorm) : 0.);
}
#endif

/**
### 4 points Gauss quadrature 

Only in 2D for now.
*/

#if dimension == 2
double norm_gauss_4p (struct RealError re)
{

  int Lnorm = re.Lnorm;

  double val = 0., volume = 0.;
  int np = 16;


   double p1 = sqrt( 3./7. - 2./7.*sqrt(6./5.) ) /2.;
   double p2 = -sqrt( 3./7. - 2./7.*sqrt(6./5.) ) /2.;
   double p3 = sqrt( 3./7. + 2./7.*sqrt(6./5.) ) /2.;
   double p4 = -sqrt( 3./7. + 2./7.*sqrt(6./5.) ) /2.;

   coord int_points[16] = {
			   {p1,p1,0.},{p1,p2,0.},{p1,p3,0.},{p1,p4,0.},
			   {p2,p1,0.},{p2,p2,0.},{p2,p3,0.},{p2,p4,0.},
			   {p3,p1,0.},{p3,p2,0.},{p3,p3,0.},{p3,p4,0.},
			   {p4,p1,0.},{p4,p2,0.},{p4,p3,0.},{p4,p4,0.}
   };

   double pds12 = (18.+sqrt(30.))/36.;
   double pds34 = (18.-sqrt(30.))/36.;
   double coef_points[16] = {
			     pds12*pds12,  pds12*pds12, pds12*pds34, pds12*pds34, 
			     pds12*pds12,  pds12*pds12, pds12*pds34, pds12*pds34, 
			     pds34*pds12,  pds34*pds12, pds34*pds34, pds34*pds34, 
			     pds34*pds12,  pds34*pds12, pds34*pds34, pds34*pds34, 
   };


   foreach(reduction(+:val) reduction(+:volume))
      if (dv() > 0.) {
	 volume += dv();
	 for (int i =0; i < np; i++) 
	    val += pow(point_error(point, re, int_points[i]), Lnorm)*coef_points[i] *dv()/4.;

      }

   return (volume ? pow(val/volume,1./Lnorm) : 0.);
}
#endif


/**
### 5 points Gauss quadrature 

in 2D and 3D.
*/

double norm_gauss_5p (struct RealError re)
{

  int Lnorm = re.Lnorm;

  double val = 0., volume = 0.;
  int np = 25;


#if dimension == 2
   double p1 = 1./3.*sqrt( 5. - 2.*sqrt(10./7.) ) /2.;
   double p2 = -1./3.*sqrt( 5. - 2.*sqrt(10./7.) ) /2.;
   double p3 = 1./3.*sqrt( 5. + 2.*sqrt(10./7.) ) /2.;
   double p4 = -1./3.*sqrt( 5. + 2.*sqrt(10./7.) ) /2.;
   double p5 = 0.;

   coord int_points[25] = {
			   {p1,p1,0.},{p1,p2,0.},{p1,p3,0.},{p1,p4,0.},{p1,p5,0.},
			   {p2,p1,0.},{p2,p2,0.},{p2,p3,0.},{p2,p4,0.},{p2,p5,0.},
			   {p3,p1,0.},{p3,p2,0.},{p3,p3,0.},{p3,p4,0.},{p3,p5,0.},
			   {p4,p1,0.},{p4,p2,0.},{p4,p3,0.},{p4,p4,0.},{p4,p5,0.},
			   {p5,p1,0.},{p5,p2,0.},{p5,p3,0.},{p5,p4,0.},{p5,p5,0.}
   };

   double pds12 = (322.+13.*sqrt(70.))/900.;
   double pds34 = (322.-13.*sqrt(70.))/900.;
   double pds5 = 128./225.;
   double coef_points[25] = {
			     pds12*pds12,  pds12*pds12, pds12*pds34, pds12*pds34, pds12*pds5,  
			     pds12*pds12,  pds12*pds12, pds12*pds34, pds12*pds34, pds12*pds5,
			     pds34*pds12,  pds34*pds12, pds34*pds34, pds34*pds34, pds34*pds5,
			     pds34*pds12,  pds34*pds12, pds34*pds34, pds34*pds34, pds34*pds5,
			     pds5*pds12,  pds5*pds12, pds5*pds34, pds5*pds34, pds5*pds5
   };
#endif
#if dimension == 3
   double p1 = 1./3.*sqrt( 5. - 2.*sqrt(10./7.) ) /2.;
   double p2 = -1./3.*sqrt( 5. - 2.*sqrt(10./7.) ) /2.;
   double p3 = 1./3.*sqrt( 5. + 2.*sqrt(10./7.) ) /2.;
   double p4 = -1./3.*sqrt( 5. + 2.*sqrt(10./7.) ) /2.;
   double p5 = 0.;

   coord int_points[125] = {
			   {p1,p1,p1},{p1,p2,p1},{p1,p3,p1},{p1,p4,p1},{p1,p5,p1},
			   {p2,p1,p1},{p2,p2,p1},{p2,p3,p1},{p2,p4,p1},{p2,p5,p1},
			   {p3,p1,p1},{p3,p2,p1},{p3,p3,p1},{p3,p4,p1},{p3,p5,p1},
			   {p4,p1,p1},{p4,p2,p1},{p4,p3,p1},{p4,p4,p1},{p4,p5,p1},
			   {p5,p1,p1},{p5,p2,p1},{p5,p3,p1},{p5,p4,p1},{p5,p5,p1},
			   {p1,p1,p2},{p1,p2,p2},{p1,p3,p2},{p1,p4,p2},{p1,p5,p2},
			   {p2,p1,p2},{p2,p2,p2},{p2,p3,p2},{p2,p4,p2},{p2,p5,p2},
			   {p3,p1,p2},{p3,p2,p2},{p3,p3,p2},{p3,p4,p2},{p3,p5,p2},
			   {p4,p1,p2},{p4,p2,p2},{p4,p3,p2},{p4,p4,p2},{p4,p5,p2},
			   {p5,p1,p2},{p5,p2,p2},{p5,p3,p2},{p5,p4,p2},{p5,p5,p2},
			   {p1,p1,p3},{p1,p2,p3},{p1,p3,p3},{p1,p4,p3},{p1,p5,p3},
			   {p2,p1,p3},{p2,p2,p3},{p2,p3,p3},{p2,p4,p3},{p2,p5,p3},
			   {p3,p1,p3},{p3,p2,p3},{p3,p3,p3},{p3,p4,p3},{p3,p5,p3},
			   {p4,p1,p3},{p4,p2,p3},{p4,p3,p3},{p4,p4,p3},{p4,p5,p3},
			   {p5,p1,p3},{p5,p2,p3},{p5,p3,p3},{p5,p4,p3},{p5,p5,p3},
			   {p1,p1,p4},{p1,p2,p4},{p1,p3,p4},{p1,p4,p4},{p1,p5,p4},
			   {p2,p1,p4},{p2,p2,p4},{p2,p3,p4},{p2,p4,p4},{p2,p5,p4},
			   {p3,p1,p4},{p3,p2,p4},{p3,p3,p4},{p3,p4,p4},{p3,p5,p4},
			   {p4,p1,p4},{p4,p2,p4},{p4,p3,p4},{p4,p4,p4},{p4,p5,p4},
			   {p5,p1,p4},{p5,p2,p4},{p5,p3,p4},{p5,p4,p4},{p5,p5,p4},
			   {p1,p1,p5},{p1,p2,p5},{p1,p3,p5},{p1,p4,p5},{p1,p5,p5},
			   {p2,p1,p5},{p2,p2,p5},{p2,p3,p5},{p2,p4,p5},{p2,p5,p5},
			   {p3,p1,p5},{p3,p2,p5},{p3,p3,p5},{p3,p4,p5},{p3,p5,p5},
			   {p4,p1,p5},{p4,p2,p5},{p4,p3,p5},{p4,p4,p5},{p4,p5,p5},
			   {p5,p1,p5},{p5,p2,p5},{p5,p3,p5},{p5,p4,p5},{p5,p5,p5}
   };

   double pds12 = (322.+13.*sqrt(70.))/900.;
   double pds34 = (322.-13.*sqrt(70.))/900.;
   double pds5 = 128./225.;
   double coef_points[125] = {
			     pds12*pds12*pds12,  pds12*pds12*pds12, pds12*pds34*pds12, pds12*pds34*pds12, pds12*pds5*pds12,  
			     pds12*pds12*pds12,  pds12*pds12*pds12, pds12*pds34*pds12, pds12*pds34*pds12, pds12*pds5*pds12,
			     pds34*pds12*pds12,  pds34*pds12*pds12, pds34*pds34*pds12, pds34*pds34*pds12, pds34*pds5*pds12,
			     pds34*pds12*pds12,  pds34*pds12*pds12, pds34*pds34*pds12, pds34*pds34*pds12, pds34*pds5*pds12,
			     pds5*pds12*pds12,  pds5*pds12*pds12, pds5*pds34*pds12, pds5*pds34*pds12, pds5*pds5*pds12,
			     pds12*pds12*pds12,  pds12*pds12*pds12, pds12*pds34*pds12, pds12*pds34*pds12, pds12*pds5*pds12,  
			     pds12*pds12*pds12,  pds12*pds12*pds12, pds12*pds34*pds12, pds12*pds34*pds12, pds12*pds5*pds12,
			     pds34*pds12*pds12,  pds34*pds12*pds12, pds34*pds34*pds12, pds34*pds34*pds12, pds34*pds5*pds12,
			     pds34*pds12*pds12,  pds34*pds12*pds12, pds34*pds34*pds12, pds34*pds34*pds12, pds34*pds5*pds12,
			     pds5*pds12*pds12,  pds5*pds12*pds12, pds5*pds34*pds12, pds5*pds34*pds12, pds5*pds5*pds12,
			     pds12*pds12*pds34,  pds12*pds12*pds34, pds12*pds34*pds34, pds12*pds34*pds34, pds12*pds5*pds34,  
			     pds12*pds12*pds34,  pds12*pds12*pds34, pds12*pds34*pds34, pds12*pds34*pds34, pds12*pds5*pds34,
			     pds34*pds12*pds34,  pds34*pds12*pds34, pds34*pds34*pds34, pds34*pds34*pds34, pds34*pds5*pds34,
			     pds34*pds12*pds34,  pds34*pds12*pds34, pds34*pds34*pds34, pds34*pds34*pds34, pds34*pds5*pds34,
			     pds5*pds12*pds34,  pds5*pds12*pds34, pds5*pds34*pds34, pds5*pds34*pds34, pds5*pds5*pds34,
			     pds12*pds12*pds34,  pds12*pds12*pds34, pds12*pds34*pds34, pds12*pds34*pds34, pds12*pds5*pds34,  
			     pds12*pds12*pds34,  pds12*pds12*pds34, pds12*pds34*pds34, pds12*pds34*pds34, pds12*pds5*pds34,
			     pds34*pds12*pds34,  pds34*pds12*pds34, pds34*pds34*pds34, pds34*pds34*pds34, pds34*pds5*pds34,
			     pds34*pds12*pds34,  pds34*pds12*pds34, pds34*pds34*pds34, pds34*pds34*pds34, pds34*pds5*pds34,
			     pds5*pds12*pds34,  pds5*pds12*pds34, pds5*pds34*pds34, pds5*pds34*pds34, pds5*pds5*pds34,
			     pds12*pds12*pds5,  pds12*pds12*pds5, pds12*pds34*pds5, pds12*pds34*pds5, pds12*pds5*pds5,  
			     pds12*pds12*pds5,  pds12*pds12*pds5, pds12*pds34*pds5, pds12*pds34*pds5, pds12*pds5*pds5,
			     pds34*pds12*pds5,  pds34*pds12*pds5, pds34*pds34*pds5, pds34*pds34*pds5, pds34*pds5*pds5,
			     pds34*pds12*pds5,  pds34*pds12*pds5, pds34*pds34*pds5, pds34*pds34*pds5, pds34*pds5*pds5,
			     pds5*pds12*pds5,  pds5*pds12*pds5, pds5*pds34*pds5, pds5*pds34*pds5, pds5*pds5*pds5
   };
#endif


   foreach(reduction(+:val) reduction(+:volume))
      if (dv() > 0.) {
	 volume += dv();

	 for (int i =0; i < np; i++) 
	   val += pow(point_error(point, re, int_points[i]), Lnorm)*coef_points[i] *dv()/pow(2,dimension);  // the sum of the weight coef_points is 2^dimension

      }

   return (volume ? pow(val/volume,1./Lnorm) : 0.);
}






/**
### 5 points Gauss quadrature whith local error registration

Only in 2D for now.

*/



#if dimension == 2
double norm_gauss_5p_local (struct RealError re, scalar local_error)
{

  int Lnorm = re.Lnorm;

  double val = 0., volume = 0.;
  int np = 25;


   double p1 = 1./3.*sqrt( 5. - 2.*sqrt(10./7.) ) /2.;
   double p2 = -1./3.*sqrt( 5. - 2.*sqrt(10./7.) ) /2.;
   double p3 = 1./3.*sqrt( 5. + 2.*sqrt(10./7.) ) /2.;
   double p4 = -1./3.*sqrt( 5. + 2.*sqrt(10./7.) ) /2.;
   double p5 = 0.;

   coord int_points[25] = {
			   {p1,p1,0.},{p1,p2,0.},{p1,p3,0.},{p1,p4,0.},{p1,p5,0.},
			   {p2,p1,0.},{p2,p2,0.},{p2,p3,0.},{p2,p4,0.},{p2,p5,0.},
			   {p3,p1,0.},{p3,p2,0.},{p3,p3,0.},{p3,p4,0.},{p3,p5,0.},
			   {p4,p1,0.},{p4,p2,0.},{p4,p3,0.},{p4,p4,0.},{p4,p5,0.},
			   {p5,p1,0.},{p5,p2,0.},{p5,p3,0.},{p5,p4,0.},{p5,p5,0.}
   };

   double pds12 = (322.+13.*sqrt(70.))/900.;
   double pds34 = (322.-13.*sqrt(70.))/900.;
   double pds5 = 128./225.;
   double coef_points[25] = {
			     pds12*pds12,  pds12*pds12, pds12*pds34, pds12*pds34, pds12*pds5,  
			     pds12*pds12,  pds12*pds12, pds12*pds34, pds12*pds34, pds12*pds5,
			     pds34*pds12,  pds34*pds12, pds34*pds34, pds34*pds34, pds34*pds5,
			     pds34*pds12,  pds34*pds12, pds34*pds34, pds34*pds34, pds34*pds5,
			     pds5*pds12,  pds5*pds12, pds5*pds34, pds5*pds34, pds5*pds5
   };

   foreach(reduction(+:val) reduction(+:volume)){
      local_error[]=0.;
      if (dv() > 0.) {
	 volume += dv();

	 for (int i =0; i < np; i++) {
           val += pow(point_error(point, re, int_points[i]), Lnorm)*coef_points[i] *dv()/4.;
           local_error[] += pow(point_error(point, re, int_points[i]), Lnorm)*coef_points[i] *dv()/4.;
        }
      }
   }

   return (volume ? pow(val/volume,1./Lnorm) : 0.);
}
#endif






