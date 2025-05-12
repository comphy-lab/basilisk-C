#ifndef _MY_VISCOSITY_H
#define _MY_VISCOSITY_H
#include "mypoisson.h"

struct Viscosity {
  vector u;
  face vector mu;
  scalar rho;
  double dt;
  int nrelax;
  scalar * res;
};

#if AXI
# define lambda ((coord){1., 1. + dt/rho[]*(mu.x[] + mu.x[1] + \
					    mu.y[] + mu.y[0,1])/2./sq(y)})
#else // not AXI
# if dimension == 1
#   define lambda ((coord){1.})
# elif dimension == 2
#   define lambda ((coord){1.,1.})
# elif dimension == 3
#   define lambda ((coord){1.,1.,1.})
#endif
#endif

#if USE_MY_SOLID
extern scalar is_solid;
extern face vector is_solid_face;
extern bool IS_SOLID_x;
extern bool IS_SOLID_y;
extern const bool is_slip_x;
extern const bool is_slip_y;
void boundarySolidVelCNoauto (vector v);

static void relaxVelX(Point point, const face vector mu, const scalar rho, 
                        const double dt, const vector r, vector u, vector w)
{
  double coef = 0.0;
  double dt_rho = dt / rho[];
  //contribution from tauxx
  double tauxx = 2. * mu.x[1] * u.x[1] + 2. * mu.x[] * u.x[-1];
  coef = coef +  2. * mu.x[1] + 2. * mu.x[]; //diagonal terms from tauxx
  double dudy = mu.y[0, 1] * u.x[0, 1] + mu.y[] * u.x[0, -1]; //contribution from dudy in tauxy
  coef = coef + mu.y[0, 1] + mu.y[]; //diagonal terms from dudy in tauxy
  coef = coef * dt / rho[] + sq(Delta)*lambda.x;
  double dvdx_yp = (u.y[1, 0] + u.y[1, 1]) / 4. - (u.y[-1, 0] + u.y[-1, 1]) / 4.;
  double dvdx_ym = (u.y[1, -1] + u.y[1, 0]) / 4. - (u.y[-1, -1] + u.y[-1, 0]) / 4.;
  double tauxy = mu.y[0, 1] * dvdx_yp - mu.y[] * dvdx_ym + dudy; //contribution from d tauxy dy
  double res = dt / rho[] * (tauxx + tauxy) + r.x[] * sq(Delta); //the second term is the contribution from advection
  w.x[] = res / coef;
}

static void relaxVelY(Point point, const face vector mu, const scalar rho, 
                        const double dt, const vector r, vector u, vector w)
{
  double coef = 0.0;
  //contribution from tauyy
  double tauyy = 2. * mu.y[0, 1] * u.y[0, 1] + 2. * mu.y[] * u.y[0, -1];
  coef = coef + 2. * mu.y[0, 1] + 2. * mu.y[];                // diagonal terms from tauyy
  double dvdx = mu.x[1] * u.y[1] + mu.x[] * u.y[-1]; //contribution from dvdx in tauxy
  coef = coef + mu.x[1] + mu.x[]; //diagonal terms from dvdx in tauxy
  coef = coef * dt / rho[] + sq(Delta)*lambda.y;
  double dudy_xp = (u.x[0, 1] + u.x[1, 1]) / 4. - (u.x[0, -1] + u.x[1, -1]) / 4.;
  double dudy_xm = (u.x[-1, 1] + u.x[0, 1]) / 4. - (u.x[-1, -1] + u.x[0, -1]) / 4.;
  double tauxy = mu.x[1] * dudy_xp - mu.x[] * dudy_xm + dvdx;          // contribution from d tauxy dx
  double res = dt / rho[] * (tauyy + tauxy) + r.y[] * sq(Delta); //the second term is the contribution from advection
  w.y[] = res / coef;
}

static void relaxVelXSolid(Point point, const face vector mu, const scalar rho, 
                        const double dt, const vector r, vector u, vector w)
{
  if((int) is_solid[] == 1)
  {
    u.x[] = 0.0;
    w.x[] = 0.0;
    return;
  }

  double coordmx = x - Delta;
  double m1x = 1.0;
  double compx = 0.0; //compensating coef
  if(IS_SOLID_x && coordmx <= 0.0)
  {
    m1x = 0.0;
    compx = 2.0 * mu.x[];
  }
  double coef = 0.0;
  //contribution from tauxx
  double tauxx = 2. * mu.x[1] * u.x[1] + 2. * mu.x[] * u.x[-1] * m1x;
  coef = coef +  2. * mu.x[1] + 2. * mu.x[] + compx; //diagonal terms from tauxx

  

  double coordmy = y - Delta;
  double m1y = 1.0;
  double compy = 0.0; //compensating coef
  if(IS_SOLID_y && coordmy <= 0.0)
  {
    m1y = 0.0;
    compy = is_slip_y ? -mu.y[] : mu.y[];
  }
  double dudy = mu.y[0, 1] * u.x[0, 1] + mu.y[] * u.x[0, -1]*m1y; //contribution from dudy in tauxy
  coef = coef + mu.y[0, 1] + mu.y[] + compy; //diagonal terms from dudy in tauxy
  coef = coef * dt / rho[] + sq(Delta)*lambda.x;
  
  double uyxmyp = (u.y[-1, 0] + u.y[-1, 1]) / 4.;
  double uyxmym = (u.y[-1, -1] + u.y[-1, 0]) / 4.;
  if(IS_SOLID_x && coordmx <= 0.0)
  {
    uyxmyp = is_slip_x ? (u.y[0, 0] + u.y[0, 1]) / 4. : -(u.y[0, 0] + u.y[0, 1]) / 4.;
    uyxmym = is_slip_x ? (u.y[0, -1] + u.y[0, 0]) / 4. : - (u.y[0, -1] + u.y[0, 0]) / 4.;
  }

  double dvdx_yp = (u.y[1, 0] + u.y[1, 1]) / 4. - uyxmyp;
  double dvdx_ym = (u.y[1, -1] + u.y[1, 0]) / 4. - uyxmym;

  if(IS_SOLID_y && coordmy <= 0.0)
  {
    dvdx_ym = 0.0; //the solid-liquid boundary (v = 0)
  }

  double tauxy = mu.y[0, 1] * dvdx_yp - mu.y[] * dvdx_ym + dudy; //contribution from d tauxy dy
  double res = dt / rho[] * (tauxx + tauxy) + r.x[] * sq(Delta); //the second term is the contribution from advection
  w.x[] = res / coef;

}

static void relaxVelYSolid(Point point, const face vector mu, const scalar rho, 
                          const double dt, const vector r, vector u, vector w)
{
  if((int) is_solid[] == 1)
  {
    u.y[] = 0.0;
    w.y[] = 0.0;
    return;
  }
  double coordmy = y - Delta;
  double m1y = 1.0;
  double compy = 0.0; //compensating coef
  if(IS_SOLID_y && coordmy <= 0.0)
  {
    m1y = 0.0;
    compy = 2.0 * mu.y[];
  }

  //this part will not be influenced by the existence of the solid
  double coef = 0.0;
  //contribution from tauyy
  double tauyy = 2. * mu.y[0, 1] * u.y[0, 1] + 2. * mu.y[] * u.y[0, -1] * m1y;
  coef = coef + 2. * mu.y[0, 1] + 2. * mu.y[] + compy;                // diagonal terms from tauyy

  double coordmx = x - Delta;
  double m1x = 1.0;
  double compx = 0.0; //compensating coef
  if(IS_SOLID_x && coordmx <= 0.0)
  {
    m1x = 0.0;
    compx = is_slip_x ? -mu.x[] : mu.x[];
  }

  double dvdx = mu.x[1] * u.y[1] + mu.x[] * u.y[-1] * m1x; //contribution from dvdx in tauxy
  coef = coef + mu.x[1] + mu.x[] + compx; //diagonal terms from dvdx in tauxy
  coef = coef * dt / rho[] + sq(Delta)*lambda.y;


  double uxymxp = (u.x[0, -1] + u.x[1, -1]) / 4.;
  double uyymxm = (u.x[-1, -1] + u.x[0, -1]) / 4.;

  if (IS_SOLID_y && coordmy <= 0.0)
  {
    uxymxp = is_slip_y ? (u.x[0, 0] + u.x[1, 0]) / 4. : -(u.x[0, 0] + u.x[1, 0]) / 4.;
    uyymxm = is_slip_y ? (u.x[-1, 0] + u.x[0, 0]) / 4. : -(u.x[-1, 0] + u.x[0, 0]) / 4.;
  }

  double dudy_xp = (u.x[0, 1] + u.x[1, 1]) / 4. - uxymxp;
  double dudy_xm = (u.x[-1, 1] + u.x[0, 1]) / 4. - uyymxm;

  if(IS_SOLID_x && coordmx <= 0.0)
  {
    dudy_xm = 0.0; //the solid-liquid boundary (u = 0)
  }

  double tauxy = mu.x[1] * dudy_xp - mu.x[] * dudy_xm + dvdx;          // contribution from d tauxy dx
  double res = dt / rho[] * (tauyy + tauxy) + r.y[] * sq(Delta); //the second term is the contribution from advection
  w.y[] = res / coef;

}

//TODO: merge those functions into one
static void relaxVelXSolidX(Point point, const face vector mu, const scalar rho, 
                        const double dt, const vector r, vector u, vector w)
{
  if((int) is_solid[] == 1)
  {
    u.x[] = 0.0;
    w.x[] = 0.0;
    return;
  }
  double coordm = x - Delta;
  double m1 = 1.0;
  double comp = 0.0; //compensating coef
  if(coordm <= 0.0)
  {
    m1 = 0.0;
    comp = 2.0 * mu.x[];
  }
  double coef = 0.0;
  //contribution from tauxx
  double tauxx = 2. * mu.x[1] * u.x[1] + 2. * mu.x[] * u.x[-1] * m1;
  coef = coef +  2. * mu.x[1] + 2. * mu.x[] + comp; //diagonal terms from tauxx

  //will not be influenced by the solid part
  double dudy = mu.y[0, 1] * u.x[0, 1] + mu.y[] * u.x[0, -1]; //contribution from dudy in tauxy
  coef = coef + mu.y[0, 1] + mu.y[]; //diagonal terms from dudy in tauxy
  coef = coef * dt / rho[] + sq(Delta)*lambda.x;
  
  double uyxmyp = (u.y[-1, 0] + u.y[-1, 1]) / 4.;
  double uyxmym = (u.y[-1, -1] + u.y[-1, 0]) / 4.;
  if(coordm <= 0.0)
  {
    uyxmyp = is_slip_x ? (u.y[0, 0] + u.y[0, 1]) / 4. : -(u.y[0, 0] + u.y[0, 1]) / 4.;
    uyxmym = is_slip_x ? (u.y[0, -1] + u.y[0, 0]) / 4. : - (u.y[0, -1] + u.y[0, 0]) / 4.;
  }

  double dvdx_yp = (u.y[1, 0] + u.y[1, 1]) / 4. - uyxmyp;
  double dvdx_ym = (u.y[1, -1] + u.y[1, 0]) / 4. - uyxmym;
  double tauxy = mu.y[0, 1] * dvdx_yp - mu.y[] * dvdx_ym + dudy; //contribution from d tauxy dy
  double res = dt / rho[] * (tauxx + tauxy) + r.x[] * sq(Delta); //the second term is the contribution from advection
  w.x[] = res / coef;
}

static void relaxVelYSolidX(Point point, const face vector mu, const scalar rho, 
                            const double dt, const vector r, vector u, vector w)
{
  if((int) is_solid[] == 1)
  {
    u.y[] = 0.0;
    w.y[] = 0.0;
    return;
  }
  //this part will not be influenced by the existence of the solid
  double coef = 0.0;
  //contribution from tauyy
  double tauyy = 2. * mu.y[0, 1] * u.y[0, 1] + 2. * mu.y[] * u.y[0, -1];
  coef = coef + 2. * mu.y[0, 1] + 2. * mu.y[];                // diagonal terms from tauyy

  double coordm = x - Delta;
  double m1 = 1.0;
  double comp = 0.0; //compensating coef
  if(coordm <= 0.0)
  {
    m1 = 0.0;
    comp = is_slip_x ? -mu.x[] : mu.x[];
  }

  double dvdx = mu.x[1] * u.y[1] + mu.x[] * u.y[-1] * m1; //contribution from dvdx in tauxy
  coef = coef + mu.x[1] + mu.x[] + comp; //diagonal terms from dvdx in tauxy
  coef = coef * dt / rho[] + sq(Delta)*lambda.y;


  double dudy_xp = (u.x[0, 1] + u.x[1, 1]) / 4. - (u.x[0, -1] + u.x[1, -1]) / 4.;
  double dudy_xm = (u.x[-1, 1] + u.x[0, 1]) / 4. - (u.x[-1, -1] + u.x[0, -1]) / 4.;

  if(coordm <= 0.0)
  {
    dudy_xm = 0.0; //the solid-liquid boundary (u = 0)
  }
  double tauxy = mu.x[1] * dudy_xp - mu.x[] * dudy_xm + dvdx;          // contribution from d tauxy dx


  double res = dt / rho[] * (tauyy + tauxy) + r.y[] * sq(Delta); //the second term is the contribution from advection
  w.y[] = res / coef;
}

static void relaxVelXSolidY(Point point, const face vector mu, const scalar rho, 
                        const double dt, const vector r, vector u, vector w)
{
  if((int) is_solid[] == 1)
  {
    u.x[] = 0.0;
    w.x[] = 0.0;
    return;
  }

  //this part will not be influenced by the existence of the solid
  double coef = 0.0;
  //contribution from tauxx
  double tauxx = 2. * mu.x[1] * u.x[1] + 2. * mu.x[] * u.x[-1];
  coef = coef +  2. * mu.x[1] + 2. * mu.x[]; //diagonal terms from tauxx

  double coordm = y - Delta;
  double m1 = 1.0;
  double comp = 0.0; //compensating coef
  if(coordm <= 0.0)
  {
    m1 = 0.0;
    comp = is_slip_y ? -mu.y[] : mu.y[];
  }

  double dudy = mu.y[0, 1] * u.x[0, 1] + mu.y[] * u.x[0, -1] * m1; //contribution from dudy in tauxy
  coef = coef + mu.y[0, 1] + mu.y[] + comp; //diagonal terms from dudy in tauxy
  coef = coef * dt / rho[] + sq(Delta)*lambda.x;


  double dvdx_yp = (u.y[1, 0] + u.y[1, 1]) / 4. - (u.y[-1, 0] + u.y[-1, 1]) / 4.;
  double dvdx_ym = (u.y[1, -1] + u.y[1, 0]) / 4. - (u.y[-1, -1] + u.y[-1, 0]) / 4.;

  if(coordm <= 0.0)
  {
    dvdx_ym = 0.0; //the solid-liquid boundary (v = 0)
  }
  double tauxy = mu.y[0, 1] * dvdx_yp - mu.y[] * dvdx_ym + dudy; //contribution from d tauxy dy
  double res = dt / rho[] * (tauxx + tauxy) + r.x[] * sq(Delta); //the second term is the contribution from advection
  w.x[] = res / coef;
}

static void relaxVelYSolidY(Point point, const face vector mu, const scalar rho, 
                            const double dt, const vector r, vector u, vector w)
{
  if((int) is_solid[] == 1)
  {
    u.y[] = 0.0;
    w.y[] = 0.0;
    return;
  }

  double coordm = y - Delta;
  double m1 = 1.0;
  double comp = 0.0; //compensating coef
  if(coordm <= 0.0)
  {
    m1 = 0.0;
    comp = 2.0 * mu.y[];
  }

  double coef = 0.0;
  //contribution from tauyy
  double tauyy = 2. * mu.y[0, 1] * u.y[0, 1] + 2. * mu.y[] * u.y[0, -1] * m1;
  coef = coef + 2. * mu.y[0, 1] + 2. * mu.y[] + comp;                // diagonal terms from tauyy

  //this part will not be influenced by the existence of the solid
  double dvdx = mu.x[1] * u.y[1] + mu.x[] * u.y[-1]; //contribution from dvdx in tauxy
  coef = coef + mu.x[1] + mu.x[]; //diagonal terms from dvdx in tauxy
  coef = coef * dt / rho[] + sq(Delta)*lambda.y;

  double uxymxp = (u.x[0, -1] + u.x[1, -1]) / 4.;
  double uyymxm = (u.x[-1, -1] + u.x[0, -1]) / 4.;

  if (coordm <= 0.0)
  {
    uxymxp = is_slip_y ? (u.x[0, 0] + u.x[1, 0]) / 4. : -(u.x[0, 0] + u.x[1, 0]) / 4.;
    uyymxm = is_slip_y ? (u.x[-1, 0] + u.x[0, 0]) / 4. : -(u.x[-1, 0] + u.x[0, 0]) / 4.;
  }

  double dudy_xp = (u.x[0, 1] + u.x[1, 1]) / 4. - uxymxp;
  double dudy_xm = (u.x[-1, 1] + u.x[0, 1]) / 4. - uyymxm;
  double tauxy = mu.x[1] * dudy_xp - mu.x[] * dudy_xm + dvdx;          // contribution from d tauxy dx
  double res = dt / rho[] * (tauyy + tauxy) + r.y[] * sq(Delta); //the second term is the contribution from advection
  w.y[] = res / coef;
}

static void relaxViscousitySolid (scalar * a, scalar * b, int l, void * data)
{
  struct Viscosity * p = (struct Viscosity *) data;
  (const) face vector mu = p->mu;
  (const) scalar rho = p->rho;
  double dt = p->dt;
  vector u = vector(a[0]), r = vector(b[0]);

#if JACOBI
  vector w[];
#else
  vector w = u;
#endif

  foreach_level_or_leaf(l)
  {
    relaxVelXSolid(point, mu, rho, dt, r, u, w);
    relaxVelYSolid(point, mu, rho, dt, r, u, w);
  }

#if JACOBI
  foreach_level_or_leaf (l)
    foreach_dimension()
      u.x[] = (u.x[] + 2.*w.x[])/3.;
#endif
  
#if TRASH
  vector u1[];
  foreach_level_or_leaf (l)
    foreach_dimension()
      u1.x[] = u.x[];
  trash ({u});
  foreach_level_or_leaf (l)
    foreach_dimension()
      u.x[] = u1.x[];
#endif
}

#endif

static void relax_viscosity (scalar * a, scalar * b, int l, void * data)
{
  struct Viscosity * p = (struct Viscosity *) data;
  (const) face vector mu = p->mu;
  (const) scalar rho = p->rho;
  double dt = p->dt;
  vector u = vector(a[0]), r = vector(b[0]);

#if JACOBI
  vector w[];
#else
  vector w = u;
#endif
  
  foreach_level_or_leaf (l) {
    foreach_dimension()
      w.x[] = (dt/rho[]*(2.*mu.x[1]*u.x[1] + 2.*mu.x[]*u.x[-1]
               #if dimension > 1
			   + mu.y[0,1]*(u.x[0,1] +
					(u.y[1,0] + u.y[1,1])/4. -
					(u.y[-1,0] + u.y[-1,1])/4.)
			   - mu.y[]*(- u.x[0,-1] +
				     (u.y[1,-1] + u.y[1,0])/4. -
				     (u.y[-1,-1] + u.y[-1,0])/4.)
               #endif
	       #if dimension > 2
			   + mu.z[0,0,1]*(u.x[0,0,1] +
					  (u.z[1,0,0] + u.z[1,0,1])/4. -
					  (u.z[-1,0,0] + u.z[-1,0,1])/4.)
			   - mu.z[]*(- u.x[0,0,-1] +
				     (u.z[1,0,-1] + u.z[1,0,0])/4. -
				     (u.z[-1,0,-1] + u.z[-1,0,0])/4.)
               #endif
			   ) + r.x[]*sq(Delta))/
    (sq(Delta)*lambda.x + dt/rho[]*(2.*mu.x[1] + 2.*mu.x[]
                                    #if dimension > 1
				      + mu.y[0,1] + mu.y[]
                                    #endif
			            #if dimension > 2
				      + mu.z[0,0,1] + mu.z[]
			            #endif
			     ));
  }

#if JACOBI
  foreach_level_or_leaf (l)
    foreach_dimension()
      u.x[] = (u.x[] + 2.*w.x[])/3.;
#endif
  
#if TRASH
  vector u1[];
  foreach_level_or_leaf (l)
    foreach_dimension()
      u1.x[] = u.x[];
  trash ({u});
  foreach_level_or_leaf (l)
    foreach_dimension()
      u.x[] = u1.x[];
#endif
}

static double residual_viscosity (scalar * a, scalar * b, scalar * resl, 
				  void * data)
{
  struct Viscosity * p = (struct Viscosity *) data;
  (const) face vector mu = p->mu;
  (const) scalar rho = p->rho;
  double dt = p->dt;
  vector u = vector(a[0]), r = vector(b[0]), res = vector(resl[0]);
  double maxres = 0.;
#if TREE
  /* conservative coarse/fine discretisation (2nd order) */

  /**
  We manually apply boundary conditions, so that all components are
  treated simultaneously. Otherwise (automatic) BCs would be applied
  component by component before each foreach_face() loop. */
  
  boundary ({u});
#if USE_MY_SOLID
  boundarySolidVelCNoauto(u);
#endif

  foreach_dimension() {
    face vector taux[];
    foreach_face(x)
      taux.x[] = 2.*mu.x[]*(u.x[] - u.x[-1])/Delta;
    #if dimension > 1
      foreach_face(y)
	taux.y[] = mu.y[]*(u.x[] - u.x[0,-1] + 
			   (u.y[1,-1] + u.y[1,0])/4. -
			   (u.y[-1,-1] + u.y[-1,0])/4.)/Delta;
    #endif
    #if dimension > 2
      foreach_face(z)
	taux.z[] = mu.z[]*(u.x[] - u.x[0,0,-1] + 
			   (u.z[1,0,-1] + u.z[1,0,0])/4. -
			   (u.z[-1,0,-1] + u.z[-1,0,0])/4.)/Delta;
    #endif
    foreach (reduction(max:maxres)) {
      double d = 0.;
      foreach_dimension()
	d += taux.x[1] - taux.x[];
      res.x[] = r.x[] - lambda.x*u.x[] + dt/rho[]*d/Delta;
#if USE_MY_SOLID
      res.x[] *= (1.0 - is_solid[]);
#endif
      if (fabs (res.x[]) > maxres)
	maxres = fabs (res.x[]);
    }
  }
#else
#if USE_MY_SOLID
  boundarySolidVelCNoauto(u);
#endif
  /* "naive" discretisation (only 1st order on trees) */
  foreach (reduction(max:maxres))
    foreach_dimension() {
      res.x[] = r.x[] - lambda.x*u.x[] +
        dt/rho[]*(2.*mu.x[1,0]*(u.x[1] - u.x[])
		  - 2.*mu.x[]*(u.x[] - u.x[-1])
        #if dimension > 1
		  + mu.y[0,1]*(u.x[0,1] - u.x[] +
			       (u.y[1,0] + u.y[1,1])/4. -
			       (u.y[-1,0] + u.y[-1,1])/4.)
		  - mu.y[]*(u.x[] - u.x[0,-1] +
			    (u.y[1,-1] + u.y[1,0])/4. -
			    (u.y[-1,-1] + u.y[-1,0])/4.)
	#endif
        #if dimension > 2
		  + mu.z[0,0,1]*(u.x[0,0,1] - u.x[] +
				 (u.z[1,0,0] + u.z[1,0,1])/4. -
				 (u.z[-1,0,0] + u.z[-1,0,1])/4.)
		  - mu.z[]*(u.x[] - u.x[0,0,-1] +
			    (u.z[1,0,-1] + u.z[1,0,0])/4. -
			    (u.z[-1,0,-1] + u.z[-1,0,0])/4.)
	#endif
		  )/sq(Delta);
#if USE_MY_SOLID
      res.x[] *= (1.0 - is_solid[]);
#endif
      if (fabs (res.x[]) > maxres)
	maxres = fabs (res.x[]);
    }
#endif
  return maxres;
}

#undef lambda

trace
mgstats viscosity (struct Viscosity p)
{
  vector u = p.u, r[];
  foreach()
    foreach_dimension()
      r.x[] = u.x[];

  face vector mu = p.mu;
  scalar rho = p.rho;
  restriction ({mu,rho});
  
#if USE_MY_SOLID
  return mg_solve ((scalar *){u}, (scalar *){r},
		   residual_viscosity, relaxViscousitySolid, &p, p.nrelax, p.res);
#else
  return mg_solve ((scalar *){u}, (scalar *){r},
		   residual_viscosity, relax_viscosity, &p, p.nrelax, p.res);
#endif
}

trace
mgstats viscosity_explicit (struct Viscosity p)
{
  vector u = p.u, r[];
  mgstats mg = {0};
  mg.resb = residual_viscosity ((scalar *){u}, (scalar *){u}, (scalar *){r}, &p);
  foreach()
    foreach_dimension()
      u.x[] += r.x[];
  return mg;
}

#endif //_MY_VISCOSITY_H