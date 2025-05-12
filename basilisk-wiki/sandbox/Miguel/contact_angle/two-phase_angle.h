/**
# Two-phase interfacial flows

This file takes vof_angle.h and heights_angle.h to implement the boundary conditions. 

This file helps setup simulations for flows of two fluids separated by
an interface (i.e. immiscible fluids). It is typically used in
combination with a [Navier--Stokes solver](navier-stokes/centered.h). 

The interface between the fluids is tracked with a Volume-Of-Fluid
method. The volume fraction in fluid 1 is $f=1$ and $f=0$ in fluid
2. The densities and dynamic viscosities for fluid 1 and 2 are *rho1*,
*mu1*, *rho2*, *mu2*, respectively. */

#include "vof_angle.h"
#include "heights_angle.h"

scalar f[], * interfaces = {f};
double rho1 = 1., mu1 = 0., rho2 = 1., mu2 = 0.;

/**
Auxilliary fields are necessary to define the (variable) specific
volume $\alpha=1/\rho$ as well as the cell-centered density. */

face vector alphav[];
scalar rhov[];

event defaults (i = 0) {
  alpha = alphav;
  rho = rhov;

  /**
  If the viscosity is non-zero, we need to allocate the face-centered
  viscosity field. */
  
  if (mu1 || mu2)
    mu = new face vector;
}

/**
The density and viscosity are defined using arithmetic averages by
default. The user can overload these definitions to use other types of
averages (i.e. harmonic). */

#ifndef rho
# define rho(f) (clamp(f,0.,1.)*(rho1 - rho2) + rho2)
#endif
#ifndef mu
# define mu(f)  (clamp(f,0.,1.)*(mu1 - mu2) + mu2)
#endif

/**
We have the option of using some "smearing" of the density/viscosity
jump. */

#ifdef FILTERED
scalar sf[];
#else
# define sf f
#endif

event properties (i++) {

  /**
  When using smearing of the density jump, we initialise *sf* with the
  vertex-average of *f*. */

#ifndef sf
#if dimension <= 2
  foreach()
    sf[] = (4.*f[] + 
	    2.*(f[0,1] + f[0,-1] + f[1,0] + f[-1,0]) +
	    f[-1,-1] + f[1,-1] + f[1,1] + f[-1,1])/16.;
#else // dimension == 3
  foreach()
    sf[] = (8.*f[] +
	    4.*(f[-1] + f[1] + f[0,1] + f[0,-1] + f[0,0,1] + f[0,0,-1]) +
	    2.*(f[-1,1] + f[-1,0,1] + f[-1,0,-1] + f[-1,-1] + 
		f[0,1,1] + f[0,1,-1] + f[0,-1,1] + f[0,-1,-1] +
		f[1,1] + f[1,0,1] + f[1,-1] + f[1,0,-1]) +
	    f[1,-1,1] + f[-1,1,1] + f[-1,1,-1] + f[1,1,1] +
	    f[1,1,-1] + f[-1,-1,-1] + f[1,-1,-1] + f[-1,-1,1])/64.;
#endif
#endif 

#if TREE
  sf.prolongation = refine_bilinear;
  boundary ({sf});

  vector h[];
  heights(f, h);
  double theta = theta_0;
  foreach_boundary(left)
	{
        h.y[0,-1] = h.y[] + (1./tan(theta));
    	if (f[] < 1. && f[] > 0.)
      		{
			if(theta <= pi/2.)
                          {
				f[-1,-1] = 1.;
                        	double contactline;
				double xghost;
				double x1;
				x1 =((y + (Delta*height(h.y[])))/Delta) +  (((cos(theta)/sin(theta))/2.));
				contactline = x1 - floor(y/Delta);                   
				xghost = contactline + ((cos(theta)/sin(theta)));
				if(xghost < 1.)
				{
				f[-1, 1] = 0.; 
				f[-1, 0] = contactline + ((cos(theta)/sin(theta))/2.);


				}
				else if (xghost >= 1.)
				{
				f[-1, 1] = ((sq(xghost-1))*tan(theta)/2.);
				f[-1, 0] = 1. - ((sq(1-contactline))*tan(theta)/2.);
				}
				else if (contactline >= 1.)
				{
				
				f[-1, 1] = ((sq(2. - contactline))*(tan(theta)/2));
				f[-1, 0] = 1.;
			        f[-1, 2] = (sq(xghost-2.))*tan(theta)/2.;
				f[0, 1]  = (sq(contactline - 1.))*tan(theta)/2.;
				

				}
                   
			}

			else if(theta > pi/2)
			{
				f[-1, 1] = 0;
                        	double contactline;
				double xghost;
				double x1;
				double phi = pi - theta;
				double x0;
				x0 = ((y + (Delta*height(h.y[])))/Delta) - floor(y/Delta);
				x1 =((y + (Delta*height(h.y[])))/Delta) -  (((cos(phi)/sin(phi))/2.));
				contactline = x1 - floor(y/Delta);                   
				xghost = contactline - ((cos(phi)/sin(phi)));
				if (xghost> 0)
				{
				f[-1, -1] = 1;
				f[-1, 0] = 1 -(contactline + xghost)/2;
				}
				else if (xghost <= 0. && contactline>0.)
				{
				f[-1, 0] = ((sq(contactline)*tan(phi))/2.);
				f[-1, -1] = 1 - ((sq(contactline)*tan(phi))/2.);
				}
				else if (contactline <= 0.)
				{
				f[-1,0] = 0;
				f[-1, -1] = ((sq(contactline)*tan(phi))/2.);
				f[0, -1] = 1 - ((sq(contactline)*tan(phi))/2.);
				f[] = f[] - ((sq(x0)*tan(phi)) + (sq(x0)*tan(phi)));
				}
			}
      		 }
	}
#endif
  
  foreach_face() {
    double ff = (sf[] + sf[-1])/2.;
    alphav.x[] = fm.x[]/rho(ff);
    if (mu1 || mu2) {
      face vector muv = mu;
      muv.x[] = fm.x[]*mu(ff);
    }
  }
  foreach()
    rhov[] = cm[]*rho(sf[]);

#if TREE  
  sf.prolongation = fraction_refine;
  boundary ({sf});

  foreach_boundary(left)
	{
        h.y[0,-1] = h.y[] + (1./tan(theta));
    	if (f[] < 1. && f[] > 0.)
      		{
			if(theta <= pi/2.)
                          {
				f[-1,-1] = 1.;
                        	double contactline;
				double xghost;
				double x1;
				x1 =((y + (Delta*height(h.y[])))/Delta) +  (((cos(theta)/sin(theta))/2.));
				contactline = x1 - floor(y/Delta);                   
				xghost = contactline + ((cos(theta)/sin(theta)));
				if(xghost < 1.)
				{
				f[-1, 1] = 0.; 
				f[-1, 0] = contactline + ((cos(theta)/sin(theta))/2.);


				}
				else if (xghost >= 1.)
				{
				f[-1, 1] = ((sq(xghost-1))*tan(theta)/2.);
				f[-1, 0] = 1. - ((sq(1-contactline))*tan(theta)/2.);
				}
				else if (contactline >= 1.)
				{
				
				f[-1, 1] = ((sq(2. - contactline))*(tan(theta)/2));
				f[-1, 0] = 1.;
			        f[-1, 2] = (sq(xghost-2.))*tan(theta)/2.;
				f[0, 1]  = (sq(contactline - 1.))*tan(theta)/2.;
				

				}
                   
			}

			else if(theta > pi/2)
			{
				f[-1, 1] = 0;
                        	double contactline;
				double xghost;
				double x1;
				double phi = pi - theta;
				double x0;
				x0 = ((y + (Delta*height(h.y[])))/Delta) - floor(y/Delta);
				x1 =((y + (Delta*height(h.y[])))/Delta) -  (((cos(phi)/sin(phi))/2.));
				contactline = x1 - floor(y/Delta);                   
				xghost = contactline - ((cos(phi)/sin(phi)));
				if (xghost> 0)
				{
				f[-1, -1] = 1;
				f[-1, 0] = 1 -(contactline + xghost)/2;
				}
				else if (xghost <= 0.) //&& contactline>0.)
				{
				f[-1, 0] = ((sq(contactline)*tan(phi))/2.);
				f[-1, -1] = 1 - ((sq(contactline)*tan(phi))/2.);
				}
				else if (contactline <= 0.)
				{
				f[-1,0] = 0;
				f[-1, -1] = ((sq(contactline)*tan(phi))/2.);
				f[0, -1] = 1 - ((sq(contactline)*tan(phi))/2.);
				f[] = f[] - ((sq(x0)*tan(phi)) + (sq(x0)*tan(phi)));
				}
			}
      		 }
	}
#endif
}