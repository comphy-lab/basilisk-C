#include "axi.h"
#include "navier-stokes/centered.h"
#include "contact.h"
#include "two-phase.h"
#include "tension.h"

#define LARGE 1e36

double max_level = 8;
double min_level = 3;
double L = 2.;
double t_out = 0.1;
double t_end = 0.11;
double theta0 = 180;
double amp = 0.e-2;
double omega = 0.277;
double dia = 5.e-3;
double gravity = 0.5;

double a0=11.02263;
double a1=-10.0165031;
double a2=-11.26938;
double a3=11.62171;
double a4=-11.39031;
double x_max=10.559;
double y_max=0.884192277431746;

double amp_pert  = 0.005;
double time_pert = 0.05;

double femax = 0.001;
double uemax = 0.001;

double gx;  //dimensionless gravity in x direction  

double percentage=0.5;
double init_amp=0.;
double rhol=1000;
double rhog=1.2;
double mul=1.e-3;
double mug=1.8e-5;
double sigma=0.07;
                                
double fcl,ycl;

/**
To set the contact angle, we need the [height-function vector
field](/src/heights.h). The contact angle boundary condition is specified based on tangential height. */
vector h[];
h.t[left] = contact_angle ( y < y_max ? 15.*pi/180. : 165.*pi/180.);

/** 
No slip boundary condition. */
u.t[left]  = dirichlet(0.);
u.n[left]  = dirichlet(0.);

int main(int argc, char * argv[])
{

  size (L);
  init_grid (64);

  f.height = h;

  /**
  The liquid phase is water, rho_l=1000 kg/m3, mu_l=1e-3 Ps s; 
  the gas phase is air, rho_g = 1.2 kg/m3, mu_g = 1.7e-5 Pa s; 
  surface tension is sigma=0.07 N/m; 
  gravity is along y direction, g = 9.8 m/s2.
  The dimensionless parameters can then be computed based on 
  rho_l, sigma, and dia. */

  rho1 = rhol/rhol, rho2 = rhog/rhol;
  mu1 = mul/sqrt(rhol*dia*sigma), mu2 = mug/sqrt(rhol*dia*sigma);
  f.sigma = 1.;

  gx  = -percentage*gravity*rhol*sq(dia)/sigma;

  run();
}

/**
The initial drop shape is first obtained by equilibrium sessile drop theory. Then the drop contour is fitted using a polynomial function. */

event init (t = 0)
{
  if (!restore (file = "dump")) {

    refine ( y < 1.5 && x < 1.5 && level < max_level);
    refine ( (y > y_max-1e-2) && (y < y_max+1e-2) && (x < 1e-2) && level < max_level);
    fraction (f, -y
            + (a0*pow(-x+x_max,0.5)
            + a1*pow(-x+x_max,1)
            + a2*pow(-x+x_max,2)
            + a3*pow(-x+x_max,3)
            + a4*pow(-x+x_max,4)));
  }
}

event logfile (i+=10)
{
  scalar posy[],posx[],posycl[];
  position (f,posx,{1,0});
  position (f,posy,{0,1});
  position (f,posycl,{0,1});
  double area=0.,vol=0.,ke=0.,ud=0.,xd=0.;

  foreach(reduction(+:area) reduction(+:vol) reduction(+:ke)
          reduction(+:ud) reduction(+:xd)) {
    if (f[] <= 1e-6 || f[] >= 1. - 1e-6)  
     {
      posx[] = nodata;
      posy[] = nodata;
    }

    if (f[] <= 1e-6 || f[] >= 1. - 1e-6 || x > L/pow(2,max_level) ){
      posycl[] = nodata;
    }

    // statistics in axisymmetric geometry 
    if (f[] > 1e-6 && f[] < 1. - 1e-6) {
      //interfacial area
      coord n = interface_normal (point, f), p;
      double alpha = plane_alpha (f[], n);
      area += pow(Delta, 1.)*plane_area_center (n, alpha, &p)*2.*pi*posy[];
    }

    double dv_axi = pow(Delta, 2.)*2.*pi*y;

    //volume
    if (f[] > 1e-6 ) {
      vol += dv_axi*f[];

      /** Set the velocity and kinetic energy near the contact line to be zero. This is important to clean up the spurious oscillations induced by the adrupt change of contact angle near the contact line location.*/
      foreach_dimension() {
         if ((y > y_max-1e-2) && (y < y_max+1e-2) && (x < 1e-2)) {
            ke += 0.;
            u.x[] = 0.;
         } else {
            ke += dv_axi*sq(u.x[]);
         }
      }

      //mean velocity
      ud += dv_axi*f[]*u.x[];

      //centroid
      xd += dv_axi*f[]*x;
    }
  }
  ke /= 2.;
  ud /= vol;
  xd /= vol;

  fprintf (ferr, "%g %g %g %g %g %g %g %g %g %g %g %ld %g %g \n", t, dt, xd, ud,
      statsf(posy).max, statsf (posx).max, statsf(f).sum,
      area, statsf(posycl).max,vol,ke, grid->tn, perf.t, perf.speed);

}

/** output drop interface (more often) and the whole snapshots (less often) */
event interface (t += t_out/10; t<= t_end) {

  char name[80];
  sprintf (name, "infc-%05.0f.dat", t*100);
  FILE * fp1 = fopen (name, "w");
  output_facets (f,fp1);

}

event snapshot (t += t_out; t <= t_end ) {
  char name[80];

  sprintf (name, "dump-%07.5f", t);
  dump (file = name);
}

/**
Adapt mesh based on the volume fraction and velocity. */
event adapt (i=10; i++) {
  adapt_wavelet ({f,u}, (double[]){femax,uemax,uemax}, minlevel = min_level, maxlevel = max_level);

#if 1
  refine ((y > y_max-1e-2) && (y < y_max+1e-2) && (x < 1e-2) && level < max_level);
  unrefine ( (y>1.5 || x>1.5) && level >= 3 );
#endif
}

/**
We add the acceleration of gravity in (-x)
direction. */

event acceleration (i++) {
  face vector av = a;
  foreach_face(x)
#if 1
  av.x[] = (t<time_pert?gx:-(percentage-0.5)*gravity*rhol*sq(dia)/sigma);
#else
  av.x[] = gx+amp*noise();
#endif
}
