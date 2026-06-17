/**
# Flow around a 2D circular cylinder with a square notch, containing a gas bubble.
This case has no direct application, and represents a proof of concept.
*/

#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "popinet/contact/contact-embed.h"
#include "view.h"

#define T 50.
#define tMovie 0.1

double R0 = 0.5;
double XC = 0.;
double YC = 0.;
double NS = 0.1; // notch size
double THETA = pi/2.; // notch tangential location


/**
Solid Geometry. A cylinder with a square notch.*/
static inline double cylinder_notch (double x, double y,
                                      double x_cylinder, double y_cylinder,
                                      double cylinder_radius,
                                      double notch_size,
                                      double notch_location)
{
  double ct = cos(notch_location), st = sin(notch_location);

  // Square center on the original circle boundary
  double xs = x_cylinder + (cylinder_radius-notch_size/2.)*ct;
  double ys = y_cylinder + (cylinder_radius-notch_size/2.)*st;

  // Local square coordinates:
  // X is radial, Y is tangential.
  double X =  (x - xs)*ct + (y - ys)*st;
  double Y = -(x - xs)*st + (y - ys)*ct;

  // Negative inside the circular solid
  double phi_circle =
    sq(x - x_cylinder) + sq(y - y_cylinder) - sq(cylinder_radius);

  // Negative inside the square notch
  double phi_square =
    max(fabs(X), fabs(Y)) - notch_size/2.;

  // Cylinder minus square
  return max(phi_circle, -phi_square);
}
/**
Gas geometry. A square notch within a cylinder.*/
static inline double notch_bubble (double x, double y,
                                   double x_cylinder, double y_cylinder,
                                   double cylinder_radius,
                                   double notch_size,
                                   double notch_location)
{
  double ct = cos(notch_location), st = sin(notch_location);

  // Same square center as in cylinder_notch()
  double xs = x_cylinder + (cylinder_radius - notch_size/2.)*ct;
  double ys = y_cylinder + (cylinder_radius - notch_size/2.)*st;

  // Local square coordinates
  double X =  (x - xs)*ct + (y - ys)*st;
  double Y = -(x - xs)*st + (y - ys)*ct;

  // Positive inside the original cylinder
  double phi_inside_circle =
    sq(cylinder_radius) - sq(x - x_cylinder) - sq(y - y_cylinder);

  // Positive inside the square notch
  double phi_inside_square =
    notch_size/2. - max(fabs(X), fabs(Y));

  // Bubble = intersection of circle and square
  return min(phi_inside_circle, phi_inside_square);
}

/**
Also, to avoid unphysical transient at startups, I implement
a simple ramp for external velocity.*/
static inline double ramp (double t, double t_ramp, double value)
{
  return value*min(t/t_ramp, 1.);
}






int maxlevel = 12;
int minlevel = 4;
double Reynolds = 40.;


int main() {
  size (32.);

  /**
  We set the origin */

  origin (-L0/2., -L0/2.);

  init_grid (1 << (maxlevel-4));


  /**
  We use a constant viscosity. */

	mu2 = 1.*R0*2./Reynolds; // outer fluid viscosity 
	/**
	We set plastron viscosity equal to 1/100 of the main fluid one.*/
  mu1 = 0.01*mu2; // plastron viscosity 
  /**
  We set the surface tension coefficient. */
  
  f.sigma = 0.01;

	/**
	Set a constant contact angle.*/
	const scalar c[] = 90.*pi/180.; // fixed contact angle...
	contact_angle = c;

	run();
}

/**
We set the boundary conditions, so to obtain a flow around a fixed cylinder. */
u.n[left]  = dirichlet(ramp(t,1.,1.));
p[left]    = neumann(0.);
pf[left]   = neumann(0.);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

/**
Must impose no-slip on embedded boundaries!*/
u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);

event init (t = 0)
{
/**
Initialize solid geometry, with automatic refinement.*/
#if TREE
  astats s;
  int iter = 0;

  do {
    iter++;

    // Recompute embedded fractions on the current grid
		solid (cs, fs, cylinder_notch (x, y, XC, YC,R0,NS,THETA));

    // Force refinement where cs varies, i.e. near the embedded boundary
    s = adapt_wavelet ({cs},
                       (double[]){1.e-30},
                       maxlevel,
                       minlevel);

  } while ((s.nf || s.nc) && iter < 100);
#endif
/**
Initialize bubble location.*/

  fraction (f,notch_bubble (x, y,XC,YC,R0,NS*1.1,THETA));
}

event adapt (i++) {
  adapt_wavelet ({cs,u,f}, (double[]){1e-3,1e-2,1e-2,1e-3}, maxlevel, minlevel);
	solid (cs, fs, cylinder_notch (x, y, XC, YC,R0,NS,THETA));
	/**
	Unrefine all uninteresting areas...such as the inlet!*/
	unrefine (x < -L0/3.);

	/**
	Purge gas that is (somehow) trapped within the solid phase.*/
	foreach()
		if(cs[]<=0)
			f[] = 0.;

}

event movie(t+=tMovie,t<=T){
  view(fov=0.4, tx = 0., ty = -R0/L0);
  draw_vof("cs", "fs",filled = -1);
	draw_vof ("f", lc = {1, 0, 0}, lw = 3);
//	squares ("u.x", linear = true);
	// Draw grid only on upper part of flow
	cells (lc = {0.7, 0.7, 0.7});

	  vectors ("u", scale = 0.01);


  save("movie.mp4");
}
/**
![Circular cylinder + single notch + gas bubble](cylinder_single_notch_bubble/movie.mp4)(width="60%")
*/

// fixme: Comparison to theory is missing (add soon) 
//
//
event compute_forces (i += 1, t<=T)
{
  coord Fp, Fmu;
  embed_force (p, u, mu, &Fp, &Fmu);
// Output forces
  static FILE * fp = NULL;
	if(pid()==0){
  	if (!fp) {
  	  fp = fopen ("force_coeff.dat", "w");
  	  fprintf (fp, "# i t dt Cx Cy\n");
  	}
		double Cx = (Fp.x+Fmu.x)*2.;
		double Cy = (Fp.y+Fmu.y)*2.;
  	fprintf (fp, "%06d %+6.5e %+6.5e %+6.5e %+6.5e\n",
  	         i, t, dt, Cx, Cy);
  	fflush (fp);
	}
}


/**
~~~pythonplot Force coefficient
import numpy as np
import matplotlib.pyplot as plt

# columns:
# i t dt Cx Cy
i,t,dt,Cx,Cy = np.loadtxt('force_coeff.dat', comments='#', unpack=True)

tTarget = 10.

m = t > tTarget
tt = t[m]
xx = Cx[m]

Cxmean = np.mean(xx)

# local maxima/minima around the mean
imax = np.where((xx[1:-1] > xx[:-2]) & (xx[1:-1] > xx[2:]))[0] + 1
imin = np.where((xx[1:-1] < xx[:-2]) & (xx[1:-1] < xx[2:]))[0] + 1

Cxplus  = np.mean(xx[imax] - Cxmean)
Cxminus = np.mean(Cxmean - xx[imin])

print('Average Cx for t > %.1f = %+6.5e' % (tTarget, Cxmean))
print('Average positive peak   = %+6.5e' % Cxplus)
print('Average negative peak   = %+6.5e' % Cxminus)
print('Cx = %+6.5e +%+6.5e -%+6.5e' % (Cxmean, Cxplus, Cxminus))

# literature reference
Cxref = 1.55
Cxerr = 0.012

plt.figure()
plt.plot(t, Cx, 'k-', label=r'$C_x$')

plt.axhline(Cxmean, color='k', linestyle='--',
            label=r'$\langle C_x\rangle$')

plt.axhline(Cxref, color='C1', linestyle='-',
            label=r'$C_{x,\mathrm{ref}}$')

plt.fill_between(t, Cxref - Cxerr, Cxref + Cxerr,
                 color='C1', alpha=0.25,
                 label=r'$C_{x,\mathrm{ref}}\pm 0.012$')

plt.xlabel(r'$t$')
plt.ylabel(r'$C_x$')

plt.text(0.97, 0.97,
         r'$t > %.1f$' '\n'
         r'$\langle C_x\rangle = %.5e$' '\n'
         r'$+\Delta C_x = %.5e$' '\n'
         r'$-\Delta C_x = %.5e$'
         % (tTarget, Cxmean, Cxplus, Cxminus),
         transform=plt.gca().transAxes,
         ha='right', va='top',
         bbox=dict(facecolor='white', alpha=0.85))

plt.ylim(0.0, 2.0)
plt.legend(loc='lower right')
plt.tight_layout()

plt.savefig('Cx_evolution.svg')
~~~
*/

/**
~~~pythonplot Force coefficient Cy
import numpy as np
import matplotlib.pyplot as plt

# columns:
# i t dt Cx Cy
i,t,dt,Cx,Cy = np.loadtxt('force_coeff.dat', comments='#', unpack=True)

tTarget = 10.

m = t > tTarget
tt = t[m]
yy = Cy[m]

Cymean = np.mean(yy)

# local maxima/minima around the mean
imax = np.where((yy[1:-1] > yy[:-2]) & (yy[1:-1] > yy[2:]))[0] + 1
imin = np.where((yy[1:-1] < yy[:-2]) & (yy[1:-1] < yy[2:]))[0] + 1

Cyplus  = np.mean(yy[imax] - Cymean)
Cyminus = np.mean(Cymean - yy[imin])

print('Average Cy for t > %.1f = %+6.5e' % (tTarget, Cymean))
print('Average positive peak   = %+6.5e' % Cyplus)
print('Average negative peak   = %+6.5e' % Cyminus)
print('Cy = %+6.5e +%+6.5e -%+6.5e' % (Cymean, Cyplus, Cyminus))

# literature reference
Cyref = 0.
Cyosc_ref = 0.333

plt.figure()
plt.plot(t, Cy, 'k-', label=r'$C_y$')

plt.axhline(Cymean, color='k', linestyle='--',
            label=r'$\langle C_y\rangle$')

plt.axhline(Cyref, color='C1', linestyle='-',
            label=r'$C_{y,\mathrm{ref}}$')

plt.fill_between(t, Cyref - Cyosc_ref, Cyref + Cyosc_ref,
                 color='C1', alpha=0.25,
                 label=r'$C_{y,\mathrm{ref}}\pm 0.333$')

plt.xlabel(r'$t$')
plt.ylabel(r'$C_y$')

plt.text(0.03, 0.97,
         r'$t > %.1f$' '\n'
         r'$\langle C_y\rangle = %.5e$' '\n'
         r'$+\Delta C_y = %.5e$' '\n'
         r'$-\Delta C_y = %.5e$'
         % (tTarget, Cymean, Cyplus, Cyminus),
         transform=plt.gca().transAxes,
         ha='left', va='top',
         bbox=dict(facecolor='white', alpha=0.85))

plt.ylim(-0.5, 0.5)
plt.legend(loc='lower left')
plt.tight_layout()

plt.savefig('Cy_evolution.svg')
~~~
*/