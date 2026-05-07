/**
# Flow around a 2D cylinder at Re = 100
Some quantitative plots, drag (Cx) and lift (Cy).
With some simple plotting options.
Largely taken from [karman.c](https://basilisk.fr/src/examples/karman.c) case.

Reference values are taken from [Legendre et al. 2009](https://arxiv.org/abs/0905.0648).
*/

/**
![Animation of the vorticity field.](cylinder_plastron/vort.mp4)(loop)
*/

//#include "grid/multigrid.h"
#include "grid/quadtree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
//#include "tension.h"
//#include "tavares/contact-embed.h"
#include "view.h"

#define R0 0.5 // solid cylinder radius
//#define R1 0.55 // plastron cylinder radius
#define xc 0.
#define yc 1. 
#define endTime 200.

double theta0, volume_vof_init;
int MAXLEVEL = 9;
int MINLEVEL = 4;

double Reynolds = 100.;

int main() {
  size (32.);

  /**
  We set the origin */

  origin (-8, -L0/2.); // at least 8 diameters beyond the inlet face
											 // so to avoid influence of confinement...

  init_grid (1 << MAXLEVEL);

  /**
  We use a constant viscosity. */

	mu2 = 1.*R0*2./Reynolds; // outer fluid viscosity 
	/**
	We set plastron viscosity equal to 1/100 of the main fluid one.*/
  mu1 = 1.*mu2; // plastron viscosity 
  /**
  We set the surface tension coefficient. */
  
//  f.sigma = 1.;

//	/**
//	Set a constant contact angle.*/
//	const scalar c[] = 10.*pi/180.; // fixed contact angle...
//	contact_angle = c;

	run();
}

/**
We set the boundary conditions, so to obtain a flow around a fixed cylinder. */
u.n[left]  = dirichlet(1.0);
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
  We define the solid cylinder (EMBED) and fluid cylinder (PLASTRON) 
  interface. */
  solid (cs, fs, (sq(x - xc) + sq(y - yc) - sq(R0)));

//  fraction (f, - (sq(x - xc) + sq(y - yc) - sq(R1)));
}

event adapt (i++) {
  adapt_wavelet ({cs,u,f}, (double[]){1e-2,1e-2,1e-2,1e-1}, MAXLEVEL, MINLEVEL);
}

//event movie(i+=10,t<=T){
//  view(fov=5, tx = 0, ty = 0);
//  draw_vof("cs", "fs",filled = -1);
//  //draw_vof ("f", filled = 1, fc = {1,0,0});
//	draw_vof ("f", lc = {1, 0, 0}, lw = 2);
//	squares ("u.x", linear = true);
//	// Draw grid only on upper part of flow
//	cells (lc = {0.7, 0.7, 0.7});
//
//
//  save("movie.mp4");
//}

event compute_forces (i += 1, t<=endTime)
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


event movies (i += 10)
{
  scalar omega[], m[];
  vorticity (u, omega);
  foreach()
    m[] = cs[] - 0.5;
  output_ppm (omega, file = "vort.mp4", box = {{-0.5,-0.5},{7.5,0.5}},
	      min = -10, max = 10, linear = true, mask = m);
}

/**
~~~pythonplot Force coefficient
import numpy as np
import matplotlib.pyplot as plt

# columns:
# i t dt Cx Cy
i,t,dt,Cx,Cy = np.loadtxt('force_coeff.dat', comments='#', unpack=True)

tTarget = 150.

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
Cxref = 1.38
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

plt.ylim(1.1, 1.5)
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

tTarget = 150.

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