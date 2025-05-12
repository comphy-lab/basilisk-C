/**
# Rising bubble dissolution

We simulate the dissolution of a single rising bubble.

![Animation of the volume fraction field.](rising_bubble/f.mp4)

![Animation of the concentration field.](rising_bubble/c.mp4)

The transfer rate compared with the analytical solution provided in [Farsoiya et al.,
2021](#farsoiya2021).

~~~pythonplot Sherwood number vs time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
	
def sherwood(filename,D,alpha):
    r0=1;
    pi=3.1416;
    rho=1;
    t1,tr,trw,vel1,vol,volo,area,dt = np.loadtxt(filename,delimiter=' ',unpack=True)
    ti=[];
    veli=[];
    voli=[];
    areai=[];
    ci=[];
    sh=[];
 
    levich=[];
    
    for i in range(0,len(t1)-1):
        ti.extend([0.5*(t1[i+1]+t1[i])]);
        veli.extend([0.5*(vel1[i+1]+vel1[i])]);
        voli.extend([0.5*(vol[i+1]+vol[i])]);
        areai.extend([0.5*(area[i+1]+area[i])]);
        ci.extend([0.5*(tr[i+1]/vol[i+1]+tr[i]/vol[i])]);
        sh.extend([abs((trw[i+1]-trw[i])/0.01/areai[i]/ci[i]/alpha)]);
        levich.extend([2.0 * np.sqrt(D*np.mean(veli[i])/pi/np.power(6.*voli[i]/pi,1./3))]);
    
    
    return ti,sh,levich


plt.rcParams['figure.figsize'] = 6, 4

ti4,sh,levich4 = sherwood('conc',0.08944,0.8)
plt.plot(ti4[0:-1:2], sh[0:-1:2],'r',label=r'Basilisk Pe = 102.4')
plt.plot(ti4, levich4,'k--',label=r'Levich')
plt.legend()
plt.ylim(0,2)
plt.ylabel(r'$k_L$',rotation=0,labelpad=15)
plt.xlabel(r'$t$')
plt.tight_layout()

plt.savefig('fig2e.png',bbox_inches='tight');
~~~

## References

~~~bib

@phdthesis {magdelaine2019hydrodynamics,
  title = {Hydrodynamique des films liquides h{\'e}t{\'e}rog{\`e}nes},
  author = {Magdelaine, Quentin},
  year = {2019},
  school = {Sorbonne universit{\'e}}
}

@article{farsoiya2021,
  title = {Bubble mediated single component gas transfer in homogeneous isotropic turbulence},
  author = {P. K. Farsoiya and Q. Magdelaine and A. Antkowiak and S. Popinet and L. Deike},
  journal = {Journal of Fluid Mechanics},
  year = {2021},
  note = {submitted}
}
~~~


We define the geometrical, temporal and resolution parameters: */

#define R0 1.
#define L 20 // size of the box

#define MIN_LEVEL 5
#define LEVEL 10
#define MAX_LEVEL 10
#define dR_refine (2.*L0/(1 << LEVEL))

#define F_ERR 1e-10

#define T_END 2.
#define DT_MAX 0.1
#define DELTA_T 0.1 // for measurements and videos

#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "phase-change.h"
#include "reduced.h"

#define D_V 0.08944
#define cinf 0
/**
The VOF tracer describing the interface, $f$, is already defined by
[two-phase.h](/src/two-phase.h). We allocate a second scalar field to describe
 concentration field and place it in the list of stracers of
[phase-change.h](../phase-change.h). */

scalar c[], *stracers = {c};

c[top]   = dirichlet(cinf);

/**
In the main function of the program, we set the domain geometry to
be twenty times larger than the bubble, the surface tension, the gravity, the
densities and the viscoties of both phases. */

int main() {
  size (L);
  origin (0., 0.);
  N = 1 << LEVEL;
  init_grid (N);

	f.sigma = 10.;
  rho1 = 1.;
  rho2 = 0.001;
  mu1 = D_V;
  mu2 = mu1/20.;

  G.x = - 7.8125;

	 
  c.inverse = false;	// false: bubble, true: droplet
  c.alpha = 0.8;	// solubility
  c.D = D_V;            // Diffusivity in outside fluid
  c.mw = 1.;            //Molecular weight
  c.tr_eq = rho2/c.mw;  //surface concentration ( constrained once we fix the
                        //density and molecular weight of the gas)

  run();
}

/**
The initial position of the interface is defined with this function: */

#define circle(x, y, R) (sq(R) - sq(x) - sq(y))

/**
Before the first step, we initialize the concentration field (after having
refined the grid around the future interface): $c_s$ in the bubble and 
$c_\infty$ in the liquid. */

event init (i = 0) {
  #if TREE
    refine (level < MAX_LEVEL && circle(x - L0*0.2 , y, (R0 - dR_refine)) < 0.
            && circle(x - L0*0.2, y, (R0 + dR_refine)) > 0.);
  #endif
  fraction (f, -circle(x- L0*0.2, y, R0));
  foreach() {
    u.x[] = 0.;
    c[] = f[]*cinf + (1. - f[])*c.tr_eq;
  }
  foreach_face()
    uf.x[] = 0.;
  boundary({ u, uf,c});
}

#if TREE
event adapt (i++) {
  adapt_wavelet ({f, c,u}, (double[]){1e-3, 1e-3, 1e-2, 1e-2,1e-2,1e-2},
		             minlevel = MIN_LEVEL, maxlevel = MAX_LEVEL);
}
#endif

/**
## Post-processings and videos

We now juste have to write  post-processing events to save tthe effective
radius of the bubble. */

event extract (t = 0; t += 0.01; t <= T_END)
{	


  double effective_radius, sum = 0;
	
  foreach()
    sum += (1. - f[])*dv();

  effective_radius = pow(3./2.*sum, 1./3.);
  fprintf (stderr, "%.17g %.17g\n", t, effective_radius);
  fflush(stderr);
  
  char name[80];
  sprintf (name, "conc");
  static FILE * fp = fopen (name, "w"); 
  
  double yb = 0., vb = 0., vbx = 0., area = 0., ci = 0., co = 0.;
  foreach (reduction(+:yb) reduction(+:vb) reduction(+:vbx)
	   reduction(+:ci) reduction(+:co)
	   reduction(+:area)) {
    double dvb = (1. - f[])*dv();
    vb += dvb;          // volume of the bubble
    yb += x*dvb;	// location of the bubble
    vbx += u.x[]*dvb;	// bubble velocity

    ci += c[]*(1. - f[])/(f[]*c.alpha + (1. - f[]))*dv();
    co += c[]*f[]*c.alpha/(f[]*c.alpha + (1. - f[]))*dv();
    
    if (f[] > 1e-6 && f[] < 1. - 1e-6) {
      coord n = interface_normal (point, f), p;      
      double alpha = plane_alpha (f[], n);
      // area of the bubble interface
      area += y*pow(Delta, dimension - 1)*plane_area_center (n, alpha, &p);      
    }
  }
	
  //~ if (i == 0)
    //~ fprintf (fp, "t ci co vbx vb vbo area dt\n");
  fprintf (fp,"%e %.12e %.12e %.12e %.12e %.12e %.12e %e\n",
	   t,
	   ci*2.*pi, co*2.*pi,
	   vbx/vb, 2.*pi*vb, 2.*pi*statsf(f).sum, 2.*pi*area,
	   dt);
  
  fflush (fp);
  
  /**
  We create a video with the concentration and volume fraction field */


  double vcs = rho2/c.mw;

 
  output_ppm (f, file = "f.mp4", box = {{0,0},{15,3}},
	      linear = true, min = 0, max = 1); 
  output_ppm (c, file = "c.mp4", box = {{0,0},{15,3}},
	      linear = true, min = 0, max = vcs);
}