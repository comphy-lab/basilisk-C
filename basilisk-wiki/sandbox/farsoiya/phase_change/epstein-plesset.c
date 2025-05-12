/**
# Epstein-Plesset test

We simulate a static bubble which shrinks due to diffusion of gas. 


The time evolution of the interface is given by 
[Epstein and Plesset, 1950](#epstein1950stability) as,
$$
\frac{d R}{d t}=-\mathcal{D} M_w\frac{\alpha c_s - c_i}{\rho}\left\{\frac{1}{R}+\frac{1}{\sqrt{\pi \mathcal{D} t}}\right\}
$$
and concentration,

$$
c = c_0 + (\alpha c_s - c_0) \frac{R(t)}{r} \operatorname{erfc}\left( \frac{r - R(t)}{\sqrt{4 \mathcal{D} t}} \right)
$$
<center>
<table>
<tr>
<td>![](epstein-plesset/p001cbt.png){ width="400px" }</td>
<td>![](epstein-plesset/p001clt.png){ width="400px" }</td>
</tr>
<tr>
<td><center>Interface evolution of the bubble</center></td> 
<td><center>Concentration field outside bubble</center></td> 
</tr>
</table>
</center>

~~~pythonplot

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 15})	
	
plt.figure()
t,rep,cep = np.loadtxt('../ref-epstein-plesset',delimiter='\t',unpack=True);
plt.plot(t/4, rep,'k',label='Epstein-Plesset')

ts, rb, cb = np.loadtxt('conc-plesset-spherical',delimiter=' ',unpack=True)
plt.plot(ts/4,rb,'k--',label='Basilisk');
# plt.ylim(0.97,1.007)

plt.legend();
plt.xlabel(r'$t\; \mathscr{D}_l/d_0^2$')
plt.ylabel(r'$r(t)/r(0)$')
plt.tight_layout()

plt.savefig('p001cbt.png')

plt.figure()

plt.plot(t/4, cep,'k',label='Analytical')

plt.plot(ts/4,cb,'k--',label='Basilisk');
# plt.ylim(0.97,1.007)

plt.legend();
plt.xlabel(r'$t\; \mathscr{D}_l/d_0^2$')
plt.ylabel(r'$c_l/c_{b0}$')
plt.tight_layout()

plt.savefig('p001clt.png')
plt.figure()

plt.savefig(' ')

~~~

## References

~~~bib

@article{epstein1950stability,
  title={On the stability of gas bubbles in liquid-gas solutions},
  author={Epstein, Paul S and Plesset, Milton S},
  journal={The Journal of Chemical Physics},
  volume={18},
  number={11},
  pages={1505--1509},
  year={1950},
  publisher={American Institute of Physics}
}

@article{farsoiya2021,
  title = {Bubble mediated single component gas transfer in homogeneous isotropic turbulence},
  author = {P. K. Farsoiya and Q. Magdelaine and A. Antkowiak and S. Popinet and L. Deike},
  journal = {Journal of Fluid Mechanics},
  year = {2021},
  note = {submitted}
}
~~~
*/

#define R0 1.
#define L 10. // size of the box

#define MIN_LEVEL 5
#define LEVEL 10
#define MAX_LEVEL 10
#define dR_refine (2.*L0/(1 << LEVEL))

#define F_ERR 1e-10

#define T_END 600
#define DT_MAX 10.
#define DELTA_T 1. // for measurements and videos


#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "phase-change.h"
#include "reduced.h"


#define D_V 1.
#define cinf 0

scalar c[], *stracers = {c};

/**
Thanks to symmetry, we only simulate quarter bubble */


c[right] = dirichlet(cinf);
c[top]   = dirichlet(cinf);


/**
In the main function of the program, we set the domain geometry to
be ten times larger than the bubble. */

int main() {
  size (L);
  origin (0., 0.);
  N = 1 << LEVEL;
  init_grid (N);

  f.sigma = 0;
  rho1 = 1.;
  rho2 = 1.;
  mu1 = 1./10.;
  mu2 = mu1/20.;

  G.x = 0.;
  Z.x = 0.;

	 
  c.inverse = false;	// false: bubble, true: drop
  c.alpha = 0.8;	// solubility
  c.D = D_V;            // Diffusivity in outside fluid
  c.mw = 0.001;         //Molecular weight
  c.tr_eq = 1.;         
	
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
    refine (level < MAX_LEVEL && circle(x, y, (R0 - dR_refine)) < 0.
            && circle(x, y, (R0 + dR_refine)) > 0.);
  #endif
  fraction (f, -circle(x, y, R0));
  //~ fraction (f, y - L0*0.5);
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
  adapt_wavelet ({f, c,u}, (double[]){1e-3, 1e-3, 1e-1, 1e-1,1e-1,1e-1},
		             minlevel = MIN_LEVEL, maxlevel = MAX_LEVEL);
}
#endif

/**
## Post-processings and videos

![Animation of the volume fraction field.](epstein-plesset/f.mp4)

![Animation of the concentration field.](epstein-plesset/c.mp4)

We now juste have to write post-processing events to save tthe effective
radius of the bubble. */


event outputs (t = 0.; t += DELTA_T; t <= T_END) { 

 double effective_radius;
	scalar fc[];
	foreach()
	  fc[] = 1. - f[];
	
  effective_radius = pow(3.*statsf(fc).sum, 1./3.);
  fflush(stderr);
  char name[80];
  sprintf (name, "conc-plesset-spherical");
  static FILE * fp = fopen (name, "w");

  fprintf (fp, "%.17g %.17g %.17g\n", t,  effective_radius, interpolate(c, (1. + 0.2)*cos(M_PI/4), (1. + 0.2)*sin(M_PI/4) , 0) );
 
  fflush (fp);
  
  output_ppm (f, file = "f.mp4", box = {{0,0},{3,3}},
	      linear = true, min = 0, max = 1); 
  output_ppm (c, file = "c.mp4", box = {{0,0},{3,3}},
	      linear = true, min = 0, max = c.tr_eq);
}
