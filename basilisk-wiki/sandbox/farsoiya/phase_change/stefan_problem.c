/**
# Stefan Problem

We simulate a static planar interface which moving due to diffusion of gas. 


The time evolution of the interface is given as [Alexiades and Solomon, 1992](#alexiades1992mathematical),

$$
Y(t) = Y_0 + 2\lambda \sqrt{D t}
$$
and concentration,

$$
c(t) = \frac{\alpha c_s - c_0}{\operatorname{erfc}(\lambda)}\operatorname{erfc}\left(\frac{y - Y_0}{2\sqrt{Dt}}\right)
$$

where lambda is solution of the equation,

$$
\lambda \operatorname{erfc}({\lambda}) e^{\lambda^2} = -\frac{\alpha c_s - c_0}{\sqrt{\pi}}\frac{M_w}{\rho_g}
$$


<center>
<table>
<tr>
<td>![](stefan_problem/p001cbt.png){ width="400px" }</td>
<td>![](stefan_problem/p001clt.png){ width="400px" }</td>
</tr>
<tr>
<td><center>Interface height</center></td> 
<td><center>Concentration at a point in liquid</center></td> 
</tr>
</table>
</center>

~~~pythonplot

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 15})	
	
plt.figure()
t,rep,cep = np.loadtxt('../ref-stefan',delimiter='\t',unpack=True);
plt.plot(t, rep,'k',label='Analytical')

ts, rb, cb = np.loadtxt('conc-stefan',delimiter=' ',unpack=True)
plt.plot(ts,rb,'k--',label='Basilisk');
# plt.ylim(0.97,1.007)

plt.legend();
plt.xlabel(r'$t$')
plt.ylabel(r'$Y(t)$')
plt.tight_layout()

plt.savefig('p001cbt.png')

plt.figure()

plt.plot(t, cep,'k',label='Analytical')

plt.plot(ts,cb,'k--',label='Basilisk');
# plt.ylim(0.97,1.007)

plt.legend();
plt.xlabel(r'$t$')
plt.ylabel(r'$c_l$')
plt.tight_layout()

plt.savefig('p001clt.png')
plt.figure()

plt.savefig(' ')

~~~

## References

~~~bib

@book{alexiades1992mathematical,
  title={Mathematical modeling of melting and freezing processes},
  author={Alexiades, Solomon},
  year={1992},
  publisher={CRC Press}
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

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "phase-change.h"
#include "reduced.h"

#define D_V 1.e-3
#define cinf 0

scalar c[], *stracers = {c};


c[top]   = dirichlet(cinf);


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

	 
  c.inverse = false;	
  c.alpha = 0.8;	
  c.tr_eq = 1.;
  c.D = D_V;
  c.mw = 1.;
  
  run();
}

event init (i = 0) {

  fraction (f, y - L0*0.5);
  foreach() {
    u.x[] = 0.;
    //~ vapor[] = f[]*cinf + (1. - f[])*vcs;
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

![Animation of the concentration field.](stefan_problem/c.mp4)

We now just have to write several post-processing events */

event outputs (t = 0.; t += DELTA_T; t <= T_END) { 

 double effective_height;
	scalar fc[];
	foreach()
	  fc[] = 1. - f[];
	
  effective_height = statsf(fc).sum/10.;

  
  char name[80];
  sprintf (name, "conc-stefan");
  static FILE * fp = fopen (name, "w");

    fprintf (fp, "%.17g %.17g %.17g\n", t,  effective_height, interpolate(c, L0*0.5, L0*0.5 + 0.2));
 
  fflush (fp);
	
  output_ppm (c, file = "c.mp4", box = {{0,0},{10,10}},
	      linear = true, min = 0, max = c.tr_eq);  
}
