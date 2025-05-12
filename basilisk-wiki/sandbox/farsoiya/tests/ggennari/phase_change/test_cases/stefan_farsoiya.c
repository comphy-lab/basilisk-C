/**
# Stefan Problem

We simulate a static planar interface which moving due to diffusion of gas. 


The time evolution of the interface is given as [Alexiades and Solomon, 1992](#alexiades1992mathematical, #farsoiya2021),

$$
Y(t) = Y_0 + 2\lambda \sqrt{D t}
$$
and concentration in liquid phase,

$$
c(t) = c_0 + \frac{\alpha c_g - c_0}{\operatorname{erfc}(\lambda)}\operatorname{erfc}\left(\frac{y - y_0}{2\sqrt{Dt}}\right)
$$

and integration of concentration in liquid phase,

$$
 \int_{Y(t)}^\infty c(t) dy= \int_{Y(t)}^\infty \left[c_0 + \frac{\alpha c_g - c_0}{\operatorname{erfc}(\lambda)}\operatorname{erfc}\left(\frac{y -y_0}{2\sqrt{Dt}}\right)\right] dy = 4\lambda \sqrt{D t} \quad \text{for} \quad c_0 = 0
$$

where lambda is solution of the equation,

$$
\lambda \operatorname{erfc}({\lambda}) e^{\lambda^2} = -\frac{\alpha c_g - c_0}{c_g(1-\alpha)\sqrt{\pi}}
$$

As discussed in Section 4.3 of [Gennari et al.,2022](#gennari2022).The analytical solution of the interface displacement is:

$$
l(t) = \frac{2}{He} \sqrt{\frac{t\, D}{\pi}}
$$

and the concentration field in the liquid domain:

$$
c(y,t) = c_\Sigma \left(1 - erf \left( \frac{y - y_\Sigma(t)}{2\sqrt{D\, t}} \right) \right)
$$

and the integration of concentration field in liquid domain,

$$
c(y,t) = \int_{y_{\Sigma(t)}}^{\infty} c_\Sigma \left(1 - erf \left( \frac{y - y_\Sigma(t)}{2\sqrt{D\, t}}  \right) \right)dy = = \frac{2}{He}c_\Sigma \sqrt{\frac{tD}{\pi}}
$$



<center>
<table>
<tr>
<td>![](stefan_farsoiya/interface_movement.png){ width="400px" }</td>
<td>![](stefan_farsoiya/conc_vs_time.png){ width="400px" }</td>
<td>![](stefan_farsoiya/mass_budget.png){ width="400px" }</td>
</tr>
<tr>
<td><center>Interface height</center></td> 
<td><center>Concentration at a point (Y(t) + 0.2) in liquid</center></td> 
<td><center>Mass budget</center></td> 
</tr>
</table>
</center>

~~~pythonplot

import numpy as np
import matplotlib.pyplot as plt
#import scipy as sp
#import mpmath as mp
from scipy.optimize import fsolve
import math
plt.rcParams.update({'font.size': 15})	
#plt.rcParams['text.usetex'] = True
mw = 1; rho = 0.1; alpha = 1; cg = 0.1; c0 =0.; D = 0.001;

def f(x):
	return x*math.erfc(x)*math.exp(x**2) + mw/rho*(alpha*cg - c0)/np.pi**0.5;

lam = fsolve(f, -2);
print(lam);

y0 = 5; domain_width = 10;
yt = lambda t: y0 + 2*lam*np.sqrt(D*t)
ct = lambda t: c0 + (alpha*cg - c0)*math.erfc((yt(t) +0.2 - y0)/np.sqrt(4*D*t))/math.erfc(lam)

ytg = lambda t: y0 - 2*np.sqrt(D*t/np.pi)
ctg = lambda t: alpha*cg*math.erfc(0.2/np.sqrt(4*D*t))

intct = lambda t: -alpha*cg*4*lam*np.sqrt(D*t)*domain_width
intctg = lambda t: alpha*cg*2*np.sqrt(D*t/np.pi)*domain_width
tlist = np.linspace(0.1,100,100)

ylist =[]
clist =[]
intclist =[]
intclistg =[]

ylistg =[]
clistg =[]

for i in range(0,len(tlist)):
	ylist.extend([yt(tlist[i])])
	clist.extend([ct(tlist[i])])
	ylistg.extend([ytg(tlist[i])])
	clistg.extend([ctg(tlist[i])])
	intclist.extend([intct(tlist[i])])
	intclistg.extend([intctg(tlist[i])])
	#print("%g %g %g %g" %(tlist[i], ylist[i], clist[i], intclist[i]))
	




plt.figure()

plt.plot(tlist[1:-1:3],ylist[1:-1:3], 'kx', label='Analytical Farsoiya et. al.')
plt.plot(tlist[1:-1:3],ylistg[1:-1:3], 'bx', label='Analytical Gennari et. al.')

ts, rb, cb,net_mass,mass_in_gas_phase_quentin,mass_in_gas_phase,mass_in_liq_phase = np.loadtxt('conc-stefan',delimiter=' ',unpack=True)

#load Gennari et. al. result
tsg, rbg, cbg,net_massg,mass_in_gas_phase_quenting,mass_in_gas_phaseg,mass_in_liq_phaseg = np.loadtxt('../stefan/conc-stefan',delimiter=' ',unpack=True)

plt.plot(ts,rb,'k--',label='Farsoiya et. al.');
plt.plot(tsg,rbg+0.3,'b--',label='Gennari et. al');
# plt.ylim(0.97,1.007)

plt.legend();
plt.xlabel(r'$t$')
plt.ylabel(r'$Y(t)$')
plt.tight_layout()

plt.savefig('interface_movement.png')

plt.figure()

#plt.plot(t, cep,'k',label='Analytical')
plt.plot(tlist[1:-1:2],clist[1:-1:2], 'kx', label='Analytical Farsoiya et. al.')
plt.plot(tlist[1:-1:3],clistg[1:-1:3], 'bx', label='Analytical Gennari et. al.')

plt.plot(ts,cb,'k--',label='Farsoiya et. al.');
plt.plot(tsg,cbg,'b--',label='Gennari et. al.');
# plt.ylim(0.97,1.007)

plt.legend();
plt.xlabel(r'$t$')
plt.ylabel(r'$c$')
plt.tight_layout()

plt.savefig('conc_vs_time.png')
plt.figure()

plt.figure()


#plt.plot(ts,net_mass,'r--',label='Net Mass');
#plt.plot(ts[1:-1:10],mass_in_gas_phase[1:-1:10],'k--x',label='Mass in gas phase');
#plt.plot(ts,mass_in_gas_phase_quentin,'b--',label='Mass in gas phase Quentin');

plt.plot(tlist[1:-1:2],intclist[1:-1:2],'kx',label=r'$\int_V c\; dv = 4\alpha c_g \lambda \sqrt{D t}$');
plt.plot(tlist[1:-1:3],intclistg[1:-1:3],'bx',label=r'$\int_V c\; dv = 2 \alpha c_g \sqrt{D t/\pi}$');

#plt.plot(tsg,net_massg+0.3,'r-',label='Net Mass Gennari');
#plt.plot(ts[1:-1:10],mass_in_gas_phase[1:-1:10],'k--x',label='Mass in gas phase');
#plt.plot(tsg,mass_in_gas_phaseg+0.3,'b-',label='Mass in gas phase Quentin Gennari');

plt.plot(ts,mass_in_liq_phase,'k--',label='Farsoiya et. al.');

plt.plot(tsg,mass_in_liq_phaseg,'b-',label='Gennari et. al.');


# plt.ylim(0.97,1.007)

plt.legend();
plt.xlabel(r'$t$')
plt.ylabel(r'total mass in liquid phase')
plt.tight_layout()

plt.savefig('mass_budget.png')
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
#define LEVEL 8
#define MAX_LEVEL 8
#define dR_refine (2.*L0/(1 << LEVEL))

#define F_ERR 1e-10

#define T_END 100
#define DT_MAX 10.
#define DELTA_T 1. // for measurements and videos

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "../../../../phase_change/phase-change.h"
#include "reduced.h"

#define D_V 1.e-3
#define cinf 0.

scalar c[], *stracers = {c};


//~ c[top]   = dirichlet(cinf);
c[top]   = neumann(0.);
c[bottom]   = neumann(0.);
c[left]   = neumann(0.);
c[right]   = neumann(0.);
u.n[top] = dirichlet(0.);
u.n[bottom] = dirichlet(0.);
u.n[left] = dirichlet(0.);
u.n[right] = dirichlet(0.);

int main() {
  size (L);
  origin (0., 0.);
  N = 1 << LEVEL;
  init_grid (N);

  f.sigma = 0;
  rho1 = 1.;
  rho2 = 0.1;
  mu1 = 1./10.;
  mu2 = mu1/20.;

  G.x = 0.;
  Z.x = 0.;

	 
  c.inverse = false;	
  c.alpha = 1.;	
 
  c.D = D_V;
  c.mw = 1.;
   c.tr_eq = rho2/c.mw;
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



//~ #if TREE
//~ event adapt (i++) {
  //~ adapt_wavelet ({f, c,u}, (double[]){1e-3, 1e-3, 1e-1, 1e-1,1e-1,1e-1},
		             //~ minlevel = MIN_LEVEL, maxlevel = MAX_LEVEL);
//~ }
//~ #endif

/**
## Post-processings and videos

![Animation of the concentration field.](stefan_problem/c.mp4)

We now just have to write several post-processing events */

event outputs (t = 0.; t += DELTA_T; t <= T_END) { 

 double effective_height;
	double net_mass = 0., mass_in_liq_phase = 0., mass_in_gas_phase = 0., mass_in_gas_phase_quentin = 0.;
	scalar fc[];
	foreach(){
	  fc[] = 1. - f[];
	  net_mass += c[]*dv();
	  mass_in_liq_phase += c[]*f[]*dv();
	  mass_in_gas_phase +=c[]*(1. - f[])*dv();
	}
	
	mass_in_gas_phase_quentin = statsf(fc).sum*c.tr_eq; // Quentin suggestion to multiply const concentration of gas phase with volume
	
  effective_height = statsf(fc).sum/10.;

  
  char name[80];
  sprintf (name, "conc-stefan");
  static FILE * fp = fopen (name, "w");

    fprintf (fp, "%.17g %.17g %.17g %.17g %.17g %.17g %.17g\n", t,  effective_height, interpolate(c, L0*0.5, effective_height + 0.2),
	net_mass, mass_in_gas_phase, mass_in_gas_phase_quentin, mass_in_liq_phase);
 
  fflush (fp);
	
  output_ppm (c, file = "c.mp4", box = {{0,0},{10,10}},
	      linear = true, min = 0, max = c.tr_eq);  
}
