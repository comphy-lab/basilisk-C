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
\lambda \operatorname{erfc}({\lambda}) e^{\lambda^2} = -\frac{\alpha c_g - c_0}{\sqrt{\pi}}\frac{M_w}{\rho_g}
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
<td>![](stefan/interface_movement.png){ width="400px" }</td>
<td>![](stefan/conc_vs_time.png){ width="400px" }</td>
<td>![](stefan/mass_budget.png){ width="400px" }</td>
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
#from scipy.optimize import fsolve #will not work in sandbox
import math
plt.rcParams.update({'font.size': 15})	
#plt.rcParams['text.usetex'] = True
mw = 1; rho = 0.1; alpha = 1; cg = 0.1; c0 =0.; D = 0.001;

def f(x):
	return x*math.erfc(x)*math.exp(x**2) + mw/rho*(alpha*cg - c0)/np.pi**0.5;

#lam = fsolve(f, -2);
#print(lam);

y0 = 5; domain_width = 10;
#yt = lambda t: y0 + 2*lam*np.sqrt(D*t)
#ct = lambda t: c0 + (alpha*cg - c0)*math.erfc((yt(t) +0.2 - y0)/np.sqrt(4*D*t))/math.erfc(lam)

ytg = lambda t: y0 - 2*np.sqrt(D*t/np.pi)
ctg = lambda t: alpha*cg*math.erfc(0.2/np.sqrt(4*D*t))

#intct = lambda t: -alpha*cg*4*lam*np.sqrt(D*t)*domain_width
intctg = lambda t: alpha*cg*2*np.sqrt(D*t/np.pi)*domain_width
tlist = np.linspace(0.1,100,100)

#ylist =[]
#clist =[]
#intclist =[]

intclistg =[]
ylistg =[]
clistg =[]

for i in range(0,len(tlist)):
	#ylist.extend([yt(tlist[i])])
	#clist.extend([ct(tlist[i])])
	ylistg.extend([ytg(tlist[i])])
	clistg.extend([ctg(tlist[i])])
	#intclist.extend([intct(tlist[i])])
	intclistg.extend([intctg(tlist[i])])
	




plt.figure()


tlist,ylist,clist,intclist = np.loadtxt('../farsoiya-analytic',delimiter=' ',unpack=True)

plt.plot(tlist[1:-1:3],ylist[1:-1:3], 'kx', label='Analytical Farsoiya et. al.')
plt.plot(tlist[1:-1:3],ylistg[1:-1:3], 'bx', label='Analytical Gennari et. al.')

ts, rb, cb,net_mass,mass_in_gas_phase_quentin,mass_in_gas_phase,mass_in_liq_phase = np.loadtxt('../farsoiya-conc-stefan',delimiter=' ',unpack=True)
tsg, rbg, cbg,net_massg,mass_in_gas_phase_quenting,mass_in_gas_phaseg,mass_in_liq_phaseg = np.loadtxt('conc-stefan',delimiter=' ',unpack=True)

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

#define F_ERR 1e-10

#include "grid/quadtree.h"
#define PHASE_CHANGE 1
#include "../navier-stokes/centered-gennari.h"
#include "../two-phase-gennari.h"
#include "../diffusion-gennari.h"
#include "../phase_change_pure_species-gennari.h"
//~ #include "view.h"
//~ #include "save_data.h"

/** 
Two-phase system properties
*/

#define RHOR 10. //density ratio
#define MUR 20. //viscosity ratio
#define D_V 1.e-3

#define Sc 1000. //Schmidt number

#define WIDTH 10
#define TEND 100.

/**
The concentration field is stored in the scalar *a_c*.
*/

scalar a_c[];
scalar * species = {a_c};
scalar m_a[];
scalar * m_list = {m_a};

int main()
{
  size (WIDTH);
  init_grid (256);

  rho1 = 1.;
  rho2 = rho1/RHOR;
  mu1 = 1./10;
  mu2 = mu1/MUR;

  f.tracers = species;

  a_c.inverse = false;
  //a_c.D = mu1/(rho1*Sc); //Diffusion coefficient
  a_c.D = D_V; //Diffusion coefficient
  a_c.mol_mass = 1.; //molar mass
  a_c.cd = rho2/a_c.mol_mass; //uniform concentration inside the bubble
  a_c.He = 1.; //Henry's coeff
  a_c.diff = true; //we turn on diffusion of the tracer
  
  tt = 0.; //we switch on volume change at t=0

  TOLERANCE = 1e-4;
  DT = 1.;
  run();
}

/**
We need outflow boundary conditions to let the liquid enter the domain as the gaseous region dissolves. 
*/

//Bcs
u.n[top] = neumann(0.);  // outflow is needed here, because poisson solver have a source term
p[top] = dirichlet(0.);
pf[top] = dirichlet(0.);
a_c[top] = neumann(0.); //normal derivate of conc is zero

/**
We initialize the two-phase domain.
*/

event init (t = 0) {
  fraction (f, y - 4.7);
}

//~ event stability (i++) {
  //~ DT = t < 1. ? 1e-4 : 1e-2; //This helps the stability
//~ }

/**
We monitor the volume of the gaseous region.
*/
//~ #if TREE
//~ event adapt (i++) {
  //~ adapt_wavelet ({f, a_c,u}, (double[]){1e-3, 1e-3, 1e-1, 1e-1,1e-1,1e-1},
		             //~ minlevel = 5, maxlevel = 7);
//~ }
//~ #endif
event log_simulation (i++) {
  double vol = 0.;
  foreach(reduction(+:vol))
    vol += (1. - f[])*dv();

  fprintf (stderr, "%g %g\n", t, vol);
  fflush (stderr);
}

event outputs (t = 0.; t += 0.1; t <= TEND) { 

 double effective_height;
	double net_mass = 0., mass_in_liq_phase = 0., mass_in_gas_phase = 0., mass_in_gas_phase_quentin = 0.;
	scalar fc[];
	foreach(){
	  fc[] = 1. - f[];
	  //~ net_mass += a_c[]*dv();
	  mass_in_liq_phase += a_c[]*f[]*dv();
	  //~ mass_in_gas_phase +=a_c[]*(1. - f[])*dv(); // Gennari told me this is not correct way to calculate gas conc in his code 
	}
	
	mass_in_gas_phase_quentin = statsf(fc).sum*a_c.cd; // Gennari recommended this. In our code this doesn't make any difference
	net_mass = mass_in_liq_phase + mass_in_gas_phase_quentin;
  effective_height = statsf(fc).sum/10.;

  
  char name[80];
  sprintf (name, "conc-stefan");
  static FILE * fp = fopen (name, "w");

    fprintf (fp, "%.17g %.17g %.17g %.17g %.17g %.17g %.17g\n", t,  effective_height, interpolate(a_c, L0*0.5, effective_height + 0.2 ),
	net_mass, mass_in_gas_phase, mass_in_gas_phase_quentin, mass_in_liq_phase);
 
  fflush (fp);
	
  //~ output_ppm (a_c, file = "c.mp4", box = {{0,0},{10,10}},
	      //~ linear = true, min = 0, max = a_c.cd);  
}
/**
We make a movie of the phase-average concentration around the interface.
*/

//~ event movie (t++) {
  //~ scalar c_c[]; //phase-average concentration
  //~ foreach()
    //~ c_c[] = f[] > F_ERR ? a_c[]/f[] : a_c.cd;
  //~ view (quat = {0.002, 0.005, 0.000, 1.000},
	//~ fov = 35, near = 0.01, far = 1000,
	//~ tx = -0.5, ty = -0.5, tz = -2.,
	//~ width = 600, height = 600);
  //~ clear();
  //~ box();
  //~ draw_vof ("f");
  //~ squares("c_c", min=0., max=0.8, linear=true);
  //~ save ("movie.mp4");
//~ }

event stop_simulation (t = TEND) {
  return 1;
}

/**
## References
~~~bib
@article{gennari2022,
title = {A phase-change model for diffusion-driven mass transfer problems in incompressible two-phase flows},
journal = {Chemical Engineering Science},
volume = {259},
pages = {117791},
year = {2022},
issn = {0009-2509},
doi = {https://doi.org/10.1016/j.ces.2022.117791},
url = {https://www.sciencedirect.com/science/article/pii/S000925092200375X},
author = {Gabriele Gennari and Richard Jefferson-Loveday and Stephen J. Pickering}
}
~~~
*/