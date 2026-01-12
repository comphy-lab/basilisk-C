/**
# Rayleigh-Bénard convection in 2D

For this example we solve the [Boussinesq
equations](https://basilisk.fr/sandbox/acastillo/convection/convection_boussinesq.h).

\newcommand{\Ra}{{\color{green}\mathsf{Ra}}}
\newcommand{\Pr}{{\color{blue}\mathsf{Pr}}}
\newcommand{\Nu}{{\color{red}\mathsf{Nu}}}

In this example, we fix $\Ra=5 \times 10^5$ and $\Pr=0.71$ which
leads to time-dependent convection patterns.

*/

#define LEVEL 6
#define RAYLEIGH 5e5
#define PRANDTL 0.71
#define MAXTIME 50.
#include "grid/multigrid.h"
#include "acastillo/convection/convection_boussinesq.h"

scalar omega[];
scalar psi[];

double tend= MAXTIME;

int main() {

  // Geometric Parameters
  L0 = 1.;
  X0 = Y0 = -0.5;

  // Physical Parameters
  Ra = RAYLEIGH; 
  Pr = PRANDTL;

  // Numerical Parameters
  DT = 0.05;
  TOLERANCE = 1e-3;
  NITERMIN = 1;

  FILE * fp = fopen("nusselts.asc", "w");
  fprintf (fp, "[1]t ");
  fprintf (fp, "[2]nutop ");
  fprintf (fp, "[3]nubot ");
  fprintf (fp, "[4]nuvol ");
  fprintf (fp, "[5]nuvisc ");
  fprintf (fp, "[6]nutemp ");
  fprintf (fp, "[7]numix \n");
  fclose (fp);

  fp = fopen("energies.asc", "w");
  fprintf (fp, "[1]t "); 
  fprintf (fp, "[2]E_kin ");
  fprintf (fp, "[3]E_pot \n");
  fclose (fp);

  N = 1 << LEVEL;
  run();
}

/**
## Boundary conditions

The top and bottom walls are iso-thermal $\theta = \mp 1/2$ at $y = \pm 1/2$,
while all other boundaries are adiabatic. In addition to no-slip walls on all
boundaries.

*/
T[top] = dirichlet(-0.5);
T[bottom] = dirichlet(0.5);

u.n[top] = dirichlet(0.);
u.n[bottom] = dirichlet(0.);
u.n[right] = dirichlet(0.);
u.n[left] = dirichlet(0.);

u.t[top] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);
u.t[right] = dirichlet(0.);
u.t[left] = dirichlet(0.);


/**
## Initial conditions
Initial conditions correspond to a linear temperature profile with small noise
and no motion.
*/

event init (t=0) {
  foreach(){
    T[] = -y + 0.05*noise();
    foreach_dimension()
      u.x[] = 0.;
  }
}

/**
## Outputs

We store the residuals from the Poisson solvers

~~~pythonplot Residuals
import numpy as np
import matplotlib.pyplot as plt

filename1 = 'log'
data = np.genfromtxt(filename1, usecols=[0,2,6,12,18], invalid_raise=False)
data = data[~np.isnan(data).any(axis=1)]
t = data[:,1]
residual_p = data[:,2]
residual_u = data[:,3]
residual_T = data[:,4]

fig, axes = plt.subplots(1,3, figsize=(4*3, 3))
ax = axes.ravel()
ax[0].semilogy(t, residual_p)
ax[1].semilogy(t, residual_u)
ax[2].semilogy(t, residual_T)
xlabels = ['Pressure', 'Velocity', 'Temperature']
for i in range(3):
  ax[i].set_xlabel('Time')
  ax[i].set_ylabel('Residual '+xlabels[i])

plt.tight_layout()
plt.savefig('plot_residuals.svg', dpi=150)
~~~

*/

void mg_print (mgstats mg)
{
  if (mg.i > 0 && mg.resa > 0.)
    fprintf (stderr, " \t - \t %d %g %g %g %d ", mg.i, mg.resb, mg.resa,
	    mg.resb > 0 ? exp (log (mg.resb/mg.resa)/mg.i) : 0.,
	    mg.nrelax);
}

event log (i++) {
  fprintf (stderr, " %d %d %g ", N, i, t);
  mg_print (mgp);
  mg_print (mgu);
  mg_print (mgT);
  fprintf (stderr, "\n");
  fflush (stderr);
}

/**
We are going to compute (and save) the heat-fluxes, the kinetic and potential
energies, and the temperature variance. 
*/

#include "acastillo/convection/global_nusselt.h"
void compute_and_store_fluxes (vector u, scalar T){
  
  double nu1 = nusselt_top(T);
  double nu2 = nusselt_bot(T);
  double nu3 = nusselt_vol(T,u);
  double nu4 = nusselt_viscous(u);
  double nu5 = nusselt_tmp(T);

  if (pid()==0){
    FILE * fp = fopen("nusselts.asc", "a");
    fprintf (fp, "%.9g ", t); 
    fprintf (fp, "%.9g ", nu1);
    fprintf (fp, "%.9g ", nu2);
    fprintf (fp, "%.9g ", nu3);
    fprintf (fp, "%.9g ", nu4);
    fprintf (fp, "%.9g \n", nu5);
    fclose (fp);
  }
}

event time_series (t += 0.05; t <= tend){

  compute_and_store_fluxes (u, T);

  double ekin = 0., epot = 0., tvar = 0.;
  foreach(reduction(+:ekin) reduction(+:epot) reduction(+:tvar)){
    epot += -dV * Pr * T[] * y;
    foreach_dimension()
      ekin += dV * 0.5 * sq(u.x[]);
    tvar += dV * sq(T[]);
  }

  if (pid()==0){
    FILE * fp = fopen("energies.asc", "a");
    fprintf (fp, "%f %.9g %.9g %.9g \n", t, ekin, epot, tvar);
    fclose (fp);
  }
}

/** 
These quantities are extremely handy, since there are some **exact** relations
that can be used to check the results (see Chapter 2 in [Castillo
(2017)](#castillo:2017)). For instance, the heat fluxes should be related to the
rates of change of the potential energy as follows:

$$
\displaystyle
\frac{\mathrm{dE}_{\text{pot}}}{\mathrm{d} t} = \frac{\Pr}{\sqrt{\Ra}} 
\left( \Nu_{\text{vol}} + \frac{\Nu_{\text{top}} + \Nu_{\text{bot}}}{2} \right)
$$

or to the kinetic energy through the viscous dissipation:

$$
\displaystyle
\frac{\mathrm{dE}_{\text{kin}}}{\mathrm{d} t} = 
\frac{\Pr}{\sqrt{\Ra}} \left( \Nu_{\text{vol}} - (\Nu_{\text{visc}} + 1) \right)
$$

or to the temperature variance through the scalar dissipation:

$$
\displaystyle
\frac{\mathrm{d E}_{\theta}}{\mathrm{d} t} = 
\frac{1}{\sqrt{\Ra}} \left( \Nu_{\theta} -  \frac{\Nu_{\text{top}} + \Nu_{\text{bot}}}{2} \right)
$$

In the statistically stable state, the rates of change are zero and the heat
fluxes are related through the classical relations from [Shraiman and 
Siggia (1990)](#shraiman:1990).

Since the viscous and scalar dissipation rates depend on higher order
derivatives, it is going to converge more slowly than quantities that depend on
the primitive variables. This is very common on convection problems.

~~~pythonplot Relation to the potential energy
import numpy as np
import matplotlib.pyplot as plt

Ra = 5e5
Pr = 0.71

filename1 = 'nusselts.asc'
filename2 = 'energies.asc'

data1 = np.genfromtxt(filename1, usecols=[0,1,2,3,4,5], skip_header=2)
data2 = np.genfromtxt(filename2, usecols=[0,1,2], skip_header=2)
t = data1[:,0]
nutop = data1[:,1]
nubot = data1[:,2]
nuvol = data1[:,3]
epot = data2[:,2]

dt = t[1] - t[0]
depot_dt = np.diff(epot)/dt

fig, axes = plt.subplots(1,2, figsize=(4*2, 3))
ax = axes.ravel()
ax[0].plot(t, epot)
ax[0].set_ylabel(r'$E_{pot}$')
ax[0].set_xlabel(r'$t$')
ax[1].plot(t[:-1], depot_dt, label='l.h.s.')
ax[1].plot(t, (Pr/np.sqrt(Ra))*(-nuvol + (nubot+nutop)/2.),label='r.h.s.', ls='--')
ax[1].set_ylabel(r'$\frac{\mathrm{d}E_{pot}}{\mathrm{d} t}$')
ax[1].set_xlabel(r'$t$')
ax[1].legend()
plt.tight_layout()
plt.savefig('plot_nusselt1.svg', dpi=150)
~~~

~~~pythonplot Relation to the kinetic energy
import numpy as np
import matplotlib.pyplot as plt

Ra = 5e5
Pr = 0.71

filename1 = 'nusselts.asc'
filename2 = 'energies.asc'

data1 = np.genfromtxt(filename1, usecols=[0,1,2,3,4,5], skip_header=2)
data2 = np.genfromtxt(filename2, usecols=[0,1,2], skip_header=2)
t = data1[:,0]
nuvol = data1[:,3]
nuvis = data1[:,4]
ekin = data2[:,1]

dt = t[1] - t[0]
dekin_dt = np.diff(ekin)/dt

fig, axes = plt.subplots(1,2, figsize=(4*2, 3))
ax = axes.ravel()
ax[0].plot(t, ekin)
ax[0].set_ylabel(r'$E_{kin}$')
ax[0].set_xlabel(r'$t$')
ax[1].plot(t[:-1], dekin_dt, label='l.h.s.')
ax[1].plot(t, (Pr/np.sqrt(Ra))*(nuvol - (nuvis + 1)),label='r.h.s.', ls='--')
ax[1].set_ylabel(r'$\frac{\mathrm{d}E_{kin}}{\mathrm{d} t}$')
ax[1].set_xlabel(r'$t$')
ax[1].legend()
plt.tight_layout()
plt.savefig('plot_nusselt2.svg', dpi=150)
~~~

~~~pythonplot Relation to the temperature variance
import numpy as np
import matplotlib.pyplot as plt

Ra = 5e5
Pr = 0.71

filename1 = 'nusselts.asc'
filename2 = 'energies.asc'

data1 = np.genfromtxt(filename1, usecols=[0,1,2,3,4,5], skip_header=2)
data2 = np.genfromtxt(filename2, usecols=[0,1,2,3], skip_header=2)
t = data1[:,0]
nutop = data1[:,1]
nubot = data1[:,2]
nutmp = data1[:,5]
tvar = 0.5*data2[:,3]

dt = t[1] - t[0]
dtvar_dt = np.diff(tvar)/dt

fig, axes = plt.subplots(1,2, figsize=(4*2, 3))
ax = axes.ravel()
ax[0].plot(t, tvar)
ax[0].set_ylabel(r'$E_{\theta}$')
ax[0].set_xlabel(r'$t$')
ax[1].plot(t[:-1], dtvar_dt, label='l.h.s.')
ax[1].plot(t, -(1/np.sqrt(Ra))*(nutmp - (nutop + nubot)/2.),label='r.h.s.', ls='--')
ax[1].set_ylabel(r'$\frac{\mathrm{d}E_{\theta}}{\mathrm{d} t}$')
ax[1].set_xlabel(r'$t$')
ax[1].legend(loc='lower left')
plt.tight_layout()
plt.savefig('plot_nusselt3.svg', dpi=150)
~~~

Based on these results one should be able to see if the spatial resolution is 
good enough and if the solver works as expected. 

We may also plot the temperature field. 

![Temperature field](temperature.mp4)

*/

void streamfunction (vector u, scalar psi) {
  vorticity (u, omega);
  poisson (psi, omega);
}

#include "view.h"
event movie (t += 0.05){
  streamfunction (u, psi);
  squares ("T", linear = false, spread = -1);
  isoline ("psi", n = 11);
  box();
  save ("temperature.mp4");
}

/**

# References

~~~bib

@hal{castillo:2017, tel-01609741}

@article{shraiman:1990,
  title={Heat transport in high-Rayleigh-number convection},
  author={Shraiman, Boris I and Siggia, Eric D},
  journal={Physical Review A},
  volume={42},
  number={6},
  pages={3650},
  year={1990},
  publisher={APS}
}


~~~
*/