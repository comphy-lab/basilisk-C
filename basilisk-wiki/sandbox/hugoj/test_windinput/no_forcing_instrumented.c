/**
# Wave energy viscous dissipation

## Goal

This code tests the ability of the multilayer solver to dissipate energy when
explicite viscosity is used.

## Analytical solution

A linear surface gravity wave will decay in a viscous fluid. The rate of this
decay is $E(t)=E_0 e^{-4\nu k^2 t}$ (Lamb 1932)

## This code

This code uses a modified version of src/layered/{hydro.h, nh.h, remapc.h} to
gather the triggering of limiters. This is just diagnostic, no modification of
the dynamic has been done. The modified files are in the folder
'modifier_layered'. Use the recipe 'prep_instrumented' from the Makefile to do
the proper linking.

## Observations








Re=40000
- N=1024 nl=15 n'est pas convergé.
- 1 mode, N=2048 nl=15 fait apparaitre des vitesses verticale très fortes ...
- 1 mode, N=4096 pareil
Instabilité à lambda / 4 ?? Pour lambda=1, L=1, instabilité à -0.25, 0, 0.25

- 2 modes, N=2048 ok !
- Diminuer Re fait augmenter les instabilités ?? + message "warning: could not conserve barotropic flux"


Cas qui fonctionne mais trop dissipatif:
k=2pi L=1 ak=0.05 g=9.81 theta_H=0.6 Re=40000 N=1024 nl=15

## References

Horace Lamb, Hydrodynamics (6th ed., 1932), Chapter XI, Article 348, "Effect of
Viscosity on Water-Waves," pp. 623–625

et Article 349 pour une dérivation plus sérieuse


*/

#include "grid/multigrid1D.h"
#if DIAG
#include "modified_layered/hydro.h"
#include "modified_layered/nh.h"
#include "modified_layered/remap.h"
#else
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
#endif // DIAG
#include "layered/perfs.h"
#include "bderembl/libs/netcdf_bas.h"


double k_, h_, ak, RE, g_ = 9.81; 

#define T0  (2.*pi/sqrt(g_*k_))
double NT0;

double p0;

int main(int argc, char * argv[])
{
  if (argc > 1) N = atoi(argv[1]); else N=1024;             // Horizontal resolution
  if (argc > 2) nl = atoi(argv[2]); else nl=60;             // Number of layers
  if (argc > 3) NT0 = atof(argv[3]); else NT0=100.;         // Number of periods targeted 
  if (argc > 4) RE = atof(argv[4]); else RE=400.;           // Reynolds number, nu=1/Re, u diffusion (isotropic)
  if (argc > 5) L0 = atof(argv[5]); else L0=100.;           // Box size (not necessarily related to peak wave number)
  if (argc > 6) k_ = 2.*pi/atof(argv[6]); else k_=2*pi/10.; // Peak wavelength 
  if (argc > 7) h_ = atof(argv[7]); else h_=5.0;            // Depth of bassin
  if (argc > 8) ak = atof(argv[8]); else ak=0.01;           // Steepness
  if (argc > 9) CFL_H = atof(argv[9]); else CFL_H=1.0;      // Barotropic CFL (advection of eta)
  if (argc > 10) CFL = atof(argv[10]); else CFL=0.5;        // Non hydrostatic CFL (advection of u, tracers)
  if (argc > 11) theta_H = atof(argv[11]); else theta_H=0.5;  // Numerical parameter to dump fast barotropic modes
  if (argc > 12) max_slope = atof(argv[12]); else max_slope=1.; // Breaking parameter
  if (argc > 13) theta = atof(argv[13]); else theta=1.3;    // minmod2: 1 => minmod, 2 => superbee
  if (argc > 14) DT = atof(argv[14]); else DT=HUGE;         // Max timestep (default is given by CFL)
  
  fprintf (stderr, "N=%d,nl=%d,NT0=%g,RE=%g,L0=%g,k_=%g,h_=%g,ak=%g,\
  CFL_H=%g, CFL=%g, theta_H=%g, max_slope=%g, DT=%g,\
  theta=%g\n",
  N,nl,NT0,RE,L0,k_,h_,ak,CFL_H,CFL,theta_H,max_slope,DT,theta);

  origin (-L0/2.);
  periodic (right);
  G = g_;
  nu = 1./RE;

  p0=1.;

  run();
}

#if DIAG
/** minmod2 but with a counter */
long lim_zero = 0, lim_calls = 0;
double logged_minmod2 (double s0, double s1, double s2) {
  lim_calls++;
  double d = minmod2 (s0, s1, s2);
  if (d == 0.) lim_zero++;
  return d;
}
#endif // DIAG

event init (i = 0)
{
  foreach() {
    zb[] = -h_;
    double H =  ak/k_*cos(k_*x) - zb[];
    double z = zb[];
    foreach_layer() {
      h[] = H/nl;
      z += h[]/2.;
      u.x[] = ak/k_*sqrt(g_*k_)*exp(k_*z)*cos(k_*x); 
      w[] = ak/k_*sqrt(g_*k_)*exp(k_*z)*sin(k_*x);  
      z += h[]/2.;
    }
  }
  // initialise netcdf file
  create_nc({zb, eta, h, u.x, w}, "out.nc");
  #if DIAG
  // initialise log file
  fprintf (stderr, "time clip_count flux_count lim_zero lim_calls remap_flat remap_clip remap_calls w_clip w_calls\n");
  #endif // DIAG
}

#if DIAG
// Replace default minmod with my instrumented version
event defaults (i = 0) {
  // Here 'theta' plays a role
  h.gradient = logged_minmod2;
  for (scalar s in tracers)
    s.gradient = logged_minmod2;
}
#endif


#if FORCING
#warning "FORCING accel branch active"
event acceleration(i++, last){
  foreach_face(x) {
    double detadx = 0.;
    detadx = (eta[] - eta[-1])/Delta;
    ha.x[0,0,nl-1] += hf.x[]*(p0 * detadx * sin(atan(detadx))); //   
  }
}
#endif

event viscous_term (i++) {
  // Note: vertical diffusion is done in diffusion.h (u.x) and nh.h (w)
  horizontal_diffusion ({u.x, w}, nu, dt);
  //implicit_horizontal_diffusion ({u.x, w}, nu, dt);
}

// Diag of wave energy
event logfile (i++; t <= NT0*T0) // target: at least 300*T0
{
  double ke = 0., gpe = 0.;
  foreach (reduction(+:ke) reduction(+:gpe)) {
    foreach_layer() {
      double norm2 = sq(w[]);
      foreach_dimension()
	norm2 += sq(u.x[]);
      ke += norm2*h[]*dv();
    }
    gpe += sq(eta[])*dv();
  }
  // fixme: divide by L0 if not =1 !
  printf ("%g %g %g\n", t/T0, ke/2., g_*gpe/2.);
}

// Writing netcdf file
event out_nc (i+=3)
{
  write_nc();
}

#if DIAG
// Log the trigger of limiters
event log_clip (i+=10) {
  fprintf (stderr, "%g ", t);
  // ** CFL clip (hydro.h)
  // This is to ensure positiveness of the layer heights
  fprintf (stderr, "%ld %ld ", clip_count, flux_count);
  clip_count = flux_count = 0;
  // ** gradient limiter
  // limiter on the gradient term (minmod but with a counter)
  // see the current file, redefinition of s.gradient for s in list(scalars)
  fprintf (stderr, "%ld %ld ", lim_zero, lim_calls);
  // maybe later add x,y position ?
  lim_zero = lim_calls = 0;
  // ** remap limiter (remapc.h)
  // this limiter is used to avoid creating new extremum when fitting a polynom
  fprintf (stderr, "%ld %ld %ld ", remap_flat, remap_clip, remap_calls);
  remap_flat = remap_calls = remap_clip = 0;
  // ** breaking limiter (nh.h and hydro.h)
  // this limiter is a hard limit on vertical velocity amplitude that is
  // restricted to be always < max_slope*sqrt(g*depth)
  fprintf (stderr, "%ld %ld ", w_clip, w_calls);
  w_clip = w_calls = 0;

  // The line in the 'log' file looks like this
  // time clip_count flux_count lim_zero lim_calls remap_flat remap_clip remap_calls w_clip w_calls
  fprintf (stderr, "\n");
}
#endif // DIAG

/**
~~~pythonplot Wave energy
import numpy as np
import matplotlib.pyplot as plt

g = 9.81
L0 = 100
ak = 0.01
k = 2*np.pi/10
lam = 2*np.pi/k*L0
Re = 400

nu = np.sqrt(g*lam**3)/Re
nu = 1/Re
T0 = (2*np.pi/np.sqrt(g*k))

def E_linwave(E0,nu,ak,k,t):
  print('\nWave decay for linear wave (theory)')
  print(f'nu={nu},ak={ak},k={k/np.pi}pi\n')
  return E0*np.exp(-4*nu*k**2*t)

data = np.loadtxt("out",skiprows=1)
time = data[:,0]
ke  = data[:,1]
gpe = data[:,2]
E = ke + gpe
E0 = E[0]
Eth = E_linwave(E0,nu,ak,k,time*T0)
#E0 = 1 

# Energy plot
fig, ax = plt.subplots(figsize=(5, 5))
ax.plot(time,2*ke/E0,color='b', label='2*ke')
ax.plot(time,2*gpe/E0,color='g', label='2*gpe')
ax.plot(time,E/E0, color='k', label='E')
ax.plot(time, Eth/E0, color='r', label=r'$E(t)=E_0 e^{-4 \nu k^2 t}$')
ax.set_xlabel("t/T0")
ax.set_ylabel("E/E0")
ax.set_ylim([0,1.0])
ax.grid()
ax.legend(loc="lower left")
plt.tight_layout()
plt.savefig("energy.png", dpi=150)


log = np.loadtxt("log", skiprows=2)
# time clip_count flux_count lim_zero lim_calls remap_flat remap_clip remap_calls w_clip w_calls
ltime = log[:,0]
clip_count = log[:,1]
flux_count = log[:,2]
lim_zero = log[:,3]
lim_calls = log[:,4]
remap_flat = log[:,5]
remap_clip = log[:,6]
remap_calls = log[:,7]
w_clip = log[:,8]
w_calls = log[:,9]

CFL_clip = np.where(flux_count==0, 0, clip_count/flux_count)
gradient_limiter = np.where(lim_calls==0, 0, lim_zero/lim_calls)
remap_clip = np.where(remap_clip==0, 0, remap_clip/remap_calls)
remap_flat = np.where(remap_flat==0, 0, remap_flat/remap_calls)
w_clip = np.where(w_calls==0, 0, w_clip/w_calls)

# limiter plot
fig, ax = plt.subplots(figsize=(6, 6))
# -> energy evolution
ax.plot(time,E/E0, color='k', label='E', ls='--')
ax.plot(time, Eth/E0, color='k', label=r'$E(t)=E_0 e^{-4 \nu k^2 t}$')
ax.set_ylabel('E/E0')
ax.set_ylim([0.,1.0])
ax.legend(loc="upper right")
# -> limiters call evolution
ax2 = ax.twinx()
ax2.plot(ltime/T0, 100*CFL_clip, c='g', label='CFL clip')
ax2.plot(ltime/T0, 100*gradient_limiter, c='b',label='gradient_limiter')
ax2.plot(ltime/T0, 100*remap_clip, c='cyan', label='remap_clip')
ax2.plot(ltime/T0, 100*remap_flat, c='yellow', label='remap_flat')
ax2.plot(ltime/T0, 100*w_clip, c='r', label='w_clip')
ax2.set_xlabel('t/T0')
ax2.set_ylabel("limiter applied / limiter call (%)")
ax2.set_title('Evolution of Energy and limiters calls (per 10 steps)')
ax2.legend(loc="lower left")
ax2.set_ylim([0,10])
ax.grid('both')
plt.savefig('limiters.png')
#plt.show()
~~~
*/
