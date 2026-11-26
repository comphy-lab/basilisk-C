/**
# Wave growth using wind forcing

## Decay of a linear wave

### Analytical solution

A linear surface gravity wave will decay in a viscous fluid. For a no stress
boundary condition at the moving surface, Lamb (1932, pp. 623–625, see also article 349) gives

\begin{aligned}
  E(t)=E_0 e^{-4\nu k^2 T_0 \tilde{t}}
\end{aligned}

with $\tilde{t} = t/T_0$.

For the case with a free slip condition at the surface, we get (Dore 1978,
Joseph and Wang 2004)

$$
\begin{aligned}
  E(t)=E_0 e^{-2\nu k^2 T_0 \tilde{t}}
\end{aligned}
$$

One of the essential ingredient for the dissipation from the model to match the
theoretical dissipation is to correctly represent currents, i.e. have enough
layers for the vertical exponential decay discretisation.

## Wind forcing

### Form drag

We use the following principle: a wind pressure is applied on positive slopes
(in the x direction) with the norm
$$
\begin{aligned}
p_s(t,x,y) = -\frac{p_0(t)}{\rho} \frac{\partial \eta}{\partial x}
\end{aligned}
$$
This pressure is added to the barotropic pressure from the deformation of the
surface (much like in the vein of the [hydro-tension](https://basilisk.fr/src/layered/hydro-tension.h) code).
The 'a_baro' macro is overloaded.

The amplitude $p_0$ can be set to a constant (growing sea) or maintained at a
specific energy level.

The energy input from this forcing in the multilayer simulation is
$$
\begin{aligned}
\frac{\partial E_{in}}{\partial t} = \int_x \int_x u \cdot a_p dx dz
\end{aligned}
$$
where $a_p$ is the barotropic acceleration 
$$
\begin{aligned}
a_p = \frac{\partial p_s}{\partial x} = \frac{pw}{\rho} \frac{\partial^2 \eta}{\partial x^2} 
\end{aligned}
$$

In 1D, the energy input is:
$$
\begin{aligned}
\frac{\partial E_{in}}{\partial t} = \int_x \int_z u \frac{p_0}{\rho}
\frac{\partial \eta}{\partial x} dx dz
\end{aligned}
$$

In the multilayer, $dz=h$, the discretised version is:
$$
\begin{aligned}
\frac{\partial E_{in}}{\partial t} = \frac{p_0}{\rho}\frac{1}{\Delta}\sum_i^N ( \sum_{k}^{nl} h_k u.x_k) (\eta[i+1] + \eta[i-1] - 2 \eta[i])
\end{aligned}
$$

This formulation was first used in the multilayer context by Rui Yang (Princeton).

### Exact forcing of to counter viscous dissipation

For a viscous dissipation $\nu$, the pressure $p_0$ of the forcing is
$$
\begin{aligned}
 p_0= 2 \rho \nu k c
\end{aligned}
$$

### Dynamic forcing to reach a target energy

See this [file](https://basilisk.fr/sandbox/hugoj/lib/wind_forcing.h).

Let's say we want to maintain a quasi-stationnary sea state. We can do this
by holding the the variance of the surface height to a constant value. At each
timestep, dissipation occurs (viscous, breaking or implicit) so energy must be
injected into the domain. We can compute an energy deficit $\Delta E$ like the
following

$$
\begin{aligned}
\Delta E = \rho g ( \overline{\eta^2}_{target} - \overline{\eta^2}(t))
\end{aligned}
$$
We relax the forcing on a $\Delta t$ timescale (chosen by the user, typically a
few period of the main wave). The pressure amplitude $p_0$ can then be inferred
by equating the energy deficit over the timescale $\Delta t$ and the energy
input by the present forcing

$$
\begin{aligned}
\frac{\partial E_{in}}{\partial t} = \frac{\Delta E }{\Delta t }
\end{aligned}
$$
Rearranging terms gives the amplitude for the current timestep
$$
\begin{aligned}
 p_0 = \frac{-\rho g}{\Delta t}( \overline{\eta^2}_{target} - \overline{\eta^2}(t))
        / \sum_{x,y} (Q \frac{\partial^2 \eta}{\partial x^2} dx)
\end{aligned}
$$

with $Q= \sum_{k}^{nl} h_k u_k$ the integrated zonal transport.

This approach is naive for two reasons:

 * Potential energy usually oscillates a lot for a real sea state, while the
    total energy is much more regular. An oscillating error mesurement can induce
    instability with the proportionnal controller presented above.

  * A proportionnal controller will never reach the target value, an integral
    component is needed.

In order to improve the controller robustness for a real sea state and to shocks
in the solution (e.g. initialisation), we propose to use a
proportionnal-integral controller (PI). The $p_0$ value at time $t$ is computed as:
$$
\begin{aligned}
 pw = rhs / \sum_{x,y} (Q \frac{\partial^2 \eta}{\partial x^2} dx)
\end{aligned}
$$

with
$$
\begin{aligned}
 rhs = K_1 \Delta E + K_2 \int_0^t \Delta E dt
\end{aligned}
$$

The values of $K_1$ and $K_2$ are set to minimise the response time of the
system while having a small overshoot. After some tuning, we use the values 
$K_1 = 2/T_0$ and $K_2 = 1/T_0^2$ with $T_0$ the period of the monochromatic wave (this file) or the main
peak period for a real seastate ([here](https://basilisk.fr/sandbox/hugoj/test_windinput/test_real_field/PM_seastate.c))

## How to run this code

~~~bash
make
~~~

## Results

### 1D linear wave

In 1D for a linear wave, the multilayer code is able to reproduce the theoretical dissipation (red vs
    pink curves). 

When setting a constant $p_0=2 \rho \nu k c$, energy is close to
constant (orange vs grey curves). 

The dynmic forcing is fast: even if the initial value of $p_0$ is far from the
theoretical value, the controller adjust fast enough the total energy so that it
always stays in +/- 5% of the target (= initial) value. We note that energy
curves computed from the simulation are not starting at exactly 1.0, this is
because of the discretisation error of the initial condition. Increasing the
number of layer reduces this error. For dynamic forcing cases, we target the
numerical initial energy, not the theoretical one, this is why the end value of
energy is slightly less than the $E/E_0=1.0$ mark.

![Energy evolution for a 1D linear wave](energy_1D_linear.png){ width="60%" }

### 2D linear wave

### 1D Stokes wave

### 2D Stokes wave

## References

~~~bib
@book{lambHydrodynamics1932,
  title = {Hydrodynamics},
  author = {Lamb, Horace},
  year = 1932,
  publisher = {Cambridge University Press},
}

@article{doreEffectsAirwaterInterface1978,
  title = {Some Effects of the Air-Water Interface on Gravity Waves},
  author = {Dore, B. D.},
  year = 1978,
  month = jan,
  journal = {Geophysical \& Astrophysical Fluid Dynamics},
  volume = {10},
  number = {1},
  pages = {215--230},
  issn = {0309-1929, 1029-0419},
  doi = {10.1080/03091927808242638},
  urldate = {2026-10-23},
  langid = {english},
}

@article{josephDissipationApproximationViscous2004,
  title = {The Dissipation Approximation and Viscous Potential Flow},
  author = {Joseph, D. D. and Wang, J.},
  year = 2004,
  month = apr,
  journal = {Journal of Fluid Mechanics},
  volume = {505},
  pages = {365--377},
  issn = {0022-1120, 1469-7645},
  doi = {10.1017/S0022112004008602},
  urldate = {2026-10-29},
  abstract = {Dissipation approximations have been used to calculate the drag on
              bubbles and drops and the decay rate of free gravity waves on
              water. In these approximations, viscous effects are calculated by
              evaluating the viscous stresses on irrotational flows. The pressure
              is not involved in the dissipation integral, but it enters into the
              power of traction integral, which equals the dissipation. A viscous
              correction of the irrotational pressure is needed to resolve the
              discrepancy between the zero-shear-stress boundary condition at a
              free surface and the non-zero irrotational shear stress. Here we
              show that the power of the pressure correction is equal to the
              power of the irrotational shear stress. The viscous pressure
              correction on the interface can be expressed by a harmonic series.
              The principal mode of this series is matched to the velocity
              potential and its coefficient is explicitly determined. The other
              modes do not enter into the expression for the drag on bubbles and
              drops. They vanish in the case of free gravity waves.},
  copyright = {https://www.cambridge.org/core/terms},
  langid = {english},
}
~~~

Ref Rui Yang
*/

// Dim [L, T, Mass, Temp]
double rho = 1000 [-3, 0, 1];
double pw = 0. [-1, -2, 1]; // Pa
/** The forcing is implemented as an overload of 'a_baro' macro */
#if HOLD_FORCING || EXACT_FORCING

/**
Wind acts in x only: qcc rewrites .x -> .y inside the
foreach_dimension()/foreach_face() passes of nh.h/hydro.h, so
multiplying the wind term by windmask.x makes it vanish in the
y-passes.
*/
coord windmask = {1, 0, 0}; 

/** Definition of the macro. 
  This has to be in 1 expression (no nested macros) and not using already used
  index such as 'i' for the macro to be compiled on GPU (CUDA). Note also that
  here I replaced gmetric(i) with 1.0 as I do not use adaptative mesh.
  this is necessary for the code to compile on CUDA.
*/

// this compiles
// #define a_baro(f,k) (1.0*(-G*(f[k]-f[k-1])/Delta - (pw/rho)*windmask.x*(f[k+1]-2.*f[k]+f[k-1])/sq(Delta)))
// #define a_baro(f,k) (1.0*(-G*(f[k]-f[k-1])/Delta - (pw/rho)*(f[k+1]-2.*f[k]+f[k-1])/sq(Delta)))

// This doesnt
//#define a_baro(f,k) (gmetric(k)*(-G*(f[k]-f[k-1])/Delta - (pw/rho)*(f[k+1]-2.*f[k]+f[k-1])/sq(Delta)))
// This doesnt either
// #define a_baro(eta,i) \
//     (1.0 * ( \
//       -G*(eta[i] - eta[i-1])/Delta \
//       - (pw/rho) * \
//         (eta[i+1] - 2.*eta[i] + eta[i-1])/(Delta*Delta) \
//   ))

// overload of the a_baro macro to add wind pressure
#define a_baro(f,k) (1.0*(-G*(f[k]-f[k-1])/Delta - (pw/rho)*(f[k+1]-2.*f[k]+f[k-1])/sq(Delta)))
//#define a_baro_p(f,k) (- (pw/rho)*(f[k+1]-2.*f[k]+f[k-1])/sq(Delta))

#endif // HOLD_FORCING || EXACT_FORCING

/** regular imports */
#include "grid/multigrid.h"
#include "layered/hydro.h"
#if NH-S
#include "layered/nh-staggered.h" // Bruno's staggered, replaced in src/layered
#else
#include "layered/nh.h"
#endif // NH-S
#include "layered/remap.h"
#include "layered/perfs.h"
#if OUTNC
#include "bderembl/libs/netcdf_bas.h"
#endif


/** ## Parameters */
double g_ = 1. [1,-2];
double k_ = 2.*pi [-1], h_ = 2. [1]; 
double RE = 20000. [0,0,0];

#if STOKES
/**
For the Stokes wave, a high steepness is used
*/
double ak = 0.35 [0,0,0]; 
#else
/**
For the linear wave, a very low steepness is used
*/
double ak = 0.01 [0,0,0];
#endif //STOKES

#include "hugoj/lib/common_waves.h"

#define T0  (2.*pi/sqrt(g_*k_)) // = sqrt(2PI) = 2.5s
double lam = 0. [1];
double NT0 = 10. [0, 0, 0]; // Duration of the experiment
double cp = 0.[1, -2];
double omega = 0. [0,-1];

// Initial state
double Ek_i = 0. [2, -2, 1];
double Ep_i = 0. [2, -2, 1];
double etam_i = 0. [1];
double etavar_i = 0. [2];

/** A function to compute the energy of the simulation */
void energy_stats (double *mean, double *variance, double *pe, double *ke)
{
  double sum = 0.;
  foreach (reduction(+:sum))
    sum += eta[];
  double m = sum/pow (N, dimension);

  double gpe = 0., ke2 = 0.;
  foreach (reduction(+:gpe) reduction(+:ke2)) {
    gpe += sq (eta[] - m)*dv();
    foreach_layer() {
      double norm2 = sq (w[]);
      foreach_dimension()
        norm2 += sq (u.x[]);
      ke2 += norm2*h[]*dv();
    }
  }
  *mean = m;
  *variance = gpe/sq(L0);  // <(eta-m)^2> since sum of dv() = L0^dimension
  *pe = 0.5*rho*G*gpe;
  *ke = 0.5*rho*ke2;
}

/** wind forcing controler is added with this [file](https://basilisk.fr/sandbox/hugoj/lib/wind_forcing.h) */
#if HOLD_FORCING
#include "hugoj/lib/wind_forcing.h"
#endif

// Note
// For now I have to use this face_fields event to add the forcing.
// When i'll be able to use a vector inside the a_baro macro I will 
// be able to use the original machinery with the a_baro macro overloaded.
#if 0
#if HOLD_FORCING || CST_FORCING
event face_fields (i++,first) {
  foreach_face (x) {
    double ax = a_baro_p (eta, 0);
    foreach_layer() {
      ha.x[] += hf.x[]*ax;
    }
  }
}
#endif // HOLD_FORCING || CST_FORCING
#endif


/** main program */

int main()
{
  L0=1.;
  origin (-L0/2.);
  periodic (right);

  N = 128;
  nl = 64;
  G = g_;
  lam = 2*pi/k_*L0;
  
  omega = sqrt(g_*k_);
  cp = omega/k_;
  nu = cp*lam/RE; 
  theta_H=0.5;         // scheme conserve energy when theta_H = 0.5
  DT=0.01;                // fixed DT to study spatial and temporal resolution separately
  NITERMIN=3;             // Forces to do more cycles to avoid any influence of the poisson solver
  
  /**
  In the case of counter balancing the viscous dissipation exactly, we set once
  the value of pw
  */
  
  pw = 2*rho*nu*k_*cp; // Forcing to exactly balance viscous diss
    
  /** 
  If the dynamic forcing is used, we use PI correction.
  */
  #if HOLD_FORCING
  ppw = 2*rho*nu*k_*cp;
  tau_f = T0/2.;  // timescale for low pass filter
  
  double Ku = 2./(T0);        // Ultimate gain
  double Tu = 2.*T0;          // Ultimate time scale
  double Td = T0/2.;
  
  K_1 = Ku;               // P
  K_2 = Ku/Tu;            // I, K2 = K_1/Ti
  K_3 = 0.*Ku*Td;         // D, K3 = K_1*Td
  #endif 

  fprintf(stderr, 
      "T0 = %f, lam=%f, g_=%f,omega=%f,nu=%g\n", T0, lam, g_, omega,nu);

  run();

  /** Study with varying initial pw */
  #if HOLD_FORCING && TESTS_P0
  ppw = 0.05;
  pw = ppw;
  dE_f = 0., sum_f = 0., sum_abs_f = 0., sum_I = 0.;
  run();

  ppw = 0.01;
  pw = ppw;
  dE_f = 0., sum_f = 0., sum_abs_f = 0., sum_I = 0.;      
  run();

  ppw = 0.;
  pw = ppw;
  dE_f = 0., sum_f = 0., sum_abs_f=0., sum_I = 0.;
  run();
  #endif
}



event init (i = 0)
{

  geometric_beta (1./3., true); // default beta is 1/nl (1.)
  
  foreach() {
    zb[] = -h_;
    #if STOKES
    eta[] = wave_stokes (x, 0., ak, k_, h_);
    #else
    eta[] = wave_monolin(0., x, ak/k_, k_);
    #endif
    double H = eta[] - zb[];
    double z = zb[];
    foreach_layer() {
      h[] = H*beta[point.l];
      
      /** 
        We use true z levels for the initialisation of the Stokes wave, but for
        the linear waves, z levels are flat. Using true z levels is adding
        energy compared to the linear theory as some z points are positive, so
        the exponential in the currents can grow fast for steeper cases.
      */

      #if STOKES
      z += H*beta[point.l]/2.; 
      u.x[] = u_x_stokes (x, z, ak, k_, h_);
      w[] = u_y_stokes (x, z, ak, k_, h_);
      z += H*beta[point.l]/2.; 
      #else
      z += h_*beta[point.l]/2.; 
      u.x[] = u_x_monolin(0., x, z, ak/k_, k_); 
      w[] = u_y_monolin(0., x, z, ak/k_, k_); 
      z += h_*beta[point.l]/2.; 
      #endif
    }

  }
  /** 
  The initial energy of the wavefield is computed.
  */
  energy_stats(&etam_i, &etavar_i, &Ep_i, &Ek_i);
  fprintf (stderr, "INITIAL eta mean = %.10f, variance = %.10f\n", etam_i, etavar_i);
  fprintf (stderr, "initial pw = %f\n", pw);
  /** If dynamic forcing is used, the target energy is the initial energy */
  #if HOLD_FORCING
  Etarget = Ep_i + Ek_i;
  fprintf(stderr, "target Et = %g\n", Etarget);
  #endif
  /** Outputs for visual inspection */
  #if OUTNC
  create_nc({zb, eta, h, u, w}, "out.nc");
  #endif

  
}
/** 
  We log potential and kinetic energy. 
 */
event logfile (i++; t <= NT0*T0) // target: at least 300*T0 
{ 
  double etam = 0., etavar_current = 0., pe = 0., ke = 0.;
  energy_stats (&etam, &etavar_current, &pe, &ke);
  #if HOLD_FORCING || EXACT_FORCING
  fprintf (stdout, "%g %g %g %g\n", t/T0, ke, pe, pw);
  #else
  fprintf (stdout, "%g %g %g\n", t/T0, ke, pe);
  #endif 
}

/** 
Isotropic viscous dissipation
*/
event viscous_term (i++) {
  // vertical diffusion is done in diffusion.h (u.x) and nh.h (w)
  horizontal_diffusion ({u, w}, nu, dt);
}

/** 
 On GPU, float 32 are used. This introduces a number of round error in the code.
 The main one is by the remapping procedure that makes the total volume of
 water drift (sign of the drift depends on the fp32 representation of nl). This
 has been fixed by enforcing that sum(h[]) = H in remap.h.

 A second one could be in the advection scheme but I still need to look at this.
*/

/**
  Optional: output fields for analysis
*/
#if OUTNC
event writenc (i++; t <= NT0*T0) 
{
  if (i%10==0){
    write_nc();
  }
}
#endif


/**
## TODO

- effect of number of layer (geometric, equidistant)
- effect of N
- understand why dissipation rate in 2D is not equal to dissipation rate in 1D
- run the 2D case with the patch that fixes BC on GPU
- use a mask to apply a_baro only along X (parser bug for now on GPU)
- add figures for 1D stokes, 2D linear, 2D stokes

## NOTES:

2/11/26

- 2D CPU (HPC) has a big energy decrease at startup (~0.3%)
- 2D CPU (HPC) is more dissipative than GPU (local)
- 2D GPU (local) is more dissipative than theory
- 2D GPU (local) with exact forcing. E grows with time whereas 2D CPU (HPC)
decay with time.

*/
/**
~~~pythonplot
import os
os.system("cd ..; python3 script.py")
~~~
*/
