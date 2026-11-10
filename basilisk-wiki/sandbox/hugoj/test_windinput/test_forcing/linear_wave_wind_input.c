/**
# Wave growth using wind forcing

## 1. Decay of a linear wave

### Analytical solution

A linear surface gravity wave will decay in a viscous fluid. For a no stress
boundary condition at the moving surface, Lamb (1932) gives

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

### Notes

## 2. Wind forcing: form drag

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
a_p = \frac{\partial p_s}{\partial x} = \frac{p0}{\rho} \frac{\partial^2 \eta}{\partial x^2} 
\end{aligned}
$$

In the multilayer, $dz=h$. The discretised energy input is (here in 1D):
$$
\begin{aligned}
\frac{\partial E_{in}}{\partial t} = \frac{p0}{\rho}\sum_i^N ( \sum_{k}^{nl} h_k u.x_k) \frac{1}{\Delta^2} (eta_[i+1} + eta[i-1] - 2 eta[i]) \frac{L}{N}
\end{aligned}
$$

This formulation was first used in the multilayer context by Rui Yang (Princeton).

## 3. Exact forcing of to counter viscous dissipation

For a viscous dissipation $\nu$, the pressure $p0$ of the forcing is
$$
\begin{aligned}
 p_0= 4 \rho \nu k c
\end{aligned}
$$

## 4. Dynamic forcing to reach a target energy

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
 p0 = \frac{-\rho g}{\Delta t}( \overline{\eta^2}_{target} - \overline{\eta^2}(t))
        / \sum_N (Q \frac{\partial^2 \eta}{\partial x^2} dx)
\end{aligned}
$$

with $Q= \sum_{k}^{nl} h_k u_k$ the integrated transport.

## References

Horace Lamb, Hydrodynamics (6th ed., 1932), Chapter XI, Article 348, "Effect of
Viscosity on Water-Waves," pp. 623–625, see also article 349 

~~bib
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

~~

Ref Rui Yang

## TODO

- better control (right final value, minimal oscillation) of dynamic forcing using integral quantities

- effect of p0 initial value
- effect of number of layer (geometric, equidistant)
- effect of N
- understand why dissipation rate in 2D is not equal to dissipation rate in 1D

*/

// Dim [L, T, Mass, Temp]

double p0 = 0.00625 [-1, -2, 1]; // Pa
double rho = 1000 [-3, 0, 1];

// Note: here I replaced gmetric(i) with 1.0 as I do not use adaptative mesh.
// this is necessary for the code to compile on CUDA.


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
index such as 'i' for the macro to be compiled on GPU (CUDA).
*/

// this compiles
// #define a_baro(f,k) (1.0*(-G*(f[k]-f[k-1])/Delta - (p0/rho)*windmask.x*(f[k+1]-2.*f[k]+f[k-1])/sq(Delta)))
#define a_baro(f,k) (1.0*(-G*(f[k]-f[k-1])/Delta - (p0/rho)*(f[k+1]-2.*f[k]+f[k-1])/sq(Delta)))

// This doesnt
//#define a_baro(f,k) (gmetric(k)*(-G*(f[k]-f[k-1])/Delta - (p0/rho)*(f[k+1]-2.*f[k]+f[k-1])/sq(Delta)))
// This doesnt either
// #define a_baro(eta,i) \
//     (1.0 * ( \
//       -G*(eta[i] - eta[i-1])/Delta \
//       - (p0/rho) * \
//         (eta[i+1] - 2.*eta[i] + eta[i-1])/(Delta*Delta) \
//   ))



#endif // HOLD_FORCING || EXACT_FORCING


#include "grid/multigrid.h"
#include "layered/hydro.h"
#include "layered/nh.h"
//#include "layered/nh-staggered.h" // Bruno's staggered, replaced in src/layered
#include "layered/remap.h"
#include "layered/perfs.h"
#if OUTNC
#include "bderembl/libs/netcdf_bas.h"
#endif

double pp0 = 0. [-1, -2, 1];

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


// Control
double K_1 = 1. [0,0,0];
double relax_dt = 0. [0,1];

double K_2 = 0. [0,0,0];
double T_f = 1.;
double sum_I = 0.;

double K_3 = 0. [0,0,0];
double last_dE=0.;

// Low pass filter
double dE_f, sum_f, tau_f, sum_abs_f;
// p0 clip
double p0_max = 100. [-1,-2,1];

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

int main()
{
  L0=1.;
  origin (-L0/2.);
  periodic (right);

  N = 64;
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
  the value of p0
  */
  pp0 = 2*rho*nu*k_*cp; // Forcing to exactly balance viscous diss
  p0 = pp0; 
    
  /** 
  If the dynamic forcing is used, we use PID correction.
  */
  dE_f = 0.;      // filtered error
  sum_f = 0.;     // filtered sum
  sum_abs_f = 0.;
  tau_f = T0/2.;  // timescale for low pass filter
  
  // double Ku = 1./(2.*T0); // Ultimate gain
  double Ku = 2./(T0);                   //
  double Tu = 2.*T0;      // Ultimate time scale
  double Td = T0/2.;
  
  K_1 = Ku;               // P
  K_2 = Ku/Tu;            // I, K2 = K_1/Ti
  K_3 = 0.*Ku*Td;         // D, K3 = K_1*Td

  // PI
  // response time = 2π·sqrt(Tu/Ku)
  // damping = sqrt(Ku·Tu)/2



  fprintf(stderr, 
      "T0 = %f, lam=%f, g_=%f,omega=%f,nu=%g\n", T0, lam, g_, omega,nu);

  // study with varying initial p0
  run();

  #if HOLD_FORCING || TESTS_P0
  pp0 /= 2.;
  p0 = pp0;
  dE_f = 0., sum_f = 0.;
  run();

  pp0 /= 2.;
  p0 = pp0;
  dE_f = 0., sum_f = 0.;      
  run();

  pp0 /= 2.;
  p0 = pp0;
  dE_f = 0., sum_f = 0.;
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
      the linear waves, z levels are flat. 

      Using true z levels is adding energy compared to the linear theory as some
      z points are positive, so the exponential in the currents can grow fast
      for steeper cases.
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
  If dynamic forcing is used, the target energy is the initial energy
  */
  energy_stats(&etam_i, &etavar_i, &Ep_i, &Ek_i);
  fprintf (stderr, "INITIAL eta mean = %.10f, variance = %.10f\n", etam_i, etavar_i);
  fprintf (stderr, "initial p0 = %f\n", p0);
  fprintf(stderr, "target Et = %g\n", Ep_i+Ek_i);
  /**
  Outputs for visual inspection 
  */
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
  fprintf (stdout, "%g %g %g %g\n", t/T0, ke, pe, p0);
  #else
  fprintf (stdout, "%g %g %g\n", t/T0, ke, pe);
  #endif 

}
/**
  ## Forcing using a controler

  IN 2D (x,z):
  input_power = -sum_xy(dp/dx Q)dx
             with p = p0.deta/dx
                  Q = sum_layer(u*h)
      
        so input_power = - p0* sum_xy(sum_layer(d2eta/dx2*u*h*dx))


  Energy deficit = E_ini - E(t)
  */
#if HOLD_FORCING
event update_p0 (i++)
{
  /** Computing the error on wave energy */
  double etam = 0., etavar_current = 0., pe = 0., ke = 0.;
  energy_stats (&etam, &etavar_current, &pe, &ke);
  double dE = (Ek_i+Ep_i) - (pe + ke);
  //double dE = 0., etam=0., etavar_current=0.;
  /**
  Computing input power, P = p0*sum
  */
  double sum = 0.; 
    foreach(reduction(+:sum)) {
    double etaxx = (eta[1] + eta[-1]- 2*eta[])/sq(Delta);
    foreach_layer()
      sum += - h[] * u.x[] * etaxx * dv(); 
  }
  /**
  Low pass filter on error and sum
  */
  dE_f += (dE  - dE_f ) * dt / tau_f;
  sum_f += (sum - sum_f) * dt / tau_f;

  /**
  safeguard: never divide by near zero sum  
  */
  sum_abs_f += (fabs(sum) - sum_abs_f) * dt / tau_f;
  double sum_den = (sum_f >= 0. ? 1. : -1.) * fmax (fabs(sum_f), 0.1*sum_abs_f);
  /**
  PID Control

  p0 = 1/sum * K1( dE + 1/Ti*int_t dE dt + Td*dE/dt ) */
  double rhs = 0.;
  
  /** P: proportional constrain */
  if (i>1){
    rhs += K_1 * dE_f;
  }

  /** I: integral constrain */
  if (fabs((rhs + K_2*sum_I)/sum_den) <= p0_max)
  { 
    // avoid to take into account spikes of energy
    sum_I += dE_f * dt;
  }
  //sum_I += dE_f * dt;
  //rhs += K_2 * sum_I;

  /**  D: derivative constrain */  
  double dedt = 0.;  
  dedt = (dE_f-last_dE)/dt;
  rhs += K_3*dedt;
  last_dE = dE_f;
  /** Finally we compute p0 */
  double dp0 = rhs/sum_den;

  p0 = pp0 + dp0;

  /** safeguard 2: p0 limiter */
  if (p0 >  p0_max) p0 =  p0_max;
  if (p0 < -p0_max) p0 = -p0_max;

  // if (i%10==0){
  //   fprintf(stderr, "i=%d etam=%.3e var=%.6e dE=%g sum=%g sum_den=%g p0=%g\n",
  //       i, etam, etavar_current, dE, sum, sum_den, p0);
  // }
  fprintf(stderr, "i=%d etam=%.3e var=%.6e dE=%g sum=%g sum_den=%g p0=%g\n",
        i, etam, etavar_current, dE, sum, sum_den, p0);
}
#endif 

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


// #if HOLD_FORCING || TESTS_P0
// event save_p0 (t>NT0*T0 - dt)
// {
//   char mycmd[100];
//   sprintf (mycmd, "mv out out_p0%.03f", pp0);
//   system (mycmd);
//   // system ("rm -f out");
//   // system ("touch out");
// }
//
// #endif

/**
NOTES:

2/11/26
- 2D CPU (HPC) has a big energy decrease at startup (~0.3%)
- 2D CPU (HPC) is more dissipative than GPU (local)
- 2D GPU (local) is more dissipative than theory
- 2D GPU (local) with exact forcing. E grows with time whereas 2D CPU (HPC)
decay with time.



*/
