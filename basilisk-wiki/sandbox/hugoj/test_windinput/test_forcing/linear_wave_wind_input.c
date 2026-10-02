/**
# Wave growth using wind forcing

## 1. Decay of a linear wave

### Analytical solution

A linear surface gravity wave will decay in a viscous fluid. The rate of this
decay is $E(t)=E_0 e^{-4\nu k^2 t}$ (Lamb 1932)

One of the essential ingredient for the dissipation from the model to match the
theoretical dissipation is to correctly represent currents, i.e. have enough
layers for the vertical exponential decay discretisation.

### Notes

## 2. Wind forcing: form drag

We use the following principle: a wind pressure is applied on positive slopes
(in the x direction) with the norm
$$
\begin{aligned}
p_s(t,x,y) = \frac{p_0(t)}{\rho} \frac{\partial \eta}{\partial x}
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

Ref Rui Yang

## TODO

- better control (right final value, minimal oscillation) of dynamic forcing using integral quantities

*/

// Faire le cas test !
// - propre
// - avec dimensions
// - simplifier

#include "grid/multigrid1D.h"
double p0 = 0.00625 [-1, -2, 1]; // Pa
double rho = 1000 [-3, 0, 1];
#if HOLD_FORCING || EXACT_FORCING
// wind pressure param
// #define wind_pressure(eta,i)  (p0/rho)*(eta[i+1] - eta[i-1])/(2*Delta)
// // total barotropic pressure
// #define p_baro(eta,i) (-G*eta[i] - wind_pressure(eta,i))
//
// #define a_baro(eta,i) (gmetric(i)*(p_baro(eta,i)-p_baro(eta,i-1))/Delta) 
#define a_baro(eta,i) \
  (gmetric(i) * ( \
      -G*(eta[i] - eta[i-1])/Delta \
      - (p0/rho) * \
        (eta[i+1] - 2.*eta[i] + eta[i-1])/(Delta*Delta) \
  ))
#endif // HOLD_FORCING || EXACT_FORCING

#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
#include "layered/perfs.h"
#if OUTNC
#include "bderembl/libs/netcdf_bas.h"
#endif

// Dim [L, T, Mass, Temp]

double g_ = 1. [1,-2];
/**
We use a linear wave mode, with very low steepness
*/
double k_ = 2.*pi [-1], h_ = 1. [1], ak = 0.01; 
double RE = 20000.;

#include "hugoj/lib/common_waves.h"

#define T0  (2.*pi/sqrt(g_*k_)) // = sqrt(2PI) = 2.5s
double lam = 0. [1];
double NT0 = 10. [0, 1];
double etam_i = 0. [1];
double etavar_i = 0. [2];
double etavar_previous =0. [2], etavar_current=0. [2];
double cp = 0.[1, -2];
double omega = 0. [0,-1];
double relax_dt = 0. [0,1];

int main()
{
  L0=1.;
  origin (-L0/2.);
  periodic (right);
  N = 64;
  nl = 100; // 60
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
  #if EXACT_FORCING 
  p0 = 2*rho*nu*k_*cp; // Forcing to exactly balance viscous diss
  #endif
  /** If the dynamic forcing is used, we set a timescale for the relaxation of
    the forcing.
    */
  relax_dt = 5*T0; 
  fprintf(stderr, "T0 = %f, lam=%f, g_=%f,omega=%f,nu=%g\n", T0, lam, g_, omega,
          nu);
  run();
}


/**
 Here a function is defined: it computes the variance and the mean of the
 surface elevation in this 1D model.
 */
void eta_stats (double *mean, double *variance)
{
  double sum = 0., var = 0.;
  foreach (reduction(+:sum))
    sum += eta[]; 
  *mean = sum/N;
  foreach (reduction(+:var))
    var += sq(eta[] - *mean);
  *variance = var/N; 
}

event init (i = 0)
{
  //geometric_beta (1./5., true); // crashes when more cells near sfx ?
  // default beta is 1/nl
  foreach() {
    zb[] = -h_;
    #if STOKES
    eta[] = wave_stokes (x, 0., ak, k_, h_);
    //eta[] = ak/k_*cos(k_*x);
    #else
    eta[] = wave_monolin(0., x, ak/k_, k_);
    #endif
    double H = eta[] - zb[];
    double z = zb[];
    foreach_layer() {
      h[] = H*beta[point.l];
      /** 
      In the linear theory, $\eta$ is almost 0 and z levels are flat. 

      Using true z levels is adding energy compared to the linear theory as some
      z points are positive, so the exponential in the currents can grow fast
      for steeper cases.
      */
      z += h_/nl/2; 
      u.x[] = u_x_monolin(0., x, z, ak/k_, k_); //ak/k_*sqrt(g_*k_)*exp(k_*z)*cos(k_*x); // 
      w[] = u_y_monolin(0., x, z, ak/k_, k_); //ak/k_*sqrt(g_*k_)*exp(k_*z)*sin(k_*x);  //
      z += h_/nl/2; 
    }

  }
  // Compute initial wave energy
  eta_stats (&etam_i, &etavar_i);
  // If dynamic forcing is used, the target energy is the initial energy
  fprintf (stderr, "INITIAL eta mean = %.10f, variance = %.10f\n", etam_i, etavar_i);
  fprintf (stderr, "initial p0 = %f\n", p0);
  #if OUTNC
  create_nc({zb, eta, h, u.x, w}, "out.nc");
  #endif
}

/* compute p0 from eta in a 'update_p0' event so that it uses eta from previous timestep,
   before applying the forcing using the macro 'a_baro'
*/

#if HOLD_FORCING
event update_p0 (i++)
{
  double etam;
  eta_stats (&etam, &etavar_current);
  double dE = G*(etavar_i - etavar_current);
  double sum = 0.;
  foreach() {
    double integrated_transport=0.;
    double etaxx = (eta[1] + eta[-1]- 2*eta[])/sq(Delta);
    foreach_layer()
      integrated_transport += h[] * u.x[]; // h * detadx * u
    sum += integrated_transport*etaxx*dv();
  }
  p0 = -rho*dE/(relax_dt*sum);
  // fprintf(stderr,
  //       "i=%d dE=%g sum=%g dt=%g p0=%g predicted=%g\n",
  //       i, dE, sum, relax_dt, p0,
  //       -p0*sum*relax_dt/rho);
}

#endif 

/** 
Isotropic viscous dissipation
*/
event viscous_term (i++) {
  // vertical diffusion is done in diffusion.h (u.x) and nh.h (w)
  horizontal_diffusion ({u.x, w}, nu, dt);
}


/** 
We log potential and kinetic energy
 */
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
  fprintf (stdout, "%g %g %g\n", t/T0, ke/2., g_*gpe/2.);
}


/**
  Optional: output fields for analysis
*/
#if OUTNC
event writenc (i+=10; t <= NT0*T0) 
{
  write_nc();
}
#endif



/**
~~~pythonplot Wave energy
import numpy as np
import matplotlib.pyplot as plt

g = 1.0
L0 = 1.0
ak = 0.01
k = 2 * np.pi / L0
lam = 2 * np.pi / k
Re = 20000
c = np.sqrt(g * k) / k
NT0 = 10

nu = c * lam / Re
T0 = 2 * np.pi / np.sqrt(g * k)

print("T0 = %f" % T0)


def E_linwave(E0, nu, ak, k, t):
    print("\nWave decay for linear wave (theory)")
    print(f"nu={nu},ak={ak},k={k / np.pi}pi\n")
    return E0 * np.exp(-2 * nu * k**2 * t)


data = {
    "current_noforcing": np.loadtxt("../no_forcing/out", skiprows=1),
    "current_exactforced": np.loadtxt("../exact_forcing/out", skiprows=1),
    "current_forced": np.loadtxt("../linear_wave_wind_input/out", skiprows=1),
}
time = {}
ke = {}
gpe = {}
E = {}

for case in data.keys():
    time[case] = data[case][:, 0]
    ke[case] = data[case][:, 1]
    gpe[case] = data[case][:, 2]
    E[case] = ke[case] + gpe[case]

E0th = 0.5 * g * (ak / k) ** 2  # E["current_noforcing"][0]
print("\ntheoretical E0 = %g m3/s2" % E0th)
print("initial energy for current sim E0 = %g m3/s2" % E["current_noforcing"][0])
print("ratio is %g \n" % (E0th / E["current_noforcing"][0]))
Eth = E_linwave(E0th, nu, ak, k, time["current_noforcing"] * T0)
E0 = 1

print("Target for E=cst for this set of parameters:")
print("nu=%g, ak=%g, k=%g" % (nu, ak, k))
print("E0 theory = %g" % E0th)
print("E0 real = %g (need more resolution to match E0th)" % E["current_noforcing"][0])
print("Initial energy missmatch = E0th - E0 = %g" % (E0th - E["current_noforcing"][0]))
print(
    "noforcing: miss match Eth-E=%g (at %d T0) "
    % (Eth[-1] - E["current_noforcing"][-1], NT0)
)
print(
    "exact forcing: miss match Eth-E=%g (at %d T0) "
    % (Eth[-1] - E["current_forced"][-1], NT0)
)

fig, ax = plt.subplots(figsize=(7, 6))

# ax.plot(
#     time["current_noforcing"],
#     2 * gpe["current_noforcing"] / E0,
#     color="g",
#     label="2Ep (current)",
#     alpha=0.5,
# )
# ax.plot(
#     time["current_noforcing"],
#     2 * ke["current_noforcing"] / E0,
#     color="b",
#     label="2Ek (current_noforcing)",
#     alpha=0.5,
# )
ax.semilogy(time["current_noforcing"], Eth / E0, color="r", label=r"$E(t)=E_0 e^{-2 \nu k^2 t}$")
ax.hlines(E0th, 0, 100, colors="gray", alpha=0.7)
ax.semilogy(
    time["current_noforcing"],
    E["current_noforcing"] / E0,
    color="pink",
    label="E (no forcing)",
    ls='--',
)
ax.semilogy(
    time["current_forced"],
    E["current_forced"] / E0,
    color="purple",
    label="E (hold forcing)",
    ls='--'
)
ax.semilogy(
    time["current_exactforced"],
    E["current_exactforced"] / E0,
    color="orange",
    label="E (exact forcing)",
    ls='--'
)

ax.set_xlabel("t/T0")
ax.set_ylabel("E/E0")
ax.set_xlim([0, 10])
ax.set_ylim([1.2e-6, 1.28e-6])
ax.legend()
ax.grid(axis="y", which="both")
plt.savefig("energy.pdf", dpi=300)
plt.show()

# plot [:2] "out" using 1:($2+$3) w l
~~~
*/
