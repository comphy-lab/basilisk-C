/**
# Held & Suarez climate model

This example illustrates how the hydrostatic [multilayer
solver](/src/layered/README) can be used to build the simple, but
realistic, climate model proposed by [Held & Suarez, 1994](#held1994)
as a benchmark for the dynamical cores of climate models. See also
[Held, 2005](#held2005) for an interesting general discussion on why
it is important to build such idealised models.

The movie below illustrates the evolution of the surface temperature
when a quasi steady state is reached. The horizontal scale is 360
degrees of longitude and the vertical $\pm$ 75 degrees of latitude.

![Surface temperature evolution. Dark red is 40$^\circ$C and dark blue
 -9$^\circ$C.](held-suarez/theta0-short.mp4)(loop autoplay width="80%")

Several realistic features of the Earth's atmosphere are reproduced
including the characteristic wavenumber of the primary baroclinic
instabilities (between 5 and 7 modes), the "dorsals" extending from
the tropics into the mid-latitudes, the "storm tracks", fronts and
curling depressions characteristic of the mid-latitudes etc.

The surface vorticity gives another view of the surface dynamics.

![Surface vorticity evolution. Dark red is 10^-4^ s^-1^ and dark blue
 -10^-4^ s^-1^.](held-suarez/omega0-short.mp4)(loop autoplay width="80%")

The three-dimensional structure is further investigated and compared
quantitatively with previous results in the [post-processing
example](held-suarez-postproc.c).
 
## Fluid isomorphism

The approach relies on exploiting the isomorphism between the
(oceanic) hydrostatic, Boussinesq equations of motion of an
incompressible fluid, with the fluid depth as vertical coordinate and
the (atmospheric) hydrostatic equations of motion of a compressible
fluid, with the pressure as vertical coordinate. Following [Marshall
et al. 2004](#marshall2004), the compressible, hydrostatic,
atmospheric equations in pressure coordinates can be written
$$
\begin{aligned}
  \frac{D}{D t} \mathbf{v}_h + f\mathbf{k} \times \mathbf{v}_h +
  \mathbf{{\nabla}}_p \Phi & = \mathcal{F}\\
  \frac{\partial \Phi}{\partial p} + \alpha & = 0\\
  \mathbf{{\nabla}}_p \cdot \mathbf{v}_h + \frac{\partial \omega}{\partial
  p} & = 0\\
  \alpha & = \alpha (\theta, p) = \frac{\partial \Pi}{\partial p} \theta\\
  \frac{D}{D t} \theta & = \frac{\mathcal{Q}_{\theta}}{\Pi}
\end{aligned}
$$
with $\mathbf{v}_h$ the horizontal velocity, $\omega=(D/Dt)p$ the
vertical velocity in pressure coordinates, $f$ the Coriolis parameter,
$\Phi = gz$ the geopotential, $\alpha$ the specific volume, $\theta =
c_pT/\Pi$ the potential temperature, with $T$ the temperature and $\Pi
= c_p(p/p_0)^\kappa$ the Exner function. We have omitted the equation
for the specific humidity, which is not modelled in the Held-Suarez
idealised case.

We will show how we can use the hydrostatic multilayer solver to solve
this system. We do this on a regular multigrid, in spherical
coordinates (longitude-latitude) and with an implicit time-integration
for the surface pressure. */

#include "grid/multigrid.h"
#include "spherical.h"
#include "layered/hydro.h"
#include "layered/implicit.h"
#include "layered/remap.h"
#include "layered/perfs.h"
#include "profiling.h"

/**
We declare a few constants and fields which will be used for time
averages. */

const double day = 86400.;

scalar mu, mtheta, vtheta;

/**
## Coriolis acceleration and linear friction

Coriolis acceleration takes its standard form. */

const double Omega = 7.292205e-5;
#define F0() (2.*Omega*sin(y*pi/180.))

/**
The pressure is obtained from the reference pressure $p_0$ and the $r$
coordinate.  */

const double P0 = 101300.;
#define Ps(r) (P0 - (r))

/**
Following Held & Suarez, a linear friction coefficient is applied,
which varies in the vertical according to
$$
k_v = k_f \max \left(0,\frac{\sigma - \sigma_b}{1 - \sigma_b}\right)
$$
*/

const double kf = 1./day, sigmab = 0.7;
#define K0() (kf*max(0., (Ps((point.l + 0.5)/nl*P0)/P0 - sigmab)/(1. - sigmab)))

#include "layered/coriolis.h"

/**
## Potential temperature and buoyancy

The buoyancy is added as a [Boussinesq buoyancy](/src/layered/dr.h) in
the isomorphic formulation (eqs. 38 and 31 of Marshall et al.) i.e.
$$
b = - g \frac{\rho(\theta)}{\rho_0} = - g \Delta\rho(\theta) =
- \alpha(\theta) = - \frac{\partial\Pi}{\partial p} \theta
$$
where we have used the notations in [dr.h](/src/layered/dr.h). The isomorphism implies that
$$
\Delta\rho(\theta) = \frac{1}{g}\frac{\partial\Pi}{\partial p} \theta
$$
*/

const double Cp = 1004., Kappa = 2./7. [0];
#define PI(p) (Cp*pow((p)/P0, Kappa))
#define dPIdp(p) (Cp*Kappa/P0*pow((p)/P0, Kappa - 1.))

/**
To compute $\partial\Pi/\partial p$ we need the pressure in (the
middle of) the layer indexed by `point.l`. We assume that the layers
are equally spaced in pressure from `P0` at the bottom of the first
layer (`point.l = 0`) to zero at the top boundary of the top layer
(`point.l = nl - 1`). */

#define drho(theta) (- dPIdp(Ps((point.l + 0.5)/nl*P0))*theta/G)

#include "layered/dr.h"

/**
## Equilibrium temperature

We define the equilibrium temperature distribution as in
Held & Suarez (p. 1826 and figure 1.b) i.e.
$$
T_{eq} = \max\left(200K, \left[315K - \Delta T_y\sin^2\phi -
\Delta T_z\log(p/p_0)\cos^2\phi\right]\frac{\Pi}{C_p}\right)
$$
The equilibrium *potential* temperature is then obtained as
$$
\theta_{eq} = C_p\frac{T_{eq}}{\Pi}
$$
*/

double Thetaeq (double Phi, double z)
{
  const double DTY = 60., DTZ = 10., T0 = 200., T1 = 315.;
  double E = PI(Ps(z));
  double Teq = (T1 - DTY*pow(sin(Phi),2) - DTZ*log(Ps(z)/P0)*pow(cos(Phi),2))*E/Cp;
  if (T0 > Teq) Teq = T0;
  return Cp*Teq/E;
}

/**
*Be careful that, from now on, $T$ designates the potential
temperature, not just the temperature.*

## Initial conditions */

const double maxlat = 75;

event init (i = 0)
{

  /**
  We can restart from a (dump) file if it exists. */
  
  if (restore ("restart"))
    event ("metric");
  else {

    /**
    Otherwise we set initial conditions. */
    
    const double DT0 = 0.001;
    foreach () {

      /**
      We bound the domain at $\pm$ maxlat using a "dry terrain". */
      
      zb[] = fabs(y) < maxlat ? 0. : 10*P0;
      if (zb[] > 0.)
        foreach_layer()
          h[] = 0., T[] = 0.;
      else {

        /**
        We set the initial potential temperature $T$ to the equilibrium
        at 45 degrees plus some high-frequency modes of 10^-3^ degrees
        amplitude. */
      
        double z = zb[];
        foreach_layer() {
          h[] = P0/nl;
          z += h[]/2.;
          T[] = Thetaeq(45.*pi/180., z) + DT0*(sin(10.*pi*x/180.) + cos(10.*pi*y/180.));
          z += h[]/2.;
        }
      }
    }
  
    /**
    We reset the fields used to store various averages/diagnostics. */
    
    reset ({mu, mtheta, vtheta}, 0.);
  }
}

/**
## Radiative forcing

Following Held & Suarez, the input of energy into the system is a
"radiative forcing" implemented as a pressure-dependent relaxation
toward the equilibrium radiative temperature i.e.
$$
\partial_t T = \dots - k_T [T - \theta_{eq}]
$$
with
$$
k_t = k_a + (k_s - k_a)\max\left(0,\frac{\sigma - \sigma_b}{1 - \sigma_b}\right)\cos^4\phi
$$
*/

event radiative_forcing (i++)
{
  const double ka = 1./40.*1./day, ks = 1./4.*1./day;
  foreach()
    if (zb[] < P0) {
      double z = zb[];
      foreach_layer() {
        z += h[]/2.;
        double sigma = Ps(z)/P0;
        double kt = ka + (ks - ka)*max(0., (sigma - sigmab)/(1. - sigmab))*pow(cos(y*pi/180.),4.);

        /**
        We use a first-order implicit time integration scheme. */
        
        T[] = (T[] + dt*kt*Thetaeq(y*pi/180., z))/(1. + dt*kt);
        z += h[]/2.;
      }
    }
}

/**
## Rigid lid correction

We use a rigid lid, which can cause small variations in the total
column height/pressure. We ensure that this height remains constant
(at P0) by expanding/contracting it slightly if needed. */

event set_eta (i++)
{
  foreach()
    if (zb[] < P0) {
      double etap = zb[];
      foreach_layer()
        etap += h[];
      etap /= P0;
      foreach_layer() {
        h[] /= etap;
        for (scalar s in tracers)
          s[] *= etap;
      }
    }
}

/**
## Main setup */

int main()
{

  /**
  The domain spans 360 x 180 degrees (bounded by maxlat = $\pm$ 75
  degrees). The radius of the earth sets the length unit (meters). */
  
  dimensions (2, 1);
  size (360);
  origin (0, - 90);
  periodic (right);  
  Radius = 6371220. [1];

  /**
  We use 256 grid points i.e. $\approx 1.4$ degrees on the equator and
  20 (equidistant) layers. */
  
  N = 256;
  nl = 20;

  /**
  We use a rigid lid. */
  
  rigid = true;
  theta_H = 1.;
  mgH.minlevel = 3; // fixme: necessary for convergence but only on GPUs?

  /**
  The timestep needs to be small enough. Not sure what controls the
  stability... */
  
  DT = 150.[0,1] * 256./N;
  TOLERANCE = 0.1;

  /**
  The solution is sensitive to limiting in the vertical remapping,
  probably due to the changes this induces in the top and bottom
  boundary conditions. Using monotonic limiting leads to higher
  stratospheric temperature variance, as well as the appearance of an
  "equatorial soliton". Results closer to the reference numerical
  results are obtained without limiting. */
  
  //  cell_lim = mono_limit;
  
  run();

  /**
  We make shorter versions of the full movies. */

  system ("ffmpeg -i theta0.mp4 -ss 00:02:00 -to 00:02:20 -crf 28 -movflags +faststart "
          "theta0-short.mp4");
  system ("ffmpeg -i omega0.mp4 -ss 00:02:00 -to 00:02:20 -crf 28 -movflags +faststart "
          "omega0-short.mp4");
}

/**
## Simple diagnostics */

event logfile (i++)
{

  /**
  The rigid-lid pressure can drift, because it is defined to within a
  constant. We substract the mean to prevent this drift. */
  
  stats s = statsf (eta_r);
  double avg = s.sum/s.volume;
  foreach()
    eta_r[] -= avg;
  fprintf (stderr, "%g %g %g %g %g %g\n",
           t, dt,
           normf (T).max, normf(u.x).max,
           s.min, s.max);
}

/**
## Time averages and variances

We allocate new fields to store the time-averaged velocities,
potential temperatures and their variance. */

event init (i = 0)
{
  mu = new scalar[nl];
  mtheta = new scalar[nl];
  vtheta = new scalar[nl];
}

/**
We compute the running mean and variance using [Welford/West weighted
incremental
algorithm](https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Weighted_incremental_algorithm). This
is important to avoid catastrophic round-off errors, particularly on
(32-bits) GPUs. */

const double tspinup = 200*day;

event average (t = tspinup; i++)
{
  double w_sum_old = t - tspinup, w_sum = w_sum_old + dt;
  foreach()
    foreach_layer() {

      /**
      The mean and variance of the potential temperature, updated
      using Welford/West's algorithm. */
         
      double mean_old = mtheta[];
      mtheta[] = mean_old + (dt/w_sum)*(T[] - mean_old);
      vtheta[] = (w_sum_old*vtheta[] + dt*(T[] - mean_old)*(T[] - mtheta[]))/w_sum;
      
      /**
      The mean zonal velocity. */
      
      mu[] = (w_sum_old*mu[] + dt*u.x[])/w_sum;
    }
}

/**
## Movies

We make movies of the averaged potential temperature and standard deviation. */

#define BOX {{X0, max(Y0, - maxlat)}, {X0 + L0, min(Y0 + L0/dimensions().x, maxlat)}}

event average_outputs (t = tspinup; t += day)
{
  output_ppm (vtheta, file = "vtheta0.mp4", spread = -1,
              n = 1024, box = BOX, linear = true);
  output_ppm (mtheta, file = "mtheta0.mp4", min = 264, max = 313,
              n = 1024, box = BOX, linear = true);
}

/**
More movies */

event movies (t += 9000)
{
  output_ppm (T, file = "theta0.mp4", min = 264, max = 313,
              n = 1024, box = BOX, linear = true);
  output_ppm (eta_r, file = "eta_r.mp4", min = -2000, max = 3000,
              n = 1024, box = BOX, linear = true);
  scalar omega[];
  vorticity (u, omega);
  output_ppm (omega, file = "omega0.mp4", min = -1e-4, max = 1e-4, map = cool_warm,
              n = 1024, box = BOX, linear = true);
}

#if _GPU && SHOW
event display (i++)
{
  output_ppm (T, fp = NULL, fps = 30,
              min = 264, max = 313, n = 1024, box = BOX, map = cool_warm, linear = false);
}
#endif

/**
## Snapshots
   
We dump regular snasphots which will be used for
[post-processing](held-suarez-postproc.c). */

event snapshot (t += 10*day) {
  dump();
}

/**
We end after 1200 days, thus computing the averages above for 1000
days, as in Held & Suarez. */

event ending (t = 1200*day);

/**
## Performances

This example runs fine on GPUs, using e.g.

~~~bash
CFLAGS=-DSHOW make held-suarez.gpu.tst
~~~

On an [RTX
4090](https://www.techpowerup.com/gpu-specs/geforce-rtx-4090.c3889)
runtime is approximately 44 minutes for a resolution of 256 x 128 (1.4
degrees), which is also about 5.8 nanoseconds per degree of freedom, or
again 107 years per day (ypd). The time per degree of freedom decreases
to 2.8 nanoseconds for 512 x 256 ($\approx$ 28 ypd) and to 2
nanoseconds for 1024 x 512 ($\approx$ 4.9 ypd).

For `N = 256` profiling gives the following:

~~~bash
   calls    total     self   % total   function
      19     0.05     0.04     30.2%   advect():/src/layered/hydro.h:404
      57     0.02     0.02     16.1%   multigrid_restriction():/src/grid/gpu/../multigrid-common.h:459
     665     0.02     0.02     15.7%   relax_hydro():/src/layered/implicit.h:104
    5092     0.01     0.01     10.1%   setup_shader():/src/grid/gpu/grid.h:1956
      19     0.01     0.01      9.2%   vertical_remapping():/src/layered/remap.h:165
      38     0.04     0.00      2.8%   mg_cycle():/src/poisson.h:92
      38     0.00     0.00      2.8%   gpu_reduction():/src/utils.h:143
      19     0.00     0.00      1.4%   gpu_reduction():/src/utils.h:171
      19     0.00     0.00      1.4%   face_fields():/src/layered/hydro.h:291
      19     0.01     0.00      1.2%   acceleration_0():/src/layered/implicit.h:206
      19     0.00     0.00      1.2%   acceleration_2():/src/layered/dr.h:76
      57     0.00     0.00      1.1%   residual_hydro():/src/layered/implicit.h:125
      19     0.05     0.00      1.0%   pressure():/src/layered/hydro.h:461  
~~~

For `N = 1024`:

~~~bash
   calls    total     self   % total   function
      19     0.14     0.14     29.1%   vertical_remapping():/src/layered/remap.h:165
      19     0.08     0.07     14.7%   advect():/src/layered/hydro.h:404
      57     0.07     0.06     13.1%   multigrid_restriction():/src/grid/gpu/../multigrid-common.h:459
    1015     0.04     0.03      7.0%   relax_hydro():/src/layered/implicit.h:104
      19     0.05     0.03      5.7%   acceleration_0():/src/layered/implicit.h:206
      19     0.10     0.02      4.8%   pressure():/src/layered/hydro.h:461
      19     0.02     0.02      4.1%   face_fields():/src/layered/hydro.h:291
      19     0.11     0.02      3.8%   pressure_0():/src/layered/implicit.h:251
      19     0.02     0.02      3.8%   acceleration_2():/src/layered/dr.h:76
    5986     0.02     0.02      3.2%   setup_shader():/src/grid/gpu/grid.h:1956
      19     0.01     0.01      2.4%   acceleration_1():/src/layered/coriolis.h:64
      19     0.01     0.01      2.4%   set_eta():held-suarez.gpu.c:263
~~~

## Todo

* Cubed sphere global coordinates would be good.
* Generalisation to a "wet" atmosphere.
* The link between surface pressure $p_s$ in [Marshall et al,
  2004](#marshall2004) and the rigid-lid pressure `eta_r` is not
  entirely clear.
* Topography can be included but will probably require a better
  understanding of the point above (link between $p_s$, geopotential
  and `eta_r`).
* The rigid lid assumption is not necessary. It would be interesting
  to study what changes when it is replaced by an implicit free
  surface.
* Using `z` as name for vertical coordinate is confusing.
* The dimension of the vertical coordinate should not be a length, as
  imposed by [/src/layered/hydro.h](/src/layered/hydro.h#175).
* Vertical remapping has been optimised but still takes a substantial
  fraction of time (twice as much as horizontal advection for `N = 1024`).

## See also

* [Post-processing for the Held-Suarez example](held-suarez-postproc.c)
  
## References

~~~bib
@article {held1994,
 author = "Isaac M. Held and Max J. Suarez",
 title = "A Proposal for the Intercomparison of the Dynamical Cores of
          Atmospheric General Circulation Models",
 journal = "Bulletin of the American Meteorological Society",
 year = "1994",
 publisher = "American Meteorological Society",
 address = "Boston MA, USA",
 volume = "75",
 number = "10",
 pages=      "1825 - 1830",
 pdf = "https://journals.ametsoc.org/downloadpdf/view/journals/bams/75/10/1520-0477_1994_075_1825_apftio_2_0_co_2.pdf"
}

@article{marshall2004,
  title={Atmosphere--ocean modeling exploiting fluid isomorphisms},
  author={Marshall, John and Adcroft, Alistair and Campin, Jean-Michel and Hill, Chris and White, Andy},
  journal={Monthly Weather Review},
  volume={132},
  number={12},
  pages={2882--2894},
  year={2004},
  doi={10.1175/MWR2835.1},
  pdf={https://journals.ametsoc.org/view/journals/mwre/132/12/mwr2835.1.pdf}
}

@article{held2005,
  title={The gap between simulation and understanding in climate modeling},
  author={Held, Isaac M},
  journal={Bulletin of the American Meteorological Society},
  volume={86},
  number={11},
  pages={1609--1614},
  year={2005},
  publisher={American Meteorological Society}
}
~~~
*/
