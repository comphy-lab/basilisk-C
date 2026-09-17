/**
# Interfacial power

Rate at which the interfacial force does work on the flow,

$$
\Psi_\sigma = \int_\Omega \mathbf{u}\cdot\phi\nabla f \, dV,
\qquad \phi = \sigma\kappa
$$

In the continuum $\Psi_\sigma = -d(\sigma A)/dt$. It measures what the discrete
force did, so it closes against the kinetic energy budget; as an estimate of
interfacial area it is one to two orders worse than `interface_area()` (see
`readme`).

Include **after** `tension.h`, which pulls in the unguarded `curvature.h`.
Needs `navier-stokes/centered.h` for `a`, `alpha`, `fm` and `u`.

## The two routes

Both integrate $\mathbf{u}\cdot\phi\nabla f$ and agree to roundoff on a uniform
grid
([test_interfacial_power_routes.c](tests_quantities/test_interfacial_power_routes.c)).
On trees `iforce.h` swaps `f`'s prolongation to the pressure's before
differencing and route 2 does not, so they differ wherever the interface sits
on a refinement boundary -- a divergence signals that a refined band no longer
contains the interface.

Use `interfacial_power_curvature()` by default. One value is independent of the
event it is read from; an accumulated $\int\Psi_\sigma\,dt$ is not, changing by
about a factor of two if sampled after `centered.h`'s `correction()`.

`interfacial_power_acceleration()` is what the solver applied, under three
constraints:

1. Call it from a `projection` event, the first hook where `a` is complete and
   `u` is still the velocity the force acted on. Same-named events run in
   reverse declaration order and `,last` does not help; check `qcc -events`.
2. `ag` is subtracted from `a`. Impossible under `REDUCED`, which folds
   buoyancy into the same $\phi$.
3. Not from `end_timestep`: `correction()` has already advanced `u`.
*/

/**
## Route 1 -- read back the applied acceleration

`alpha = fm/rho` on faces, so `a*fm/alpha` is $\phi\nabla f$. `ag` is the
uniform body acceleration already in `a` (`{0,0,-1./At}` for gravity along
$-z$, `{0}` if none), removed with the same face weights: with `sigma = 0`
this returns machine zero.

`foreach_dimension()` does not rotate the members of a plain `coord`, hence
the per-component subtraction.
*/

double interfacial_power_acceleration (coord ag)
{
  face vector av = a;
  double psi = 0.;
  foreach(reduction(+:psi)) {
    double d = 0.;
    foreach_dimension()
      d += u.x[]*(av.x[]*fm.x[]/alpha.x[] + av.x[1]*fm.x[1]/alpha.x[1])/2.;

    d -= u.x[]*ag.x*(fm.x[]/alpha.x[] + fm.x[1]/alpha.x[1])/2.;
    d -= u.y[]*ag.y*(fm.y[]/alpha.y[] + fm.y[0,1]/alpha.y[0,1])/2.;
#if dimension == 3
    d -= u.z[]*ag.z*(fm.z[]/alpha.z[] + fm.z[0,0,1]/alpha.z[0,0,1])/2.;
#endif
    psi += dv()*d;
  }
  return psi;
}

/**
## Route 2 -- recompute the potential

Rebuilds $\phi = \sigma\kappa$ and forms $\phi\nabla f$ with `iforce.h`'s
stencil, including its `nodata` fallbacks (both neighbours defined: average;
one: take it; neither: zero). Costs one `curvature()` call.
*/

/**
The force itself, $\mathbf{F}_\sigma = \sigma\kappa\nabla f$, per cell --
`interfacial_power_curvature()`'s potential $\phi\nabla f$ kept as a vector
rather than immediately dotted with $\mathbf{u}$, so it can be correlated
against a field other than the solver's own `u` (a spectral flux term, say)
instead of only integrated against it.
*/

void surface_tension_force (scalar c, vector Fs)
{
  scalar phi[];
  curvature (c, phi, c.sigma, add = false);
  foreach()
    foreach_dimension() {
      double pl = (phi[]   < nodata && phi[-1] < nodata) ? (phi[] + phi[-1])/2. :
                   phi[]   < nodata ? phi[]   :
                   phi[-1] < nodata ? phi[-1] : 0.;
      double pr = (phi[1]  < nodata && phi[]   < nodata) ? (phi[1] + phi[])/2. :
                   phi[1]  < nodata ? phi[1]  :
                   phi[]   < nodata ? phi[]   : 0.;
      Fs.x[] = (pl*(c[] - c[-1]) + pr*(c[1] - c[]))/(2.*Delta);
    }
}

/**
$\Psi_\sigma$ and its density, as `u`.$\mathbf{F}_\sigma$ built from
`surface_tension_force()`. `d` returns the per-cell density, the return value
its volume integral.
*/

double interfacial_power_curvature (scalar c, scalar d)
{
  vector Fs[];
  surface_tension_force (c, Fs);
  double psi = 0.;
  foreach(reduction(+:psi)) {
    double di = 0.;
    foreach_dimension()
      di += u.x[]*Fs.x[];
    d[] = di;
    psi += dv()*di;
  }
  return psi;
}
