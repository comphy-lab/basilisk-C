/**
# Capillary waves on a thin liquid jet

The motion and axisymmetric perturbations riding on a thin straight
liquid jet can conveniently be described within the frame of a
one-dimensional approximation. Disregarding viscous effects, this
approximation reads:
$$
\frac{\partial h}{\partial t} + v \frac{\partial h}{\partial z}  = 
-\frac{1}{2} h \frac{\partial v}{\partial z} 
$$
$$
\frac{\partial v}{\partial t} + v \frac{\partial v}{\partial z} = 
- \frac{\gamma}{\rho} \frac{\partial p}{\partial z} 
$$ 
where $h(z,t)$ and $v(z,t)$ stand respectively for the radius of the
jet and the local velocity in a thin slice.  The ``pressure'' $p$
actually only amounts to the capillary contribution of pressure, and
is fully determined by the geometry of the jet:
$$
p(h) = \frac{1}{h \left(1 + h_z^2\right)^{1/2}} - \frac{h_{zz}}{\left(1+h_z^2\right)^{3/2}}.
$$
Note that here, the full expression for the curvature is used (see e.g
[Eggers and Dupont, J. Fluid Mech., 262
(1994)](http://dx.doi.org/10.1017/S0022112094000480) for a
discussion.)

## Conservative form

This set of equations can be rewritten in the following conservative form:
$$
\frac{\partial h^2}{\partial t} + \frac{\partial h^2 v}{\partial z}  = 0
$$
$$
\frac{\partial h^2 v}{\partial t} + \frac{\partial h^2 v^2}{\partial z} = 
\frac{\gamma}{\rho} \frac{\partial}{\partial z}\left(h^2 K\right) 
$$
with $K$ denoting the following quantity:
$$
K = \frac{h_{zz}}{\left(1+h_z^2\right)^{3/2}} + \frac{1}{h \left(1 + h_z^2\right)^{1/2}}.
$$
(yes, this is now a `+' sign -- see [Li and Fontelos, Phys. Fluids,
15(4) (2003)](http://dx.doi.org/10.1063/1.1556291).)

This formulation reveals the mass per unit length $\rho h^2$ and the
momentum per unit length $\rho h^2 v$ of a thin slice as natural
variables of the problem. 

## Plateau-Rayleigh instability

We will now solve the governing equations in the unstable regime, extract
the instability growth rate and compare it with the theoretical prediction.

We use a 1D Cartesian grid and the generic loop. */

#include "grid/cartesian1D.h"
#include "run.h"

/**
## Variables

We define the conservative variables as a *scalar* mass per unit length
and a *face vector* momentum per unit length. */

scalar h2[];
face vector h2u[];

double dt, wavenumber, growth_rate, pulsation;
int counter=0;
FILE * fp;

/**
## Boundary conditions
*/

h2u.n[left]  = dirichlet(0);
h2[left]     = neumann(0);
h2u.n[right] = dirichlet(0);
h2[right]    = neumann(0);

/**
## Scanning the unstable region

We explore the unstable region corresponding to a dimensionless
wavenumber $k R < 1$. We set the numerical domain to two wavelengths.
*/

int main() {
  for (wavenumber = 0.02; wavenumber <= 1.; wavenumber += 0.02) { 
    counter++;
    growth_rate = sqrt(0.5 * (pow(wavenumber,2.) - pow(wavenumber,4.)));
    L0 = 2. * 2. * pi / wavenumber; // domain = two wavelengths
    N = 100;
    DT = 5e-5;
    run();
  }
}

/**
## Size of the perturbation

The size of the initial disturbance $\epsilon$ is set to:
*/

#define EPS 0.001

/**
## Initialization

At initial time a placeholder for the monitoring of the amplitude
is created. The structure of the initial disturbance matches that
of un unstable mode (both in shape and velocity structure).
*/
event init (t = 0) {

  char name[80];
  sprintf (name, "amp-%i.dat", counter);
  fp = fopen(name,"w");
  fprintf(fp, "Wavenumber %g \n", wavenumber);
  
  foreach()
    h2[] = sq(1. - EPS*EPS / 4. + EPS*cos(wavenumber*x)); // Note the O(EPS^2) term to ensure mass conservation
  boundary ({h2});
  if (growth_rate >= 0) {
    foreach_face()
      h2u.x[] = -2. * growth_rate * EPS / wavenumber * sin(wavenumber*x) * h2[];
  }
  boundary_flux ({h2u});
}

/**
## End of the simulation
*/
event end (t = 5) {

}

/**
## Monitoring the growth

Every 100 timesteps, the relative amplitude of the disturbance is monitored in a file.
*/
event growth (i += 100) {
  stats s = statsf (h2);
  fprintf (fp,"%d %g %g %g\n", i, t, 
	  (1.-sqrt(s.min)) / (EPS), 
	  (sqrt(s.max)-1.) / (EPS));
}

/**
## Explicit resolution

The integration is fully explicit. Note that the size of timestep is adjusted by 
computing the CFL and the distance to the next event. The $K$ term is computed by composing
stencil operations, and the convective flux is approximated in the Harlow-Welch style.*/
event integration (i++) {
  double dt = DT;
  scalar u[];
  scalar h2uflux[];
  foreach_face() {
    u[] = 2.*h2u.x[]/(h2[] + h2[-1]);
    h2uflux[] = sq(h2u.x[1] + h2u.x[])/(4.*h2[]);
    double un = fabs(u[]);
    if (un > 0. && CFL*Delta/un < dt)
      dt = CFL*Delta/un;
  }
  dt = dtnext (t, dt);

  boundary({h2uflux});

  scalar h[];
  foreach()
    h[] = sqrt(h2[]);
  boundary ({h});

  scalar K[];
  foreach() {
    double hz = (h[1] - h[-1])/(2.*Delta);
    double hzz = (h[1] + h[-1] - 2.*h[])/sq(Delta);
    double ds = sqrt(1. + sq(hz));
    K[] = hzz/pow(ds,3) + 1./(h[]*ds);
  }
  boundary ({K});

  foreach_face()
    h2u.x[] += - dt*(h2uflux[] - h2uflux[-1])/Delta +
    dt*(h2[]*K[] - h2[-1]*K[-1])/Delta;

  scalar dh2 = K;
  foreach()
    dh2[] = - (h2u.x[1] - h2u.x[])/Delta;
  foreach()
    h2[] += dt*dh2[];

  boundary ({h2});
  boundary_flux ({h2u});
}
/**
Then compile and run.

The postprocessing of the output files can be done with this gnuplot script

~~~bash
reset
set terminal pngcairo size 1024,768 enhanced font 'Arial,18'
set xlabel "Wavenumber"
set ylabel "Growth rate"

do for [i=1:50] {
  stats 'amp-'.i.'.dat' every ::::0 using 2 nooutput
  k = STATS_min
  f(x) = a + b*x
  fit f(x) 'amp-'.i.'.dat' u 2:(log($4)) every ::1 via a,b
  set object circle at first k,b radius char 0.5 \
  fillstyle empty border lc rgb '#aa1100' lw 2
}
set samples 1000
set xrange [0:1.2]
plot real(sqrt(0.5*(x**2-x**4))) t 'Theory'
~~~ 

![rp](/sandbox/liquid_jet/pics/rp.png)

*/