/**
# Mesh dependent dynamic contact angle

Moving contact lines are often simulated by combining a static contact angle
model with no-slip boundary conditions on the solid surface. However, this
combination leads to a singularity in the viscous stresses at the contact line,
which is typically regularized by introducing a small slip length (although
other regularization mechanisms exist).

In VOF simulations, the interface is advected by the face velocity, located half
a grid cell away from the solid boundary. This introduces an implicit
(numerical) slip length that removes the contact line singularity. However, this
also results in mesh-dependent solutions.

In this file, we implement the models proposed by [Afkhami et al.
(2009)](afkhami2009) and [Legendre and Maglio (2015)](#legendre2015). These
models exploit a [Cox–Voinov](https://youtu.be/K3gcJXDVjZU?si=K-dlcYATp7g-xa3l)-type
law to dynamically adjust the contact angle based on the local capillary number.

The two expressions for the dynamic contact angle considered here are:

$$
g(\theta_{d,1}) = g(\theta_s) + Ca_{cl} \log \left(\dfrac{L}{\lambda}\right)
\quad \text{with} \quad
g(\theta) = \int_0^\theta \dfrac{x - \sin{x}\cos{x}}{2\sin{x}} \, dx
$$

which simplifies to ([Legendre and Maglio, 2015](#legendre2015)):

$$
\theta_{d,1}^3 = \theta_s^3 + Ca_{cl} \log\left(\dfrac{L}{\lambda}\right)
$$

and ([Afkhami et al., 2018](#afkhami2009)):

$$
\cos{\theta_{d,2}} = \cos{\theta_{s}} + 5.63 \, Ca_{cl} \log\left(\dfrac{L}{\lambda}\right)
$$

where the capillary number at the contact line is defined as
$Ca_{cl} = \mu u_{cl} / \sigma$.

The inner and outer length scales are defined differently in the two models.
They are summarized in the following table (adapted from [Legendre and Maglio,
2015](#legendre2015)).

| Name  | Contact angle | Outer length $L$ | Inner length $\lambda$ | Slip length $\lambda_N$ |
|-------|---------------------------|------------------|------------------|-------------------|
| Stat1 | $\theta_d = \theta_s$     | —                | —                | $0$               |
| Stat2 | $\theta_d = \theta_s$     | —                | —                | $\Delta/2$        |
| Stat3 | $\theta_d = \theta_s$     | —                | —                | $\Delta_{32}/2$   |
| Dyn1  | $\theta_d = \theta_{d2}$  | $10^{-6}$        | $10^{-9}$        | $0$               |
| Dyn2  | $\theta_d = \theta_{d2}$  | $\Delta/2$       | $10^{-9}$        | $0$               |
| Dyn3  | $\theta_d = \theta_{d2}$  | $\Delta/2$       | $10^{-9}$        | $\Delta/2$        |
| Dyn4  | $\theta_d = \theta_{d1}$  | $K$              | $\Delta/2$       | $\Delta/2$        |

## Withdrawing Plate

The first configuration used to test the convergence issue consists in a square
cavity with an imposed left wall velocity, which transports the interface
upwards from a flat configuration. We consider the impact of different contact
angle models on the time evolution of the interface height.

The combination of static contact angle and no-slip condition does not converge
in the contact line height (and therefore the interface shape).

~~~gnuplot (Stat1) Interface shape
set size square
set xr[0:1]
set yr[0:1]
set grid
set arrow from -0.1,0.4 to -0.1,0.6 lw 1.5
set label "U_0" at -0.15,0.5

plot "facets-5" u 1:2 w l lw 1.2 lc -1 dt 4 t "LEVEL 5", \
     "facets-6" u 1:2 w l lw 1.2 lc -1 dt 3 t "LEVEL 6", \
     "facets-7" u 1:2 w l lw 1.2 lc -1 dt 2 t "LEVEL 7", \
     "facets-8" u 1:2 w l lw 1.2 lc -1 dt 1 t "LEVEL 8", \
     "facets-9" u 1:2 w l lw 1.2 lc  7 dt 1 t "LEVEL 9"
set term pop
~~~

~~~gnuplot (Stat1) Contact line height
reset
set yr[0.4:0.8]
set xlabel "tau [-]"
set ylabel "Contact line height"
set grid
set key top left

plot "<grep 'level 5' log" u 3:4 w lp lw 1.2 lc -1 dt 4 t "LEVEL 5", \
     "<grep 'level 6' log" u 3:4 w lp lw 1.2 lc -1 dt 3 t "LEVEL 6", \
     "<grep 'level 7' log" u 3:4 w lp lw 1.2 lc -1 dt 2 t "LEVEL 7", \
     "<grep 'level 8' log" u 3:4 w lp lw 1.2 lc -1 dt 1 t "LEVEL 8", \
     "<grep 'level 9' log" u 3:4 w lp lw 1.2 lc  7 dt 1 t "LEVEL 9"
~~~

Introducing a Navier slip boundary condition with a slip length smaller than
half the grid size (i.e. under-resolved) is effectively equivalent to imposing a
no-slip condition in terms of contact line convergence. However, the presence of
slip still affects the dynamics, slowing down the interface evolution at all
refinement levels compared to the no-slip case.

~~~gnuplot (Stat2) Interface shape
reset
set xr[0:1]
set yr[0:1]
set size square
set grid
set arrow from -0.1,0.4 to -0.1,0.6 lw 1.5
set label "U_0" at -0.15,0.5

plot "../plate-stat2/facets-5" u 1:2 w l lw 1.2 lc -1 dt 4 t "LEVEL 5", \
     "../plate-stat2/facets-6" u 1:2 w l lw 1.2 lc -1 dt 3 t "LEVEL 6", \
     "../plate-stat2/facets-7" u 1:2 w l lw 1.2 lc -1 dt 2 t "LEVEL 7", \
     "../plate-stat2/facets-8" u 1:2 w l lw 1.2 lc -1 dt 1 t "LEVEL 8", \
     "../plate-stat2/facets-9" u 1:2 w l lw 1.2 lc  7 dt 1 t "LEVEL 9"
set term pop
~~~

~~~gnuplot (Stat2) Contact line height
reset
set yr[0.4:0.8]
set xlabel "tau [-]"
set ylabel "Contact line height"
set grid
set key top left

plot "<grep 'level 5' ../plate-stat2/log" u 3:4 w lp lw 1.2 lc -1 dt 4 t "LEVEL 5", \
     "<grep 'level 6' ../plate-stat2/log" u 3:4 w lp lw 1.2 lc -1 dt 3 t "LEVEL 6", \
     "<grep 'level 7' ../plate-stat2/log" u 3:4 w lp lw 1.2 lc -1 dt 2 t "LEVEL 7", \
     "<grep 'level 8' ../plate-stat2/log" u 3:4 w lp lw 1.2 lc -1 dt 1 t "LEVEL 8", \
     "<grep 'level 9' ../plate-stat2/log" u 3:4 w lp lw 1.2 lc  7 dt 1 t "LEVEL 9"
~~~

We can achieve convergence on the interface position by using a slip length
which is larger than the grid size. In this case, it is equal to half the size
of the coarsest grid. However, if the viable grid refinements are much larger
than the physical slip length, this method is not adequate since it introduces
too much slip. For example, in the following plots the interface position does
converge, but at a much lower elevation compared to the no-slip case.

~~~gnuplot (Stat3) Interface shape
reset
set xr[0:1]
set yr[0:1]
set size square
set grid
set arrow from -0.1,0.4 to -0.1,0.6 lw 1.5
set label "U_0" at -0.15,0.5

plot "../plate-stat3/facets-5" u 1:2 w l lw 1.2 lc -1 dt 4 t "LEVEL 5", \
     "../plate-stat3/facets-6" u 1:2 w l lw 1.2 lc -1 dt 3 t "LEVEL 6", \
     "../plate-stat3/facets-7" u 1:2 w l lw 1.2 lc -1 dt 2 t "LEVEL 7", \
     "../plate-stat3/facets-8" u 1:2 w l lw 1.2 lc -1 dt 1 t "LEVEL 8", \
     "../plate-stat3/facets-9" u 1:2 w l lw 1.2 lc  7 dt 1 t "LEVEL 9"
set term pop
~~~

~~~gnuplot (Stat3) Contact line height
reset
set yr[0.4:0.8]
set xlabel "tau [-]"
set ylabel "Contact line height"
set grid
set key top left

plot "<grep 'level 5' ../plate-stat3/log" u 3:4 w lp lw 1.2 lc -1 dt 4 t "LEVEL 5", \
     "<grep 'level 6' ../plate-stat3/log" u 3:4 w lp lw 1.2 lc -1 dt 3 t "LEVEL 6", \
     "<grep 'level 7' ../plate-stat3/log" u 3:4 w lp lw 1.2 lc -1 dt 2 t "LEVEL 7", \
     "<grep 'level 8' ../plate-stat3/log" u 3:4 w lp lw 1.2 lc -1 dt 1 t "LEVEL 8", \
     "<grep 'level 9' ../plate-stat3/log" u 3:4 w lp lw 1.2 lc  7 dt 1 t "LEVEL 9"
~~~

This configuration was used by [Afkhami et al. (2009)](#afkhami2009mesh) to tune
the outer length scale of their dynamic contact angle model, which was then
applied to the following sessile droplet configuration.

## Sessile droplet

We can see that for each level of refinement, the droplet correctly relaxes towards a
steady configuration (the scheme is well-balanced). However, the spreading dynamics
changes depending on the resolution: as we increase the refinement we tend toward the
no-slip condition, which justifies the slower evolution of the higher refinements.

~~~gnuplot (Stat1) Interface shape
reset
set xr[0:1]
set yr[0:1]
set size square
set grid

plot "../drop-stat1/facets-5" u 1:2 w l lw 1.2 lc -1 dt 4 t "LEVEL 5", \
     "../drop-stat1/facets-6" u 1:2 w l lw 1.2 lc -1 dt 3 t "LEVEL 6", \
     "../drop-stat1/facets-7" u 1:2 w l lw 1.2 lc -1 dt 2 t "LEVEL 7", \
     "../drop-stat1/facets-8" u 1:2 w l lw 1.2 lc -1 dt 1 t "LEVEL 8", \
     "../drop-stat1/facets-9" u 1:2 w l lw 1.2 lc  7 dt 1 t "LEVEL 9"
set term pop
~~~

~~~gnuplot (Stat1) Normalized radius
reset
set xlabel "tau [-]"
set ylabel "Normalized radius"
set grid
set key top left

plot "<grep 'level 5' ../drop-stat1/log" u 3:4 w lp lw 1.2 lc -1 dt 4 t "LEVEL 5", \
     "<grep 'level 6' ../drop-stat1/log" u 3:4 w lp lw 1.2 lc -1 dt 3 t "LEVEL 6", \
     "<grep 'level 7' ../drop-stat1/log" u 3:4 w lp lw 1.2 lc -1 dt 2 t "LEVEL 7", \
     "<grep 'level 8' ../drop-stat1/log" u 3:4 w lp lw 1.2 lc -1 dt 1 t "LEVEL 8", \
     "<grep 'level 9' ../drop-stat1/log" u 3:4 w lp lw 1.2 lc  7 dt 1 t "LEVEL 9"
~~~

~~~gnuplot (Stat1) Capillary number
reset
set xlabel "tau [-]"
set ylabel "Capillary number [-]"
set grid
set key top right

plot "<grep 'level 5' ../drop-stat1/log" u 3:5 w lp lw 1.2 lc -1 dt 4 t "LEVEL 5", \
     "<grep 'level 6' ../drop-stat1/log" u 3:5 w lp lw 1.2 lc -1 dt 3 t "LEVEL 6", \
     "<grep 'level 7' ../drop-stat1/log" u 3:5 w lp lw 1.2 lc -1 dt 2 t "LEVEL 7", \
     "<grep 'level 8' ../drop-stat1/log" u 3:5 w lp lw 1.2 lc -1 dt 1 t "LEVEL 8", \
     "<grep 'level 9' ../drop-stat1/log" u 3:5 w lp lw 1.2 lc  7 dt 1 t "LEVEL 9"
~~~

~~~gnuplot (Stat1) Dynamic angle
reset
set xlabel "tau [-]"
set ylabel "Dynamic angle [deg]"
set grid
set key top left

plot "<grep 'level 5' ../drop-stat1/log" u 3:6 w lp lw 1.2 lc -1 dt 4 t "LEVEL 5", \
     "<grep 'level 6' ../drop-stat1/log" u 3:6 w lp lw 1.2 lc -1 dt 3 t "LEVEL 6", \
     "<grep 'level 7' ../drop-stat1/log" u 3:6 w lp lw 1.2 lc -1 dt 2 t "LEVEL 7", \
     "<grep 'level 8' ../drop-stat1/log" u 3:6 w lp lw 1.2 lc -1 dt 1 t "LEVEL 8", \
     "<grep 'level 9' ../drop-stat1/log" u 3:6 w lp lw 1.2 lc  7 dt 1 t "LEVEL 9"
~~~

Introducing the model by [Legendre and Maglio (2015)](#legendre2015) we can
enforce the convergence of the spreading dynamics. The model is really sensitive
to the way the capillary number and the contact line is computed, and to how the
interface position is reconstructed.

~~~gnuplot (Dyn4) Interface shape
reset
set xr[0:1]
set yr[0:1]
set size square
set grid

plot "../drop-dyn4/facets-5" u 1:2 w l lw 1.2 lc -1 dt 4 t "LEVEL 5", \
     "../drop-dyn4/facets-6" u 1:2 w l lw 1.2 lc -1 dt 3 t "LEVEL 6", \
     "../drop-dyn4/facets-7" u 1:2 w l lw 1.2 lc -1 dt 2 t "LEVEL 7", \
     "../drop-dyn4/facets-8" u 1:2 w l lw 1.2 lc -1 dt 1 t "LEVEL 8", \
     "../drop-dyn4/facets-9" u 1:2 w l lw 1.2 lc  7 dt 1 t "LEVEL 9"
set term pop
~~~

~~~gnuplot (Dyn4) Normalized radius
reset
#set xr[0:0.6]
#set yr[0.95:1.45]
set xlabel "tau [-]"
set ylabel "Normalized radius"
set grid
set key top left

plot "<grep 'level 5' ../drop-dyn4/log" u 3:4 w lp lw 1.2 lc -1 dt 4 t "LEVEL 5", \
     "<grep 'level 6' ../drop-dyn4/log" u 3:4 w lp lw 1.2 lc -1 dt 3 t "LEVEL 6", \
     "<grep 'level 7' ../drop-dyn4/log" u 3:4 w lp lw 1.2 lc -1 dt 2 t "LEVEL 7", \
     "<grep 'level 8' ../drop-dyn4/log" u 3:4 w lp lw 1.2 lc -1 dt 1 t "LEVEL 8", \
     "<grep 'level 9' ../drop-dyn4/log" u 3:4 w lp lw 1.2 lc  7 dt 1 t "LEVEL 9"
~~~

~~~gnuplot (Dyn4) Capillary number
reset
set xlabel "tau [-]"
set ylabel "Capillary number [-]"
set grid
set key top right

plot "<grep 'level 5' ../drop-dyn4/log" u 3:5 w lp lw 1.2 lc -1 dt 4 t "LEVEL 5", \
     "<grep 'level 6' ../drop-dyn4/log" u 3:5 w lp lw 1.2 lc -1 dt 3 t "LEVEL 6", \
     "<grep 'level 7' ../drop-dyn4/log" u 3:5 w lp lw 1.2 lc -1 dt 2 t "LEVEL 7", \
     "<grep 'level 8' ../drop-dyn4/log" u 3:5 w lp lw 1.2 lc -1 dt 1 t "LEVEL 8", \
     "<grep 'level 9' ../drop-dyn4/log" u 3:5 w lp lw 1.2 lc  7 dt 1 t "LEVEL 9"
~~~

~~~gnuplot (Dyn4) Dynamic angle
reset
set xlabel "tau [-]"
set ylabel "Dynamic angle [deg]"
set grid
set key top right

plot "<grep 'level 5' ../drop-dyn4/log" u 3:6 w lp lw 1.2 lc -1 dt 4 t "LEVEL 5", \
     "<grep 'level 6' ../drop-dyn4/log" u 3:6 w lp lw 1.2 lc -1 dt 3 t "LEVEL 6", \
     "<grep 'level 7' ../drop-dyn4/log" u 3:6 w lp lw 1.2 lc -1 dt 2 t "LEVEL 7", \
     "<grep 'level 8' ../drop-dyn4/log" u 3:6 w lp lw 1.2 lc -1 dt 1 t "LEVEL 8", \
     "<grep 'level 9' ../drop-dyn4/log" u 3:6 w lp lw 1.2 lc  7 dt 1 t "LEVEL 9"
~~~

## Conclusions

Among the dynamic models tested, only the formulation by [Afkhami et al.
(2018)](#afkhami2009) produces converging results, while the models by [Legendre
and Maglio (2015)](#legendre2015) do not match the expected behavior.

Interestingly, the latter study reports the opposite conclusion, suggesting that
*it may be dependent on the method used to solve the system of equations*. This
is consistent with the present results, as the simulations in [Afkhami et al.
(2018)](#afkhami2009) were performed using Gerris.

It is also important to note that a straightforward implementation of these
models does not guarantee correct behavior. In particular:

1. The interface position should be computed using height functions to ensure
   accurate contact line dynamics.
2. At high levels of refinement, small oscillations in the contact line velocity
   can induce large variations in the dynamic contact angle, potentially
   jeopardizing the simulation. To mitigate this issue, we introduce a
   relaxation factor that limits the variation of the contact angle between
   consecutive time steps.
*/

/**
### Compiler macros

We can run all the cases in the aforementioned papers by varying the following
parameters at compile-time. */

#ifndef DYNAMIC_ANGLE
# define DYNAMIC_ANGLE 0
#endif

#ifndef SLIP_LENGTH
# define SLIP_LENGTH 0
#endif

#ifndef THETA_LEGENDRE
# define THETA_LEGENDRE 0
#endif

#ifndef LOUT
# define LOUT 0
#endif

#ifndef LIN
# define LIN 0
#endif

#ifndef MAXLEVEL
# define MAXLEVEL 8
#endif

/**
### Simulation setup

We resolve the Navier--Stokes equations for a two-phase system with surface
tension and contact angle. Gravity is introduced only for the withdrawing plate
case. */

//#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "contact.h"
#include "two-phase.h"
#include "tension.h"
#include "reduced.h"

#if SESSILE
const double theta0 = 60.;
double thetad = theta0;
double u0 = 0., slip = 0.;
#else
const double theta0 = 90.;
double thetad = theta0;
double u0 = 1., slip = 0.;
#endif

vector h[];
h.t[left] = contact_angle (thetad*pi/180.);

double HALF_GRID_SIZE = 0., HALF_COARSE_SIZE = 0.;
u.n[left] = dirichlet (0.);
u.t[left] = navier (u0, SLIP_LENGTH);
p[left] = neumann (0.);

int maxlevel;
double R0 = 0.5, tend = 1.4;

int main (void) {
  rho1 = rho2 = 1.;
  mu1 = mu2 = 0.25;

  f.sigma = 7.5;
  f.height = h;

#if !SESSILE
  G.y = -9.81;
#endif

  for (maxlevel = 5; maxlevel <= MAXLEVEL; maxlevel++) {
    HALF_COARSE_SIZE = 0.5*L0/(1 << 5);
    HALF_GRID_SIZE = 0.5*L0/(1 << maxlevel);

    init_grid (1 << maxlevel);
    run();
  }
}

#define circle(x,y,R) (sq(R) - sq(x) - sq(y))

event init (i = 0) {
#if SESSILE
  fraction (f, circle (x, y, R0));
#else
  fraction (f, y - 0.4);
#endif
}

/**
### Dynamic angle implementation

The different dynamic contact angle models. */

#define THETA_MIN 30
#define THETA_MAX 150

#if THETA_LEGENDRE
static inline double gfun (double x) {
  x = (x < 0) ? 0 : x;
  return pow (x, 3.)/9. - 0.00183985*pow (x, 4.5) +
    1.845823*1e-6*pow (x, 12.258487);
}

static inline double ginv (double x) {
  x = (x < 0) ? 0 : x;
  return pow (9.*x, 1./3.) + 0.0727387*x -
    0.0515388*sq (x) + 0.00341336*pow (x, 3.);
}

double dynamic_angle (double Ca) {
  double thetad3 = pow (theta0*pi/180., 3.) - Ca*log (2.*LOUT/LIN);
  double thetad = pow (thetad3, 1./3.)*180./pi;
  return clamp (thetad, THETA_MIN, THETA_MAX);

  // I don't know why but the "general" model below is too unstable
  // Perhaps resolving the integral rather than using the polynomial is better
  //double thetad = ginv (gfun (theta0*pi/180.) - Ca*log (LOUT/LIN))*180./pi;
  //return clamp (thetad, THETA_MIN, THETA_MAX);
}
#elif THETA_AFKHAMI
double dynamic_angle (double Ca) {
  double costhetad = cos (theta0*pi/180.) + 5.63*Ca*log (LOUT/LIN);
  double thetad = acos (clamp (costhetad, -1, 1))*180./pi;
  return clamp (thetad, THETA_MIN, THETA_MAX);
}
#else
double dynamic_angle (double Ca) {
  return theta0;
}
#endif

/**
The following function checks if an interfacial cell along the
`left` boundary contains the contact line.

fixme: these functions should be made more general considering
any boundary orientation, not only the `left` boundary. */

static inline
int is_contact (Point point, scalar f) {
  if (f[] == 0. || f[] == 1.)
    return false;
  else {
    for (int i = -1; i <= 1; i += 2)
      if (theta0 < 90. && f[0,i] == 0. && is_boundary (neighbor(-1)))
        return true;
      else if (theta0 >= 90. && f[0,i] == 1. && is_boundary (neighbor(-1)))
        return true;
  }
  return false;
}

static inline
int is_adjacent (Point point, scalar f) {
  if (f[] == 0. || f[] == 1.)
    return false;
  else {
    for (int i = -1; i <= 1; i += 2)
      if (is_contact (neighborp(0,i), f))
        return true;
  }
  return false;
}

/**
The best estimation of the interface position is obtained using the height
functions. First we try to use the height function values at the face of the
cell. If the ghost value is undefined we use only the value in the internal
cell, if it is undefined we employ the PLIC reconstruction. */

double interface_position (Point point, scalar f) {
  if (f.height.y.i) {
    vector h = f.height;
    if (h.y[] != nodata && h.y[-1] != nodata)
      return y + Delta*0.5*(height (h.y[]) + height (h.y[-1]));
    else if (h.y[] != nodata)
      return y + Delta*height (h.y[]);
    else {
      coord m = interface_normal (point, f), p;
      double alpha = plane_alpha (f[], m);
      plane_area_center (m, alpha, &p);
      return y + Delta*p.y;
    }
  }
  return nodata;
}

/**
The capillary number is computed from the contact line velocity, which is
obtained as the tangential velocity interpolated on the x coordinate
corresponding to the cell face, and the y coordinate being the barycenter of the
interface segment. Note that this model is very sensitive to the way the contact
line velocity is computed. */

scalar iscontact[];

double capillary (void) {
  foreach() {
    if (is_contact (point, f))
      iscontact[] = 1.;
    else
      iscontact[] = 0.;
  }

  double ucl = 0.;
  foreach (serial)
    if (is_contact (point, f))
      ucl = interpolate (u.y, x, y + Delta*height (h.y[]));
  return ucl*mu1/f.sigma;
}

/**
If we change the contact angle based on the Capillary number at each time step,
the system shows large oscillations, especially at higher refinement levels,
which eventually propagate in an unrecoverable manner. To overcome this issue,
we introduce a relaxation parameter. This is not specified in the aforementioned
papers, but I noticed that it helps a lot. Try it yourself setting the
relaxation factor to 1. */

#if DYNAMIC_ANGLE
event vof (i++) {
  double relax = 0.02;
  double Ca = capillary();
  thetad = relax*dynamic_angle (Ca) + thetad*(1. - relax);
}
#endif

#if TREE
event adapt (i++) {
  adapt_wavelet ({f,u.x,u.y}, {1e-2,1e-2,1e-2}, maxlevel);
}
#endif

/**
### Post-processing

Using the height functions for computing the interface position gives more
accurate results compared to the mycs reconstruction. */

event logger (t += 0.05) {
  double xcl = HUGE;
  foreach_boundary (left, serial)
    if (is_contact (point, f))
      xcl = min (xcl, interface_position (point, f));
  double Ca = capillary();

#if SESSILE
  fprintf (stderr, "level %d %g %g %g %g\n", maxlevel,
      t, xcl/R0, Ca, thetad);
#else
  fprintf (stderr, "level %d %g %g %g %g\n", maxlevel,
      t, xcl, Ca, thetad);
#endif
}

#if MOVIE
#include "view.h"

event movie (t += 0.01) {
  clear();
  view (tx = -0.5, ty = -0.5);
  draw_vof ("f", lw = 2.);
  squares ("u.y", spread = -1);
  box();
  save ("movie.mp4");
}
#endif

event end (t = tend) {
  char name[80];
  sprintf (name, "facets-%d", maxlevel);
  FILE * fp = fopen (name, "w");
  output_facets (f, fp);
  fclose (fp);
}

/**
## References

~~~bib
@article{afkhami2009,
  title={A mesh-dependent model for applying dynamic contact angles to VOF simulations},
  author={Afkhami, Shahriar and Zaleski, Stephane and Bussmann, Markus},
  journal={Journal of computational physics},
  volume={228},
  number={15},
  pages={5370--5389},
  year={2009},
  publisher={Elsevier},
  pdf={https://web.njit.edu/~shahriar/download/AZB_JCP09.pdf}
}

@hal{legendre2015, hal-01340390}
~~~
*/
