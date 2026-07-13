/**
# Mesh-dependent dynamic contact angle

A static contact angle combined with a no-slip condition is singular at the
contact line. In VOF simulations this singularity is regularized implicitly,
since the interface is advected by the face velocity, located half a grid cell
away from the wall; but the resulting effective slip length -- and therefore the
apparent contact angle -- then depends on the mesh resolution.

## Usage

~~~literatec
#include "contact.h"
#include "dynamic.h"

vector h[];
h.t[left] = contact_angle (thetad*pi/180.);

int main() {
  ...
  HALF_GRID_SIZE = L0/(1 << maxlevel);
  HALF_COARSE_SIZE = L0/(1 << minlevel);
  SLIP_LENGTH = 0.;
  LNUM = HALF_GRID_SIZE, LAPP = 1e-9;
  dynamic_angle = dynamic_polynomial;
  ...
}
~~~

This header adjusts the contact angle dynamically as a function of the local
capillary number, following [Afkhami et al. (2009)](#afkhami2009) and [Legendre
and Maglio (2015)](#legendre2015), so that the apparent angle converges as the
mesh is refined. See [dynamic.c](/sandbox/ecipriano/run/dynamic.c) for the
derivation of the two laws, the Stat1-3/Dyn1-4 table of configurations, and
validation cases. */

extern double theta0, mu1;
extern scalar f;

double thetad = 90.;          // dynamic contact angle [deg]
double DYN_THETA_MIN = 2.;    // lower clamp for thetad [deg]
double DYN_THETA_MAX = 198.;  // upper clamp for thetad [deg]
double DYN_RELAX = 0.02;      // relaxation factor for the update of thetad
double SLIP_LENGTH = 0.;      // slip length
double HALF_GRID_SIZE = 0.;   // Delta/2 at the current refinement level
double HALF_COARSE_SIZE = 0.; // Delta/2 at the coarsest refinement level
double LNUM = 0.;             // numerical contact angle length scale
double LAPP = 0.;             // equilibrium/apparent contact angle scale

/**
The clamp keeps `thetad` away from the 0/180 degrees singularities of the laws
below, while still allowing some overshoot past the physical range so that the
relaxation in the `vof` event (see below) is not clipped every step during
transients. */

typedef enum {
  STAT1 = 1, // static angle, no slip
  STAT2,     // static angle, slip length = Delta/2
  STAT3,     // static angle, slip length = coarsest Delta/2
  DYN1,      // dynamic angle (polynomial law), L = 1e-6, lambda = 1e-9
  DYN2,      // dynamic angle (polynomial law), L = Delta/2, lambda = 1e-9
  DYN3,      // dynamic angle (polynomial law), L = Delta/2, lambda = 1e-9, with slip
  DYN4       // dynamic angle (Afkhami law), L = Delta/2, lambda = 0.04*R
} DYN_SETUP;

DYN_SETUP dyn_setup = STAT1;

/**
## Implementation

Three dynamic angle laws are provided, corresponding to the expressions derived
in [dynamic.c](/sandbox/ecipriano/run/dynamic.c): `dynamic_polynomial` and
`dynamic_legendre` both implement the $\theta_{d,1}$ law of [Legendre and Maglio
(2015)](#legendre2015) -- the former using the full polynomial fit of $g$ and
its inverse, the latter using the closed-form cubic approximation valid for
$\theta_d < 3\pi/4$ -- while `dynamic_afkhami` implements the $\theta_{d,2}$ law
of [Afkhami et al. (2009)](#afkhami2009).

A simulation selects a law by pointing `dynamic_angle` at one of the functions
below (or leaves it `NULL` for a static angle equal to `theta0`). A user law
with the same signature can be assigned in the same way. */

static inline double gfun (double x) {
  x = (x < 0) ? 0 : x;
  return pow (x, 3.)/9. - 0.00183985*pow (x, 4.5) +
    1.845823*1e-6*pow (x, 12.258487);
}

static inline double ginv (double x) {
  x = (x < 0) ? 0 : x;
  return cbrt (9.*x) + 0.0727387*x -
    0.0515388*sq (x) + 0.00341336*pow (x, 3.);
}

double dynamic_polynomial (double Ca) {
  double thetad = ginv (gfun (theta0*pi/180.) + Ca*log (LNUM/LAPP))*180./pi;
  return clamp (thetad, DYN_THETA_MIN, DYN_THETA_MAX);
}

double dynamic_legendre (double Ca) {
  double thetad3 = pow (theta0*pi/180., 3.) + 9.*Ca*log (LNUM/LAPP);
  double thetad = cbrt (thetad3)*180./pi;
  return clamp (thetad, DYN_THETA_MIN, DYN_THETA_MAX);
}

double dynamic_afkhami (double Ca) {
  double costhetad = cos (theta0*pi/180.) - 5.63*Ca*log (LNUM/LAPP);
  double thetad = acos (clamp (costhetad, -1, 1))*180./pi;
  return clamp (thetad, DYN_THETA_MIN, DYN_THETA_MAX);
}

double (* dynamic_angle) (double Ca) = NULL;

/**
## Legendre & Maglio's table

Each function below configures `dyn_setup`, `SLIP_LENGTH` and (for the dynamic
cases) `dynamic_angle`, `LNUM` and `LAPP` to match one row of the Stat1-3/Dyn1-4
table reproduced in [dynamic.c](/sandbox/ecipriano/run/dynamic.c). `HALF_GRID_SIZE`/`HALF_COARSE_SIZE` must be
set beforehand, since several rows are expressed relative to them. */

void setup_stat1 (void) {
  dyn_setup = STAT1;
  SLIP_LENGTH = 0.;
  dynamic_angle = NULL;
}

void setup_stat2 (void) {
  dyn_setup = STAT2;
  SLIP_LENGTH = HALF_GRID_SIZE;
  dynamic_angle = NULL;
}

void setup_stat3 (void) {
  dyn_setup = STAT3;
  SLIP_LENGTH = HALF_COARSE_SIZE;
  dynamic_angle = NULL;
}

void setup_dyn1 (void) {
  dyn_setup = DYN1;
  SLIP_LENGTH = 0.;
  LNUM = 1e-6, LAPP = 1e-9;
  dynamic_angle = dynamic_polynomial;
}

void setup_dyn2 (void) {
  dyn_setup = DYN2;
  SLIP_LENGTH = 0.;
  LNUM = HALF_GRID_SIZE, LAPP = 1e-9;
  dynamic_angle = dynamic_polynomial;
}

void setup_dyn3 (void) {
  dyn_setup = DYN3;
  SLIP_LENGTH = HALF_GRID_SIZE;
  LNUM = HALF_GRID_SIZE, LAPP = 1e-9;
  dynamic_angle = dynamic_polynomial;
}

void setup_dyn4 (void) {
  dyn_setup = DYN4;
  SLIP_LENGTH = 0.;
  LNUM = HALF_GRID_SIZE, LAPP = 0.02;
  dynamic_angle = dynamic_afkhami;
}

/**
## Interface position

Best estimate of the interface position from the height functions, which is
more accurate than the PLIC (`mycs`) reconstruction used as a fallback where
the height is undefined. Reusable for post-processing the contact-line
location. */

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
## Capillary number at the contact line

Computes $Ca_{cl} = \mu u_{cl} / \sigma$, where $u_{cl}$ is the tangential
velocity at the contact line, sampled half a cell away from the wall (the
"numerical slip" location which advects the VOF interface). The contact line
is detected with `is_contact_x()` from [contact-small.h](/src/contact-small.h).

The $y$ coordinate of the sampling point depends on the interface
orientation. At moderate angles the interface crosses the wall-tangential
columns close to the contact line and `h.y` locates it accurately. At small
(or large) angles the interface becomes nearly parallel to the wall: `h.y` is
undefined -- or worse, defined but crossing the interface
$\cot(\theta_d)/2$ cells away from the contact line -- so we revert to the
intersection between the PLIC reconstruction and the wall face. The selection
is based on the dominant component of the interface normal, consistently with
the validity region of the height functions.

`sign` orients the velocity regardless of which side of the contact line
($y+1$ or $y-1$) the liquid phase occupies. The dynamic angle laws are quite
sensitive to how this velocity is evaluated.

fixme: only works on the left boundary; generalizing `is_contact_x` and the
face lookup to an arbitrary boundary is left for future work. */

double capillary (void) {
  if (f.height.x.i) {
    vector h = f.height;
    double ucl = 0., sign = 0.;
    foreach (serial)
      if (is_contact_x (point, f, theta0*pi/180.)) {
        coord m = interface_normal (point, f);
        if (fabs (m.y) > fabs (m.x) && h.y[] != nodata)
          ucl = interpolate (u.y, x, y + Delta*height (h.y[]));
        else {
          double alpha = plane_alpha (f[], m);
          double ycl = fabs (m.y) > 1e-10 ?
            clamp ((alpha + 0.5*m.x)/m.y, -0.5, 0.5) : 0.;
          ucl = interpolate (u.y, x, y + Delta*ycl);
        }
        sign = (f[0,1] < f[0,-1]) ? 1. : -1.;
      }
    return sign*ucl*mu1/f.sigma;
  }
  else
    return 0.;
}

/**
## Update

The dynamic angle is refreshed after VOF advection. Updating it every step
directly from the law produces large oscillations at high resolution, so we
relax towards the target value; `relax = 1` disables this. */

event init (i = 0) {
  thetad = theta0;
}

event vof (i++) {
  if (dynamic_angle) {
    double relax = DYN_RELAX;
    double Ca = capillary();
    thetad = relax*dynamic_angle (Ca) + thetad*(1. - relax);
  }
}

/**
## References

~~~bib
@article{afkhami2009,
  title={A mesh-dependent model for applying dynamic contact angles to VOF simulations},
  author={Afkhami, Shahriar and Zaleski, Stephane and Bussmann, Markus},
  journal={Journal of computational physics},
  volume={228}, number={15}, pages={5370--5389}, year={2009},
  publisher={Elsevier},
  pdf={https://web.njit.edu/~shahriar/download/AZB_JCP09.pdf}
}

@hal{legendre2015, hal-01340390}
~~~
*/
