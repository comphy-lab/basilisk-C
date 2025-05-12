/**
# A warm axi-symmetric blob

This page examplifies the axi-symmetric flow solver and adds gravity in
[Boussinesq's
limit](https://en.wikipedia.org/wiki/Boussinesq_approximation_(buoyancy)). 

## Axis-symmetry

Nature can be laden with symmetries. An example is a perfect vortex
ring when using a proper [cylindircal coordinate
system](https://en.wikipedia.org/wiki/Cylindrical_coordinate_system). A
3D Cartesian flow field $\mathbf{u}(x,y,z)$ can be transformed into
Cylindrical coordinates ($\mathbf{u}(r,\phi,z)$). For an axi-symmetric
flow there is no dependence on the azimuthal coordinate
$\frac{\partial\mathbf{u}}{\partial\phi}=0$. The thee-dimensional
equations for fluid flow can then be re-written in terms of only two
dimensions; the radial distance $r$ and the axial coordinate $z$. With
Basilisk's axi-symmetric flow solver, $\mathtt{x} = z$ and $\mathtt{y}
= r$. It is switched-on as follows:

Do not `#include "grid/octree.h"` 

but,
 */
#include "axi.h"
#include "navier-stokes/centered.h"
/**
## Buoyancy

In order to include the most prominent effects of density differences
in the atmospheric boundary layer due to temperature variations, it
makes sense to employ the aforementioned Boussinesq approximation. The
advantages include that we do not need to worry about the proper
thermodynamics on a solver level and the discription naturally
includes the acceleration of gravity ($g$).

If we consider small temperature variations in an ideal gas with
respect to an absolute reference temperature (e.g. $T_{ref} = 280K$
with variations of $10K$ above a canopy), we can introduce a buoyancy
parameter ($b$) to describe the
[body-force](https://en.wikipedia.org/wiki/Body_force) acting on a
fluid parcel of temperature $T$ due to gravity ($g$).

$$b = g\frac{T - T_{ref}}{T_{ref}}.$$

Note 1: The units of $b$ are units of acceleration ($ms^{-2}$).  
Note 2: $T$ is formally the [virtual potential
temperature](http://glossary.ametsoc.org/wiki/Virtual_potential_temperature).

As such, we define a buoyancy field, which will be advected and
diffused. Further, the acceleration if coupled to the flow solver via
an acceleration face-vector field.
 */
#include "tracer.h"
#include "diffusion.h"

scalar b[], * tracers = {b};

face vector av[];

/**
## Numerical setup

We consider a fluid with a diffusivity of momentum ($\nu$) equal to
the diffusivity of heat $\kappa$ (i.e. $\frac{\nu}{\kappa} = 1$). In
which, a warm *spherical* blob of radius $R$ is placed. The atmospheric
buoyancy profile corresponds to a linear stable stratification:

$$b(z) = N^2 z,$$

with $N$ the [Brunt-Vaisala
frequency](https://en.wikipedia.org/wiki/Brunt%E2%80%93V%C3%A4is%C3%A4l%C3%A4_frequency). The
warm spherical blob of size $R$ is placed at $z = 0$ with $b_{blob} =
b_0$.

### Quick Dimensional analysis

From the remaining parameters ($\nu, R, N, b_o$ (not $\kappa$, $T_{ref}$ nor $g$)) we can
construct two dimensionless numbers:

$$Re = \frac{b_oR}{N\nu},$$

$$\Pi = R\frac{N^2}{b_0}.$$

Here we set $Re = 1000$ and $\Pi = 2$.
 */
double R = 1, b0 = 1;    // Normalized values;
double Re = 1000, Pi = 0.2; // Dimensionless numbers
double nu, Nbv;            // Dependent variables

b[left]  = dirichlet (sq(Nbv)*x);
b[right] = dirichlet (sq(Nbv)*x);

int maxlevel = 9;

int main() {
  L0  = 10*R;
  Nbv = sqrt(Pi*b0/R);
  nu  = b0*R/(Nbv*Re); 
  X0  = -L0/3;  // Axial coordinate
  a   = av;     // Link Gravity
  const face vector muc[] = {nu, nu};
  mu = muc;
  run();
}

#define RADIUS (sqrt(sq(x) + sq(y)))
event init (t = 0) {
  refine (RADIUS > R - R/5  && RADIUS < R + R/5  && level < maxlevel - 1); 
  refine (RADIUS > R - R/20 && RADIUS < R + R/20 && level < maxlevel);
  foreach() {
    b[] = sq(Nbv)*x;
    if (RADIUS < R)
      b[] = b0;
  }
  boundary ({b});
  DT = 0.1;
}
/**
## Gravity

Gravity can now be added by overloading the `acceleration` event
following the syntax dictated by the flow solver.
*/

event acceleration (i++) {
  boundary ({b});
  foreach_face(x) // x is the axial coordinate "z"
    av.x[] = face_value(b,0);
}
/**
## Advection and diffusion of the $b$ field

The advection is automatically taken care of by `tracer.h`, for
diffusion we overload the `tracer_diffusion` event. 
 */

event tracer_diffusion (i++) {
  diffusion (b, dt, mu);
}

/**
## Output

We output movies of the usual suspects

![The "Vertical" velocity field](axiblob/ux.mp4)

![The Buoyancy field](axiblob/b.mp4)

![The level of refinement](axiblob/l.mp4)

The blob rises to a level where it is neutrally buoyant and flows a
bit outward.

Note that there is no reason to assume an axi-symetric evolution for
this case. 
 */

event movies (t += 0.1) {
  scalar lev[];
  foreach()
    lev[] = level;
  output_ppm (u.x, file = "ux.mp4", n = 300, min = -1, max = 1);
  output_ppm (lev, file = "l.mp4", n = 300, max = maxlevel);
  output_ppm (b, file = "b.mp4", n = 300, min = -2*b0, max = 2*b0);
}

event adapt (i++) {
  adapt_wavelet ({b, u}, (double[]){b0/20, 0.01, 0.01}, maxlevel);
}

event stop (t = 20);

/**
## To do

* Remove the warm blob and the related code.  
* Add a vortex ring injector at a suitable boundary.  
* Construct dimenless-number group for the new case. 
* Add code for quantitative analysis.  
* Do the relevant analysis ;)
*/
