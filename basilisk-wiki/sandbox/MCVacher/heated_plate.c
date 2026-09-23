/**
# Natural convection along a heated vertical plate

Numerical CFD simulation (Basilisk) of the natural convection boundary layer
along a vertical isothermal plate ($x = 0$), maintained at $T_0 > T_\infty$.

## Non-dimensionalisation
The Navier-Stokes equations are non-dimensionalised by the thermal diffusivity
$\alpha$ and the domain length $L_0$:

$$x^* = \frac{x}{L_0}, \quad t^* = \frac{\alpha \, t}{L_0^2}, \quad
\mathbf{u}^* = \frac{L_0}{\alpha}\,\mathbf{u}, \quad
\theta = \frac{T - T_\infty}{T_0 - T_\infty}$$

The non-dimensional momentum equation reads:

$$\frac{\partial \mathbf{u}^*}{\partial t^*} + (\mathbf{u}^* \cdot \nabla)\mathbf{u}^* = -\nabla p^* + Pr \, \nabla^2 \mathbf{u}^* + Ra \cdot Pr \; \theta \, \mathbf{e}_y$$

and the energy equation:

$$\frac{\partial \theta}{\partial t^*} + \mathbf{u}^* \cdot \nabla \theta = \nabla^2 \theta$$

with $Ra = \dfrac{g \beta (T_0 - T_\infty) L_0^3}{\nu \alpha}$ and $Pr = \dfrac{\nu}{\alpha}$.
The buoyancy term `Ra * Pr * θ` appears directly in the `acceleration` event.

## Simulation parameters
 *   Ra = 1e3    (Rayleigh number)
 *   Pr = 0.71   (air)
 */

#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"

#define Ra 1e3
#define Pr 0.71

double tmax     = 0.5;
double t_output = 0.005;

scalar T[];
scalar * tracers = {T};

mgstats mgT;

face vector D[];

int main()
{
  L0 = 1.;
  origin (0., 0.);

  N = 256;

  DT = 1e-4;
  TOLERANCE = 1e-4;
  NITERMIN = 4;

  // Left: hot wall
  u.n[left] = dirichlet(0.);
  u.t[left] = dirichlet(0.);

  // Right: open boundary
  u.n[right] = neumann(0.);
  u.t[right] = neumann(0.);

  // Bottom: open boundary
  u.n[bottom] = neumann(0.);
  u.t[bottom] = neumann(0.);

  // Top: open boundary
  u.n[top] = neumann(0.);
  u.t[top] = neumann(0.);

  // Temperature
  T[left]   = dirichlet(1.);
  T[right]  = dirichlet(0.);
  T[bottom] = dirichlet(0.);
  T[top]    = neumann(0.);

  run();
}

event init (i = 0)
{
  a = new face vector;
  mu = new face vector;

  face vector muc = mu;

  foreach_face() {
    muc.x[] = fm.x[] * Pr;
    D.x[] = fm.x[];
  }

  foreach() {
    T[] = 0;
    u.x[] = 0.;
    u.y[] = 0.;
  }
}


/**
Thermal diffusion
*/

event tracer_diffusion (i++)
{
  mgT = diffusion (T, dt, D);
}


/**
Buoyancy force
*/

event acceleration (i++)
{
  face vector av = a;

  foreach_face(y)
    av.y[] += Ra * Pr * (T[] + T[0,-1]) / 2.;
}

event logfile (i++)
{
  double umax = 0.;

  foreach()
    umax = max(umax, norm(u));

  fprintf(stderr,
          "i=%d t=%g umax=%g\n",
          i, t, umax);
}


/**
Save final solution
*/

event res_save (t = end)
{
  FILE * fp = fopen("res_final.txt", "w");

  foreach()
    fprintf(fp, "%g %g %g %g %g %g\n",
            x, y, u.x[], u.y[], p[], T[]);

  fclose(fp);
}


/**
Movies
*/

event movie (t = 0; t += t_output; t <= tmax)
{
  output_ppm(T,
             file = "T.mp4",
             n = 512,
             min = 0.,
             max = 1.,
             linear = true,
             map = jet);

  output_ppm(u.x,
             file = "ux.mp4",
             n = 512,
             spread = -1,
             linear = true,
             map = cool_warm);

  output_ppm(u.y,
             file = "uy.mp4",
             n = 512,
             spread = -1,
             linear = true,
             map = cool_warm);
}


event end (t = tmax)
{
  printf("-----END-----\n");
}

/**
# Outputs
![Temperature](heated_plate/T.mp4)
![Vertical velocity](heated_plate/uy.mp4)
![Horizontal velocity](heated_plate/ux.mp4)
*/

/**
# Theory: self-similar solution

The non-dimensional temperature is $\theta = \dfrac{T - T_\infty}{T_0 - T_\infty} \in [0,1]$.

The similarity variable is built from the thermal boundary layer scale
$\delta_T \sim y \, Ra_y^{-1/4}$:

$$\eta = \frac{x}{\delta_T} = \frac{x}{y} Ra_y^{1/4}, \qquad Ra_y = \frac{g \beta \Delta T \, y^3}{\nu \alpha}$$

The stream function is introduced as:

$$\psi = \alpha \, Ra_y^{1/4} \, F(\eta, Pr)$$

Non-dimensional equations ($Pr > 1$):

$$\boxed{\begin{cases} F''' + \dfrac{1}{Pr}\!\left(\dfrac{F'^2}{2} - \dfrac{3}{4} F F''\right) = \theta \\[8pt] \dfrac{3}{4} F \theta' = \theta'' \end{cases}}$$

Boundary conditions:

$$\eta = 0 : \quad F = 0, \quad F' = 0, \quad \theta = 1$$
$$\eta \to \infty : \quad F' \to 0, \quad \theta \to 0$$
*/