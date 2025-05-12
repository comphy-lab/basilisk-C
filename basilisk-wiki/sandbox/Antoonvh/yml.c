/**
![Yuk Man Lau's setup](yml/mov.mp4)

# Particle through an interface

This setup was proposed by Yuk Man Lau: 

A cylindrical particle with radius $R$ is forced through a
fluid-fluid interface with a velocity $U$. The two imissible fluids
both have a viscosity $\mu$ and density $\rho$. There exists a surface
tension assoiated with the fluid-fluid interface of $\sigma$.

Given the *five* system parameters $R, U, \mu, \rho$ and $\sigma$,
with *three* base units of lenght, time and mass, we may identify two
dimensionless parameters that describe the system:

$$Re = \frac{\rho UR}{\mu},$$
$$Ca = \frac{\mu U}{\sigma}.$$

We choose $R=1$, $U=1$ and $\rho = 1$. Note that the simulation is run
in the frame-of-reference that co-moves with the particle.
 */
double Re = 10.;
double Ca = 1.;

#include "embed.h" 
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "view.h"
/**
In order to include the effect of surface tension the code under
`tension.h` is *mimiced*, except for the zero divisions...
 */
#include "iforce.h"
#include "curvature.h"
attribute {
  double sigma;
}
/**
We also set the in and out flow boundary conditions, 
 */
f[left] = 1.;
u.n[left] = dirichlet(1.);
u.n[right] = neumann(0.);
p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);

/**
and a no-slip condition on the cylinder.
 */
u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);
f[embed] = 0.;
/**
The domain size `L0`$=40R$. 
 */
face vector muc[];
int main(){
  L0 = 40.;
  X0 = -L0/2.;
  stokes = true; //Re is small
  mu = muc;
  f.sigma = 1./(Re*Ca);
  init_grid (256);
  run();
}
/**
The viscosity is defined in a special event such that it considers the
embedded boundary.
*/
event properties(i++){
  foreach_face()
    muc.x[] = fm.x[]/Re;
  boundary ((scalar*){muc});
}

event init (t = 0){
  CFL = 0.5;
  vertex scalar phi[];
  foreach_vertex()
    phi[] = sq(x + 5.) + sq(y) - sq(1.);
  fractions (phi, cs, fs);
  foreach()
    u.x[] = (cm[] > 0);
  DT = 0.01;
}
/**
   Grid adaptation is carried out.
*/

event adapt(i++)
  adapt_wavelet ({cs, f, u.x, u.y},
		 (double[]){0.005, 0.0005, 0.02, 0.02}, 9);

/**
A movie is generated to display the results.
*/
#if 1
event mov (t += 0.1){
  scalar omg[];
  vorticity(u, omg);
  view(fov = 15, width = 800, height = 800);
  draw_vof("f");
  draw_vof("f", filled = 1, fc = {0.7, 0.7, 0.7});
  draw_vof("cs", "fs");
  draw_vof("cs", "fs", filled = -1, fc = {0.3, 0.3, 0.3});
  squares("omg", min = -1, max = 1, map = cool_warm);
  mirror({0,1}){
    draw_vof("cs", "fs");
    draw_vof("cs", "fs", filled = -1, fc = {0.3, 0.3, 0.3});
    cells();
  }
  save("mov.mp4");
}
#endif

event the_end (t = 35);

/**
   The surface tension event:
 */
event acceleration (i++)
{
  for (scalar f in interfaces)
    if (f.sigma) {
      /**
	 If $\phi$ is already allocated, we add $\sigma\kappa$,
	 otherwise we allocate a new field and set it to
	 $\sigma\kappa$. */
      scalar phi = f.phi;
      if (phi.i)
	curvature (f, phi, f.sigma, add = true);
      else {
	phi = new scalar;
	curvature (f, phi, f.sigma, add = false);
	f.phi = phi;
      }
    }
}
