/**
# Impact of a hyperelastic ball on a solid

It has been pointed by Snoijer et al. (2020)
that an (hyper)elastic behavior of a material can obtained 
from the Olroyd-B model in the limit of
$$
\lambda \rightarrow \infty \quad \text{with} \quad \mu = \eta_p/\lambda \quad \text{finite}
$$
where $\lambda$ is the relaxation parameter, $\mu$ is analogue to the elastic shear modulus and $\eta_p$  is the polymeric viscosity.

We show that the Basilisk's viscoelastic solver can capture this behaviour by modifying [fall.c](/src/test/fall.c)
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "log-conform.h"
#include "curvature.h"


#define RHO_r 0.001
#define FR 2.26
#define LEVEL 7

#define LAMBDA 100000
#define MUE 1.7

scalar lambdav[], mupv[];

/**
The ball comes from the right. We allow the fluid to get through that
boundary. */

u.n[right] = neumann(0);
p[right]   = dirichlet(0);

/**
The wall is at the left side. We apply a no-slip boundary condition and a 
non-wetting condition for the VOF tracer. */

u.t[left] = dirichlet(0);
tau_qq[left] = dirichlet(0);
f[left]   = 0.;

int main() {

  /**
  The domain spans $[0:2.6]\times[0:2.6]$. */

  size (2.6);
  init_grid (1 << LEVEL);
  
  /**
  The densities and viscosities are defined by the parameters above. */

  rho1 = 1.;
  rho2 = RHO_r;
  mu1 = 0.01;
  mu2 = 0.001;

  /**
  The viscoelastic fields will be set below. */
  
  mup = mupv;
  lambda = lambdav;

  /**
  We set a maximum timestep. This is necessary for proper temporal
  resolution of the viscoelastic stresses. */
  
  DT = 1e-3;
  run();
}

event init (t = 0) {

  /**
  At a wall of normal $\mathbf{n}$ the component of the viscoelastic
  stress tensor $tau_p_{nn}$ is zero. Since the left boundary is a wall, we
  set $tau_p_{xx}$ equal to zero at that boundary. */
  
  scalar s = tau_p.x.x;
  s[left] = dirichlet(0.);

  /**
  The drop is centered on (2,0) and has a radius of 0.5. */

  fraction (f, - sq(x - 2.) - sq(y) + sq(0.5));

  /**
  The initial velocity of the droplet is -1. */

  foreach()
    u.x[] = - f[];
}

/**
We add the acceleration of gravity. */

event acceleration (i++) {
  face vector av = a;
  foreach_face(x)
    av.x[] -= 1./sq(FR);
}

/**
We update the viscoelastic properties according to the limit. */

event properties (i++) {
  foreach() {
    mupv[] = MUE*LAMBDA*clamp(f[],0,1);
    lambdav[] = LAMBDA*clamp(f[],0,1);
  }
}

/**
We adapt the solution at every timestep based on the interface and
velocity errors. */

#if TREE
event adapt (i++) {
  adapt_wavelet ({f, u.x, u.y}, (double[]){1e-2, 5e-3, 5e-3},
		 maxlevel = LEVEL, minlevel = LEVEL - 2);
}
#endif

/**
We track the spreading diameter of the droplet. */

event logfile (i += 20; t <= 5) {
  scalar pos[];
  position (f, pos, {0,1});
  fprintf (stderr, "%g %g\n", t, 2.*statsf(pos).max);
}

/**
We generate a movie of the interface shape. */

#include "view.h"

event viewing (i += 10) {
  view (width = 400, height = 400, fov = 20, ty = -0.5,
	quat = {0, 0, -0.707, 0.707});

  clear();
  draw_vof ("f", lw = 2);
  squares ("u.x", linear = true);
  box (notics = true);
  mirror ({0,1}) {
    draw_vof ("f", lw = 2);
    squares ("u.y", linear = true);
    box (notics = true);
  }
  save ("movie.mp4");
} 

/**
## Results

![Animation of the interface shape. The color field on the right-hand-side 
(resp. l.h.s.) is the radial (resp. axial) velocity component.](ballhyper/movie.mp4)

The ball rebounds like a rabbit! *Ni tan mal* for the viscoelastic model...!

~~~gnuplot Time evolution of the maximum diameter
reset
set ylabel 'Maximum diameter'
set xlabel 't'
plot 'log' w l lw 2 t 'Basilisk'
~~~

## References

~~~bib
@article{Snoijer2020,
  title={The relationship between viscoelasticity and elasticity},
  author={Snoeijer, J H and Pandey, A and Herrada, M A and Eggers, J},
  journal={Proc. R. Soc. A},
  volume={476},
  pages={20200419},
  year={2020}
}
~~~
*/
