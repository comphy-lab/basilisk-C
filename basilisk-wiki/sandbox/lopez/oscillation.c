/**
# Oscillating walls in a two-layer sandwich of fluid and solid. 

This is the first test case reported in 
[Sugiyama et al. 2011](#sugiyama2011) who
proposed a numerical scheme for fluid-structure coupled problems
involving interaction of (incompressible) hyperlastic materials with
(incompressible) fluids by modifying slightly available classic
numerical schemes for incompressible fluid dynamics.*/

#include "navier-stokes/centered.h"
#include "src/hyperelastic.h"
#include "vof.h"

#define VW 1    //Amplitude of the wall velocity
#define OM M_PI //Frequency of the oscillation
#define C1 2.5  //Solid is a Neo-Hookean of G = 5
#define MU_F 1  //Fluid viscosity
#define LS 0.5  //Width of the solid layer
#define LF 0.5  //Width of the fluid layer 

face vector visc[];
scalar beta1v[];
scalar f[];
scalar * interfaces = {f};

u.t[top] = dirichlet (VW*sin(OM*t));

/** 
We will consider the upper half of the sandwich so the bottom is at
rest.*/

u.t[bottom] = dirichlet (0.);

int MAXLEVEL = 6;

int main() {
  L0 = LS + LF ;
  origin (-L0/2, 0);
  periodic (right);
  DT = .005;
  mu = visc;
  beta1 = beta1v;
  stokes = true;
  init_grid (1 << MAXLEVEL);
  //  vector v = B.x;
  // v.n[bottom] = dirichlet(0.);
  run();
}

event init (i = 0) {
  vertex scalar psi[];

  /**
  If a planar interface coincides with cell faces, i.e. the VOF scalar
  jumps from 0 to 1 in adjacent cells across the interfaces,
  *fractions()* complains. Displacing a little bit the interface avoid
  the problem. */

  foreach_vertex ()
    psi[] = -y + (LS + 0.0001); 
  fractions (psi, f);

   /** 
   We follow the idea of Sugiyama el al (2011). In their work, it is
   defined $\hat{\mathbf{B}} = \phi^q \mathbf{B}$ being $\phi$ the solid
   volume fraction and $q$ an elegible exponent. Note that
   $\hat{\mathbf{B}}$ is, in contrast to $\mathbf{B}$, defined in all the
   computational domain. The equation governing its temporal evolution is
   still UCD($\hat{\mathbf{B}}$) = 0. */

  foreach() {
    u.x[] = 0.;
    foreach_dimension()
      B.x.x[] = sqrt(clamp(f[],0, 1.));
  }
}

event properties (i++) {
  foreach_face() {
    double T = (f[]+f[-1])/2.;
    visc.x[] = MU_F*(1-T);
  }

  /**
   As in Sugiyama et al (2011) we will assume that the hyperelastic
   stresses in half-full cells will be proportional to the volume
   fraction of solid. So $\beta_1$ will be proportional to
   $\sqrt{f}$. */

  foreach()
    beta1v[] = 2.*C1*sqrt(clamp(f[],0,1.));
}

/**
In an imcompressible solid, the determinant of the Cauchy left
deformation is an invariant, $(det(\mathbf{B}))_t = 0$. */

event logfile (i++) {
  scalar invariant[];
  foreach () 
    invariant[] = B.x.x[]*B.y.y[]-B.y.x[]*B.x.y[];

  printf("t = %g sum = %g\n", t, statsf(invariant).sum);
}

/**
Velocity profiles at some particular instants are saved. */


event vprofile (t = {2, 2.8, 4, 4.8}) {
  char name[80];
  sprintf (name, "vprof-%g", t);
  FILE * fp = fopen (name, "w");
  foreach () 
    fprintf(fp,"%g %g \n", y, u.x[]);
  fclose(fp);
}

#if 0
#if QUADTREE
event gfsview (i += 4) {
  static FILE * fpg = popen("gfsview2D -s vista.gfv","w");
  output_gfs (fpg, t= t);
}
#endif
#endif

/**
~~~gnuplot Velocity field distribution at instant $t=0$ and $t=0.8$. 
set xlabel 'y'
set ylabel 'v_x'
unset key
plot 'oscillation.sugiyama' u 1:2 w l, 'oscillation.sugiyama' u 1:3 w l, \
     'vprof-2' u 1:2, 'vprof-4' u 1:2,					\
     'vprof-2.8' u 1:2, 'vprof-4.8' u 1:2
~~~

## References

~~~bib
@article{sugiyama2011,
  title={A full Eulerian finite difference approach for solving fluid-structure coupling problems},
  author={Sugiyama, Kazuyasu and Ii, Satoshi and Takeuchi, Shintaro and Takagi, Shu and Matsumoto, Yoichiro},
  journal={Journal of Computational Physics},
  volume={230},
  number={3},
  pages={596-627},
  year={2011},
  publisher={Elsevier},
  pdf={https://arxiv.org/pdf/1009.3609}
}
~~~
*/
