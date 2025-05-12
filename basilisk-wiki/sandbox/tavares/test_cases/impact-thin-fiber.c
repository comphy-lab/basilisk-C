
/** 
# Droplet impact  on a thin fiber 

 We study the behavior of o a droplet impacting a thin fiber following the experimental study of 
 [Lorenceau and. al. 2004](#Lorenceau_2004) and a related numerical application [Wang and. al. 2018](#Wang_2018).  
 */

#include "grid/multigrid3D.h"
//#include "grid/octree.h" 
#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "../contact-embed.h"
#include "view.h"

/** The density ratio is 1000 and the dynamic viscosity ratio 50 */
#define RHOR 1000
#define MUR 50
#define T 4.
#define Rd 0.0009 //Physical radius of the droplet
#define Rf 0.000375*1.1 //Physical radius of the fiber
#define RATIO Rf/Rd // Ratio between the two physical radius
#define Df Dd*RATIO // Diameter of the fiber
#define Yi 0.9
#define shift 0.035*RATIO
#define yc 0.5

/** We define the Reynolds *Re_d*, the Weber numbers *We_d* and the Froude *Fr_d* based on the droplet diameter */

#if CASE2 // Second case (the droplet detaches from the fiber)
  const double Red=830; 
  const double Wed=37;
  const double Frd=10;
#else // First case (the droplet is captured by the fiber)
  const double Red=175;
  const double Wed=1.6;
  const double Frd=2;
#endif

/** We choose a unit length for the diameter of the droplet */
const double Dd=1.;
/** We consider an hydrophilic sphere */
double theta0=15; //instead of $10^\circ$ in the reference papers

int MAXLEVEL=5;

int main() {

  /**
  The domain is $4^3$ */
  size (4);
  origin (-L0/2, -L0/2, -L0/2);
  N = 1 << MAXLEVEL;
  init_grid (N);

  /**
  We set the properties  giving the dimensionless number */
  rho1 = 1. [0];  // 
  rho2 = rho1/RHOR;
  mu1 = 1./Red;
  mu2 = 1./(MUR*Red);
  f.sigma = 1./Wed;

 // /*
 //  We reduce the tolerance on the Poisson and viscous solvers to
 //  improve the accuracy. */
 //  TOLERANCE = 1e-4 [*];;


  /** We set the contact_angle value. */
  const scalar c[] = theta0*pi/180.;
  contact_angle = c;
  run();
}


event init (t = 0) {
  /**
  We define the unit solid sphere. */
  solid (cs, fs, sq(x) + sq(y-yc) - sq(Df/2.));
  fractions_cleanup (cs, fs);

  #if ADAPT 
    refine (sq(x) + sq(y) + sq(z) < sq(1.1*Dd/2.) && level < MAXLEVEL);
  #endif

  /** and the droplet */
  fraction (f, - (sq(x-shift) + sq(y - (yc + Yi)) + sq(z) - sq(Dd/2.)) );

  foreach()
    u.y[] = -1*f[];
}

/**
No-slip boundary condition is applied at the top and the bottom wall.
*/

u.t[bottom] = dirichlet(0);
u.t[top]  = neumann(0);

/**
We make sure there is no flow through the left and right boundary*/
uf.n[right] = 0.;
uf.n[left] = 0.;

/**
 No slip boundary condition on the fiber
 */
// u.n[embed] = dirichlet(0.);
// u.t[embed] = dirichlet(0.);
// u.r[embed] = dirichlet(0.);

/**
We add the acceleration of gravity. */

event acceleration (i++) {
  face vector av = a;
  foreach_face(y)
    av.y[] -= 1/sq(Frd);
}

/**
We save the vertical velocity field statistics. */

 event logfile (i++; t<=T) {
  norm n = normf(u.y);
  fprintf(stdout, "%g %g %g %g\n", t, dt, n.avg, n.max);

 }


#if 0
event cleanfghost (i++, t<=T, last) { // clean fluid ghost cells
  foreach()
    if (cs[]==0) f[]=0.;
}
#endif


event movie (i+=10; t <= T) {
view (quat = {-0.058, 0.184, 0.021, 0.981},
      fov = 30, near = 0.01, far = 1000,
      tx = -0.021, ty = 0.094, tz = -3.041,
      width = 1140, height = 678, bg = {1,1,1});
  draw_vof (c = "cs", s = "fs", filled = -1, fc = {0.83,0.83,0.83});
  draw_vof (c = "f", fc = {0.647,0.114,0.176}, lw = 2);
  save ("movie.mp4");
}

#if TREE
event adapt (i++) {
fractions_cleanup(cs,fs,smin = 1.e-10);
#if 1
  scalar f1[];
  foreach()
    f1[] = f[];
  adapt_wavelet ({f1,cs}, (double[]){1e-3,1e-2}, minlevel = 3, maxlevel = MAXLEVEL);
#else
  adapt_wavelet ({f}, (double[]){1e-4}, minlevel = 3, maxlevel = MAXLEVEL);
#endif
}
#endif

/**

![Droplet impact on a thin fiber, the droplet is captured by the fiber in the first case.](impact-thin-fiber/movie.mp4)
![Droplet impact on a thin fiber, the droplet detaches from the fiber in the second case.](impact-thin-fiber2/movie.mp4)


## References

~~~bib
@article{Wang_2018,
  author    = {Wang, Sheng and Desjardins, Olivier},
  journal   = {Chemical Engineering Science},
  title     = {Numerical study of the critical drop size on a thin horizontal fiber: Effect of fiber shape and contact angle},
  year      = {2018},
  issn      = {0009-2509},
  month     = sep,
  pages     = {127--133},
  volume    = {187},
  doi       = {DOI: 10.1016/j.ces.2018.04.040},
  publisher = {Elsevier BV},
}

@Article{Lorenceau_2004,
  author    = {Lorenceau, Élise and Clanet, Christophe and Quéré, David},
  journal   = {Journal of Colloid and Interface Science},
  title     = {Capturing drops with a thin fiber},
  year      = {2004},
  issn      = {0021-9797},
  month     = nov,
  number    = {1},
  pages     = {192--197},
  volume    = {279},
  doi       = {10.1016/j.jcis.2004.06.054},
  publisher = {Elsevier BV},
}

~~~
*/








