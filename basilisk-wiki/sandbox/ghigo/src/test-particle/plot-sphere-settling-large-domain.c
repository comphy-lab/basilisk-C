/**
# Settling sphere in a large container */

#include "grid/octree.h"
#include "../myembed.h"
#include "../mycentered.h"
#include "../myembed-particle.h"
#include "lambda2.h"
#include "view.h"

/**
## Reference solution */

#define d    (2./12.)

/**
We also define the shape of the domain. */

#define sphere(x,y,z) (sq ((x)) + sq ((y)) + sq ((z)) - sq ((d)/2.))

void p_shape (scalar c, face vector f, coord p)
{
  vertex scalar phi[];
  foreach_vertex() {
    phi[] = HUGE;
    for (double xp = -(L0); xp <= (L0); xp += (L0))
      for (double yp = -(L0); yp <= (L0); yp += (L0))
	for (double zp = -(L0); zp <= (L0); zp += (L0))
	  phi[] = intersection (phi[],
				(sphere ((x + xp - p.x),
					 (y + yp - p.y),
					 (z + zp - p.z))));
  }
  boundary ({phi});
  fractions (phi, c, f);
  fractions_cleanup (c, f,
		     smin = 1.e-13, cmin = 1.e-13);
}

/**
We finally define the particle parameters. */

const double p_r = (2.56); // Ratio of solid and fluid density
const coord  p_i = {(p_moment_inertia_sphere ((d), 2.56)),
		    (p_moment_inertia_sphere ((d), 2.56)),
		    (p_moment_inertia_sphere ((d), 2.56))}; // Particle moment of interia
const double p_v = (p_volume_sphere ((d))); // Particle volume
const coord  p_g = {0., -(9.81), 0.};

int main ()
{
  for (int i = 1; i <= 2; i++) {
    for (int it = 20; it <= 40; it += 10) {

      /**
      The domain is $32^3$. It needs to be sufficiently big to allow
      for long settling times without the particle reaching the
      bottom. */

      L0 = 32.;
      size (L0);
      origin (-L0/2., 0., -L0/2.);

      char name[80];
      sprintf (name, "case-%d/lmax-13/dump-level-13-t-%d", i, it);
      bool restart = restore (name);
      assert (restart == true);
	
      sprintf (name, "case-%d/lmax-13/p_dump-level-13-t-%d", i, it);
      particle pp_restore;
      bool p_restart = p_restore (name, &pp_restore);
      assert (p_restart == true);

      /**
      We then need to initialize the face fraction *fs* since it is
      not dumped. */
  
      p_shape (cs, fs, p_p);

      /**
      Snapshots. */
  
      scalar l2[];
      lambda2 (u, l2);
      boundary ({l2});

      scalar omega[];
      vorticity (u, omega); // Vorticity in xy plane
      boundary ({omega});

      // phi
      double phi = pi/2. - fabs (acos ((p_u.z)/(sqrt(sq (p_u.x) + sq (p_u.y) + sq (p_u.z)))));

      // theta
      for (int itheta = 0; itheta <= 12; itheta++) {

	double theta = M_PI/6.*((double) itheta);

	double lambda2 = -0.05;
	  
	clear(); clear();
	view (fov = 7, theta = theta, relative = false,
	      tx = -(p_p.x)/(L0), ty = -(p_p.y + 25.*(d))/(L0),
	      bg = {1,1,1},
	      width = 400, height = 800);

	draw_vof ("cs", "fs");
	isosurface ("l2", lambda2, fc = {0.3,0.4,0.6}, lc = {0.3,0.4,0.6});
	/* cells (n = {sin (theta), 0, cos (theta)}, alpha = p_p.x*sin(theta) + p_p.z*cos(theta)); */
	sprintf (name, "case-%d/lmax-13/l2-level-13-t-%d-theta-%d-pis6.png", i, it, itheta);
	save (name);

	clear(); clear();
	view (fov = 2, theta = theta, relative = false,
	      tx = -(p_p.x)/(L0), ty = -(p_p.y + 7.*(d))/(L0),
	      bg = {1,1,1},
	      width = 400, height = 800);

	draw_vof ("cs", "fs");
	isosurface ("l2", lambda2, fc = {0.3,0.4,0.6}, lc = {0.3,0.4,0.6});
	cells (n = {sin (theta), 0, cos (theta)}, alpha = p_p.x*sin(theta) + p_p.z*cos(theta));
	sprintf (name, "case-%d/lmax-13/l2-zoom-level-13-t-%d-theta-%d-pis6.png", i, it, itheta);
	save (name);
      }
      
      free_grid();
    }
  }

  return 0;
}
