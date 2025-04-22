/**
# Hydrostatic balance with refined embedded boundaries in 3D

This is very close to [hydrostatic2.c]() but in 3D. */

#include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "view.h"

/**
We use a similar porous medium as in [/src/examples/porous3D.c](). */

void porous (scalar cs, face vector fs)
{
  int ns = 60; // 160, 80
  coord pc[ns];
  double R[ns];
  srand (0);
  double a = 0.5, b = 0.04;
  for (int i = 0; i < ns; i++) {
    foreach_dimension()
      pc[i].x = a*noise();
    R[i] = 2.*(0.02 + b*fabs(noise()));
  }
    
  vertex scalar phi[];
  foreach_vertex() {
    phi[] = HUGE;
    for (double xp = -L0; xp <= L0; xp += L0)
      for (double yp = -L0; yp <= L0; yp += L0)
	for (double zp = -L0; zp <= L0; zp += L0)
	  for (int i = 0; i < ns; i++)
	    phi[] = intersection (phi[], (sq(x + xp - pc[i].x)
					  + sq(y + yp - pc[i].y)
					  + sq(z + zp - pc[i].z)
					  - sq(R[i])));
    phi[] = - phi[];
  }
  fractions (phi, cs, fs);
  fractions_cleanup (cs, fs);
}

const face vector G[] = {1.,2.,3.};

void restriction_exact (Point point, scalar s)
{
  s[] = (cs[] != 0.)*(G.x[]*x + G.y[]*y + G.z[]*z); // exact pressure
}

void refine_exact (Point point, scalar s)
{
  foreach_child()
    s[] = (cs[] != 0.)*(G.x[]*x + G.y[]*y + G.z[]*z); // exact pressure
}

int main()
{
  init_grid (1 << 5);
  origin (-0.5, -0.5, -0.5);

  /**
  The events of the Navier-Stokes solver are called "by hand". */
  
  event ("metric");
  event ("defaults");

  porous (cs, fs);

#if 0
  refine (level < 7 && cs[] > 0 && cs[] < 1);
  porous (cs, fs);
#else
  adapt_wavelet ({cs}, (double[]){0.01}, 7);
  porous (cs, fs);
  adapt_wavelet ({cs}, (double[]){0.01}, 7);
  porous (cs, fs);
#endif

  a = G;

  TOLERANCE = 1e-9;
  NITERMAX = 10;
  alpha = fm;
  dt = 1.;

  foreach()
    p[] = (cs[] != 0.)*(G.x[]*x + G.y[]*y + G.z[]*z); // exact pressure
#if 0
  p.restriction = restriction_exact;
  p.refine = p.prolongation = refine_exact;
  p.dirty = true;
#endif
  
  event ("acceleration");
#if 1
  event ("projection");
#else
  foreach_face()
    uf.x[] -= alpha.x[] ? dt*alpha.x[]*face_gradient_x (p, 0) : 0.;

  face vector gf[];
  foreach_face()
    gf.x[] = fm.x[] ? fm.x[]*a.x[] - alpha.x[]*(p[] - p[-1])/Delta : 0.;
  
  trash ({g});
  foreach()
    foreach_dimension()
      g.x[] = (gf.x[] + gf.x[1])/(fm.x[] + fm.x[1] + SEPS);

  correction (dt);
#endif

  /**
  We check the convergence rate and the norms of the velocity field
  (which should be negligible). */
  
  fprintf (stderr, "mgp %.10f %.10f %d %d %d\n",
	   mgp.resb, mgp.resa, mgp.i, mgp.minlevel, mgp.nrelax);
  fprintf (stderr, "umax %.12f %.12f %.12f %.12f %.12f\n",
	   normf(u.x).max, normf(u.y).max,
	   normf(uf.x).max, normf(uf.y).max, normf(uf.z).max);

  /**
  The pressure is hydrostatic, in each of the pores. 

  ![Cross-section of the pressure field.](hydrostatic3/p.png)
  */

  view (fov = 19, width = 400, height = 400);
  squares ("p", spread = -1);
  save ("p.png");

  p.nodump = false;
  scalar ep[];
  foreach()
    ep[] = p[] - (cs[] != 0.)*(G.x[]*x + G.y[]*y);
  dump ();

#if 0  
  p[bottom] = dirichlet (y);
  foreach() {
    if (cs[] > 0.)
      printf ("%g %g %g %g %g %g %g %g %g %g %g\n",
	      x, y, z, p[], u.x[], u.y[], u.z[], g.x[], g.y[], g.z[], cs[]);
    p[] = G.x[]*x + G.y[]*y;
  }
   
  foreach_face (x)
    printf ("fx %g %g %g %g\n", x, y, z, face_gradient_x (p, 0));
  foreach_face (y)
    printf ("fy %g %g %g %g\n", x, y, z, face_gradient_y (p, 0));
  foreach_face (z) {
    int i = 0;
    if (fs.z[0,0,i] > 0. && fs.z[0,0,i] < 1. &&
	embed_face_gradient_z (point, p, i) != 0.) {
      printf ("fz %g %g %g %g %g\n", x, y, z, (p[0,0,i] - p[0,0,i-1])/Delta,
	      embed_face_gradient_z (point, p, i));
      scalar a = p;
      {
	coord p = embed_face_barycentre_z (point, i);
	int j = sign(p.x), k = sign(p.y);
	printf ("bary %g %g %g %d %d\n", p.x, p.y, p.z, j, k);
	p.x = fabs(p.x), p.y = fabs(p.y);
	printf ("valz %g\n",
		(((a[0,0,i] - a[0,0,i-1])*(1. - p.x) +
		  (a[j,0,i] - a[j,0,i-1])*p.x)*(1. - p.y) + 
		 ((a[0,k,i] - a[0,k,i-1])*(1. - p.x) +
		  (a[j,k,i] - a[j,k,i-1])*p.x)*p.y)/Delta);
	for (int l = -1; l <= 1; l++)
	  for (int m = -1; m <= 1; m++)
	    printf ("point %d %d %g %g %g %g\n", l, m, fs.z[l,m,i],
		    (a[l,m,i] - a[l,m,i-1])/Delta, a[l,m,i], a[l,m,i-1]);
	exit (1);
      }
    }
  }
#endif
}
