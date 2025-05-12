/**
# Vortex shedding behind a cuboid at Reynolds = 300

This case is adapted from [examples/sphere.c] and change the embeded gemetry to cuboid. 
*/

#include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h"
#include "view.h"
#include "lambda2.h"

/** 
we define the function used for calculation of phi with the condition of inside cuboid or outside of cuboid. If the coord is inside, the distance is minimum to each wall. If the coord is outside, the distance is calculated based on dx, dy, dz. One side note is this function is not applicable when the cuboid rotate. 
*/
#define max3(a,b,c) max(max(a,b), c)

double Point2Rec (const coord *P,
                  double * min_x,
                  double * max_x,
                  double * min_y,
                  double * max_y,
                  double * min_z,
                  double * max_z)
{
  double d2;
  if ((*P).x > *min_x && (*P).x < *max_x && (*P).y > *min_y && (*P).y < *max_y && (*P).z > *min_z && (*P).z < *max_z){
    d2 = max3(max(*min_x-(*P).x, (*P).x-*max_x), max(*min_y-(*P).y, (*P).y-*max_y), max(*min_z-(*P).z, (*P).z-*max_z));
  }
  else{
    double dx = max3(*min_x-(*P).x, 0, (*P).x-*max_x);
    double dy = max3(*min_y-(*P).y, 0, (*P).y-*max_y);
    double dz = max3(*min_z-(*P).z, 0, (*P).z-*max_z);
    d2 = sqrt(sq(dx)+sq(dy)+sq(dz));
  }
  return d2;
}

/** the geometry of the cuboid */
double xw = -10, W = 1, H = 10, Lw = 10; 

int maxlevel = 8;

face vector muv[];

int main(){
  init_grid(128);
  size(50.);
  X0 = -L0/2;
  mu = muv;
  run();
}

/**
The viscosity  */

event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]/300.;
}

/**
The boundary conditions are inflow with unit velocity on the
left-hand-side and outflow on the right-hand-side. */

u.n[left]  = dirichlet(1.);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

/**
The boundary condition is no slip on the embedded boundary. */

u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);
u.r[embed] = dirichlet(0.);

event init (t = 0) {

  /**
  We initially refine only in a sphere, slightly larger than the solid
  sphere. */

  // refine (sq(x+9.5)  < sq(1) && y  < 11 && level < maxlevel);
  refine (sq(x+xw)  < sq(W) && y  < H+1 && sq(z-L0/2)  < sq(Lw/2) && level < maxlevel);
  vertex scalar phi[];
  double min_x = xw, max_x = xw + W, min_y = Y0, max_y = Y0 + H, min_z = L0/2-Lw/2, max_z = L0/2+Lw/2; //L0/2-5
  foreach_vertex() {
    coord cc = {x, y, z};
    phi[] = Point2Rec(&cc, &min_x, &max_x, &min_y, &max_y, &min_z, &max_z);
  }
  boundary({phi});
  fractions (phi, cs, fs);
  /**
  We set the initially horizontal velocity to unity everywhere
  (outside the sphere). */
  
  foreach()
    u.x[] = cs[] ? 1. : 0.;
}

/** log the number of iterations */
event logfile (i++)
  fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);

/**
We use view function to create the animated isosurface of $\lambda_2$
for $1 <= t <= 60$. */

event movies (t = 1; t += 0.25; t <= 60)
{
  #if (dimension == 2)
    scalar omega[];
    vorticity (u, omega);
    boundary ({omega});
    draw_vof ("cs", "fs", filled = -1, fc = {0.5,0.1,0.2});
    draw_vof ("cs", "fs");
    squares ("omega", linear = true, map = cool_warm);
    mirror ({0,-1})
    cells();
  #elif (dimension == 3)
  /**
  Here we compute two new fields, $\lambda_2$ and the vorticity
  component in the $y-z$ plane. */

  scalar l2[], vyz[];
  foreach()
    vyz[] = ((u.y[0,0,1] - u.y[0,0,-1]) - (u.z[0,1] - u.z[0,-1]))/(2.*Delta);
  boundary ({vyz});
  lambda2 (u, l2);

  view (fov = 15, theta = -0.8, phi = 0.4, 
	tx = -0.30, ty = 0.1, bg = {65./256,157./256,217./256},
	width = 1080, height = 1080);

  draw_vof ("cs", "fs", fc = {0.5,0.1,0.2});
  draw_vof ("cs", "fs");
  isosurface ("l2", -0.01, color = "vyz", min = -1, max = 1,
	      linear = true, map = cool_warm);
  #endif
  save ("mov_cub.mp4");

}
/** save the dump file */

event dumper (t += 10) {
  
  p.nodump = false;
  char str[99];
  sprintf (str, "dump3D%g", t);
  dump(str);
}

/**
We set an adaptation criterion with an error threshold of 0.02 on all
velocity components and $10^{-2}$ on the geometry. */

event adapt (i++) {
  astats s = adapt_wavelet ({cs,u}, (double[]){1e-2,0.02,0.02,0.02},
			    maxlevel, 4);
  fprintf (ferr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
}


  
  
  

