/**
# A rapidly growing tree

This page tests the combination of the centered solver with an
increasing number of embedded tree branches. Using this, the new tree
is automatically initilaized with a nearly consistent complex flow
field.

![It seems to work well in 2D](add_branches/mov.mp4)
 */

#include "embed.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/double-projection.h"
#include "view.h"
#include "treegen.h"
#define BVIEW 1
#include "../particles.h"

Branch * trees;
double Re = 500;
int maxlevel = 10, minlevel = 4;
int nt;

void prolongate_ratio (Point point, scalar s) {
  foreach_child() {
    if (s[] != nodata)
      s[] += s[]*Delta;
  }
}

u.t[embed] = dirichlet (0.);
u.n[embed] = dirichlet (0.);
#if dimension == 3
u.r[embed] = dirichlet (0.);
#endif

scalar J[], res[];
face vector muc[];

int main () {
  periodic (left);
  mu  = muc;
  trees = tree_skeleton();
#if _MPI
  MPI_Bcast (trees, sizeof(Branch)*nb, MPI_BYTE, 0, MPI_COMM_WORLD);
#endif
nt = nb;
  L0 = 55;
  X0 = Z0 = -L0/2;
  nb = 1;
  srand(0); //Reproducible tree
  res.prolongation = prolongate_ratio;
  N = (1 << 9);
  run();
}

event init0 (t = 0) {
  foreach()
    u.x[] = cs[];
  boundary (all);
}

event init (t = 0) {
  if (n_part < 1)
    init_particles_2D_square_grid (10, -19., 15., 15.);
  for (scalar s in all) //Naive refinement to circumvent assertions 
    if (s.refine == refine_embed_linear || s.prolongation == refine_embed_linear) 
      s.refine = s.prolongation = refine_bilinear;
  do {
    tree_interface (trees, cs, fs, J);
    foreach() {
      res[] = nodata;
      if (cs[] > 0 && cs[] < 1) {
	res[] = 1/sqrt(trees[(int)(J[] + 0.5)].R/Re);
      }
    }
    boundary ({res, u});
  } while (adapt_wavelet ({res, u}, (double[]){2, 0.1, 0.1, 0.1},
			  maxlevel, minlevel).nf > 10);
  // Correct for the naive refinement
  foreach() {
    foreach_dimension()
      u.x[] *= cs[];
  }
  foreach_face() // uf[embed] = 0, this fixes assertions in timestep.h
    uf.x[] *= fs.x[];
  DT = 0.1;
}
/**
A branch is added every-so-often by increasing `nb` and calling all
`init` events.
 */
event add_branch (t = 15; t += 15) {
  nb = min(nb++, nt);
  event ("init");
}

event properties (i++) {
  foreach_face()
    muc.x[] = fs.x[]/Re;
  boundary ((scalar*){muc});
}

event mov (t += 0.25) {
  scalar omg[];
  vorticity (u, omg);
  scatter (loc, pc = {sin(pid()), cos(pid()), sin(2.3*pid())});
  draw_vof ("cs", "fs", filled = -1, fc= {98./128., 78./128, 44./128.});
  squares ("omg", min = -5, max = 5, map = cool_warm);
  mirror ({0,-1})
    cells();
  save ("mov.mp4");
}

event adapt (i++) {
  foreach() {
    res[] = nodata;
    if (cs[] > 0 && cs[] < 1) {
      res[] = 1/sqrt(trees[(int)(J[] + 0.5)].R/Re);
    }
  }
  boundary ({res});
  adapt_wavelet ({res, u}, (double[]){2, 0.1, 0.1, 0.1}, maxlevel, minlevel);
}

event stop (t = 200) {
  free (trees);
  return 1;
}

