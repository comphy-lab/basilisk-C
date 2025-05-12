/**
# Potential flow past a digitally-generated tree, 

on a digitally-generated surface:

![Do not mind the blue sky](colonizedtree/tree.mp4)

![Very realistic.... Image via
 [Pixabay](https:/pixabay.com)](https://cdn.pixabay.com/photo/2018/03/14/09/49/tree-3224754__340.jpg)
 */

#include "grid/octree.h"
#include "embed.h"
#include "poisson.h"
#include "run.h"
#include "colonization.h"
#include "../ABL/perlin.h"
#include "../tracer-particles.h"
#include "view.h"
#include "../scatter2.h"

scalar pot[], laplace[]; 
Particles P;

pot[front] = dirichlet (1);
pot[back] = dirichlet (0);
pot[embed] = neumann (0.);
vector u[];
u.n[front] = neumann(0);
u.n[back] = neumann(0);
int ja = 1000; double R1 = 3, H = 1.5, Hg = 1;

double sy = 1.4;
int main () {
  L0 = 10;
  X0 = Z0 = -L0/2;
  N = 64;
  run();
}

/**
## Initialization

In the initialization stage we first find the tree nodes via the space
colonization algorithm and convert those to Branches. Next, the tree
interface can be reconstructed in terms of tree-cell-volume
fractions. A periodic Perlin topography is added before an equation
for the flow potential (`pot`) is solved. This is then converted into
the velocity field. Finally this field is visualized with tracer
particles.
 */

face vector fp[];
scalar cp[];
event init (t = 0) {
  //Seed atraction points in a prolate semi spheriod
  srand (0);
  coord atr[ja];
  int j = 0;
  while (j < ja) {
    double R = R1*pow(fabs(noise()), 1./3.1);
    coord vec = {noise(), noise(), noise()};
    normalize (&vec);
    if (vec.y > 0) {
      foreach_dimension()
	atr[j].x = R*vec.x;
      atr[j].y *= sy;
      atr[j].y += H + Hg;
      j++;
    }
  }

  //Get tree nodes
  int tn = 1;
  Tnode * nodes = calloc (tn, sizeof(Tnode));
  nodes[0].x = 0;
  nodes[0].y = Hg;
  nodes[0].z = 0;
  tn = colonize (atr, ja, &nodes, 1, 2, 0.2, 0.3);

 // Convert nodes to branches
  Branch * trees = calloc (499, sizeof(Branch)); //?
  printf ("# nodes: %d\n", tn);
  nb = nodestotree (nodes, tn, &trees);
  printf ("# branches: %d\n", nb);
  fflush (stdout);

  //Convert branches to volume fractions
  int it = 0;
  do {
    tree_interface (trees, cs, fs);
    boundary({cs});
    printf ("# Tree Reconstruction Iteration %d\n", ++it);
  } while (adapt_wavelet ({cs}, (double[]){0.1}, 9, 4).nf > (grid->tn/100));
  tree_interface (trees, cs, fs);
  boundary({cs, fs});
  fractions_cleanup (cs, fs);

  //Obtain surface volume and face fractions
  cp.prolongation = cs.prolongation;
  vertex scalar phi[];
  int nx = 5;
  int ny = 3;
  init_perlin (nx, ny);
  it = 0;
  do {
    foreach_vertex()
      phi[] =  y - Hg - (perlin (x, z, nx, ny));
    boundary ({phi});
    fractions (phi, cp, fp);
    boundary ({cp});
    printf ("# Surface Reconstruction Iteration %d\n", ++it);
  }  while (adapt_wavelet ({cs, cp}, (double[]){0.1, 0.01}, 8, 4).nf > (grid->tn/1000));
  foreach_vertex()
    phi[] =  y - Hg - (perlin (x, z, nx, ny));
  boundary ({phi});
  fractions (phi, cp, fp);
 
  //Combine tree and surface
  foreach()
    cs[] = min (cs[], cp[]);
  foreach_face()
    fs.x[] = min(fs.x[], fp.x[]);
  boundary ({cs, fs});
  fractions_cleanup (cs, fs);
  
  //Obain potential
  printf ("# Start solving for the Potential\n");
  NITERMIN = 10;
  mgstats mg = poisson (pot, laplace, alpha = fs,
			embed_flux = embed_flux, tolerance = HUGE);
  printf ("# Solved with a max. residual %g\n", mg.resa);

  //Get potential flow field
  foreach() {
    foreach_dimension()
      u.x[] = center_gradient(pot) * cs[];
  }
  boundary ((scalar*){u});

  //Place tracer particles
  P = init_tp_square (n = 30, ym = H + Hg, l = 2*R1, zp = -L0/4);
  DT = 0.1;

  //free stuff!
  free (trees);
  free (nodes);
  free (gradp);
}

event set_dtmax (i++) {
  dt = dtnext(DT);
}

double th = 0;
event movie (t += 1.) {
  view (bg = {0.3, 0.3, 0.9}, theta = th, phi = 0.2);
  view (ty = -.3);
  draw_vof ("cs", "fs", fc = {0.73,0.4,0.2});
  translate (y = 0.01)
    draw_vof ("cp", "fp", fc = {0.2,0.9,0.2});
  scatter (P);
  box();
  save ("tree.mp4");
  th += 0.006;
}

event stop (t = 500);
