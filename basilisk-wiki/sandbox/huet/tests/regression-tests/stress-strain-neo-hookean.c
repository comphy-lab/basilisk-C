/**
# Uniaxial stretch of a flat membrane

In this case, we define a rectangular flat membrane composed of just two
triangles, stretch it in and plot the non-zero principal tension $T_1$ with
respect to its corresponding the principal stretches
$e_1 = \frac{1}{2}(\lambda_1^2 -1)$:
$$
T_1 = \frac{1}{\lambda_2}\frac{\partial W}{\partial \lambda_1},
$$
with $\lambda_i = \frac{l}{l_0}$ the stretch ratio (ratio of current over
initial length) and $W(\lambda_1, \lambda_2)$ the strain energy function of the
considered elastic law.

Following [Barthès-Biesel](barthes2002effect), we ensure the deformation of the
membrane satisfies $T_2 = 0$, allowing to express $\lambda_2$ as a function of
$\lambda_1$.

## Setup
We define the relevant quantities: domain size, number of capsules (membranes),
number of Lagrangian points (nodes), octree level, characteristic size of the
rectangular membrane, etc.
*/
#define L0 1.
#define NCAPS 1
#define LEVEL 1
#define H0 (L0/8.)
/** Change the SKALAK macro to 0 to consider a neo-Hookean membrane. */
#define SKALAK 0

/**
Technically we don't need the octree grid nor the Navier-Stokes solver, but we
keep them because the stress computations are already conveniently called using
the event framework.
And this case is super fast so we don't really care :)
*/
#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "lagrangian_caps/capsule-ft.h"
#if SKALAK
  #include "lagrangian_caps/skalak-ft.h"
#else
  #include "lagrangian_caps/neo-hookean-ft.h"
#endif
#include "lagrangian_caps/view-ft.h"

int main(int argc, char* argv[]) {
  origin(-.5*L0, -.5*L0, -.5*L0);
  N = 1 << LEVEL;
  run();
}

/** In the event below we manually create a rectangular membrane composed of
four nodes (the vertices of the rectangles), five edges (its sides plus one
diagonal) and two triangles. */
event init (i = 0) {
  /** Initialize two triangles */
  lagMesh* mesh = &CAPS(0);
  mesh->nln = 4;
  mesh->nle = 5;
  mesh->nlt = 2;
  mesh->nodes = malloc(mesh->nln*sizeof(lagNode));
  mesh->edges = malloc(mesh->nle*sizeof(Edge));
  mesh->triangles = malloc(mesh->nlt*sizeof(Triangle));
  /** Create nodes, edges and triangles */
  mesh->nodes[0].pos.x = -H0;
  mesh->nodes[0].pos.y = -H0/2;
  mesh->nodes[0].pos.z = 0.;
  mesh->nodes[1].pos.x = -H0;
  mesh->nodes[1].pos.y = H0/2;
  mesh->nodes[1].pos.z = 0.;
  mesh->nodes[2].pos.x = H0;
  mesh->nodes[2].pos.y = H0/2;
  mesh->nodes[2].pos.z = 0.;
  mesh->nodes[3].pos.x = H0;
  mesh->nodes[3].pos.y = -H0/2;
  mesh->nodes[3].pos.z = 0.;
  mesh->edges[0].node_ids[0] = 0;
  mesh->edges[0].node_ids[1] = 1;
  mesh->edges[1].node_ids[0] = 1;
  mesh->edges[1].node_ids[1] = 2;
  mesh->edges[2].node_ids[0] = 2;
  mesh->edges[2].node_ids[1] = 3;
  mesh->edges[3].node_ids[0] = 3;
  mesh->edges[3].node_ids[1] = 0;
  mesh->edges[4].node_ids[0] = 1;
  mesh->edges[4].node_ids[1] = 3;
  mesh->triangles[0].node_ids[0] = 0;
  mesh->triangles[0].node_ids[1] = 1;
  mesh->triangles[0].node_ids[2] = 3;
  mesh->triangles[0].edge_ids[0] = 0;
  mesh->triangles[0].edge_ids[1] = 3;
  mesh->triangles[0].edge_ids[2] = 4;
  mesh->triangles[1].node_ids[0] = 1;
  mesh->triangles[1].node_ids[1] = 2;
  mesh->triangles[1].node_ids[2] = 3;
  mesh->triangles[1].edge_ids[0] = 1;
  mesh->triangles[1].edge_ids[1] = 2;
  mesh->triangles[1].edge_ids[2] = 4;
  /** Create neighbors informations */
  mesh->nodes[0].nb_neighbors = 2;
  mesh->nodes[0].nb_triangles = 1;
  mesh->nodes[0].neighbor_ids[0] = 1;
  mesh->nodes[0].neighbor_ids[1] = 3;
  mesh->nodes[0].edge_ids[0] = 0;
  mesh->nodes[0].edge_ids[1] = 3;
  mesh->nodes[0].triangle_ids[0] = 0;
  mesh->nodes[1].nb_neighbors = 3;
  mesh->nodes[1].nb_triangles = 2;
  mesh->nodes[1].neighbor_ids[0] = 0;
  mesh->nodes[1].neighbor_ids[1] = 2;
  mesh->nodes[1].neighbor_ids[2] = 3;
  mesh->nodes[1].edge_ids[0] = 0;
  mesh->nodes[1].edge_ids[1] = 1;
  mesh->nodes[1].edge_ids[2] = 4;
  mesh->nodes[1].triangle_ids[0] = 0;
  mesh->nodes[1].triangle_ids[1] = 1;
  mesh->nodes[2].nb_neighbors = 2;
  mesh->nodes[2].nb_triangles = 1;
  mesh->nodes[2].neighbor_ids[0] = 1;
  mesh->nodes[2].neighbor_ids[1] = 0;
  mesh->nodes[2].edge_ids[0] = 1;
  mesh->nodes[2].edge_ids[1] = 2;
  mesh->nodes[2].triangle_ids[0] = 0;
  mesh->nodes[3].nb_neighbors = 3;
  mesh->nodes[3].nb_triangles = 2;
  mesh->nodes[3].neighbor_ids[0] = 2;
  mesh->nodes[3].neighbor_ids[1] = 0;
  mesh->nodes[3].neighbor_ids[2] = 1;
  mesh->nodes[3].edge_ids[0] = 2;
  mesh->nodes[3].edge_ids[1] = 3;
  mesh->nodes[3].edge_ids[2] = 4;
  mesh->nodes[3].triangle_ids[0] = 1;
  mesh->nodes[3].triangle_ids[1] = 0;

  store_initial_configuration(mesh);
}

/**
Below we impose the position of the membrane in order to span a range of
principal strains, while ensuring the tension $T_2$ is always zero. This results in a non-zero stretch ratio in the y-direction:

$\lambda_2^{NH} = 1/\sqrt{\lambda_1}$ for the neo-Hookean law, and

$$\lambda_2^{SK} = \sqrt{\frac{1 + C \lambda_1^2}{1 + C \lambda_1^4}}$$
for the Skalak law.
*/
event acceleration (i++) {
  lagMesh* mesh = &CAPS(0);
  double strain_comp = 0.01*i;
  double stretch1 = sqrt(2*strain_comp +1);
  #if SKALAK
    double stretch2 = sqrt((1 + AREA_DILATATION_MODULUS*sq(stretch1))/
      (1 + AREA_DILATATION_MODULUS*sq(sq(stretch1))));
  #else // neo-hookean membrane
    double stretch2 = 1/sqrt(stretch1);
  #endif
  mesh->nodes[0].pos.x = -stretch1*H0;
  mesh->nodes[1].pos.x = -stretch1*H0;
  mesh->nodes[2].pos.x = stretch1*H0;
  mesh->nodes[3].pos.x = stretch1*H0;
  mesh->nodes[0].pos.y = -.5*H0*stretch2;
  mesh->nodes[1].pos.y = .5*H0*stretch2;
  mesh->nodes[2].pos.y = .5*H0*stretch2;
  mesh->nodes[3].pos.y = -.5*H0*stretch2;
  mesh->nodes[0].pos.z = 0.;
  mesh->nodes[1].pos.z = 0.;
  mesh->nodes[2].pos.z = 0.;
  mesh->nodes[3].pos.z = 0.;
  for(int j=0; j<4; j++) foreach_dimension() mesh->nodes[j].lagForce.x = 0.;
}

event logfile (i++) {
  comp_elastic_stress(&CAPS(0));
  for(int i=0; i<CAPS(0).nlt; i++) {
    fprintf(stderr, "%d, %g, %g, %g, %g\n", i,
    .5*(sq(CAPS(0).triangles[i].stretch[0]) - 1.),
    CAPS(0).triangles[i].tension[0],
    .5*(sq(CAPS(0).triangles[i].stretch[1]) - 1.),
    CAPS(0).triangles[i].tension[1]);
  }
}

event end (i = 100) {
  lagMesh* mesh = &CAPS(0);
  free(mesh->nodes);
  free(mesh->edges);
  free(mesh->triangles);
  return 0;
}

/**
## Results
~~~gnuplot
set title "Stress-strain response of flat elastic neo-Hookean (bottom) \nand Skalak (top) membranes under an uniaxial elongation"
set grid
set xlabel "Principal strain"
set ylabel "Principal stress"
set key top left reverse Left

set label at 0.68,1.2 "Neo-Hookean law"

l(x) = sqrt(2*x + 1)
nh(x) = (l(x)**3 - 1)/l(x)**1.5
sk(x) = l(x)*(l(x)**2 - 1)*sqrt((1 + l(x)**2)/(1 + l(x)**4))*((1 + l(x)**4)/(1 + l(x)**2) + 1/(1 + l(x)**4))

set samples 20
set pointsize 1.5

plot 'log' using 4:(3*$5) every 10 w p pt 1 lc rgb "red" title "This study", \
nh(x) w l lc -1 title "Exact"


~~~

## References
~~~bib
@Article{barthes2002effect,
  author    = {Barthes-Biesel, Dominique and Diaz, Anna and Dhenin, Emmanuelle},
  journal   = {Journal of Fluid Mechanics},
  title     = {Effect of constitutive laws for two-dimensional membranes on flow-induced capsule deformation},
  year      = {2002},
  pages     = {211--222},
  volume    = {460},
  file      = {:effect-of-constitutive-laws-for-two-dimensional-membranes-on-flow-induced-capsule-deformation.pdf:PDF},
  groups    = {Biological flows},
  publisher = {Cambridge University Press},
}
~~~
*/
