/**
# The dynamics of a blob of vorticity

## Initial conditions and simulation parameters

The evolution of a blob of vorticity is simulated as in [../2d/blob.c]()
but in 3D.  The initial conditions are read from a gnuplot-compatible binary
file, obtained from `vorticity_18deg.mat`. We preprocess this field using
`blob.m`, which essentially extends the domain, centers the vorticity patches
and applies a Hanning window to ensure that vorticity vanishes smoothly away
from the center. The resulting vorticity field looks something like this

![Initial vorticity field](blob/omega.mp4)

The goal is to see, wether the initial vortiity distribution results in
unexpected growth of the vortex core size. One of the vortex is destabilized,
but doesn't translate in a thicker core, and it doesn't trigger a merging.
*/

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#define MINLEVEL 4
#define RAD (sqrt(sq(x) + sq(y)))

int main(){
  L0 = 20;
  X0 = Y0 = Z0 = -L0/2;
  init_grid (1 << MINLEVEL);
  periodic(top);
  periodic(left);
  periodic(front);
  run();
}

event adapt (i++){
  adapt_wavelet ((scalar*){u}, (double[]){1e-4,1e-4,1e-4}, MAXLEVEL, MINLEVEL);
}

#include "../../input_fields/auxiliar_input.h"
event init (t = 0){

  refine  (RAD < 10 && level < MAXLEVEL-1);
  refine  (RAD < 5 && level < MAXLEVEL);

  vector psi[], omega[];
  foreach(){
    foreach_dimension(){
      u.x[] = 0.;
      omega.x[] = 0.;
      psi.x[] = 0.;
    }
  }

  FILE * fp = fopen("./omega.bin", "r");
  if (!fp) printf("Binary file not found");
  input_matrix(omega.z,fp,N,X0,Y0,L0);
  fclose (fp);

  fp = fopen("./axial_velocity.bin", "r");
  if (!fp) printf("Binary file not found");
  input_matrix(u.z,fp,N,X0,Y0,L0);
  fclose (fp);

  boundary ((scalar*){psi,omega});
  poisson (psi.x, omega.x);
  poisson (psi.y, omega.y);
  poisson (psi.z, omega.z);

  foreach(){
    u.x[] -= ((psi.z[0,1,0] - psi.z[0,-1,0]) - (psi.y[0,0,1] - psi.y[0,0,-1]))/(2.*Delta);
    u.y[] -= ((psi.x[0,0,1] - psi.x[0,0,-1]) - (psi.z[1,0,0] - psi.z[-1,0,0]))/(2.*Delta);
    u.z[] -= ((psi.y[1,0,0] - psi.y[-1,0,0]) - (psi.x[0,1,0] - psi.x[0,-1,0]))/(2.*Delta);
  }
  boundary ((scalar *){u});

  fp = fopen("vortex_z0.asc", "w");
  fputs ("[1]t\t [2]tag\t [3]Gamma\t [4]mu_x\t [5]mu_y\t [6]mu_z\t [7]M20\t [8]M02\t [9]M11\t [10]a\t [11]b\t [12]c\t [13]e\t [14]maxvor \n", fp);
  fclose(fp);

  fp = fopen("blob3d.asc", "w");
  fputs ("[1]t\t [2]ekin\t [3]enstrophy\t [4]helicity \n", fp);
  fclose (fp);
}

/**
## Ellipticity of the vortex dipole
The goal is to evaluate several quantities that are pertinent to describe the
flow field such as circulation, sizes, and ellipticity. When a vorticity field
possesses an elliptical shape, these quantities are sufficient to define its
geometry.
*/


#include "ellipticity.h"
event logfile (t += 0.05) {
  vector omega[];
  vorticity3d(u, omega);
  stats s = statsf (omega.z);

  scalar m[];
  foreach()
    m[] = abs(omega.z[]) > s.stddev/100;

  FILE * fp = fopen("vortex_z0.asc", "a");
  vorticity_moments_plane(omega.z, m, fp, (coord){0,0,1}, 0.);
  fclose(fp);
}

/**
# Additional outputs
*/

#include "lambda2.h"
#include "view.h"
event movie (t += 0.05) {
  scalar l2[];
  lambda2 (u, l2);
  stats f = statsf (l2);

  vector omega[];
  vorticity3d(u, omega);
  stats s = statsf (omega.z);

  view (theta = pi/6, phi = pi/6, fov=35);
  isosurface ("l2", -f.stddev);
  squares ("l2", linear = false, alpha=-L0/2);
  box();
  save ("lambda2.mp4");

  isosurface ("omega.z",  s.stddev);
  isosurface ("omega.z", -s.stddev);
  squares ("omega.z", linear = false, alpha=-L0/2);
  box();
  save ("omega.mp4");
}

#include "../../output_fields/output_vtu_foreach.h"
event snapshots (t += 2.5) {

  scalar l2[];
  lambda2 (u, l2);
  stats f = statsf (l2);

  vector omega[];
  vorticity3d(u, omega);
  stats s = statsf (omega.z);

  static int nf = 0;
  char name[80];
  sprintf(name, "blob_%3.3d", nf);
  output_vtu ((scalar *) {l2}, (vector *) {u, omega}, name);
  nf++;
}

event slices (t += 0.25) {

  scalar l2[];
  lambda2 (u, l2);
  stats f = statsf (l2);

  vector omega[];
  vorticity3d(u, omega);
  stats s = statsf (omega.z);

  static int ns = 0;
  char name[80];
  sprintf(name, "blob_x0_%3.3d", ns);
  output_vtu_plane ((scalar *) {l2}, (vector *) {u, omega}, name, (coord){1,0,0}, 0.0);

  sprintf(name, "blob_y0_%3.3d", ns);
  output_vtu_plane ((scalar *) {l2}, (vector *) {u, omega}, name, (coord){0,1,0}, 0.0);

  sprintf(name, "blob_z0_%3.3d", ns);
  output_vtu_plane ((scalar *) {l2}, (vector *) {u, omega}, name, (coord){0,0,1}, 0.0);

  view (theta = pi/6, phi = pi/6, fov=35);
  isosurface ("l2", -f.stddev);
  squares ("omega.z", linear = false, alpha=-L0/2);
  box();
  sprintf(name, "omegaz_%3.3d.png", ns);
  save (name);
  ns++;
}

event output (t = 2.5) {
  dump();
}
