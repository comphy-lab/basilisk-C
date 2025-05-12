#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "PointTriangle.h"
#include "view.h"
#include "lambda2.h"

#define MINLEVEL 4
#define RAD (sqrt(sq(x) + sq(y)))

int n_seg;
double as = 0.1, Hs = 1.0, Rs = 1.0, n_turns;
double u_frame=0.534953599887698;

scalar l2[];
vector omega[], U[];
face vector Uf[];


int main() {
  L0 = 8;
  X0 = Y0 = Z0 = -L0/2;
  DT = 0.1;
  N = 1<<MINLEVEL;

  n_turns = L0/Hs;
  n_seg = 64 * n_turns;

  double reynolds= 8000;
  const face vector muc[] = {1./reynolds,1./reynolds,1./reynolds};
  mu = muc;

  periodic(back);
  run();
}

double delta_init = 1e-4, k_init = 0.5;
#include "filaments.c"
#include "3d/ellipticity.h"

coord P={0,0,0}, tvec={1,0,0}, nvec={0,1,0}, bvec={0,0,1};
double _alpha=0;

event init (t = 0){
  // 1. Load base flow and disable non-linear terms (stokes=True)
  // Here, the base flow is stored at the center U and faces Uf
  if (pid()==0)
    printf("Loading base flow \n");

  restore("../helical/profiles/l10z/dump_relaxed");
  stokes = true;
  scalar L2[];
  lambda2(u, L2);

  foreach()
    foreach_dimension()
      U.x[] = u.x[];
  boundary((scalar *){U});

  foreach_face()
    Uf.x[] = fm.x[]*face_value (U.x, 0);
  boundary((scalar *){Uf});
  mgpf = project (Uf, pf, alpha, dt/2., mgpf.nrelax);

  // 2. Generate the velocity field for a perturbed filament, then substract
  // from the base flow to obtain an initial perturbation
  helical_filament(u);
  foreach()
    u.z[] -= u_frame;
  boundary ((scalar *){u});
  lambda2 (u, l2);

  view(camera="iso", fov=4*L0);
  isosurface ("l2", 0);
  isosurface ("L2", 0, fc = {0.,0.7,0.7});
  box();
  save ("lambda2.png");

  // 3. Print the base flow and initial perturbation
  scalar * list = {p, u, U};
  for (scalar s in list){
    char name[80];
    view(camera="iso", fov=4*L0);
    squares (s.name, linear = false, n = {1,0,0});
    box();
    sprintf(name, "init_linearised_x0_%s.png", s.name);
    save (name);
  }

  // 4. Substract from the base flow to obtain an initial perturbation
  foreach()
    foreach_dimension()
      u.x[] -= U.x[];
  boundary ((scalar *){u});

  list = {p, u};
  for (scalar s in list){
    char name[80];
    view(camera="iso", fov=4*L0);
    squares (s.name, linear = false, n = {1,0,0});
    box();
    sprintf(name, "init_perturbation_x0_%s.png", s.name);
    save (name);
  }
  dump("dump_init");

  // restore("dump_init");
  // for (scalar s in list){
  //   char name[80];
  //   view(camera="iso", fov=4*L0);
  //   squares (s.name, linear = false, n = {1,0,0});
  //   box();
  //   sprintf(name, "init_linearised_x0_%s.png", s.name);
  //   save (name);
  // }

  FILE * fp = fopen("ekin.asc", "w");
  fclose(fp);
}

event advection_term (i++) {
  prediction();
  mgpf = project (uf, pf, alpha, dt/2., mgpf.nrelax);
  advection ((scalar *){U}, uf, dt, (scalar *){g});
  advection ((scalar *){u}, Uf, dt);
}

//## Outputs

double tsample1=0.01, tsample2=0.01, tsample3=0.01, tsample4=0.5;
event end_timestep (t += tsample1, last){
  lambda2 (u, l2);
  vorticity3d(u, omega);
}

event logfile (i++) {
  double ek = 0.;
  foreach(reduction(+:ek))
    foreach_dimension()
      ek += sq(u.x[])*dv();

  FILE * fp;
  if (pid()==0) {
    fp = fopen("ekin.asc", "a");
    fprintf (fp, "%.5f %.15g \n", t, ek);
    fclose(fp);
  }
}

#include "helical_export.h"
event movie (t += tsample1){
  char name[80];
  scalar * list = {p, u, omega};
  for (scalar s in list){
    view(camera="iso", fov=4*L0);
    squares (s.name, linear = false, n = {1,0,0});
    box();
    sprintf(name, "linearised_x0_%s.mp4", s.name);
    save (name);
  }

  stats f = statsf (u.z);

  view(camera="iso", fov=4*L0);
  isosurface ("u.z", f.stddev);
  isosurface ("u.z", -f.stddev);
  box();
  save ("linearised.mp4");
}

event stop (t = 3.0)
  dump();
