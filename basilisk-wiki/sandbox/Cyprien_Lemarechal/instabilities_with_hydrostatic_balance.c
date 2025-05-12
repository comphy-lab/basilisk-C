#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"
// #include "../output_vtu_foreach.h" //paraview

scalar b[];
scalar * tracers = {b};
face vector av[];

int maxlevel = 9, minlevel = 5;

double nu = 1e-6;

double timestep_ini =1.;

int end_simu =400;

b[bottom]   = dirichlet (b[]); 
b[top]  = neumann (0.);

pf[top]   = dirichlet (0.);
p[bottom] = neumann (-b[]); 

int main() {
  periodic (left);
  periodic (front);

  const face vector muc[] = {nu, nu, nu};

  L0 = 3.;
  mu = muc;
  a = av;
  
  foreach_dimension()
    u.x.refine = refine_linear;
  p.refine = p.prolongation = refine_linear;
  
  init_grid (1 << minlevel);
  
  run();
}

event init (t = 0) {
  DT = 0.05;
  
  // pressure with an analytical relation
  foreach()
    p[] = y*y/2. - L0*y + (L0 * L0)/2.;
    // p[] =  -y +3.; //homogeneous fluid
  boundary({p}); 	
    
  foreach()
    b[] = y - L0;
    // b[] = -1; //homogeneous fluid
  boundary({b});

  //Another way to initialize pressure
  /*
  foreach()
    p[] = 0.;
  boundary({p});

  foreach()
    p[] = p[0, 1, 0]  -  Delta * (b[] + b[0, 1, 0])/2.;
  boundary({p});
  */
}

event acceleration (i++) {
  foreach_face(y)
    av.y[] = (b[] + b[0,-1])/2.;

  boundary (all);
}

event tracer_diffusion (i++)
  diffusion (b, dt, mu);

event adapt (i++) {
  //trying to give time for p to be balanced
  DT = 0.05;
  
  double ue = 0.1;
  if (i > 50){
    DT = timestep_ini;
    adapt_wavelet ((scalar *){b, u}, (double[]){ue, ue, ue, ue}, maxlevel, minlevel);
  }
}

event print(i++){
  double avg_uy = 0., vol = 0.;
  foreach(reduction(+:avg_uy) reduction(+:vol)) {
    vol += dv();
    avg_uy += dv() * sqrt(sq(u.y[]));
  }
  avg_uy /= vol;
    
  if (pid() == 0){
  printf("iteration %d\n", i);
  printf("dt :%f\n", dt);
  printf("avg_uy :%f\n", avg_uy);
  printf("*******");
  }
}

// event export_to_paraview (i++){
  // char name[350], namebis[80];
  // sprintf(name, "wip%d",i);
  // sprintf(namebis, "wip%d",i);

  // output_vtu((scalar *) {av, u , b, p}, (vector *) {0}, name);
 // }

event end (i= end_simu) {
}