#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "navier-stokes/perfs.h"

scalar s[];
scalar * tracers = {s};

double U0;
double R_d=0.01;
double Re=40;

int P_U=1;
int maxlevel;
double MAXTIME = 500;

face vector muv[];

FILE * fpmax; 

int main(int argc, char * argv[]) {
  
  L0  = 1;      // Domain size
  U0  = 1;      // Jet velocity
  maxlevel = 9;  

  if (argc > 1)
  {
    maxlevel = atoi(argv[1]);    //maximum refinement level (2^maxlevel)
  	MAXTIME = atoi(argv[2]); //maximum runtime
  	Re = atof(argv[3]);      // Re number 
	n = atof(argv[4]);       // refinement of the mesh 
    P_U = atof(argv[4]);         // Inlet condition: Uniform (0) Poiseuille (1) 
   }
  
  N = 2 << maxlevel; 

  TOLERANCE = 1e-4 [*]; 

  // Boundary conditions
  if (P_U == 1){
    u.n[bottom] = dirichlet (U0*(x > -R_d)); 
    u.t[bottom] = dirichlet(0.);
    s[bottom]   = dirichlet (U0*(x > -R_d));
  }else{
    u.n[bottom] = (x > -R_d) ? dirichlet(3/2*R_d*(1-pow(x/R_d,2)):dirichlet(0.); // normalized such that we have the same flow rate 
    u.t[bottom] = dirichlet(0.);
    s[bottom]   = (x > -R_d) ? dirichlet(3/2*R_d*(1-pow(x/R_d,2)):dirichlet(0.);
  }

  u.n[top] = dirichlet(0);
  p[top]   = neumann(0);

  u.n[left] = y < R_d ? dirichlet(-U0) : dirichlet(0.);
  u.t[left] = y < R_d ? neumann(0.) : dirichlet(0.);

  u.n[right] = dirichlet(0.);
  p[right]   = neumann(0.);

  //Initialize grid
  origin(-L0, 0);
  init_grid(N);
  
  //Save domain parameters to file
  char param_dim[80];
  if (pid() ==0){
    sprintf(param_dim, "param_dim.txt");
    FILE * fparam = fopen(param_dim, "w");
    fprintf(fparam, "%g %g %g %d \n", L0, R_d, U0, N);
      fclose(fparam);
   }

  mu = muv;

  fpmax = fopen("log.dat", "a"); 

  run();
}

event init (t = 0) {
  if (!restore("restart")) {
    fprintf(stderr, "Starting new simulation.\n");
  } else {
    fprintf(stderr, "Restarting from previous dump. \n");
  }
}

event properties (i++) {
  foreach_face()
    muv.x[] = fm.x[] * U0 * R_d / Re;
}

double sum_u = 0.;
double sum_u_t = 0.;
double res_u = 0.;
double sum_p = 0.;
double sum_p_t = 0.;
double res_p = 0.;

event logfile (i++) { 
  sum_u_t = 0.;
  sum_p_t = 0.;

  foreach_face(reduction(+:sum_u_t)) {
    sum_u_t += sqrt(sq(u.x[]) + sq(u.y[]));
  }
  foreach(reduction(+:sum_p_t)) {
    sum_p_t += p[];
  }

  res_u = fabs(sum_u_t - sum_u);
  res_p = fabs(sum_p_t - sum_p);

  sum_u = sum_u_t;
  sum_p = sum_p_t;
  if (pid() == 0){
        fprintf(stderr, "%d %g %g %g\n", i, t, res_u, res_p); 
        fprintf(fpmax, "%d %g %g %g\n", i, t, res_u, res_p);
    }
}

event profile (t = end) {
  char name[80];
  sprintf(name, "res_end.txt");
  FILE * fpres = fopen(name, "w");
  foreach(serial)
    fprintf(fpres, "%g %g %g %g %g %g \n", x, y, u.x[], u.y[], p[], s[]);
  fclose(fpres);
  
  printf("-----END-----\n");
}

event movie_frames (t = 0; t += 0.5; t <= 500) {
  char ufile[80], sfile[80];
  sprintf(ufile, "uY.mp4");
  sprintf(sfile, "s.mp4");
  output_ppm(u.y, file = ufile, n = 512, min = -U0, max = +U0, linear = true);
  output_ppm(s,   file = sfile, n = 512, min = 0.,  max = U0,  linear = true);
}

#if TREE
event adapt (i++) {
  double uemax = 4e-3;
  adapt_wavelet ((scalar *){u}, (double[]){uemax,uemax}, maxlevel, maxlevel -3);
}
#endif