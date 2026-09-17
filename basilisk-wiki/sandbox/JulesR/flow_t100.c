#include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "tracer.h"//for buoyancy field, not to be confused with tracer particles
#include "diffusion.h"
#include "view.h"

// Geometry
void geometry (double H, double L, double Ry) {
  vertex scalar phi[];
  foreach_vertex() {
    phi[] = H/2.- fabs(y);//channel height
    phi[] = max(phi[], fabs(x) - L/2.);//channel length
    phi[] = min(phi[], Ry/2. - fabs(y));//top and bottom of tanks
    phi[] = min(phi[], H/2. - fabs(z));//expand to 3D
  }
  fractions (phi, cs, fs);
}

double H = 2., L = 60., Rx = 4.1, Ry = 6.1, incline = 4*pi/180;
double forc_b = 20., forc_u = 0.8;

//Mesh
int minlevel = 4, maxlevel = 12;
double btolerance = 0.03, vtolerance = 0.015;

//buoyancy
scalar b[], * tracers = {b};//buoyancy field
#define STRAT(x) (-b0*tanh(k*x)*(cs[] > 0))
double k = 5, b0 = 1;

//flow
double Re = 800., Sc = 7.;
face vector muc[], av[], D[];
scalar rhof[];
double rhomean = 0.9985945, deltarho = 3.7e-5;
int t_end = 100;

u.n[embed] = dirichlet (0);
u.t[embed] = dirichlet (0);
u.n[left] = dirichlet (0);
u.t[left] = dirichlet (0);
u.n[right] = dirichlet (0);
u.t[right] = dirichlet (0);
#if dimension == 3
u.r[embed]=dirichlet(0);
u.r[left] = dirichlet (0);
u.r[right] = dirichlet (0);
#endif

int main() {
  L0 = L + 2*Rx;
  X0 = Y0 = Z0 = -L0/2;
  N = 1 << 7;
  mu = muc;
  a = av;
  DT = 0.1;
  run();
}

event init (t = 0) {
  astats as = {.nf = 9999, .nc = 9999};//nf=#refined cells, nc=#coarsened cells
  while (as.nf > 10 || as.nc > 10) {//continue refining the mesh until only few cells (10) get changed 
    geometry(H, L, Ry);
    foreach(){
      b[] = STRAT(x);
    }
    as = adapt_wavelet({cs, b, u}, (double[]){1e-2, btolerance, vtolerance, vtolerance, vtolerance}, maxlevel, minlevel);//{csmax, bmax, uxmax, uymax}
  }
  b.gradient = minmod2;//Slope limiter
  srand(time(NULL));//seed for random number generator to define particle position
}

event log_perf(i++){
  if (i == 0){
    fprintf (stdout, "t dt mgp.i mgpf.i mgu.i grid->tn perf.t perf.speed\n");
  }
  fprintf (stdout, "%.2f %g %d %d %d %ld %g %g\n", t, dt, mgp.i, mgpf.i, mgu.i, grid->tn, perf.t, perf.speed);
  fflush(stdout);
}

event log_gradb(i++){
  if (i == 0){
    fprintf (stderr, "t vel gradb\n");
  }
  vector gradb[];
  gradients({b}, {gradb});
  double u_tot_sum = 0, grad_tot_sum = 0, N_sum = 0;
  double Nx = L*pow(2, maxlevel)/L0, Nyz = H*pow(2, maxlevel)/L0;//interpolate with resolution based on mesh max level
  coord p;
  coord box[2] = {{-L/2., -H/2., -H/2.}, {L/2., H/2., H/2.}};
  coord n = {round(Nx), round(Nyz), round(Nyz)};
  foreach_region(p, box, n, reduction(+:u_tot_sum) reduction(+:grad_tot_sum) reduction(+:N_sum)){
    u_tot_sum += sqrt(sq(u.x[]) + sq(u.y[]) + sq(u.z[]));
    grad_tot_sum += sqrt(sq(gradb.x[]) + sq(gradb.y[]) + sq(gradb.z[]));
    N_sum += 1;
  }
  fprintf(stderr, "%f %f %f\n", t, u_tot_sum/N_sum, grad_tot_sum/N_sum);
}

event movie (t += 0.2) {
  char s[80];
  sprintf (s, "t = %.1f", t);
  double a = 80./pow(2., 12.);
  char m[80];

  view(fov = 5);//fov = 5 only center of channel visible
  draw_string (s, pos = 1, size = 80, lc = {1, 1, 1}, lw = 2);//depict time in movie
  squares("level", min = minlevel, max = maxlevel, alpha = -a);
  colorbar(min = minlevel, max = maxlevel);
  sprintf(m, "mov_level_ml%i.mp4", maxlevel);
  save(m);
  if (fabs(remainder(t, 5)) < 0.05){//save snapshot of b every 5 timesteps
    sprintf(m, "snapshots/snapshot_level_ml%i_t%.0f.png", maxlevel, t);
    save(m);  
  }

  view(fov = 17);//fov = 17 channel and small part of tanks visible
  draw_string (s, pos = 1, size = 80, lc = {1, 1, 1}, lw = 2);
  squares("b", min = - 1.1, max = + 1.1, map = jet, alpha = -a);
  sprintf(m, "mov_b_ml%i.mp4", maxlevel);
  save(m);
  if (fabs(remainder(t, 5)) < 0.05){
    sprintf(m, "snapshots/snapshot_b_ml%i_t%.0f.png", maxlevel, t);
    save(m);  
  }
}

event forcing (i++, last) {
  double x_forcing = (L0 + L)/4.;
  foreach() {
    if (fabs(x) > x_forcing){
      b[] -= forc_b*(fabs(x)-x_forcing)/(L0/2 - x_forcing)*dt*(b[] - STRAT(x));//forcing gets stronger near the end of the reservoir
      foreach_dimension(){
        u.x[] -= forc_u*(fabs(x)-x_forcing)/(L0/2 - x_forcing)*dt*u.x[];
      }
    } 
  }
}

event properties (i++) {
  foreach_face(){
    muc.x[] = fm.x[] * 1./Re;//visocisity has to me defined each timestep as mesh is adaptive (every timestep new cells can be made)
  }
}

event acceleration (i++) { //Buoyancy
  coord dir = {sin(incline), cos(incline), 0};
  double Gprime = 0.25;//dimensionless reduced gravity (contains density difference), equal to Richardson number
  foreach_face()
    av.x[] = Gprime * dir.x*(b[] + b[-1])/2.;//b is cell-centered, acceleration is defined on cell faces.
}

event tracer_diffusion (i++) {
  foreach_face()
    D.x[] = fm.x[] * 1./(Re*Sc);
  diffusion (b, dt, D);
}

event adapt (i++) {
  adapt_wavelet({cs, b, u}, (double[]){1e-2, btolerance, vtolerance, vtolerance, vtolerance}, maxlevel, minlevel);
}

event dump (t = 50; t += 10){
  char name[80];
  output_field(all);
  dump(file = "restart_th4_Re800");
}

event stop (t = t_end);
