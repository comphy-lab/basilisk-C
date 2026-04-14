/**
# Mesh study using particles

Use both flow and particle statistics to study the mesh. Particles are deleted if they hit the channel wall or reach the tanks.
*/

#include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#define ADD_PART_MEM coord u; coord u2; long unsigned int tag; coord start; double start_b; double start_t; double end_b; int del; int cat;
#include "stokes-particles.h"
#include "tracer.h"//for buoyancy field, not to be confused with tracer particles
#include "diffusion.h"
#include "view.h"
#include "scatter2.h"

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
  //output_ppm (phi, file = "phi.png");
}
double H = 2., L = 60., Ry = 6.1, Rx = 9.1, incline = 2*pi/180;

//Mesh
int minlevel = 4, maxlevel = 10;
double btolerance = 6e-2, vtolerance = 3e-2;

//buoyancy
scalar b[], * tracers = {b};//buoyancy field
#define STRAT(x) (-b0*tanh(k*x)*(cs[] > 0))
double k = 5, b0 = 1;

//flow
double Re = 650., Sc = 7.;
face vector muc[], av[], D[];
scalar rhof[];
double rhomean = 1.0125, rhodelta = 0.0125;

//particles
double radius = 2e-2;
Particles center;
int Nparticles = 500;
int t_end = 150;

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
  rho = rhof;
  G.x = -9.81*sin(incline);
  G.y = -9.81*cos(incline);
  G.z = 0.;
  a = av;
  DT = 0.1;
  run();
}

event init (t = 0) {
  center = new_inertial_particles(0);
  astats as = {.nf = 9999, .nc = 9999};//nf=#refined cells, nc=#coarsened cells
  while (as.nf > 10 || as.nc > 10) {//continue refining the mesh until only few cells (10) get changed 
    geometry(H, L, Ry);
    foreach(){
      b[] = STRAT(x);
      rhof[] = rhomean + rhodelta * - b[];
    }
    as = adapt_wavelet({cs, b, u}, (double[]){1e-2, btolerance, vtolerance, vtolerance, vtolerance}, maxlevel, minlevel);//{csmax, bmax, uxmax, uymax, uzmax}
  }
  b.gradient = minmod2;//Slope limiter
  fprintf(stderr, "%s %i %s %i %s %.3f %s %.3f %s %.0f %s %.0f %s %.1f %.1f %s %i\n", "#dimension is", dimension, "Nparticles", Nparticles, "particle radius is", radius, "incline", incline, "Re", Re, "Sc", Sc, "tank Ry x Rx", Ry, Rx, "maxlevel", maxlevel);
  srand(time(NULL));//seed for random number generator to define particle radius
}
/**
Add particles as soon as exchange flow has reached the tanks at t=70
*/
event particles_add(t = 70){//
  for (int tag = 1; tag <= Nparticles; tag += 1){
    particle p;
    foreach_dimension(){
      p.x = H * (((double)rand())/((double)RAND_MAX) -0.5);
      p.u.x = 0.;//set particle velocity initially at zero
    }
    p.u2.x = 1.0125;//mass density particle = average density fluid
    p.u2.y = radius;
    p.del = 0;
    p.tag = tag;
    add_particle(p, center);
  }
}
/**
Check if particles hit the channel wall or reached the tanks.
*/
event particles_boundaries(i++){
  foreach_particle(){
    foreach_point(p().x, p().y, p().z, serial){//check if particles crossed or are in contact with channel wall
      if (cs[] > 0. && cs[] < 1.){
        coord b;//center of boundary relative to cell center
        coord n;//normal vector of boundary (points from fluid to solid)
        embed_geometry(point, &b, &n);
        double distancey = fabs(y + b.y*Delta) - fabs(p().y);
        double distancez = fabs(z + b.z*Delta) - fabs(p().z);
        if (distancey <= p().u2.y || (dimension == 3 && distancez <= p().u2.y)){//if distance is smaller than particle radius (also works if distance < 0 bc particle is outside boundary)
          p().del = 1;
          p().cat = 1;
        }
      }
    }

    if (fabs(p().x) > L/2.){//check if particles reached the tank
      p().del = 1;
      p().cat = 2;
    }
  }
}

event particles_delete(i++){
  remove_particles(center, p().del > 0);
}

event log_perf(i++){
  if (i == 0){
    fprintf (stdout, "t dt mgp.i mgpf.i mgu.i grid->tn perf.t perf.speed\n");
  }
  fprintf (stdout, "%.2f %g %d %d %d %ld %g %g\n", t, dt, mgp.i, mgpf.i, mgu.i, grid->tn, perf.t, perf.speed);
  fflush(stdout);
}
/**
Log particle statistics: average absolute particle buoyancy (buoyancy of the fluid at the position of the particle)
*/
event log_particles(i++){
  if (i == 0){
    fprintf(stderr, "%s\n", "t, average absolute buoyancy, Nparticles, mag u, mag grad b, N");
  }
  double absb_sum = 0, N_sum = 0;
  foreach_particle(){
    foreach_point(p().x, p().y, p().z, reduction(+:absb_sum) reduction(+:N_sum)){
      double part_b = interpolate_linear(point, b, p().x, p().y, p().z);//calculate particle buoyancy
      absb_sum += fabs(part_b);
      N_sum += 1;
    }
  }
  fprintf(stderr, "%.3f %f %.0f ", t, absb_sum/N_sum, N_sum);
}
/**
Log flow statistics: average velocity inside the duct and the average magnitude of the gradient of b
*/
event log_gradb(i++){
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
  fprintf(stderr, "%f %f %.0f\n", u_tot_sum/N_sum, grad_tot_sum/N_sum, N_sum);
}
/**
Make zoomed-in movie of mesh level, zoomed-out movies of only the flow and flow with particles plotted.
*/
event movie (t += 0.2) {
  view(fov = 5);
  char s[80];
  sprintf (s, "t = %.1f", t);
  draw_string (s, pos = 1, size = 80, lc = {1, 1, 1}, lw = 2);//depict time in movie
  #if dimension == 2
  draw_vof("cs", "fs", filled = -1, lw = 2);
  #endif
  double a = 80./pow(2., 12.);
  squares("level", min = minlevel, max = maxlevel, alpha = a);
  colorbar(min = minlevel, max = maxlevel);
  char m[80];
  sprintf(m, "snapshot_level_ml%i.png", maxlevel);
  save(m);
  sprintf(m, "mov_level_ml%i.mp4", maxlevel);
  save(m);

  view(fov = 15);
  draw_string (s, pos = 1, size = 80, lc = {1, 1, 1}, lw = 2);//depict time in movie
  #if dimension == 2
  draw_vof("cs", "fs", filled = -1, lw = 2);
  #endif
  squares("b", min = - 1.1, max = + 1.1, map = jet, alpha = a);
  sprintf(m, "snapshot_flow_ml%i.png", maxlevel);
  save(m);
  sprintf(m, "mov_flow_ml%i.mp4", maxlevel);
  save(m);

  view(fov = 15);
  draw_string (s, pos = 1, size = 80, lc = {1, 1, 1}, lw = 2);//depict time in movie
  #if dimension == 2
  draw_vof("cs", "fs", filled = -1, lw = 2);
  #endif
  squares("b", min = - 1.1, max = + 1.1, map = jet, alpha = a);
  scatter(center, s = 3, pc = {1, 1, 1});
  sprintf(m, "snapshot_particles_ml%i.png", maxlevel);
  save(m);
  sprintf(m, "mov_particles_ml%i.mp4", maxlevel);
  save(m);
}

/**
Forcing is applied to be able to use smaller tanks while keeping the density inside the tanks constant. Inspired by work done by Zhu et all^[Zhu, Lu, Amir Atoufi, Adrien Lefauve, et al. “Stratified Inclined Duct: Direct Numerical Simulations.” Journal of Fluid Mechanics 969 (August 2023): A20. https://doi.org/10.1017/jfm.2023.502.]

*/
event forcing (i++, last) {
  double x_forcing = (L0 + L)/4.;
  foreach() {
    if (fabs(x) > x_forcing){
      b[] -= 10*(fabs(x)-x_forcing)/(L0/2 - x_forcing)*dt*(b[] - STRAT(x));//forcing gets stronger near the end of the reservoir
      foreach_dimension(){
        u.x[] -= 0.4*(fabs(x)-x_forcing)/(L0/2 - x_forcing)*dt*u.x[];
      }
    } 
  }
}
/**
The dimensional density of the fluid is changed based on the (dimensionless) buoyancy field. This is because the particles have a dimensional density.
*/
event properties (i++) {
  foreach(){
    rhof[] = rhomean + rhodelta * - b[];
  }
  foreach_face(){
    muc.x[] = fm.x[]*1./Re;//visocisity has to me defined each timestep as mesh is adaptive (every timestep new cells can be made)
  }
}

event acceleration (i++) { //Buoyancy
  coord dir = {sin(incline), cos(incline), 0};
  double Gprime = 0.25;//dimensionless reduced gravity, equal to Richardson number
  foreach_face()
    av.x[] = Gprime * dir.x*(b[] + b[-1])/2.;//b is cell-centered, acceleration is defined on cell faces.
}

event tracer_diffusion (i++) {
  foreach_face()
    D.x[] = fm.x[]*1./(Re*Sc);
  diffusion (b, dt, D);
}

event adapt (i++) {
  adapt_wavelet({cs, b, u}, (double[]){1e-2, btolerance, vtolerance, vtolerance, vtolerance}, maxlevel, minlevel);
}

event stop (t = t_end);