/**
# Particles in Stratified Inclined Duct (SID) experiment

Add Stokes particles to SID using two particle densities and two release locations per density. Boussinesq approximation is used to compute flow velocity. Navier-Stokes solver uses field rho while particles use field rhof for the local fluid densities. This is to ensure that the particles can still `feel' the local fluid density without breaking the Boussinesq approximation. The adapted particle code can be found here [here](particle_rhof.h)
*/

#include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#define ADD_PART_MEM coord u; coord u2; long unsigned int tag; long unsigned int group_tag; coord start; int del;
#include "stokes-particles.h"
#include "tracer.h"
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
}

double H = 2., L = 60., Rx = 4.1, Ry = 6.1, incline = 4*pi/180;
double forc_b = 20., forc_u = 0.8;

//Mesh
int minlevel = 4, maxlevel = 12;
double btolerance = 0.02, vtolerance = 0.01;

//buoyancy
scalar b[], * tracers = {b};//buoyancy field
#define STRAT(x) (-b0*tanh(k*x)*(cs[] > 0))
double k = 5, b0 = 1;

//flow
double Re = 800., Sc = 7.;
face vector muc[], av[], D[];
scalar rhof[];
double rhomean = 0.9985945, deltarho = 3.7e-5;

//particles
double radius = 2e-2;
Particles neutral_center, neutral_tanks, negative_center, negative_wall;
int Nparticles = 10000, t_release = 101, t_measure = 130, t_end = 140;//start measuring neutral_tanks particles at t=130

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
  G.x = - rhomean/(2 * deltarho) * sin(incline);
  G.y = - rhomean/(2 * deltarho) * cos(incline);
  G.z = 0.;
  a = av;
  DT = 0.1;
  run();
}

event init (t = 0) {
  neutral_center = new_inertial_particles(0);
  neutral_tanks = new_inertial_particles(0);
  negative_center = new_inertial_particles(0);
  negative_wall = new_inertial_particles(0);
  if (!restore (file = "flow_t100/restart_th4_Re800")){
    astats as = {.nf = 9999, .nc = 9999};//nf=#refined cells, nc=#coarsened cells
    while (as.nf > 10 || as.nc > 10) {//continue refining the mesh until only few cells (10) get changed 
      geometry(H, L, Ry);
      foreach(){
        b[] = STRAT(x);
        rhof[] = 1 - deltarho / (2 * rhomean) * b[];
      }
      as = adapt_wavelet({cs, b, u}, (double[]){1e-2, btolerance, vtolerance, vtolerance, vtolerance}, maxlevel, minlevel);//{csmax, bmax, uxmax, uymax}
    }
    b.gradient = minmod2;//Slope limiter
    srand(time(NULL));//seed for random number generator to define particle position
  }
  else {//restart
    geometry(H, L, Ry);//face fraction field is not included in dumpfile
  }
}

event add_particles(t = t_release){
  int group_tag = 0;
  foreach_P_in_list(inertial_particles){
    for (int tag = 0; tag < Nparticles; tag += 1){
      particle p;
      p.u2.y = radius;
      p.del = 0;
      p.tag = group_tag*Nparticles + tag;//combine group tag and tag to get unique value
      p.group_tag = group_tag;
      p.z = (((double)rand())/((double)RAND_MAX) -0.5) * 0.9 * H;//starting depth random between -0.9*H/2 and 0.9*H/2 so over most of the channel depth
      if (group_tag == 0){//neutral center
        p.u2.x = 1;//mass density particle = average density fluid
        p.x = 0.;//center
        p.y = 0.;
        foreach_dimension(){
          p.u.x = 0.;//set particle velocity initially at zero
          p.start.x = p.x;//save starting position
        }
      }
      else if (group_tag == 1){//neutral tanks
        p.u2.x = 1;//mass density particle = average density fluid
        double x_sign;
        if (tag <= Nparticles/2){//plot half of them in left and half in right tank
          x_sign = + 1.;
        } else{
          x_sign = -1.;
        }
        //initialize particles in the tanks, while keeping distance from the tank walls
        p.x = x_sign * (L/2. + Rx/2);
        p.y = - x_sign * 0.9 * H/2.;
        foreach_dimension(){
          p.u.x = 0.;//set particle velocity initially at zero
          p.start.x = p.x;//save starting position
        }
      }
      else if (group_tag == 2){//negative center
        p.u2.x = 1 + deltarho / (2 * rhomean);//mass density particle = density heavier fluid
        p.x = 0.;
        p.y = 0.;
        foreach_dimension(){
          p.u.x = 0.;//set particle velocity initially at zero
          p.start.x = p.x;//save starting position
        }
      }
      else {
        p.u2.x = 1 + deltarho / (2 * rhomean);//mass density particle = density heavier fluid
        p.x = -L/4.;//starting x on 1/4 of channel length
        p.y = (H/2. - 4*radius) ;//starting y a small distance (two times particle diameter) from channel wall
        foreach_dimension(){
          p.u.x = 0.;//set particle velocity initially at zero
          p.start.x = p.x;//save starting position
        }
      }
      add_particle(p, P);
    }
    group_tag += 1;
  }
}

event log_perf(i++){
  if (i == 0){
    fprintf (stdout, "t dt mgp.i mgpf.i mgu.i grid->tn perf.t perf.speed\n");
  }
  fprintf (stdout, "#%.3f %g %d %d %d %ld %g %g\n", t, dt, mgp.i, mgpf.i, mgu.i, grid->tn, perf.t, perf.speed);
  fflush(stdout);
}

event log_particles_pos(t = t_release, i += 2){
  foreach_particle(){
    fprintf(stdout, "%.3f %li %li %f %f %f\n", t, p().group_tag, p().tag, p().x, p().y, p().z);
  }
}

event log_particles_MSD(t = t_release, i += 2){
  if (t == t_release){
    fprintf(stderr, "%s %i %s %i %s %.3f %s %.3f %s %.0f %s %.0f %s %.1f %.1f %s %i\n", "#dimension is", dimension, "Nparticles", Nparticles, "particle radius is", radius, "incline", incline, "Re", Re, "Sc", Sc, "tank Ry x Rx", Ry, Rx, "maxlevel", maxlevel);
    fprintf(stderr, "%s\n", "t, (RMSD for x, y, z, Nparticles) for all 4 particle groups (neutral_center, neutral_tanks, negative_center, negative_wall)");
  }
  fprintf(stderr, "%f ", t);//print time
  double COM[4][4];//save center of mass for all particle lists for x, y, z
  memset (COM, 0, sizeof(COM));
  int index = 0;
  foreach_P_in_list(inertial_particles){
    foreach_particle_in(P, reduction(+:COM), reduction(+:index)){
      if (p().del == 0){
        if (index == 1){
          COM[index][0] += fabs(p().x);//group 1 (neutral tanks) has two symmetric clusters
          COM[index][1] += fabs(p().y);//so take absolute value in x and y direction
        } else {
          COM[index][0] += p().x;
          COM[index][1] += p().y;
        }
        COM[index][2] += p().z;
        COM[index][3] += 1;
      }
    }
    index += 1;
  }
  //compute mean squared displacemet (MSD) relative to COM
  double RMSD[4][4];
  memset (RMSD, 0, sizeof(RMSD));
  index = 0;
  foreach_P_in_list(inertial_particles){
    foreach_particle_in(P, reduction(+:RMSD) reduction(+:COM) reduction(+:index)){
      if (p().del == 0){
        if (index == 1){
          RMSD[index][0] += sq(fabs(p().x) - COM[index][0]/COM[index][3]);
          RMSD[index][1] += sq(fabs(p().y) - COM[index][1]/COM[index][3]);
        } else {
          RMSD[index][0] += sq(p().x - COM[index][0]/COM[index][3]);
          RMSD[index][1] += sq(p().y - COM[index][1]/COM[index][3]);
        }
        RMSD[index][2] += sq(p().z - COM[index][2]/COM[index][3]);
        RMSD[index][3] += 1;
      }
    }
    index += 1;
  }
  //print results
  for (int index_group = 0; index_group < 4; index_group ++){
    for (int d = 0; d < 3; d ++){
      fprintf(stderr, "%f ", COM[index_group][d]/COM[index_group][3]);//print center of mass in x, y, z by taking average
    }
    for (int d = 0; d < 3; d ++){
      fprintf(stderr, "%f ", RMSD[index_group][d]/RMSD[index_group][3]);//print RMSD in x, y, z by taking average
    }
    fprintf(stderr, "%.0f ", RMSD[index_group][3]);//print number of particles per group
  }
}

event log_gradb(i += 2){//print average vel and grad b in channel to check mesh convergence
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
  if (t < t_release){
    fprintf(stderr, "%f ", t);//print time
  }
  fprintf(stderr, "%f %f\n", u_tot_sum/N_sum, grad_tot_sum/N_sum);
}

event delete_before_measuring(t = t_measure){
  foreach_particle_in(neutral_tanks){//delete particles that have not reached the center region of the channel before the measurement starts
    if (fabs(p().x) > L*3/8.){
      p().del = 1;
    }
  }
}

event particles_boundaries(i++){
  foreach_particle(){
    if (fabs(p().z) > H/2 || (fabs(p().x) < L/2. && fabs(p().y) > H/2.) || (fabs(p().x) > L/2. && fabs(p().y) > Ry/2.)){//if particle is outside domain
      p().del = 1;
      fprintf(stderr, "%s %.3f %li %f %f %f\n", "#outside", t, p().group_tag, p().x, p().y, p().z);
    }

    foreach_point(p().x, p().y, p().z, serial){//if particles are in contact with channel wall
      if (cs[] < 0.){
        p().del = 1;
        fprintf(stderr, "%s %.3f %li %f %f %f\n", "#cs<0", t, p().group_tag, p().x, p().y, p().z);
      }      
      else if (cs[] > 0. && cs[] < 1.){
        coord b;//center of boundary relative to cell center
        coord n;//normal vector of boundary (points from fluid to solid)
        embed_geometry(point, &b, &n);
        double distancey = fabs(y + b.y*Delta) - fabs(p().y);
        double distancez = fabs(z + b.z*Delta) - fabs(p().z);
        if (distancey <= p().u2.y || (dimension == 3 && distancez <= p().u2.y)){//if distance is smaller than particle radius (also works if distance < 0 bc particle is outside boundary)
          p().del = 1;
          fprintf(stderr, "%s %.3f %li %li %f %f %f\n", "#wall", t, p().group_tag, p().tag, p().x, p().y, p().z);
        }
      }
    }

    if (fabs(p().x) > L/2. && p().x * p().start.x < 0){//if particles reached the (opposite) tank
      p().del = 1;
      fprintf(stderr, "%s %.3f %li %f %f %f\n", "#tank", t, p().group_tag, p().x, p().y, p().z);
    }
  }
}

event particles_delete(i++){
  remove_particles(neutral_center, p().del > 0);
  remove_particles(neutral_tanks, p().del > 0);
  remove_particles(negative_center, p().del > 0);
  remove_particles(negative_wall, p().del > 0);
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
  sprintf(m, "videos/mov_level_ml%i.mp4", maxlevel);
  save(m);
  if (fabs(remainder(t, 5)) < 0.05){//save snapshot of b every 5 timesteps
    sprintf(m, "snapshots/snapshot_level_zoom_ml%i_t%.0f.png", maxlevel, t);
    save(m);  
  }

  view(fov = 5);//fov = 17 channel and small part of tanks visible
  draw_string (s, pos = 1, size = 80, lc = {1, 1, 1}, lw = 2);
  squares("b", min = - 1.1, max = + 1.1, map = jet, alpha = -a);
  sprintf(m, "videos/mov_b_zoom_ml%i.mp4", maxlevel);
  save(m);
  if (fabs(remainder(t, 5)) < 0.05){
    sprintf(m, "snapshots/snapshot_b_zoom_ml%i_t%.0f.png", maxlevel, t);
    save(m);  
  }

  view(fov = 17);//fov = 17 channel and small part of tanks visible
  draw_string (s, pos = 1, size = 80, lc = {1, 1, 1}, lw = 2);
  squares("b", min = - 1.1, max = + 1.1, map = jet, alpha = -a);
  sprintf(m, "videos/mov_b_ml%i.mp4", maxlevel);
  save(m);
  if (fabs(remainder(t, 5)) < 0.05){
    sprintf(m, "snapshots/snapshot_b_ml%i_t%.0f.png", maxlevel, t);
    save(m);  
  }

  int group_tag = 0;
  foreach_P_in_list(inertial_particles){
    view(fov = 17);
    draw_string (s, pos = 1, size = 80, lc = {1, 1, 1}, lw = 2);
    squares("b", min = - 1.1, max = + 1.1, map = jet, n = {0, 0, 1}, alpha = -a);
    scatter(P, s = 1.5, pc = {0, 0, 0});
    sprintf(m, "group_%i/mov_particles_ml%i.mp4", group_tag, maxlevel);
    save(m);
    if (fabs(remainder(t, 5)) < 0.05){
      sprintf(m, "group_%i/snapshot_particles_ml%i_t%.0f.png", group_tag, maxlevel, t);
      save(m);  
    }

    view(fov = 5);
    draw_string (s, pos = 1, size = 80, lc = {1, 1, 1}, lw = 2);
    squares("b", min = - 1.1, max = + 1.1, map = jet, n = {0, 0, 1}, alpha = -a);
    scatter(P, s = 2, pc = {0, 0, 0});
    sprintf(m, "group_%i/mov_particles_zoom_ml%i.mp4", group_tag, maxlevel);
    save(m);
    if (fabs(remainder(t, 5)) < 0.05){
      sprintf(m, "group_%i/snapshot_particles_zoom_ml%i_t%.0f.png", group_tag, maxlevel, t);
      save(m);  
    }
    group_tag += 1;
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
  foreach(){
    rhof[] = 1 - deltarho / (2 * rhomean) * b[];
  }
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

event dump (t = 110, t += 5){
  char name[80];
  output_field(all);
  dump(file = "restart_flow");
  pdump(fname = "restart_particles");
}

event stop (t = t_end);
