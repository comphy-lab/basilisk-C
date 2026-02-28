#include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "tracer.h"//for buoyancy field, not to be confused with tracer particles
#include "diffusion.h"
#include "view.h"

// Geometry
void geometry (double H, double L, double Hr) {
  vertex scalar phi[];
  foreach_vertex() {
    phi[] = H/2.- fabs(y);//channel height
    phi[] = max(phi[], fabs(x) - L/2.);//channel length
    phi[] = min(phi[], Hr/2. - fabs(y));//top and bottom of tanks
    phi[] = min(phi[], H/2. - fabs(z));//expand to 3D
  }
  fractions(phi, cs, fs);
  //output_ppm(phi, file = "phi.png");
}

double H = 2., L = 60., Hr = 16.1, Lr = 10.1, incline = 2*pi/180;

//Mesh
int minlevel = 4, maxlevel = 12;
double btolerance = 4e-2, vtolerance = 2e-2;

//buoyancy
scalar b[], * tracers = {b};//buoyancy field
#define STRAT(x) (-b0*tanh(k*x)*(cs[] > 0))
double k = 5, b0 = 1.;

//flow
double Re = 650., Sc = 7.;
face vector muc[], av[], D[];

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
  L0 = L + 2*Lr;
  X0 = Y0 = Z0 = -L0/2;
  N = 1 << 7;
  mu = muc;
  a = av;
  DT = 0.1;
  run();
}

event init (t = 0) {
  astats as = {.nf = 9999, .nc = 9999};//nf = #refined cells, nc = #coarsened cells
  while (as.nf > 20 || as.nc > 20) {//continue refining the mesh until only few cells (20) get changed 
    geometry(H, L, Hr);
    foreach(){
      b[] = STRAT(x);
    }
    as = adapt_wavelet({cs, b, u}, (double[]){1e-2, btolerance, vtolerance, vtolerance}, maxlevel, minlevel);//{csmax, bmax, uxmax, uymax}
  }
  
  b.gradient=minmod2;//Slope limiter
  fprintf(stderr, "%s %i %s %.3f %s %.0f %s %.0f %s %.1f %.1f %s %i\n", "dimension is", dimension, "incline", incline, "Re", Re, "Sc", Sc, "tank Hr x Lr", Hr, Lr, "maxlevel", maxlevel);
  fprintf(stderr, "%s %s %s %s\n", "x", "time", "froude", "eta");
  //fprintf(stderr, "%s\n", "time, mean(abs(u.y)), N_sum");
}

event log_froude(t += 0.1){
  double dx = 0.1, dyz = 0.01;
  for (double i = -L/2.; i <= L/2.; i += dx){//iterate in x-direction (horizontal)
    double froude = 0.;
    double u1 = 0., u2 = 0.;
    double b1 = 0., b2 = 0.;
    double eta = 0.;
    double count1 = 0., count2 = 0., count_eta = 0.;
    //compute interface heigth
    for (double j = -H/2.; j <= H/2.; j += dyz){//iterate in y-direction (into plane)
      double eta_y = copysign(H/2., i);//starting value eta is +H/2 for x<0 and -H/2 for x>0 
      for (double k = -H/2.; k <= H/2.; k += dyz){//iterate in z-direction (vertical)
        foreach_point(i, k, j, serial){
          if (i < 0 && b[] < 0.){
            eta_y = max(eta_y, k);
          }
          if (i > 0 && b[] > 0.){
            eta_y = min(eta_y, k);//save interface heigth
          }
        }
      }
      eta += eta_y;//save y-dependent eta value
      count_eta += 1.;
    }
    //compute layer-averaged velocity and buoyancy
    for (double j = -H/2.; j <= H/2.; j += dyz){//iterate in y-direction (into plane)
      for (double k = -H/2.; k <= H/2.; k += dyz){//iterate in z-direction (vertical)
        foreach_point(i, k, j, serial){
          if (k >= eta/count_eta){
            u1 += u.x[];
            b1 += b[];
            count1 += 1;
          }
          else {
            u2 += u.x[];
            b2 += b[];
            count2 += 1;
          }
        }
      }
    }
    if (count1 == 0){
      count1 += 1;
    }
    if (count2 == 0){
      count2 += 1;
    }
    double deltarho = fabs(b1/count1 - b2/count2);
    double Ri = 0.25;
    froude = sq(u2/count2 - u1/count1)/(deltarho*Ri*cos(incline)*H);
    fprintf(stderr, "%.1f %.1f %f %f\n", i, t, froude, eta/count_eta);
  }
}

event log_vel(i++){
  double v_sum = 0., count = 0.;
  foreach(serial){
    v_sum += fabs(u.y[]);
    count += 1;
  }
  //fprintf(stderr, "%.1f %f %.0f\n", t, v_sum/count, count);
}

event movie (t += 0.2) {
  view(fov = 5);//fov = 15 channel and small part of tanks visible with Hr, Lr = 8, 10, fov = 12 when Hr, Lr = 16, 20
  char s[80];
  sprintf (s, "t = %.2f", t);
  draw_string (s, pos = 1, size = 80, lc = {1, 1, 1}, lw = 2);//depict time in movie
  #if dimension == 2
  draw_vof("cs", "fs", filled = -1, lw = 2);
  #endif
  double a = 80./pow(2., 12.);
  squares("level", min = minlevel, max = maxlevel, alpha = a);
  colorbar(min = minlevel, max = maxlevel);
  save("mov_level.mp4");
  save("snapshot_level.png");
  view(fov = 15);//fov = 15 channel and small part of tanks visible with Hr, Lr = 8, 10, fov = 12 when Hr, Lr = 16, 20
  draw_string (s, pos = 1, size = 80, lc = {1, 1, 1}, lw = 2);//depict time in movie
  #if dimension == 2
  draw_vof("cs", "fs", filled = -1, lw = 2);
  #endif
  squares("b", min = - 1.1, max = + 1.1, map = jet, alpha = a);
  save("mov.mp4");
  save("snapshot.png");
}

event dampen (i++, last) {
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

event properties (i++) {
  foreach_face(){
    muc.x[] = fm.x[]*1./Re;//visocisity has to me defined each timestep as mesh is adaptive (every timestep new cells can be made)
  }
}

event acceleration (i++) { //Buoyancy
  coord dir = {sin(incline), cos(incline), 0};
  double Gprime = 0.25;//reduced gravity (contains density difference), equal to Richardson number
  foreach_face(){
    av.x[] = Gprime * dir.x * (b[] + b[-1])/2.;//b is cell-centered, acceleration is defined on cell faces.
  }
}

event tracer_diffusion (i++) {
  foreach_face()
    D.x[]=fm.x[]*1./(Re*Sc);
  diffusion (b, dt, D);
}

event adapt (t = 0, i++) {
  adapt_wavelet({cs, b, u}, (double[]){1e-2, btolerance, vtolerance, vtolerance}, maxlevel, minlevel);
}

event dump (t = 100; t <= t_end; t += 50){
  char name[80];
  sprintf (name, "dump_files/dump-%.1f", (double) t);
  dump(file = name);
}

event stop (t = t_end);
