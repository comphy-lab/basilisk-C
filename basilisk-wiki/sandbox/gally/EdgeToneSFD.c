/**
# Edgetone's baseflow calculated with the Selective Frequency method.
Original code from [Lourenco](https://basilisk.fr/sandbox/Llourenco/)
*/

#include "embed.h"
#include "navier-stokes/centered.h"
#include "view.h"
#include "SFD.h"

double Reynolds = 150; 
int maxlevel = 10;
double b = 1.;   // jet diameter
double U0 = 1.;  // for top hat
double W = 5.;   // edge-distance
double freq = 0.067;
double OMEG ;
double eps = 1e-3;
double delta_erf;
double h, R;
double Tramp;
double tau;
double ramp;

double freq_SFD = 0.1;

bool SFD_toggle;
event adapt_toggle (i++)
  SFD_toggle = (t >= 0.);

scalar omega[];
scalar base[];
scalar un[];
face vector muv[];
double cf = 0.; 
double du;
double Tend = 2000;

FILE * fp;

int main(int argc, char * argv[]) {

  display_control(Reynolds,0,1000);
  display_control(cf,0,1);
  display_control(maxlevel,0,15);

  cf = 0.;
  L0 = 50; 
  origin (0., -L0/2.);
  N =512; 
  mu = muv;
  delta_erf = b/10.;
  OMEG = 2*M_PI*freq;
  Tramp = 50;
  tau = 0.5;
 
  init_grid (N);

  fp = fopen("data.dat","w");
    
  dt = 0.1;
  DT = 0.1;
  
  run();
}

event properties (i++) {
  foreach_face()
    muv.x[] = fm.x[]*(b*U0)/Reynolds;
  vorticity (u, omega);
  refine(x > 0.8*W && fabs(y) < h && level < maxlevel);
}

u.n[left] = dirichlet( (min(t/Tramp,1.))*(cf + (1-cf)*0.5*(erf((y+b/2.)/delta_erf)-erf((y-b/2.)/delta_erf))));
u.t[left] = dirichlet(min(t/Tramp,1.));

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

u.n[top] = dirichlet(0.);
u.t[top] = dirichlet(0.);

u.n[bottom] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);

u.n[embed] = dirichlet(0.); 
u.t[embed] = dirichlet(0.);

event init (t = 0) {
  h=0.1;
  R=h/2;

  refine(x > 0.8*W && fabs(y) < h && level < maxlevel);
  solid (cs, fs,
         min(
             max(max(W - x, x - L0),
                 max(y - h/2., -y - h/2.)),

             max(sq(x - W) + sq(y) - sq(R),
                 x - W)
         ));
}

event logfile (i++, t<=Tend) {
  vorticity(u,omega);

  fprintf (stderr, "%d %g %g %g %g %g %g %g %g\n", 
           i, t, dt, 
           interpolate(u.y, 0.1*W, 0.), 
           interpolate(u.y, 0.9*W, 0.),
           interpolate(u.x, 0.05*W, 0.),
           interpolate(p, 0.9*W, W/2),
           interpolate(p, 0.9*W, -W/2),
           interpolate(omega, 0.9*W, 0.));  

  fprintf (fp, "%d %g %g %g %g %g %g %g %g\n", 
           i, t, dt, 
           interpolate(u.y, 0.1*W, 0.), 
           interpolate(u.y, 0.9*W, 0.),
           interpolate(u.x, 0.05*W, 0.),
           interpolate(p, 0.9*W, W/2),
           interpolate(p,0.9*W, -W/2),
           interpolate(omega, 0.9*W, 0.));
}

event movie (t += 2; t <= Tend) {
  scalar m[];

  foreach()
    m[] = cs[];

  vorticity (u, omega);

  output_ppm (omega,
              file = "vorticity.mp4",
              n = N,
              min = -0.5, max = 0.5,
              linear = true,
              mask = m);
}

event field_export (t = Tend) {

  static FILE * fp = NULL;
  if (!fp)
    fp = fopen("edgetonefield.dat", "w");

  vorticity (u, omega);

  for (double xp = 0; xp <= 50; xp += 0.05)
    for (double yp = -10; yp <= 10; yp += 0.05)
      fprintf(fp, "%g %g %g\n", xp, yp, interpolate(omega, xp, yp));

  fclose(fp);
}