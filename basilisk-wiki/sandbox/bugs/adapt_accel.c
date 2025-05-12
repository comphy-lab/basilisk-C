/**
# Boussinesq Rayleigh-Taylor instability example */

#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"

#define temp 4
#define MAXLEVEL 9

scalar T[];
scalar * tracers = {T};
mgstats mgT;
face vector av[];

double Ra, Pr; 

int main()
{
  L0 = 2.;
  X0 = Y0 = -1;
  DT = 0.2;
  init_grid (1 << MAXLEVEL);
  Ra = 1e9;
  Pr = 1;
  const face vector muc[] = {Pr/sqrt(Ra),Pr/sqrt(Ra)};
  mu = muc;
  a = av;
  run();
}

event init (t = 0)
{	
  mask (y >  0.5 ? top :
        y < -0.5 ? bottom : none);
  foreach()
    T[] = (y < 0.0) + 0.01*noise();
}

/**
Apply Boussinesq "gravity". */

event acceleration (i++)
{	
  const face vector D[] = {1./sqrt(Ra), 1./sqrt(Ra)};
  mgT = diffusion (T, dt, D);
  coord delta = {0,1};
  foreach_face()
    av.x[] = delta.x*(T[] + T[-1])/2.;
}

event logfile (t <= temp; t += 0.1) 
{
  scalar omg[];
  vorticity (u, omg);
  stats s = statsf (omg);
  fprintf (ferr, "%g %d %g %g %g\n", t, i, dt, s.sum, s.max);
}

/**
We generate animations of the
[vorticity](adapt_accel/RT_vorticity.mpg), [level of
refinement](adapt_accel/RT_level.mpg),
[temperature](adapt_accel/RT_T.mpg). */

event movie (t += 0.005; t <= temp) 
{
  static FILE * fp = popen ("ppm2mpeg > RT_vorticity.mpg", "w");
  scalar omg[];	
  vorticity (u, omg); 
  output_ppm (omg, fp, min = -10, max = 10);
  static FILE * fp2 = popen ("ppm2mpeg > RT_level.mpg", "w");
  foreach()
    omg[] = level;
  output_ppm (omg, fp2, min = 5, max = MAXLEVEL);
  static FILE * fp1 = popen ("ppm2mpeg > RT_T.mpg", "w");
  output_ppm (T, fp1, min = 0, max = 1);
}

event adapt (i++) {
  adapt_wavelet ((scalar *){T,u}, (double[]){0.01,0.01,0.01}, MAXLEVEL);
}
