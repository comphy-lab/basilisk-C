
/**
# Shear stress around a sphere*/

//Geometric parameters of the drop
#define D (0.1)
#define RADIUS (D/2.)

// Physical parameters
#include "advection.h"
#include "vof.h"
#include "view.h"

scalar f[];
scalar * interfaces = {f}, * tracers = NULL;

#define MAXTIME (1.)
int maxlevel = 12; 
double uemax = 0.01;

int main (int argc, char * argv[])
{
  if (argc > 1)
    maxlevel = atoi (argv[1]);
  if (argc > 2)
    uemax = atof (argv[2]);

  size(1.);
  origin (-L0/2, 0, -L0/2);
  init_grid (128);
 
  scalar poly[];
  foreach() {
    poly[] = 0.5*(x+y);   
    u.x[] = poly[];
    u.y[] = poly[];
  }
  boundary ((scalar *){u});

  refine (sq(x) + sq(y-0.5)+ sq(z) < 2.*sq(RADIUS) && level < maxlevel);
  fraction (f, sq(RADIUS) - sq(x) - sq(y-0.5) - sq(z));
  boundary ({f});

  fprintf (ferr, "%g %g %ld %g %g %d\n", 
      t, dt, grid->tn, perf.t, perf.speed, npe());

  char name[80];
  sprintf (name, "snapshot-%g", t); 
  scalar omega[], omegabis[], omegaexa[];
  foreach(){
    if (f[] > 1e-6 && f[] < 1. - 1e-6) {
      coord n = mycs (point, f); 
      normalize (&n);
#if dimension == 2
      omega[] = sqrt(sq(((u.x[0,1] - u.x[0,0]) + (u.y[1,0] - u.y[0,0]) - (u.y[1,0] - u.y[0,0])*n.x*n.x - (u.x[0,1] - u.x[0,0])*n.x*n.x)*n.y/(Delta))
          + sq(( - (u.y[1,0] - u.y[0,0])*n.x*n.y -  (u.x[0,1] - u.x[0,0])*n.x*n.y)*n.y/(Delta)));
      omegabis[] = sqrt(sq(((u.x[0,1] - u.x[0,-1]) + (u.y[1,0] - u.y[-1,0]) - (u.y[1,0] - u.y[-1,0])*n.x*n.x - (u.x[0,1] - u.x[0,-1])*n.x*n.x)*n.y/(2.*Delta))
          + sq(( - (u.y[1,0] - u.y[-1,0])*n.x*n.y -  (u.x[0,1] - u.x[0,-1])*n.x*n.y)*n.y/(2.*Delta)));
      omegaexa[] = sqrt(sq(0.5*(2+x+y)-0.5*(2+x+y)*n.x*n.x*n.y) + sq(-0.5*(2+x+y)*n.x*n.y*n.y));
#else
      omega[] = sqrt(sq(((u.x[0,1] - u.x[0,0]) + (u.y[1,0] - u.y[0,0]) - (u.y[1,0] - u.y[0,0])*n.x*n.x - (u.x[0,1] - u.x[0,0])*n.x*n.x)*n.y/(Delta))
          + sq(( - (u.y[1,0] - u.y[0,0])*n.x*n.y -  (u.x[0,1] - u.x[0,0])*n.x*n.y)*n.y/(Delta)) + sq(( - (u.y[1,0] - u.y[0,0])*n.x*n.z -  (u.x[0,1] - u.x[0,0])*n.x*n.z)*n.y/(Delta)));
      omegabis[] = sqrt(sq(((u.x[0,1] - u.x[0,-1]) + (u.y[1,0] - u.y[-1,0]) - (u.y[1,0] - u.y[-1,0])*n.x*n.x - (u.x[0,1] - u.x[0,-1])*n.x*n.x)*n.y/(2.*Delta))
          + sq(( - (u.y[1,0] - u.y[-1,0])*n.x*n.y -  (u.x[0,1] - u.x[0,-1])*n.x*n.y)*n.y/(2.*Delta)) + sq(( - (u.y[1,0] - u.y[-1,0])*n.x*n.z -  (u.x[0,1] - u.x[0,-1])*n.x*n.z)*n.y/(2.*Delta)));
      omegaexa[] = sqrt(sq((0.5*(2. + y + x) - 0.5*(2. + x + y)*n.x*n.x)*n.y) + sq( -0.5*(2. + x + y)*n.x*n.y*n.y) + sq( -0.5*(2. + x + y)*n.x*n.z*n.y));
#endif
    }
    else {
      omega[] = 0.;
      omegabis[] = 0.;
      omegaexa[] = 0.;
    }
  }
  boundary ({omega,omegabis,omegaexa});
  dump (name);
  norm b = normf (omega);
  norm c = normf (omegabis);
  norm d = normf (omegaexa);
  fprintf (stdout, "%g %g %g %g %g %g %g %g %g %g\n", t, b.avg, b.rms, b.max, c.avg, c.rms, c.max, d.avg, d.rms, d.max);
  fflush (fout);
  run();
}

/**

## Results

~~~gnuplot Evolution of omega in function of time
set xlabel 'Time'
set ylabel 'omega2'
set logscale y
plot 'out' u 1:2 w p t 'omega avg', \
'out' u 1:5 w p t 'omegabis avg', \
'out' u 1:8 w p t 'omegaexa avg'
~~~

*/
