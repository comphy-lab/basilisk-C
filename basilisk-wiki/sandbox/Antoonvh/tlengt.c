/**
# Length of a loop of particles

`plength` is more accurate than simple distance.

The length of two curves is estimated from the particles placed on those
curves.

~~~gnuplot
set key left outside
set size ratio -1
plot 'log' w l lw 2 t 'curve 1', 'circ' w l lw 2 t 'curve 2'
~~~

The estimates are compared against Wolfram-Alpha's reference value

~~~gnuplot Error Curve 1
reset
set size square
set logscale xy
set xr [50: 25600]
set grid
set xlabel 'Particles'
set ylabel 'Error in length'
plot 'out' , 'out' u 1:3, 50000000*x**(-4), 1000*x**(-2)
~~~

The higher-order method is less robust for unresolved features.

~~~gnuplot Error Curve 2
reset
set size square
set logscale xy
set xr [1.5: 256]
set grid
set xlabel 'Particles'
set ylabel 'Error in length'
plot 'out' u 4:5 , 'out' u 4:6, 50*x**(-4), 10*x**(-2)
~~~

But only needs a few point to resolve a circle: The 6-point estimate is less than 1% less off.
 */
#define X(T) (sin(T) + 0.5*erf(2*cos(3*T)))
#define Y(T) (cos(T) + 0.4*exp(2*sin(4*T)))
#define ANS (27.6374980235184900345111)

#include "particle.h"
#include "run.h"

// Simple length:
double length (Particles p) {
  double lenp = 0;
  coord p1 = {0,0}, pp = {0};
  foreach_particle_in(p) {
    if (p1.x == 0 && p1.y == 0)
      p1 = (coord){x, y};
    else {
      double d = 0;
      foreach_dimension()
	d += sq(p().x - pp.x);
      lenp += sqrt(d);
    }
    pp = (coord){x, y};
  }
  double d = 0;
  foreach_dimension()
    d += sq(p1.x - pp.x);
  lenp += sqrt(d);
  return lenp;
}

int np;
int main() {
  for (np = 100; np <= 12800; np*= 2) 
    run();
}

event init (t = 0) {
  Particles part1 = new_particles(np);
  Particles part2 = new_particles(np/33);
  foreach_particle_in(part1) {
    double T = _j_particle*2.*pi/(np);
    p().x = X(T);
    p().y = Y(T);
    
  }
  foreach_particle_in(part2) {
    double Theta = _j_particle*2.*pi/(np/33);
    p().x = cos(Theta);
    p().y = sin(Theta);
  }
  printf ("%d %g %g %d %g %g\n", np, fabs(ANS-plength(part1)), 
          fabs(ANS-length(part1)), np/33, fabs(2*pi-plength(part2)),
          fabs(2*pi-length(part2)));
  if (np == 3200) {
    foreach_particle_in(part1)
      fprintf (stderr, "%g %g\n", x, y);
    FILE * fp = fopen("circ", "w");
    foreach_particle_in(part2)
      fprintf (fp, "%g %g\n", x, y);
  }
  
}
