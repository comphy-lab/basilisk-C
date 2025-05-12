/**
# A Mach-3 wind tunnel with a step

This is the classical test case first proposed by [Emery,
1968](#emery1968) (see also [Woodward & Colella, 1984](#woodward1984)).

We use the compressible flow solver, the initial density, pressure and
velocity are given. */

#include "compressible.h"

double rho0 = 1.4, u0 = 3., p0 = 1.;
int maxlevel = 9;

/**
The left boundary is prescribed constant inflow and the right boundary
is simple outflow. */

w.n[left]  = dirichlet (rho0*u0);
w.n[right] = neumann(0);

/**
The domain is 3 units long. We use a 256x256 initial grid. */

int main() {
  L0 = 3.;
  N = 256;
  run();
}

event init (t = 0) {

  /**
  Masking is used to define the top of the tunnel (1 unit wide) and
  the step. */
  
  mask (y > 1 ? top : x > 0.6 && y < 0.2 ? bottom : none);

  /**
  These are the initial conditions for density, momentum and
  energy. */
  
  foreach() {
    rho[] = rho0;
    w.x[] = rho0*u0;
    E[] = p0/(gammao - 1.) + rho0*sq(u0)/2.;
  }
  boundary ({rho, w.x, E});
}

/**
We log some statistics. */

event logfile (i++) {
  if (i == 0)
    fprintf (ferr,
	     "t dt grid->tn perf.t perf.speed\n");
  fprintf (ferr, "%g %g %ld %g %g\n", 
	   t, dt, grid->tn, perf.t, perf.speed);
}

/**
We make movies of the density and level of refinement... */

event movie (t += 0.01) {
  static FILE * fp = popen ("ppm2mp4 rho.mp4", "w");
  output_ppm (rho, fp, box = {{0,0},{3,1}}, linear = true, n = 512);

  static FILE * fp1 = popen ("ppm2mp4 level.mp4", "w");
  scalar l[];
  foreach()
    l[] = level;
  output_ppm (l, fp1, box = {{0,0},{3,1}}, min = 0, max = maxlevel, n = 512);
}

/**
... and dump some snapshots. */

event snapshots (t += 0.5; t <= 3) {
  char name[80];
  sprintf (name, "dump-%g", t);
  dump (file = name);

  sprintf (name, "rho-%g", t);
  FILE * fp = fopen (name, "w");
  output_matrix (rho, fp, 1 << maxlevel, true);
  fclose (fp);
}

/**
Adaptation is only on density, down to maxlevel. */

event adapt (i++) {
  adapt_wavelet ({rho}, (double[]){1e-2}, maxlevel);
}

/**
The results look OK.

<p><center>
<video width="512" height="170" controls>
<source src="step/rho.mp4" type="video/mp4">
Your browser does not support the video tag.
</video> 
</center></p>

<p><center>
<video width="512" height="170" controls>
<source src="step/level.mp4" type="video/mp4">
Your browser does not support the video tag.
</video>
</center></p>

The figure below can be compared to that in [Woodward & Colella,
1984](#woodward1984), last frame of figure 3, with the same choice of
contour levels.

~~~gnuplot Isolines of density at $t=3$.
set contour base
unset surface
set cntrparam levels incremental 0.2673, (6.383 - 0.2673)/29., 6.383
set table 'contours.txt'
splot [0:3][0:0.99]'rho-3' binary u 2:1:3 w l
unset table
set size ratio -1
set grid
plot [0:3][0:1]'contours.txt' w l lt -1 t ''
~~~

## References

~~~bib
@article{emery1968,
  title={An evaluation of several differencing methods for 
  inviscid fluid flow problems},
  author={Emery, Ashley F},
  journal={Journal of Computational Physics},
  volume={2},
  number={3},
  pages={306--331},
  year={1968},
  publisher={Elsevier}
}

@article{woodward1984,
  title={The numerical simulation of two-dimensional 
  fluid flow with strong shocks},
  author={Woodward, Paul and Colella, Phillip},
  journal={Journal of computational physics},
  volume={54},
  number={1},
  pages={115--173},
  year={1984},
  publisher={Elsevier},
  url={http://seesar.lbl.gov/anag/publications/colella/A3_1984.pdf}
}
~~~
*/
