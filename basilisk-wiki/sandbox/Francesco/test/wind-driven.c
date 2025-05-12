/**
## Wind-driven lake

This is a simple test case of a wind-driven lake where we can compare
results with an analytical solution. For the bottom of the domain we
impose *no-slip* condition (that is the default), for the top we
impose a Neumann condition (see [viscous friction between
layers](http://www.basilisk.fr/src/multilayer.h#viscous-friction-between-layers)
for details). */

#include "grid/cartesian1D.h"
#include "saint-venant.h"

/**
Analytical solution for the vertical profile of velocity at the center of the domain.*/

#define uan(z)  z/4.*(3.*z-2.)

int main() {
  L0 = 1.;
  G = 100.;
  N = 64;
  nu = 1.;
  dut = unity;

  /**
  We vary the number of layers. */
  
  for (nl = 4; nl <= 32; nl *= 2)
    run();
}

/**
We set the initial water level to 1 and we allocate a scalar field
*uc* to check the convergence of the velocity in the first layer.*/

scalar uc[];

event init (i = 0) {
  vector u0 = ul[0];
  foreach() {
    h[] = 1.;
    uc[] = u0.x[];
  }
}

/**
We check for convergence. */

event logfile (t += 0.1; i <= 100000) {
  vector u0 = ul[0];
  double du = change (u0.x, uc);
  if (i > 0 && du < 1e-9)
    return 1;
}

/**
We compute the error between the numerical solution and the 
analitycal solution. */

event error (t = end) {
  static FILE * fp = fopen("error","w");
  double e[nl];
  double emax = 0.;
  int i = 0;
  foreach() {
    if (i++ == N/2) {
      int l = 0;
      double z = 0.;
      double sumh = 0.;
      for (vector u in ul) {
	sumh += h[]*layer[l];
	z = sumh - h[]*layer[l]*0.5;
	e[l] = fabs(u.x[] - uan (z));
	if (e[l] > emax) 
	  emax = e[l];
	l++;
      }
      fprintf (fp, "%d %g\n", nl, emax);
    }
  }
}

/**
We save the horizontal velocity profile at the center of the domain
and the two components of the velocity field for the case with 32
layers.*/

event output (t = end) {
  char name[80];
  sprintf (name, "uprof-%d", nl);
  FILE * fp = fopen (name, "w");
  int i = 0;
  foreach() {
    if (i++ == N/2) {
      int l = 0;
      double z = zb[] + h[]*layer[l]/2.;
      for (vector u in ul)
	fprintf (fp, "%g %g\n", z, u.x[]), z += h[]*layer[l];
    }
    if (nl == 32) {
      double z = zb[];
      int l = 0;
      scalar w;
      vector u;
      for (w,u in wl,ul) {
	printf ("%g %g %g %g\n", x, z, u.x[], w[]);
	z += layer[l++]*h[];
      }
      printf ("\n");
    }
  }
  fclose (fp);
}

/**
## Results

~~~gnuplot Numerical and analytical velocity profiles at the center of the lake.
set xr [0:1]
set xl 'z'
set yl 'u'
set key left top
plot [0:1]x/4.*(3.*x-2.) t 'analytical', \
          'uprof-4' t '4 layers', \
          'uprof-8' t '8 layers', \
          'uprof-16' t '16 layers', \
          'uprof-32' t '32 layers'
~~~

~~~gnuplot Convergence of the error between the numerical and analytical solution with the number of layers.
reset
set logscale
set xlabel 'nl'
set ylabel 'max|e|'
set xtics 4,2,32
set grid
fit a*x+b 'error' u (log($1)):(log($2)) via a,b
plot [3:36]'error' u 1:($2) pt 7 t '', \
     exp(b)*x**a t sprintf("%.0f/N^{%4.2f}", exp(b), -a)
~~~

~~~gnuplot Horizontal velocity field (32 layers).
reset
set pm3d
set pm3d map interpolate 10,1
unset key
set xlabel 'x'
set ylabel 'z'
# jet colormap
set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25 0 0.5647 1,	\
0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392, 0.625 1 0.9333 0, 0.75 1 0.4392 0,	\
0.875 0.9333 0 0, 1 0.498 0 0 )

splot 'out' u 1:2:3
~~~

~~~gnuplot Vertical velocity field (32 layers).
set pm3d
set pm3d map interpolate 10,1
unset key
set ylabel 'z'
splot 'out' u 1:2:4
~~~

## See also

* [Same test with Gerris](http://gerris.dalembert.upmc.fr/gerris/tests/tests/lake.html#river)
*/
