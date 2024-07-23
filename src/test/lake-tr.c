/**
# Lake flowing into itself

We consider a periodic domain, 10 km long, with a 50-metres-deep
central lake. The domain is tilted with a slope of 1/1000 and a
Darcy-Weissbach friction is imposed on the bottom.

~~~gnuplot Topography
set xlabel 'x (m)'
set ylabel 'z (m)'
zb(x) = - 50.*exp(-(x/1000.)**2) - 0.001*x 
plot [-5000:5000] zb(x) w l t ''
~~~

We use either the explicit or the implicit Saint-Venant solver. */

#include "grid/multigrid1D.h"
#if EXPLICIT
# include "saint-venant.h"
#elif ML
# include "layered/hydro.h"
# include "layered/nh.h"
#else
# include "saint-venant-implicit.h"
#endif

/**
## Domain setup and initial conditions */

#define LEVEL 10

double slope = 0.001, eta0 = 0.5, f = 0.025;

int main() {
  L0 = 10000.;
  X0 = -L0/2.;
  G = 9.81;
  N = 1 << LEVEL;
  periodic (right);

#if ML
  /**
  We add some damping by off-centering the implicit scheme. */
  
  theta_H = 0.6;
#endif

  DT = 10;  
  run();
}

event init (i = 0)
{
  foreach() {
    zb[] = - 50.*exp(-sq(x/1000.));
    h[] = eta0 - zb[];
  }
}

/**
## Source terms

We add the slope explicitly and the Darcy--Weissbach friction
implicitly. */

event source (i++) {
#if EXPLICIT || ML
  foreach()
    u.x[] = (u.x[] + G*slope*dt)/(1. + dt*f/8.*u.x[]/h[]);
#else
  foreach()
    q.x[] = (q.x[] + h[]*G*slope*dt)/(1. + dt*f/8.*q.x[]/sq(h[]));
#endif
}

/**
## Outputs

We check for steady-state. */

scalar herr[], uerr[];
event init_herr(i = 0) {
  foreach()
    herr[] = uerr[] = 0;
}

event logfile (i++) {
  double dh = change (h, herr),
#if EXPLICIT || ML
    du = change (u.x, uerr);
#else
    du = change (q.x, uerr);
#endif
  fprintf (stderr, "%g %g %g\n", t, dh, du);
#if 0
  if (i > 100 && dh < 1e-7 && du < 1e-7) {
    return 1;
  }
#endif
}

/**
We output the free-surface, Froude number etc.. */

event printprofile (t = {600, 3600, 7200})
{
  char name[80];
  sprintf (name, "prof-%g", t);
  FILE * fp = fopen (name, "w");
  foreach() {
#if EXPLICIT || ML
    fprintf (fp, "%g %g %g %g %g %g\n", x, h[], u.x[], zb[],
	     u.x[]/sqrt(G*h[]), u.x[]*h[]);
#else
    fprintf (fp, "%g %g %g %g %g %g\n", x, h[], q.x[]/h[], zb[],
	     q.x[]/(h[]*sqrt(G*h[])), q.x[]);
#endif
  }
  fclose (fp);
}

#if 0
event gnuplot (i += 10)
{
  static FILE * fp = popen ("gnuplot", "w");
  if (i == 0)
    fputs ("set term x11\nset yrange [-0.2:2]\n", fp);
  fprintf (fp, "plot '-' u 1:2 w l\n");
  foreach()
    fprintf (fp, "%g %g\n", x, u.x[]/sqrt(G*h[]));
  fprintf (fp, "e\n");
  fflush (fp);
}
#endif

/**
## Results

After two hours a steady-state is reached. Both schemes give
comparable results, even for the transient solution at 10 minutes.

~~~gnuplot Evolution of the free surface
set xlabel 'x (m)'
set ylabel 'z (m)'
set key top right
plot [][-6:]								\
     'prof-600' u 1:(zb($1)+$2) w l t 't = 10 min (implicit)',		\
     '../lake-tr-ml/prof-600' u 1:(zb($1)+$2) w l t			\
     't = 10 min (implicit ML)',					\
     '../lake-tr-explicit/prof-600' u 1:(zb($1)+$2) w l t		\
     't = 10 min (explicit)',						\
     'prof-7200' u 1:(zb($1)+$2) w l t 't = 2 hours (implicit)',	\
     '../lake-tr-ml/prof-7200' u 1:(zb($1)+$2) w l t			\
     't = 2 hours (implicit ML)',					\
     '../lake-tr-explicit/prof-7200' u 1:(zb($1)+$2) w l t		\
     't = 2 hours (explicit)',						\
     zb(x) t 'topography'
~~~

The flow is supercritical at the entrance to the lake.
   
~~~gnuplot Evolution of the Froude number
set ylabel 'Froude'
set key top right
plot 'prof-600' u 1:5 w l t 't = 10 min (implicit)', \
     '../lake-tr-ml/prof-600' u 1:5 w l t 't = 10 min (implicit ML)', \
     '../lake-tr-explicit/prof-600' u 1:5 w l t 't = 10 min (explicit)', \
     'prof-7200' u 1:5 w l t 't = 2 hours (implicit)',			\
     '../lake-tr-ml/prof-7200' u 1:5 w l t 't = 2 hours (implicit ML)', \
     '../lake-tr-explicit/prof-7200' u 1:5 w l t 't = 2 hours (explicit)'
~~~

We can check the performance for all schemes.

~~~bash
cat lake-tr/out
cat lake-tr-explicit/out
cat lake-tr-ml/out
~~~

On my computer, this gives for the implicit scheme:

~~~bash 
# Multigrid, 4460 steps, 1.42008 CPU, 1.42 real, 3.22e+06 points.step/s, 14 var
~~~ 

for the multilayer implicit scheme:

~~~bash
# Multigrid, 4524 steps, 1.46437 CPU, 1.464 real, 3.16e+06 points.step/s, 14 var
~~~

and for the explicit scheme :

~~~bash
# Multigrid, 32519 steps, 8.59387 CPU, 8.685 real, 3.83e+06 points.step/s, 16 var
~~~

The gain in number of timesteps is a factor of ~7.3 and the gain in
runtime a factor of ~5.6, which reflects the fact that the implicit
scheme is slightly slower (per grid point) than the explicit scheme. */
