/**
# Post-processing for the Held-Suarez example

This example illustrates how to post-process the outputs of the [Held
& Suarez example](held-suarez.c) to obtain classical diagnostics.

The time-averaged zonal velocity compares closely with Figure 2 of
[Held & Suarez, 1994](held-suarez.c#held1994). The easterly trade
winds at the surface reach -8 m/s at $\pm$ 15 degrees latitude, the
mid-latitude westerlies 8 m/s at $\pm$ 45 degrees. The jet streams in
the upper troposphere reach 30 m/s.

~~~gnuplot Time-averaged zonal velocity {width="640px"}
set term png enhanced size 1280,800 font ",18"
set pm3d map interpolate 4,4 at b
unset key
# jet colormap
set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25 0 0.5647 1,	\
0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392, 0.625 1 0.9333 0, 0.75 1 0.4392 0,	\
0.875 0.9333 0 0, 1 0.498 0 0 )
set contour base
set cntrparam levels incremental -12,4,32
set cbrange [-12:32]
set cntrlabel onecolor
set cntrlabel font ",14"
set xtics -90,30,90
set xlabel 'Latitude'
set ylabel 'Altitude'
set output 'zonal_velocity.png'
splot [-90:90][0:1]'out' index 'zonal velocity' u 1:($2/19.):3 lt 3 lc rgb "#000000", \
     '' index 'zonal velocity' u 1:($2/19.):3 w labels
~~~

The variance of the zonal temperature also compares well with Figure 3
of [Held & Suarez, 1994](held-suarez.c#held1994).

Note that the early calculations done with Gerris and reported in
[Popinet, 2012](#popinet2012), were correct but that the graph
presented in slide 27 (left) is the *potential* temperature, not the
temperature, which explains the difference in the stratosphere
compared to the graph on slide 27 (right) from Ringler et al, 2000.

~~~gnuplot Time-averaged zonal temperature variance {width="640px"}
set cbrange [*:*]
set cntrparam levels incremental 5,5,40
set cbrange [0:40]
nl = 20
Kappa = 2./7. # R/cp
P0 = 101300.
Ps(z) = P0 - (z)
PI(i) = (Ps(P0*(i + 0.5)/nl)/P0)**Kappa
set output 'zonal_variance.png'
splot [-90:90][0:1]'out' index 'zonal variance' u 1:($2/19.):($3*PI($2)**2) lt 3 lc rgb "#000000", \
     '' index 'zonal variance' u 1:($2/19.):($3*PI($2)**2) w labels
~~~

So do the time-averaged zonal temperature and zonal potential temperature.

~~~gnuplot Time-averaged zonal temperature {width="640px"}
set cbrange [*:*]
set cntrparam levels incremental 190,5,305
set cbrange [190:305]
set output 'zonal_temperature.png'
splot [-90:90][0:1]'out' index 'zonal temperature' u 1:($2/19.):($3*PI($2)) lt 3 lc rgb "#000000", \
     '' index 'zonal temperature' u 1:($2/19.):($3*PI($2)) w labels
~~~

~~~gnuplot Time-averaged zonal potential temperature
set term svg enhanced size 640,400 font ",9"
set cntrparam levels incremental 270,5,350
set cbrange [*:*]
unset surface
unset colorbox
set cntrlabel font ",7"
splot [-90:90][0:1]'out' index 'zonal temperature' u 1:($2/19.):3 lt 3 lc rgb "#000000", \
     '' index 'zonal temperature' u 1:($2/19.):3 w labels
~~~
*/

#include "grid/multigrid.h"
#include "utils.h"

/**
This macro computes the zonal average of `expr` at level `l`. */

macro zonal_average (double expr, int l)
{{
  const double maxlat = 75;
  double delta = L0/N;
  for (double y = - 90 + delta/2.; y <= maxlat; y += delta)
    if (y > - maxlat) {
      double sum = 0.;
      coord p;
      coord box[2] = {{X0, y - delta/2.}, {X0 + L0, y + delta/2.}};
      coord n = {N, 1};
      foreach_region (p, box, n, reduction(+:sum))
        sum += expr;
      printf ("%g %d %g\n", y, l, sum/n.x);
    }
}}

int main ()
{
  if (!restore (file = "../held-suarez/dump", list = all)) {
    fprintf (stderr, "could not restore\n");
    exit (1);
  }
  restriction (all);
  fields_stats();

  printf ("# zonal velocity\n");
  for (int l = 0; l < 20; l++) {
    char name[80] = "mu";
    if (l > 0)
      sprintf (name, "mu%d", l);
    scalar mu = lookup_field (name);
    zonal_average (mu[], l);
    printf ("\n");
  }

  printf ("\n# zonal variance\n");
  for (int l = 0; l < 20; l++) {
    char name[80] = "vtheta";
    if (l > 0)
      sprintf (name, "vtheta%d", l);
    scalar vtheta = lookup_field (name);
    zonal_average (vtheta[], l);
    printf ("\n");
  }
  
  printf ("\n# zonal temperature\n");
  for (int l = 0; l < 20; l++) {
    char name[80] = "mtheta";
    if (l > 0)
      sprintf (name, "mtheta%d", l);
    scalar mtheta = lookup_field (name);
    zonal_average (mtheta[], l);
    printf ("\n");
  }
}

/**
## Todo

* The "vertically averaged zonal spectra of the eddy variance of zonal
  wind" (Figure 4 of H&S, 1994) would be nice too.

## References

~~~bib
@InProceedings{popinet2012,
  author =	 {S. Popinet},
  title =	 {Quadtree-adaptive global atmospheric modelling on
                  parallel systems (invited talk)},
  year =	 2012,
  booktitle =	 {Weather and Climate Prediction on Next Generation
                  Supercomputers, Exeter, UK, 22-25 October},
  pdf =          {/sandbox/popinet/newton.pdf}
}
~~~
*/
