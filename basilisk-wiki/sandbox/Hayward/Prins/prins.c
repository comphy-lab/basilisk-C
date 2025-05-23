/**
# Dam-break Wave Generation

[Prins, 1958](#prins1958) conducted a dam break experiment to investigate the
characteristics of the resultant waves. The 2D lab-scale physical model was
used to investigate a range of elevated and depressed water columns
instantaneously released, and compare the results to the wave theories of the
time.

In a flume of depth $H$, a water column is held behind the origin of height difference $Q$ and of horizontal extent $L$. These parameters are varied throughout the experimental series. The column is released from rest (in the experiment by means of quickly pulling a slide) and the resultant waves are measured down flume.

It is compared here to results produced by the
[non-hydrostatic multilayer](/src/layered/nh.h) and
[Navier-Stokes/VOF](/src/navier-stokes/centered) solvers. Comparisons are
also made to the [Green-Naghdi](/src/green-naghdi.h) and
[Saint-Venant equations](/src/saint-venant.h). */


#include <sys/stat.h>
#include "grid/multigrid1D.h"
#if GN
# include "green-naghdi.h"
#elif NH
# include "layered/hydro.h"
# include "layered/nh.h"
# include "layered/remap.h"
# include "layered/perfs.h"
#else // if SV
# include "saint-venant.h"
#endif

/**
Being research undertaken in California ***many*** years ago, units were of the
Imperial standard (feet, inches, pounds, etc.)

In the physical experiment, the end of the flume contained a simple 'wave absorber'. The domain for the numerical experiment is extended to prevent reflections. */

#define GRAVITY 9.81
#define LENGTH 30.
#define ENDTIME 30.
#define L 0.09144  // 0.3 ft in m
#define H 0.15240  // 0.5 ft in m
#define Q 0.09144  // +0.3 ft in m


/**
This example has parameters set for one of the two large time-series
illustrations in the publication.

The multilayer run has 30 layers imposed. */


int main()
{
  struct stat st = {0};
  if (stat("./profiles", &st) == -1) {
    mkdir("./profiles", 0755);
  }
  X0 = -L;
  L0 = LENGTH + L;
  G = GRAVITY;
  N = 8192; // 2**13
#if NH
  nl = 30; // number of layers
  breaking = 0.07;
#endif
  run();
}


/**
The initial conditions are set to instantiate the water column at rest.
This resembles a dam-break situation. */


event init (i = 0)
{
  foreach() {
#if NH
    foreach_layer() {
      h[] = x < 0. ? (Q+H)/nl : H/nl; 
      u.x[] = 0.;
    }
#else  // if SV or GN
    h[] = x < 0. ? Q+H : H;
    u.x[] = 0.;
#endif
  }
}

/**
Numerical guages are positioned at identical locations as the physical. */

Gauge gauges[] = {
  {"x05",  1.524},
  {"x15",  4.572},
  {"x25",  7.62 },
  {"x35", 10.668},
  {"x45", 13.716},
  {NULL}
};

/**
Fool-proof check and progress is printed. */

event output (i++) {
  if (i == 0) {
#if GN
    fprintf (ferr, "Serre-Green-Naghdi\n");
#elif NH
    fprintf (ferr, "Multilayer\n");
#else // if SV
    fprintf (ferr, "Saint Venant\n");
#endif
    fprintf (ferr, "t dt\n");
  }
  fprintf (ferr, "%g %g\n", t, dt);

  // Gauges
  output_gauges (gauges, {eta});
}

/**
Additionally output are free-surface profiles early in the simulation to compare between models, and an example of the multilayer distribution. */

event profile_layers_output (t = 0.35) {
#if NH
  static FILE * fp_layer = fopen ("profile_layers_t035.dat", "w");
  foreach() {
    int count = 0;
    foreach_layer() {
      fprintf (fp_layer, "%d %g %g\n", count, x, h[]);
      count++;
    }
  }
#endif
}

event profiles (t <= 1.5; t += 0.05) {
  char name[60];
  sprintf(name, "profiles/profile_t_%g.dat", t);
  FILE * fp = fopen (name, "w");
  foreach()
    fprintf (fp, "%g %g\n", x, eta[]);
  fclose(fp); 
}

event end (t = ENDTIME) {
}



/**
## Results

*(Must replace with GNUplot)*

The multilayer scheme has very good agreement with the physical experiment and the N-S solution, and accurately predicts the phase arrivals and heights. Green-Naghdi approaches some consensus, but never accurately models the characteristics of the first arrival.

<img src="https://i.ibb.co/Z1h0cQq/Screenshot-2020-12-01-at-22-01-34.png" alt="" width=800px/>




Leading crest height at the first gauge ( $x=5 ft$ ). Red dashed line from physical experiment. Key shows variation of $L$.

<img src="https://i.ibb.co/FqyCrJB/Screenshot-2020-12-01-at-22-02-19.png" alt="" width=487px/> 




Profile comparisons of the different models at $t=0.35$.

<img src="https://i.ibb.co/X4Zqg4S/Screenshot-2020-12-01-at-22-03-09.png" alt="" width=640px/>




Time progressing profiles of the multilayer solution compared with the Navier-Stokes VOF solution.

<img src="https://i.ibb.co/BZ8Kh5f/NH-time.jpg" alt="" width=500px/><img src="https://i.ibb.co/QMKPfFD/facets-lvl15.jpg" alt="" width=507px/>




*/



/**
## To-do

 - Actually use GNUplot (for once) and replace the plots.
 - Upload equivalent N-S case.
 - Reformat to use consistant units...!!

*/




/**
## References

~~~bib
@article{prins1958,
  title = {Characteristics of waves generated by a local disturbance},
  journal = {Eos, Transactions American Geophysical Union},
  volume = {39},
  number = {5},
  pages = {865-874},
  doi = {10.1029/TR039i005p00865},
  url = {https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/TR039i005p00865},
  pdf = {https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1029/TR039i005p00865},
  year = {1958}
}
~~~
*/

