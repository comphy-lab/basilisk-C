/**
# Performance monitoring (for the layered solver)

This logs simple statistics available for the [layered
solver](/src/layered/hydro.h). */

event perfs (i += 1) {
  static FILE * fp = fopen ("perfs", "w");
  if (i == 0)
    fprintf (fp,
	     "t dt mgp.i mgp.nrelax grid->tn perf.t perf.speed npe perf.ispeed\n");
  static double start = 0.;
  if (i > 10 && perf.t - start < 1.) return 0;
  fprintf (fp, "%g %g %d %d %ld %g %g %d %g\n", 
	   t, dt,
#if NH
	   mgp.i, mgp.nrelax,
#elif IMPLICIT_H
	   mgH.i, mgH.nrelax,
#else
	   0, 0,
#endif
	   grid->tn*nl, perf.t, perf.speed*nl, npe(), perf.ispeed*nl);
  fflush (fp);
  start = perf.t;
}

/**
If we have a display (and gnuplot works), a graph of the statistics is
displayed and updated at regular intervals (10 seconds as defined in
[perfs.plot]()). */

event perf_plot (i = 10) {
  if (getenv ("DISPLAY"))
    popen ("gnuplot -e 'set term x11 noraise title perfs' "
	   "$BASILISK/layered/perfs.plot "
	   "& read dummy; kill $!", "w");
}

/**
## See also

* [Performance monitoring (for the Navier-Stokes
   solvers)](/src/navier-stokes/perfs.h).
* [Continuous profiling monitoring](/src/profiling.h).
*/
