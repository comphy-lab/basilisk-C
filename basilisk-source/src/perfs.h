/**
# Performance monitoring */

event perfs (i += 1) {
  static FILE * fp = fopen ("perfs", "w");
  if (i == 0)
    fprintf (fp,
	     "t dt grid->tn perf.t perf.speed npe perf.ispeed maxrss\n");
  static double start = 0.;
  if (i > 10 && perf.t - start < 1.) return 0;
  fprintf (fp, "%g %g %ld %g %g %d %g ",
	   t, dt,
	   grid->tn*nl, perf.t, perf.speed*nl, npe(), perf.ispeed*nl);
@if _GNU_SOURCE
  struct rusage usage;
  getrusage (RUSAGE_SELF, &usage);
  fprintf (fp, "%ld\n", usage.ru_maxrss);
@else
  fputs ("0\n", fp);
@endif
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
	   "$BASILISK/perfs.plot "
	   "& read dummy; kill $!", "w");
}

/**
## See also

* [Performance monitoring (for the layered
  solver)](/src/layered/perfs.h).
* [Performance monitoring (for the Navier-Stokes
   solvers)](/src/navier-stokes/perfs.h).
* [Continuous profiling monitoring](/src/profiling.h).
*/
