

/* Function for outputting the performances of the dlmfd solver */
void output_dlmfd_perf (const timing d, const int i, particle * p) {

double mpitimings[npe()];
  
  timing s = d;
  static FILE * dlmfdperf = fopen ("dlmfd-perf", "w");

  /* global timer/timings */
  timing gns = timer_timing (perf.gt, i, perf.tnc, mpitimings);
  
  /* Compute the speed with the total number of constrained cells by the dlmfd problem (that number is multiplied by the number of iterations) */
  s.speed = p[0].tcells/s.real;

  /* Display and write the statistics of the dlmfd solver */
  fprintf (fout,"\n# dlmfd stats: " GRIDNAME
	   ", %ld average constrained cells, %ld average Lagrange multipliers, %ld total cells\n", p[0].tcells/i, p[0].tmultipliers/i, perf.tnc/i);

  fprintf (dlmfdperf,"\n# dlmfd stats: " GRIDNAME
	   ", %ld average constrained cells, %ld average Lagrange multipliers, %ld total cells\n", p[0].tcells/i, p[0].tmultipliers/i, perf.tnc/i);
  
  fprintf (fout,"\n# dlmfd stats: " GRIDNAME 
	   ", %d steps, %g CPU, %.4g real, %.3g points.step/s\n",
	   i, s.cpu, s.real, s.speed);
  fprintf (dlmfdperf,"\n# dlmfd stats: " GRIDNAME 
	   ", %d steps, %g CPU, %.4g real, %.3g points.step/s\n",
	   i, s.cpu, s.real, s.speed);
#if _MPI
  fprintf (fout,
	   "# dlmfd MPI stats: %d procs, MPI: min %.2g (%.2g%%) "
	   "avg %.2g (%.2g%%) max %.2g (%.2g%%)\n",
	   npe(),
	   s.min, 100.*s.min/s.real,
	   s.avg, 100.*s.avg/s.real,
	   s.max, 100.*s.max/s.real);

  fprintf (dlmfdperf,
	   "# dlmfd MPI stats: %d procs, MPI: min %.2g (%.2g%%) "
	   "avg %.2g (%.2g%%) max %.2g (%.2g%%)\n",
	   npe(),
	   s.min, 100.*s.min/s.real,
	   s.avg, 100.*s.avg/s.real,
 	   s.max, 100.*s.max/s.real);
#endif
  
  /* Display and write the ratio dlmfd/total  */
  fprintf (fout, "\n# ratio dlmfd/total: " GRIDNAME
	   ", %d steps, %g CPU, %.4g real, \n", i, s.cpu/gns.cpu, s.real/gns.real);
  fprintf (dlmfdperf, "\n# ratio dlmfd/total: " GRIDNAME
	   ", %d steps, %g CPU, %.4g real, \n", i, s.cpu/gns.cpu, s.real/gns.real);
  
#if _MPI
  
  fprintf (fout,
	   "# ratio MPI dlmfd/total: %d procs, MPI: min %.2g "
	   "avg %.2g max %.2g\n",
	   npe(),
	   s.min/gns.min,
	   s.avg/gns.avg, 
	   s.max/gns.max);
 fprintf (dlmfdperf,
	  "# ratio MPI dlmfd/total: %d procs, MPI: min %.2g "
	  "avg %.2g max %.2g\n",
	  npe(),
	  s.min/gns.min,
	  s.avg/gns.avg, 
	  s.max/gns.max);
#endif
 
 /* Write in the file with pointer dlmfdperf the total statistics */
 s = gns;
 fprintf (dlmfdperf,
	  "\n# " GRIDNAME 
	  ", %d steps, %g CPU, %.4g real, %.3g points.step/s, %d var\n",
	  i, s.cpu, s.real, s.speed, (int) (datasize/sizeof(double)));
#if _MPI
 fprintf (dlmfdperf,
	  "# %d procs, MPI: min %.2g (%.2g%%) "
	  "avg %.2g (%.2g%%) max %.2g (%.2g%%)\n",
	  npe(),
	  s.min, 100.*s.min/s.real,
	  s.avg, 100.*s.avg/s.real,
	  s.max, 100.*s.max/s.real);
#endif
 fflush (dlmfdperf);
 
}
