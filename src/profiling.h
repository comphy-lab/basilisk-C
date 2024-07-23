/**
# Continuous profiling monitoring

Including this file will display (if the DISPLAY variable is set and
gnuplot works), continuous [profiling
information](README.trace#graphing-tools) about the running solver. */

#if TRACE > 1
event profiling (i += 19) {
  static FILE * fp = fopen ("profiling", "w");
  trace_print (fp, 1);
}

event profiling_plot (i = 40) {
  if (getenv ("DISPLAY"))
    popen ("gnuplot -e 'set term x11 noraise noenhanced title profiling' "
	   "$BASILISK/profiling.plot "
	   "& read dummy; kill $!", "w");
}
#endif
