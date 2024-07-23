event gnuplot (i += 20)
{
  static FILE * fp = popen ("gnuplot 2> /dev/null", "w");
  if (i == 0)
    fprintf (fp, "set term x11\n");
  FILE * fp1 = fopen ("gnuplot", "w");
  foreach_leaf() {
    double H = 0.;
    foreach_layer()
      H += h[];
    fprintf (fp1, "%g %g %g", x, zb[] + H, zb[]);
    double z = zb[];
    foreach_layer() {
      fprintf (fp1, " %g", z);
      z += h[];
    }
    fprintf (fp1, "\n");
  }
  fclose (fp1);
  fprintf (fp,
	   "set title 'nl = %d, t = %.2f'\n"
	   "p [%g:%g][0:]'gnuplot' u 1:3:2 w filledcu lc 3 t '',"
	   " '' u 1:(-1):3 t '' w filledcu lc -1", nl, t,
	   X0, X0 + L0);
  int i = 4;
  foreach_layer()
    fprintf (fp, ", '' u 1:%d w l lw 2 t ''", i++);
  fprintf (fp, "\n");
  fflush (fp);

  usleep(100000);
}
