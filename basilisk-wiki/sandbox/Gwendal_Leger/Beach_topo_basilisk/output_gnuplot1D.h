



void snapshot_fluid (char filename[], scalar eta, scalar zb, double t, double z0, double z1) { // Takes a snapshot of the fluid (displays eta and zb)
  FILE * fp = fopen ("plot_data", "w");
  foreach(serial)
    fprintf (fp, "%g %g %g\n", x, eta[], zb[]);
  fclose (fp);

  double zb_min = statsf(zb).min;
  if ((z0 == 0.) && (z1 == 0.)) { // By "default", we compute the min and max values on the z axis
    z0 = zb_min;
    z1 = max (statsf(zb).max, statsf(eta).max);
  }
  fp = fopen ("plot.plot", "w");
  fprintf (fp,
	   "set term pdfcairo\n"
	   "set encoding utf8\n"
	   "set output '%s.pdf'\n"
	   "set title 't=%.6g'\n"
	   "set xlabel 'x'\n"
	   "set ylabel 'elevation'\n"
	   "set yrange [%g:%g]\n"
	   "file = 'plot_data'\n"
	   "plot \\\n"
	   "  file using 1:2:3 with filledcurves lc 3 notitle, \\\n"
	   "  file using 1:2 with lines lc rgb 'blue' title 'η', \\\n"
	   "  file using 1:3 with filledcurves above y1=%.3g lc rgb 'black' notitle",
	   filename, t, z0, z1, zb_min);
  fclose (fp);
  
  system ("gnuplot plot.plot");
}
