/**
# Interpolation on halos

~~~gnuplot
r(x,y)=sqrt(x*x+y*y)
set xrange [0:0.7]
set xlabel 'Radial coordinate'
set ylabel 'Error on halos'
plot '< awk "(\$3==10){print \$0}" out' u (r($1,$2)):($6) t 'level 10', \
     '< awk "(\$3==9){print \$0}" out' u (r($1,$2)):($6) t 'level 9', \
     '< awk "(\$3==8){print \$0}" out' u (r($1,$2)):($6) t 'level 8', \
     '< awk "(\$3==7){print \$0}" out' u (r($1,$2)):($6) t 'level 7', \
     '< awk "(\$3==6){print \$0}" out' u (r($1,$2)):($6) t 'level 6', \
     '< awk "(\$3==5){print \$0}" out' u (r($1,$2)):($6) t 'level 5', \
     '< awk "(\$3==4){print \$0}" out' u (r($1,$2)):($6) t 'level 4'
~~~
*/

scalar h[];

int main (int argc, char ** argv)
{
  int n = 2048;
  init_grid (n);

  double R0 = 0.1;
  foreach()
    h[] = exp(-(x*x + y*y)/(R0*R0));
  
  /* initial coarsening (see halo.c) */
  double tolerance = 1e-4;
  while (adapt_wavelet ({h}, &tolerance, 11).nc);

#if 0
  // we reinitialise h just to make sure that trash() does its job
  trash ({h});
  foreach()
    h[] = exp(-(x*x + y*y)/(R0*R0));
#endif

  double max = 0.;
  for (int l = 0; l < depth(); l++)
    foreach_halo (prolongation, l) 
      foreach_child() {
        double e = exp(-(x*x+y*y)/(R0*R0)) - h[];
	printf ("%g %g %d %d %g %g\n", x, y, level, cell.neighbors, h[], e);
	if (fabs(e) > max)
	  max = fabs(e);
      }

  fprintf (stderr, "maximum error on halos: %g\n", max);

  return (max > tolerance);
}
