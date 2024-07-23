/**
# Tangential interpolation on face vector fields

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
face vector u[];

double R0 = 0.1;
u.n[right] = dirichlet(exp(-(x*x + y*y)/(R0*R0)));
u.n[top] = dirichlet(exp(-(x*x + y*y)/(R0*R0)));

int main (int argc, char ** argv)
{
  int n = 2048;
  init_grid (n);

  origin (-1, -1);
  foreach()
    h[] = exp(-(x*x + y*y)/(R0*R0));
  
  /* initial coarsening (see halo.c) */
  double tolerance = 2e-4;
  while (adapt_wavelet ({h}, &tolerance, 11, list = {h}).nc)
    ;

  foreach_face(x) u.x[] = exp(-(x*x + y*y)/(R0*R0));
  foreach_face(y) u.y[] = exp(-(x*x + y*y)/(R0*R0));
  //  output_cells (stdout);

  double max = 0;
  foreach_face (x)
    for (int i = -1; i <= 1; i++) {
      double xu = x, yu = y + i*Delta;
      double e = exp(-(xu*xu+yu*yu)/(R0*R0)) - u.x[0,i];
      if (fabs(e) > max)
	max = fabs(e);
    }

  double maxv = 0;
  foreach_face (y)
    for (int i = -1; i <= 1; i++) {
      double xv = x + i*Delta, yv = y;
      double e = exp(-(xv*xv+yv*yv)/(R0*R0)) - u.y[i,0];
      if (fabs(e) > maxv)
	maxv = fabs(e);
      printf ("%g %g %d %d %g %g\n", 
	      xv, yv, level, cell.neighbors, u.y[i,0], e);
    }

  fprintf (stderr, "maximum error on halos: %g %g\n", max, maxv);

  return (max > tolerance || maxv > tolerance || max != maxv);
}
