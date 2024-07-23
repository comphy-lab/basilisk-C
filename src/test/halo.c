/**
# Definition of halo cells after coarsening 

~~~gnuplot
set size ratio -1
unset grid
unset xtics
unset ytics
unset border
set key outside
plot			  \
     'out' w l t "" lc 0, \
     '< grep halo log | awk "(\$3==4){print \$0}"' pt 7 lc 1 t 'halo 4', \
     '< grep halo log | awk "(\$3==3){print \$0}"' pt 7 lc 2 t 'halo 3'
~~~
*/

scalar h[];

int main (int argc, char ** argv)
{
  init_grid (32);

  origin (-0.5, -0.5);
  double R0 = 0.1;
  foreach()
    h[] = exp(-(x*x + y*y)/(R0*R0));
  
  /* initial coarsening */

  while (adapt_wavelet ({h}, (double []){1e-2}, 5).nc);

  for (int l = 0; l < depth(); l++)
    foreach_halo (prolongation, l)
      fprintf (stderr, "%g %g %d %d halo\n", x, y, level, cell.neighbors);
  output_cells (stdout);
}
