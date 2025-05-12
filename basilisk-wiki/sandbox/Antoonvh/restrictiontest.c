/**
Does restriction work on trees?
*/
#include "grid/bitree.h" //May be a quadtree or octree.
scalar s[];

int main() {
  init_grid (16);
  unrefine (x < 0.5);
  foreach()
    s[] = 1.;
  restriction ({s});
  foreach_cell()
    printf ("%g %d %g\n",x, level, s[]);
  
  FILE * fp = fopen("boundary", "w");
  boundary ({s});  
  foreach_cell()
    fprintf (fp, "%g %d %g\n",x, level, s[]);
}
/**
~~~gnuplot No, it does not
reset
set yr [-0.5:4.5]
set xr [-0.1:1.1]
set key off
set size square
set xlabel 'x'
set ylabel 'level'
set size ratio 1
plot 'out' using 1:2 , \
     '' using 1:2:(sprintf("    %.0f",$3)) with labels
~~~

~~~gnuplot boundary behaves as intended
reset
set yr [-0.5:4.5]
set xr [-0.1:1.1]
set xlabel 'x'
set ylabel 'level'
set key off
set size square
set size ratio 1
plot 'boundary' using 1:2 , \
     '' using 1:2:(sprintf("    %.0f",$3)) with labels
~~~
*/