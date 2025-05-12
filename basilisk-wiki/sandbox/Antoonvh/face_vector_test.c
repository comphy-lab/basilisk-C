/**
# Is it an issue to initialize as 3D `const face vector` in 2D?
*/

int main() {
  init_grid (2);
  const face vector f[] = {1., 2., 3.};
  foreach_face()
    fprintf (stderr, "%g %g %g\n", x - 0.03, y, f.x[]);
  output_cells (stdout);
}
/**
The result ...
~~~gnuplot ... looks good
blue="#5082DC"
set terminal svg enhanced size 400,400 font "times,30"
set output 'boundary.svg'
unset key 
unset border
unset tics
unset xlabel
unset ylabel
set size ratio -1

plot 'out' w l lw 2 lc rgb "#7F7F7F", \
     'log' u 1:2:3 with labels tc rgb blue
~~~

No need for pre-processor commands, 

~~~c
#if dimension < 2
const face vector f[] = {1.};
#elif dimension < 3
const face vector f[] = {1., 2.};
#else
const face vector f[] = {1., 2., 3.};
#endif
This code was *not* compiled ...
~~~
*/
