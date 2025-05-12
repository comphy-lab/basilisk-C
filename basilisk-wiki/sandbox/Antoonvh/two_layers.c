/**
# Define two layers of boundary cells

This is a rather ugly solution. The boundary number idintifier `b` is
out of scope, so each dimension gets its own formulation.

Note that this does not work for Cartesian grids, use MG or Tree instead.
*/
#define BGHOSTS 2

@define n_x neighbor.i
@define n_y neighbor.j
@define n_z neighbor.k
foreach_dimension() {
  @define layer_nr_x (n_x < GHOSTS ? (GHOSTS - n_x) : n_x - (1 << level) - GHOSTS + 1) 
  @define dirichlet_x(a) (val(_s,0,0,0) + layer_nr_x*2.*(a - val(_s,0,0,0)))
  @define neumann_x(a)   (val(_s,0,0,0) + layer_nr_x*Delta*a)
}

int main() {
  init_grid(2);
  scalar s[];
 /**
 The price to pay is that you must use `_x` etc. eventough the boundary's name already implies this info.  
 */
  s[left]   = neumann_x   (1.);
  s[right]  = dirichlet_x (1.);
  s[bottom] = dirichlet_y (2.);
  s[top]    = neumann_y   (1.);
  
  foreach()
    s[] = 1.;
  boundary({s});
  /**
  We print the cells and their values:
  */
  output_cells (stdout);
  
  foreach()
    fprintf (stderr, "%g %g %g\n", x, y, s[]);
  
  foreach_boundary(left)
    for (int i = -BGHOSTS; i < 0; i++)
      fprintf (stderr, "%g %g %g\n", x + (i + 0.5)*Delta, y, s[i]);
  
  foreach_boundary(right)
    for (int i = 1; i <= BGHOSTS; i++)
      fprintf (stderr, "%g %g %g\n", x + (i - 0.5)*Delta, y, s[i]);
    
  foreach_boundary(bottom)
    for (int i = -BGHOSTS; i < 0; i++)
      fprintf (stderr, "%g %g %g\n", x , y + (i + 0.5)*Delta, s[0,i]);
    
  foreach_boundary(top)
    for (int i = 1; i <= BGHOSTS; i++)
      fprintf (stderr, "%g %g %g\n", x, y + (i - 0.5)*Delta, s[0,i]);
  
}

/**
~~~gnuplot Second row of ghost cells seems OK
blue="#5082DC"
set terminal @PNG enhanced size 500,500 font "times,20"
set output 'boundary.png'
unset key 
unset border
unset tics
unset xlabel
unset ylabel
set size ratio -1

plot 'out' w l lw 2 lc rgb "#7F7F7F", \
     'log' u 1:2:3 with labels tc rgb blue
~~~
*/
