/**
# Test My Quadtree implementation

*/
#include "myQT.h" 
scalar s[];
int main() {
  L0 = 64;
  X0 = Y0 = -L0/2;
  init_grid (4);
  scalar a[];
  foreach()
    a[] = 0;
  refine (sq(x) + sq(y) < sq(L0/5.) && level < 4);
  output_cells();

  // Iterate faces:
  FILE * fp1 = fopen ("faces", "w");
  foreach_face(x) {
    fprintf (fp1, "%g %g %g 0\n", x, y, L0/20);
  }
  foreach_face(y)
    fprintf (fp1, "%g %g 0 %g\n", x, y, L0/20);
  fclose (fp1);
  /**
The face iterator should traverse over te coarse faces at resolution
boundaries.

~~~gnuplot it works
set key off
set size square
plot 'out' w l lw 2, 'faces' u ($1-$3/2):($2-$4/2):3:4 with vectors head filled
~~~
   */
  // All faces

  for (int l = depth(); l >=0; l--) {
    char fname[99];
    sprintf (fname, "Level%d", l);
    FILE * fp5 = fopen (fname, "w");
    foreach_face_level(l) {
      coord cc = {x, y, z};
      coord dir = {0};
      foreach_dimension() { // It does not rotate itself!  
	is_face_x() {
	  cc.x -= Delta/2;
	  dir.x = Delta/2.5;
	}
      }
      fprintf (fp5, "%g %g %g %g\n", cc.x, cc.y, dir.x, dir.y);
    }
    fclose (fp5);
  }
  /**
Active faces on on levels

~~~gnuplot They are there
set term svg size 800, 800
set key top right
plot 'out' w l lw 2, \
 'Level4' u ($1-$3/2):($2-$4/2):3:4 w vectors head filled lw 2 t 'L4',\
 'Level3' u ($1-$3/2):($2-$4/2):3:4 w vectors head filled lw 3 t 'L3',\
 'Level2' u ($1-$3/2):($2-$4/2):3:4 w vectors head filled lw 4 t 'L2',\
 'Level1' u ($1-$3/2):($2-$4/2):3:4 w vectors head filled lw 5 t 'L1',\
 'Level0' u ($1-$3/2):($2-$4/2):3:4 w vectors head filled lw 6 t 'L0'
~~~

   */
  
  // Count faces
  FILE * fp3 = fopen ("nrs", "w");
  int nf = 0;
  foreach_face()
    fprintf (fp3, "%g %g %d\n", x , y, nf++);
  fclose (fp3);

  /**
Counting the "leaf faces":

~~~gnuplot No faces should be iterated twice
set key off
plot 'out' w l lw 2, 'nrs' u 1:2:3 w labels
~~~
   */  
#define FUNC ((x - 10) + 2*y)  
  // boundary_flux
  face vector fv[];
  fv.t[bottom] = dirichlet (FUNC);
  fv.t[top] = dirichlet (FUNC);
  fv.t[left] = dirichlet (FUNC);
  fv.t[right] = dirichlet (FUNC);
  foreach_face()
    fv.x[] = FUNC;
  flux_prolongate({fv});
  FILE * fp = fopen ("dataxy", "w");
  foreach()  {
    fprintf (fp, "%g %g %.3g %.3g\n", x, y, (fv.x[0] + fv.x[1])/2  , FUNC);
    fprintf (fp, "%g %g %.3g %.3g\n", x, y, (fv.y[0] + fv.y[0,1])/2, FUNC);
  }
  fclose (fp);
  /**
  `flux_plongate()` should prolongate the values;

~~~gnuplot The four numbers in each cells should be equal
set key off
set size square
plot 'dataxy' u 1:($2):3 w labels, '' u 1:($2):4 w labels, 'out' w l lw 2
~~~

When consistent box boundaries are set, it works.
  */
  

  FILE * fp2 = fopen ("grads", "w");

  s[left]   = dirichlet (FUNC);
  s[right]  = dirichlet (FUNC);
  s[bottom] = dirichlet (FUNC);
  s[top]    = dirichlet (FUNC);
  
  foreach()
    s[] = FUNC;
  vector ds[];
  my_gradients ({s}, {ds});
  foreach()
    fprintf (fp2, "%g %g %g %g\n", x, y, ds.x[], ds.y[]);
  fclose (fp2);
  /**
  We also test the my_gradients function

~~~gnuplot 1s over 2s are good
plot 'grads' u 1:($2+1.5):3 w labels, '' u 1:($2-1.5):4 w labels, 'out' w l lw 2
~~~
  */
}
