/**
# Test of vertex stencils 

~~~gnuplot reset field
set size ratio -1
set xrange [-0.5:1.5]
set yrange [-0.5:1.5]
set grid
unset key
plot for[i = 0:3] "log-".i  u 1:2:3 i 0 w labels
~~~

~~~gnuplot foreach_vertex
set size ratio -1
set xrange [-0.5:1.5]
set yrange [-0.5:1.5]
set grid
unset key
plot for[i = 0:3] "log-".i  u 1:2:3 i 1 w labels
~~~

~~~gnuplot foreach_vertex (different values on each proc)
reset
set size ratio -1
set xrange [-0.5:1.5]
set yrange [-0.5:1.5]
set grid
unset key
plot for[i = 0:3] "log-".i  u 1:2:3 i 2 w labels
~~~

~~~gnuplot Restriction
reset
set size ratio -1
set xrange [-0.5:1.5]
set yrange [-0.5:1.5]
set grid
unset key
plot for[i = 0:3] "log-".i  u 1:2:3 i 3 w labels
~~~

~~~gnuplot Prolongation
reset
set size ratio -1
set xrange [-0.5:1.5]
set yrange [-0.5:1.5]
set grid
unset key
plot for[i = 0:3] "log-".i  u 1:2:3 i 4 w labels
~~~

~~~gnuplot use of point
reset
set size ratio -1
set xrange [-0.5:1.5]
set yrange [-0.5:1.5]
set grid
unset key
plot for[i = 0:3] "log-".i  u 1:2:3 i 5 w labels
~~~

*/

#include "grid/multigrid.h"
#include "bderembl/libs/inner-vertex.h"
#include "Antoonvh/my_vertex.h"


/**
   Ideally one should use this function to print. See problem in last test
*/

static inline void print_scalar_stencil (Point point, scalar s){

    fprintf (qerr, "%g %g %g\n", x, y, s[]);
    fprintf (qerr, "%g %g %g\n", x + Delta, y, s[1]);
    fprintf (qerr, "%g %g %g\n", x - Delta, y, s[-1]);
    fprintf (qerr, "%g %g %g\n", x, y + Delta, s[0,1]);
    fprintf (qerr, "%g %g %g\n", x, y - Delta, s[0,-1]);

}

int main()
{
  init_grid (4);

  vertex scalar s[];

  s[left]   = 10.;
  s[right]  = 10.;
  s[top]    = 10.;
  s[bottom] = 10.;
  
  s.restriction = restriction_coarsen_vert;
  s.prolongation = refine_vert;

/**
   Test 1: reset field and boundary
*/

  reset ({s}, 5.);
  
  foreach_vertex() {
    fprintf (qerr, "%g %g %g\n", x, y, s[]);
    fprintf (qerr, "%g %g %g\n", x + Delta, y, s[1]);
    fprintf (qerr, "%g %g %g\n", x - Delta, y, s[-1]);
    fprintf (qerr, "%g %g %g\n", x, y + Delta, s[0,1]);
    fprintf (qerr, "%g %g %g\n", x, y - Delta, s[0,-1]);
  }   

/**
   Test 2: foreach_vertex
*/
  fprintf(qerr,"\n\n");
  reset ({s}, 0.);

  foreach_vertex()
    s[] = x + y;
  
  foreach_vertex() {
    fprintf (qerr, "%g %g %g\n", x, y, s[]);
    fprintf (qerr, "%g %g %g\n", x + Delta, y, s[1]);
    fprintf (qerr, "%g %g %g\n", x - Delta, y, s[-1]);
    fprintf (qerr, "%g %g %g\n", x, y + Delta, s[0,1]);
    fprintf (qerr, "%g %g %g\n", x, y - Delta, s[0,-1]);
  }   


/**
   Test 3: foreach_vertex with different values at edges of domains
*/

  fprintf(qerr,"\n\n");
  reset ({s}, 0.);

  foreach_vertex()
    s[] = pid();

  foreach_vertex() {
    fprintf (qerr, "%g %g %g\n", x, y, s[]);
    fprintf (qerr, "%g %g %g\n", x + Delta, y, s[1]);
    fprintf (qerr, "%g %g %g\n", x - Delta, y, s[-1]);
    fprintf (qerr, "%g %g %g\n", x, y + Delta, s[0,1]);
    fprintf (qerr, "%g %g %g\n", x, y - Delta, s[0,-1]);
  }   


/**
   Test 4: Restriction. Vertex restriction is not preserving area integral
*/

  fprintf(qerr,"\n\n");
  reset ({s}, 2.);

  restriction({s});

  boundary_level({s}, 0); // needed here

  foreach_vertex_level(0) {
    fprintf (qerr, "%g %g %g\n", x, y, s[]);
    fprintf (qerr, "%g %g %g\n", x + Delta, y, s[1]);
    fprintf (qerr, "%g %g %g\n", x - Delta, y, s[-1]);
    fprintf (qerr, "%g %g %g\n", x, y + Delta, s[0,1]);
    fprintf (qerr, "%g %g %g\n", x, y - Delta, s[0,-1]);
  }   


/**
   Test 5: Prolongation. Integral is preserved with refine_vert
*/
  fprintf(qerr,"\n\n");
  reset ({s}, 1.);

/**
   At level 0, a vertex field has 1 interior point.
*/
  boundary_level({s}, 0);      


  for (int l = 0; l < depth(); l++) {
    foreach_vertex_level(l)
      refine_vert (point, s);
    boundary_level({s}, l+1);
  }

  foreach_vertex() {
    fprintf (qerr, "%g %g %g\n", x, y, s[]);
    fprintf (qerr, "%g %g %g\n", x + Delta, y, s[1]);
    fprintf (qerr, "%g %g %g\n", x - Delta, y, s[-1]);
    fprintf (qerr, "%g %g %g\n", x, y + Delta, s[0,1]);
    fprintf (qerr, "%g %g %g\n", x, y - Delta, s[0,-1]);
  }

/**
   Bonus test: use of point
   This test is producing strange results:

- BC are not applied as it is done in test 1
- print occurs at cell center

*/
  fprintf(qerr,"\n\n");
  reset ({s}, 5.);
  
  foreach_vertex() {
    print_scalar_stencil(point, s);
  }   


}
