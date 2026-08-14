/** 
# Testing `output_vtkhdf()` and `output_vtkhdf_slice()`

Tests VTKHDF output functions

~~~pythonplot
import sys, os
if os.path.exists('test_output_vtkhdf.py'):
  sys.path.append('.')
else:
  sys.path.append('../')
import test_output_vtkhdf
data = test_output_vtkhdf.print_vtkhdf_info('domain.hdf')
~~~

*/

#include "acastillo/output_fields/vtkhdf/output_vtkhdf.h"
#define MAXLEVEL 4
#define ASPECTRATIO 8
#define r2 (sq(x) + sq(y))

int main(){

  #if !TREE
    #if dimension == 2
      dimensions (nx = 1, ny = ASPECTRATIO);
    #else
      dimensions (nx = 1, ny = ASPECTRATIO, nz = 1);
    #endif
  #endif

  L0 = 1.0;
  X0 = Y0 = Z0 = -L0 / 2;
  N = 1 << (MAXLEVEL-1);
  init_grid(N);

  #if TREE
    double outer_radius = 0.25;
    double inner_radius = 0.1 ;
    refine(((r2 < sq(outer_radius)) && (r2 > sq(inner_radius))) && level < MAXLEVEL);
  #endif

  // Create test fields 
  scalar f[], p[];
  vector u[];
  foreach(){
    f[] = point.level;
    p[] = pid();
    u.x[] = x;
    u.y[] = y;
    #if dimension == 3
      u.z[] = z;
    #endif
  }

  // And write a vtkhdf file
  output_vtkhdf({f,p}, {u}, "domain.hdf");

#if dimension > 2
  output_vtkhdf_slice({f,p}, {u}, "slice_x.hdf", (coord){1,0,0}, 0);
  output_vtkhdf_slice({f,p}, {u}, "slice_y.hdf", (coord){0,1,0}, 0);
  output_vtkhdf_slice({f,p}, {u}, "slice_z.hdf", (coord){0,0,1}, 0);
#endif

  /**
  Verify the file that was just written against the analytical fields set
  above: `u` must reproduce the centroid of each cell recomputed from the
  connectivity, and `f` must match the refinement level. The report goes to
  stderr, i.e. to the `log` that Basilisk diffs against
  `test_output_vtkhdf.ref`, so this is what makes the test pass or fail.

  The reference was recorded with 4 MPI ranks, which is what `run_tests.py`
  uses. It cannot be made rank-independent: the file holds one partition per
  rank, so the number of points depends on the decomposition. */
  if (pid() == 0) {
    system ("python3 ../test_output_vtkhdf.py --check domain.hdf 1>&2");
#if dimension > 2
    /** The slices are checked against the plane they were asked for, in
    addition to the analytical fields. */
    system ("python3 ../test_output_vtkhdf.py --check --plane=x=0 slice_x.hdf 1>&2");
    system ("python3 ../test_output_vtkhdf.py --check --plane=y=0 slice_y.hdf 1>&2");
    system ("python3 ../test_output_vtkhdf.py --check --plane=z=0 slice_z.hdf 1>&2");
#endif
  }
}
