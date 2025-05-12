/**
Checks boundary conditions for vertex fields in parallel (in 3D). */

#include "grid/octree.h"

#define shape(x,y,z) (sq(x) + sq(y) + sq(z))

int main()
{
  size (16);
  init_grid (4);
  origin (-L0/2., -L0/2., -L0/2.);
  refine (sq(x) + sq(y) + sq(z) < sq(L0/4.) && level < 3);

  scalar s[];
  foreach_dimension() {
    s[left] = dirichlet(shape(x,y,z));
    s[right] = dirichlet(shape(x,y,z));
  }
  foreach()
    s[] = shape(x,y,z);
  
  vertex scalar phi[];
  foreach_vertex()
    phi[] = (s[] + s[-1] + s[0,-1] + s[-1,-1] +
	     s[0,0,-1] + s[-1,0,-1] + s[0,-1,-1] + s[-1,-1,-1])/8.;

  foreach() {
    fprintf (qerr, "%g %g %g %g\n", x - Delta/3., y - Delta/3., z - Delta/3.,
	     phi[]);
    fprintf (qerr, "%g %g %g %g\n", x - Delta/3., y + Delta/3., z - Delta/3.,
	     phi[0,1]);
    fprintf (qerr, "%g %g %g %g\n", x + Delta/3., y - Delta/3., z - Delta/3.,
	     phi[1,0]);
    fprintf (qerr, "%g %g %g %g\n", x + Delta/3., y + Delta/3., z - Delta/3.,
	     phi[1,1]);
    fprintf (qerr, "%g %g %g %g\n", x - Delta/3., y - Delta/3., z + Delta/3.,
	     phi[0,0,1]);
    fprintf (qerr, "%g %g %g %g\n", x - Delta/3., y + Delta/3., z + Delta/3.,
	     phi[0,1,1]);
    fprintf (qerr, "%g %g %g %g\n", x + Delta/3., y - Delta/3., z + Delta/3.,
	     phi[1,0,1]);
    fprintf (qerr, "%g %g %g %g\n", x + Delta/3., y + Delta/3., z + Delta/3.,
	     phi[1,1,1]);
  }

#if _MPI
  debug_mpi (NULL);
#else
  output_cells();
#endif  
}
