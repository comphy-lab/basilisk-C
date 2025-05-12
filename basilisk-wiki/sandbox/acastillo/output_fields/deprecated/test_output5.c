#include "grid/octree.h"
#include "FEDataStructures.h"
#include "CatalystAdaptor.h"

#define r2 (sq(x) + sq(y))
int main(int argc, char* argv[]){
	L0 = 1 ;
	X0 = Y0 = Z0 = -L0/2;
  N = 1<<4;
  init_grid (N);

#if TREE
 refine(((r2 < sq(0.25)) && (r2 > sq(0.1))) && level < 6);
#endif

  scalar s[];
  vector u[];
  foreach ()
  {
    s[] = pid();
    foreach_dimension()
    {
      u.x[] = sq(x) + sq(y) + sq(z);
    }
  }
  boundary({s, u});

  CAttributes cattributes;
  CGrid cgrid = (CGrid){ .no_points = 0, .points = 0, .no_cells = 0, .cells = 0};

  InitializeGrid(&cgrid);
  InitializeAttributes(&cattributes, &cgrid);

#ifdef USE_CATALYST
  do_catalyst_initialization(argc, argv);
#endif

  double time = 0.;
#ifdef USE_CATALYST
  for (int i = 0; i < 10; i++) {
    time = i * 0.1;
    foreach ()
      s[] = 1/(sq(x) + sq(y) + sq(z)) + 0.1*noise();
    UpdateFields(&cattributes, time, s, u);
    do_catalyst_execute(i, time, &cgrid, &cattributes);
  }
#endif


#ifdef USE_CATALYST
  do_catalyst_finalization();
#endif
}
