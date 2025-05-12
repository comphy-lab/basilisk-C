
/**

Run on clusters          : CFLAGS=-DCLUSTER=1 make clean testDumpStl.tst
For a restart            : CFLAGS="-DCLUSTER=1 -DRESTART=1" make clean testDumpStl.tst
Compilation for Parallel : The file Dump should be created with a sequential run
                           CC='mpicc -D_MPI=8' CFLAGS='-DCLUSTER=1 -DRESTART=1' make clean testDumpStl.tst
Run on workstation       : make clean mousse.tst
                           then use the navigator to open the html file to run it
*/

#include "grid/octree.h"
#include "embed.h"

#include "distance.h"
#include "curvature.h"
#include "view.h"
#include "run.h"

#if !defined(CLUSTER)
#include "display.h"
#endif

int maxlevel  = 8;

int main(){
  N=1<<7; /* default discretization of the domain */
  periodic(right);
  run();
}

event init (i = 0) {
  scalar d[];
  #if !defined(RESTART)
    printf("We load the stl file\n");
    coord * p = input_stl (fopen ("../outside_rot.stl", "r"));
    coord min, max;
    bounding_box (p, &min, &max);  
    double maxl = -HUGE;
    foreach_dimension()
      if (max.x - min.x > maxl)
    	maxl = max.x - min.x;
    size (maxl);

    X0 = (max.x + min.x)/2. - L0/2;
    Z0 = (max.z + min.z)/2. - L0/2;
    Y0 = (max.y + min.y)/2. - L0/2;

    distance (d, p);
  
    boundary ({d});
    vertex scalar phi[];  
    foreach_vertex()
      phi[] = -(d[] + d[-1] + d[0,-1] + d[-1,-1] + d[0,0,-1] + d[-1,0,-1] + d[0,-1,-1] + d[-1,-1,-1])/8.;
    
    boundary ({phi});
    fractions (phi, cs, fs);
    fractions_cleanup (cs, fs);

    dump(file = "dump", list={phi}); /* we save only phi */
  #endif
  
  #if defined(RESTART)
    vertex scalar phi[];
    fprintf(stdout,"In restore\n");
    if (!restore (file = "dump")){
      fprintf(stderr,"dump file not present\n");
      fprintf(stderr,"Simulation aborded\n");
      fprintf(stderr,"\n");
      fflush(stderr);fflush(stdout);
      return 1;
    }  
    
    boundary ({phi});
    fractions (phi, cs, fs);
    fractions_cleanup (cs, fs);       
  #endif

  while (adapt_wavelet ({cs}, (double[]){1e-6}, maxlevel).nf);

  scalar kappa[];
  curvature (cs, kappa);
  view (theta = 0.5, phi = 0.5);
  draw_vof("cs", "fs", color = "kappa");
  box();
  save ("kelvin.png");

  clear();
  view ();
  box();
  cells();
  save ("plan.png");
}
