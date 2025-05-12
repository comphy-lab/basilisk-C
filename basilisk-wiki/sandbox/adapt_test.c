#include "terrain.h"
#include "saint-venant.h"
#define MAXLEVEL 11
#define MINLEVEL 4
#define ETAE 5.e-3

scalar slide[];
double t_end=0;

struct Layers {
  int (* iterate) (void);
};

void Init_Layers (struct Layers p) {
    terrain(slide,"/data/EMILY_BASILISK/terrain/slide/gers_refined",NULL);
  do { 
    foreach() {
      h[]= ( val(slide.n,0,0)<1.e-10 ? 0. : slide[]); 
      if (val(slide.n,0,0)<1.e-10) {
      };
      eta[] = zb[] + h[];
    }
  } while (p.iterate && p.iterate());
}

int main()
{

  #if QUADTREE
  // 32^2 grid points to start with 
  init_grid(1 << MINLEVEL);
  #else // Cartesian
  // 1024^2 grid points
  init_grid(1 << MAXLEVEL);
  #endif

// the domain is
  size (30000.);
  origin(662000,191000);

  run();
}

int adapt() {
#if QUADTREE
  scalar eta[];
  foreach()
    eta[] = h[] > dry ? h[] + zb[] : 0;
  boundary ({eta});
  
  astats s = adapt_wavelet ({eta}, (double[]){ETAE},
			    MAXLEVEL, MINLEVEL);
  fprintf (stderr, "%% refined %d cells, coarsened %d cells\n", s.nf, s.nc);
  return s.nf;
#else // Cartesian
  return 0;
#endif
}
event init(i=0)
{
/**
The initial still water surface is at $z=0$ so that the water depth 
$h$ is... */
  conserve_elevation();
  foreach() {
   zb[] = 0.;
   eta[]=0.;
  }
  Init_Layers( iterate = adapt );
  boundary ({h,zb,eta}); 
}


/**
## Adaptivity
   
We apply our `adapt()` function at every timestep. */
event do_adapt (i++) adapt();

