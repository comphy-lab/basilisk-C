/**
##fraction demostration
Here we shoe some exemples can do with fraction 
*/
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "view.h"  


#define eq1(x, y) ( sq(3.) - sq(y) - sq(x)) // r large
#define eq2(x, y) ( sq(x) + sq(y) - sq(2.)) // rsmall

scalar f1[],f2[],f3[];

int main()
{
  origin (0 , -L0/2.);
  N = 256;
  init_grid(N);

  fraction (f3, eq1(x,y)*eq2(x,y));
  clear();
  box();
  draw_vof ("f3");
  squares ("f3", min = 0, max = 1);
  save ("f3.png");
  run();
}



/**
##Some little trick for get a less mesh number for a complexe VOF with little volume
*/
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "view.h"  
scalar f0[];
int main()
{
  init_grid(4);
  char name[10]; 

  for (int ii=5; ii<=10; ii++){	
    fraction (f0, 0.25 - sq(y)-sq(x));
    adapt_wavelet ({f0}, (double[]){0.00001}, ii);	
    //refine (f0[] > 0. && f0[]<1. && level < ii);
    sprintf (name, "lv_%d.png",ii);  	
    clear();
    view( tx = -0.5, ty = -0.5);
    box();
    draw_vof ("f0");
    squares ("f0", min = 0, max = 1);
    cells();
    save (name);
  }
  run();
}

int main()
{
  // origin (0 , -L0/2.);
  init_grid(4);
  //refine (level < 7);
  char name[10];
  char nameold[10];

  for (int ii=5; ii<=10; ii++){	
    fraction (f0, 0.25 - sq(y)-sq(x));
    adapt_wavelet ({f0}, (double[]){0.00001}, ii);	
    //refine (f0[] > 0. && f0[]<1. && level < ii);
    sprintf (name, "lv_%d.png",ii);  	
    clear();
    view( tx = -0.5, ty = -0.5);
    box();
    draw_vof ("f0");
    squares ("f0", min = 0, max = 1);
    cells();
    save (name);
  }
  run();
}
