/**
  This example shows that the minlevel argument (say minlevel=MINLEVEL) in adapt_wavelet does not behave as expected if the grid was previously at lower level than MINLEVEL. We show this by initializing the grid to some low level INILEVEL and then applying adapt_wavelet with minlevel=MINLEVEL where MINLEVEL > INILEVEL. */  

#include "fractions.h"
#include "view.h"

/**
   Define maximum, minimum, and initial grid levels. */
#define MAXLEVEL 7
#define MINLEVEL 6
#define INILEVEL 4

/**
   Set the scalar field... */
scalar f[];

int main()
{
/**
  ...whose attributes we set to the usual vof requirements. */
  f.prolongation = f.refine = fraction_refine;
  L0 = 1.;
  origin (-L0/2., -L0/2.);
/**
  Initialize the grid to the initial level, which is lower than the "minimum" level we will use in the following adaptivity statement. */
  init_grid (1 << INILEVEL);
  view (width = 700, height = 700);
/**
  Set up a circle for field f. */
  fraction (f, sq(x) + sq(y) - sq(0.25));
  boundary ({f});
/**
  Refine to a maximum level of MAXLEVEL and minimum of MINLEVEL. */
  adapt_wavelet ({f}, (double[]){1e-3}, MAXLEVEL, MINLEVEL);
  draw_vof ("f", lc = {1., 0., 1.}, lw=3);
  cells ();
  squares ("level", min = INILEVEL, max = MINLEVEL);
  save ("output.png");
}
/**
  The colouring makes it clear that MINLEVEL is not uniformly achieved in the grid as would be expected. In fact MINLEVEL is nowhere achieved.
  ![Compare](mingrid/output.png)
    */

