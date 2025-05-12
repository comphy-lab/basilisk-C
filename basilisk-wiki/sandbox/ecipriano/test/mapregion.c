/**
# Map Region

We want to test the function *mapregion()*, which calculates
an Heaviside function that maps regions of the domain close
to the interface. This test case will not change the world,
but it will remind me how the *mapregion()* function behaves. */

#include "utils.h"
#include "mapregion.h"
#include "view.h"

double mEvapVal = 0.;

scalar H[], f[];

void write_picture (char* name) {
  clear();
  draw_vof ("f", lw = 1.5);
  squares ("H", spread = -1);
  cells();
  save (name);
}

int main (void) {
  origin (-0.5*L0, -0.5*L0);
  init_grid (1 << 5);

  /**
  **Case 1**: Inverse is false, Narrow is false, 0 Layers.
  */

  fraction (f, sq(0.2) - sq(x) - sq(y));
  mapregion (H, f, nl=0, inverse=false, narrow=false);
  write_picture ("case1.png");
  fprintf (stderr, "case1 = %g\n", statsf(H).sum);

  /**
  **Case 2**: Inverse is true, Narrow is false, 0 Layers.
  */

  fraction (f, sq(0.2) - sq(x) - sq(y));
  mapregion (H, f, nl=0, inverse=true, narrow=false);
  write_picture ("case2.png");
  fprintf (stderr, "case2 = %g\n", statsf(H).sum);

  /**
  **Case 3**: Inverse is true, Narrow is false, 1 Layers.
  */

  fraction (f, sq(0.2) - sq(x) - sq(y));
  mapregion (H, f, nl=1, inverse=true, narrow=false);
  write_picture ("case3.png");
  fprintf (stderr, "case3 = %g\n", statsf(H).sum);

  /**
  **Case 4**: Inverse is false, Narrow is true, 2 Layers.
  */

  fraction (f, sq(0.2) - sq(x) - sq(y));
  mapregion (H, f, nl=2, inverse=false, narrow=true);
  write_picture ("case4.png");
  fprintf (stderr, "case4 = %g\n", statsf(H).sum);

  /**
  **Case 5**: Inverse is true, Narrow is true, 2 Layers.
  */

  fraction (f, sq(0.2) - sq(x) - sq(y));
  mapregion (H, f, nl=2, inverse=true, narrow=true);
  write_picture ("case5.png");
  fprintf (stderr, "case5 = %g\n", statsf(H).sum);

  /**
  **Case 6**: Inverse is true, Narrow is false, 1 Layers, Nointerface
  is true.
  */

  fraction (f, sq(0.2) - sq(x) - sq(y));
  mapregion (H, f, nl=1, inverse=true, narrow=false,
      nointerface=true);
  write_picture ("case6.png");
  fprintf (stderr, "case6 = %g\n", statsf(H).sum);
}

/**
## Results

![Case 1: Inverse is false, Narrow is false, 0 Layers](mapregion/case1.png)

![Case 2: Inverse is true, Narrow is false, 0 Layers](mapregion/case2.png)

![Case 3: Inverse is true, Narrow is false, 1 Layers](mapregion/case3.png)

![Case 4: Inverse is false, Narrow is true, 2 Layers](mapregion/case4.png)

![Case 5: Inverse is true, Narrow is true, 2 Layers](mapregion/case5.png)

![Case 6: Inverse is true, Narrow is true, 2 Layers, Nointerface is true](mapregion/case6.png)
*/
