/**
# Unreasonably long compilation times

This was fixed by [this patch](/src/?changes=20231101171232), but is
still quite long on the server (approx. 25 seconds and only 2.5
seconds on my i7 laptop).

The compilation of the small program below takes very long, and this
behaviour depends on some unexpected aspects.

1) The number of loop iterations (`imax`) changes the time, altough it saturates.
2) The repeated calls to `my_fun()` are also important
3) Setting the `_FORTIFY_SOURCE=2` flag makes this issue (much) worse

See also [here](/sandbox/Antoonvh/tcbar.c).
*/

#include "utils.h"

void my_fun (Colormap map) {
  int imax = 1000; 
  double cmap [NCMAP][3];
  map (cmap);
  for (int i = 0; i < imax; i++)
    colormap_color (cmap, (float)i/(imax - 1), 0, 1);
}

int main() {
  int imax = 1000;
  for (int i = 0; i < imax; i++) {
    // Multiple calls
    my_fun (jet);
    my_fun (jet);
    my_fun (jet);
    my_fun (jet);
    my_fun (jet);
  }
}
  
