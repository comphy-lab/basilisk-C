/**
# The interactive bview does not recognize most keys to `save()`.

It does not run on the website. The terminal reads:

~~~
# t = 0, fields = { s }
#         name:          min          avg       stddev          max
#            s:            0            0            0            0
unknown key 'sort'
~~~

These images are only generated when `bview2D` is available.

![The result](bv/nja.png)

![Is sorting taken care of?](bv/nja.svg)
*/
#include "utils.h"

scalar s[];

int main() {
  FILE * fp = fopen ("vsd.bv", "w");
  fprintf (fp,
	   "view (width = 200, height = 200, theta = 0.6, tx = -0.6, ty = -0.5);\n"
	   "begin_translate (z = -0.05)\n"
	   "cells(lw = 3);\n"
	   "end_translate()\n"
	   "squares (\"s\", map = cool_warm);\n"
	   "save (\"nja.svg\", sort = GL2PS_BSP_SORT)\n"
	   "save (\"nja.png\")\n"
           "stop()\n" 
	   ); //The last comment does not exist
  fflush (fp);
  init_grid(2);
  dump("simpledump");
  system ("bview2D simpledump vsd.bv");
}
