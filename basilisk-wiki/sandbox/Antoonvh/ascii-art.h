/**
# Live `Ascii`-art visualization

Plot grey-scale Ascii art in the termnial based on scalar field data.
*/
#include "utils.h"
#include <sys/ioctl.h>
#define SCALE " .:-=+*#%@"
#define NSCALE (9)
#define clampd(x,min,max) (x < min ? min: x > max ? max : x)
#define IND(v, min, max) (int)(NSCALE*(clampd(v, min, max) - min)/(max - min) + 0.5)

struct Aart {
  scalar s;
  FILE * fp;
  int n;
  double min, max, spread, z;
  bool linear, clear; //interp and clear screen
};

void ascii_art (struct Aart p) {
  scalar s = p.s;
  if (p.min == 0 && p.max == 0) {
    stats s = statsf (p.s);
    if (p.spread < 0.)
      p.min = s.min, p.max = s.max;
    else {
      double avg = s.sum/s.volume, spread = (p.spread ? p.spread : 5.)*s.stddev;
      p.min = avg - spread; p.max = avg + spread;
    }
  }
  if (!p.fp)
    p.fp = stdout;
  if (!p.n) {
    if (p.fp == stdout || p.fp == stderr) {
      struct winsize w;
      ioctl(0, TIOCGWINSZ, &w); 
      p.n = min(2*(w.ws_row) - 2, w.ws_col);
    } else
      p.n = N;
  }
  double Delt = L0/p.n;
  int j = 0;
  if (p.clear) {
    while (j++ < p.n/2) {
      fprintf(p.fp, "\033[A");
    }
  } else
    puts("\x1b[2J");
  for (double yp = Y0 + L0 - Delt/2.; yp > Y0; yp -= 2.*Delt) {
    for (double xp = X0 + Delt/2; xp < X0 + L0; xp += Delt) {
      Point point = locate (xp, yp, p.z);
      if (point.level > 0) {
	double val = p.linear ? interpolate(s, xp, yp, p.z) : s[];
	fputc(SCALE[IND(val,p.min,p.max)], p.fp);
      }
#if _MPI
      MPI_Barrier(MPI_COMM_WORLD);
#endif
    }
    fputc('\n', p.fp);
  }
  fflush (p.fp);
}
/**
## Example

* [A dipole-wall collision](taa.c)
 */
