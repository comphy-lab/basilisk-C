/**
# Counting connected regions of each phase

Tags the connected regions of one phase of a VOF field `f`, then computes for
each region its volume, centre of mass, volume-weighted average velocity and
interfacial area, and appends one line per region to a file.

Columns written, in 3D: `i  t  j  volume  bx  by  bz  ux  uy  uz  area`,
and in 2D: `i  t  j  volume  bx  by  ux  uy  area`

- `i`, `t`       : time-step index and time
- `j`            : region index (0-based)
- `volume`       : $\int f\,dv$ over the region
- `b`            : centre of mass
- `u`            : volume-weighted average velocity
- `area`         : interfacial area (perimeter in 2D), from the reconstructed
                   VOF planes

Uses the globals `f` and `u`. */

#include "tag.h"

/**
  - `m`         : pre-allocated tag scalar (recycled between calls to save memory)
  - `threshold` : VOF threshold separating the two phases
  - `above`     : 1 → count regions where f > threshold (droplets)
                  0 → count regions where f < threshold (bubbles)
  - `filename`  : output file (opened in append mode)
*/
void count_phase(scalar m, double threshold, int above, const char *filename, int i, double t) {

  foreach()
    m[] = above ? (f[] > threshold) : (f[] < 1. - threshold);
  int n = tag(m);

  double volume[n], b_x[n], b_y[n], vel_x[n], vel_y[n], area[n];
  memset(volume, 0, n * sizeof *volume);
  memset(b_x,    0, n * sizeof *b_x);
  memset(b_y,    0, n * sizeof *b_y);
  memset(vel_x,  0, n * sizeof *vel_x);
  memset(vel_y,  0, n * sizeof *vel_y);
  memset(area,   0, n * sizeof *area);
#if dimension == 3
  double b_z[n], vel_z[n];
  memset(b_z,    0, n * sizeof *b_z);
  memset(vel_z,  0, n * sizeof *vel_z);
#endif

/**
The accumulators all reduce the same way, so the clause is named once here
rather than spelled out on the `foreach` below -- which also keeps the only
2D/3D difference in the loop to a single place, since the `z` accumulators do
not exist in 2D. */

#if dimension == 3
#define COUNT_PHASE_REDUCTIONS						\
  reduction(+:volume[:n])						\
  reduction(+:b_x[:n])   reduction(+:b_y[:n])   reduction(+:b_z[:n])	\
  reduction(+:vel_x[:n]) reduction(+:vel_y[:n]) reduction(+:vel_z[:n])	\
  reduction(+:area[:n])
#else
#define COUNT_PHASE_REDUCTIONS						\
  reduction(+:volume[:n])						\
  reduction(+:b_x[:n])   reduction(+:b_y[:n])				\
  reduction(+:vel_x[:n]) reduction(+:vel_y[:n])				\
  reduction(+:area[:n])
#endif

  foreach (COUNT_PHASE_REDUCTIONS) {
    if (m[] > 0) {
      int j = m[] - 1;
      double fval = above ? f[] : (1. - f[]);
      volume[j] += dv() * fval;

      coord p = {x, y, z};
      foreach_dimension(){
        b_x[j] += dv() * fval * p.x;
        vel_x[j] += dv() * fval * u.x[];
      }

      if (fval > threshold && fval < 1. - threshold) {
        coord n = interface_normal (point, f), p;
        double alpha = plane_alpha (fval, n);
        area[j] += pow(Delta, dimension - 1) * plane_area_center (n, alpha, &p);
      }
    }
  }

  if (pid() == 0) {
    
    /**
    The file is appended to across calls, so the column header is written only
    when it is still empty, i.e. on the first call. */

    FILE *fp = fopen(filename, "a");
    if (ftell(fp) == 0)
#if dimension == 3
      fprintf(fp, "# i t j volume bx by bz ux uy uz area\n");
#else
      fprintf(fp, "# i t j volume bx by ux uy area\n");
#endif

    for (int j = 0; j < n; j++){
      fprintf(fp, "%d %g %d %g ", i, t, j, volume[j]);
#if dimension == 3
      fprintf(fp, "%g %g %g ",   b_x[j]/volume[j],   b_y[j]/volume[j],   b_z[j]/volume[j]);
      fprintf(fp, "%g %g %g ", vel_x[j]/volume[j], vel_y[j]/volume[j], vel_z[j]/volume[j]);
#else
      fprintf(fp, "%g %g ",   b_x[j]/volume[j],   b_y[j]/volume[j]);
      fprintf(fp, "%g %g ", vel_x[j]/volume[j], vel_y[j]/volume[j]);
#endif
      fprintf(fp, "%g\n", area[j]);
    }
    fflush(fp);
    fclose(fp);
  }
}
