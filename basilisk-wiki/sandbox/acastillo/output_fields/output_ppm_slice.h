#ifndef OUTPUT_PPM_SLICE_H
#define OUTPUT_PPM_SLICE_H

#include "output.h"

/**
## *output_ppm_slice()*: 2D Portable PixMap (PPM/PNG/MP4) slice of a 3D field

Exports a colormapped 2D image or animation frame of scalar *f* on an
axis-aligned *x*, *y*, or *z* plane defined by normal *n* and intercept
*alpha* (`n.x*x + n.y*y + n.z*z = alpha`), matching the convention of
`output_xmf_slice()` and `view.h`'s `squares()`:

* `n = {1, 0, 0}, alpha = x0` — `(y, z)` slice (`y` horizontal, `z` vertical)
* `n = {0, 1, 0}, alpha = y0` — `(x, z)` slice (`x` horizontal, `z` vertical)
* `n = {0, 0, 1}, alpha = z0` — `(x, y)` slice (`x` horizontal, `y` vertical)

Works identically on CPU, CPU+MPI, single-GPU CUDA, and Multi-GPU CUDA+MPI.
When `res <= 0` (the default), horizontal pixel resolution defaults to the
native grid cell count along the horizontal axis of the selected plane,
preserving physical aspect ratio across all three orientations.
*/

trace
void output_ppm_slice (scalar f,
		       char * file = NULL,
		       coord n = {0, 0, 1},
		       double alpha = 0,
		       int res = 0,
		       double min = 0, double max = 0, double spread = 5,
		       bool linear = false,
		       coord box[2] = {{0, 0, 0}, {0, 0, 0}},
		       scalar mask = {-1},
		       Colormap map = jet,
		       char * opt = NULL,
		       FILE * fp = stdout,
		       FILE * checksum = NULL)
{
  if (!min && !max) {
    stats s = statsf (f);
    if (spread < 0.)
      min = s.min, max = s.max;
    else {
      double avg = s.sum/s.volume;
      min = avg - spread*s.stddev;
      max = avg + spread*s.stddev;
    }
  }

  coord dom[2] = {
    {X0, Y0, Z0},
    {X0 + L0,
     Y0 + L0*Dimensions.y/Dimensions.x,
#if dimension > 2
     Z0 + L0*Dimensions.z/Dimensions.x
#else
     Z0
#endif
    }
  };

  bool has_box = (box[1].x > box[0].x) || (box[1].y > box[0].y)
#if dimension > 2
              || (box[1].z > box[0].z)
#endif
  ;
  coord b[2];
  if (has_box) {
    b[0].x = box[0].x, b[0].y = box[0].y, b[0].z = box[0].z;
    b[1].x = box[1].x, b[1].y = box[1].y, b[1].z = box[1].z;
  }
  else {
    b[0].x = dom[0].x, b[0].y = dom[0].y, b[0].z = dom[0].z;
    b[1].x = dom[1].x, b[1].y = dom[1].y, b[1].z = dom[1].z;
  }

  int axis = 2; // 0: x-normal (y,z), 1: y-normal (x,z), 2: z-normal (x,y)
#if dimension > 2
  double ax = fabs(n.x), ay = fabs(n.y), az = fabs(n.z);
  if (ax >= ay && ax >= az)
    axis = 0;
  else if (ay >= ax && ay >= az)
    axis = 1;
  else
    axis = 2;

  double eps = 1e-9*L0;
  if (axis == 0) {
    double pos = (n.x != 0.) ? alpha/n.x : alpha;
    if (pos < dom[0].x + eps) pos = dom[0].x + eps;
    if (pos > dom[1].x - eps) pos = dom[1].x - eps;
    b[0].x = b[1].x = pos;
  }
  else if (axis == 1) {
    double pos = (n.y != 0.) ? alpha/n.y : alpha;
    if (pos < dom[0].y + eps) pos = dom[0].y + eps;
    if (pos > dom[1].y - eps) pos = dom[1].y - eps;
    b[0].y = b[1].y = pos;
  }
  else {
    double pos = (n.z != 0.) ? alpha/n.z : alpha;
    if (pos < dom[0].z + eps) pos = dom[0].z + eps;
    if (pos > dom[1].z - eps) pos = dom[1].z - eps;
    b[0].z = b[1].z = pos;
  }
#endif

  int base_cells = N/Dimensions.x;
  int default_nu = (axis == 0) ? base_cells*Dimensions.y
			       : base_cells*Dimensions.x;
  int nu = (res > 0) ? res : default_nu;
  if (nu < 2) nu = 2;
  if (nu % 2) nu++;

  double u0 = (axis == 0) ? b[0].y : b[0].x;
  double u1 = (axis == 0) ? b[1].y : b[1].x;
  double v0 = (axis == 2) ? b[0].y : b[0].z;
  double v1 = (axis == 2) ? b[1].y : b[1].z;

  double delta = (u1 - u0)/nu;
  int nv = (int)((v1 - v0)/delta + 0.5);
  if (nv < 2) nv = 2;
  if (nv % 2) nv++;

  coord cn = {1, 1, 1}, p;
  if (axis == 0)
    cn.x = 1, cn.y = nu, cn.z = nv;
  else if (axis == 1)
    cn.x = nu, cn.y = 1, cn.z = nv;
  else
    cn.x = nu, cn.y = nv, cn.z = 1;

  Color ** ppm = (Color **) matrix_new (nv, nu, sizeof(Color));
  unsigned char * ppm0 = &ppm[0][0].r;
  int len = 3*nu*nv;
  memset (ppm0, 0, len*sizeof (unsigned char));
  double cmap[NCMAP][3];
  (* map) (cmap);

  if (linear) {
    if (mask.i >= 0)
      boundary ({f, mask});
    else
      boundary ({f});
  }

#if _MPI
  foreach_region (p, b, cn, reduction(max:ppm0[:len]))
#else
  foreach_region (p, b, cn, cpu)
#endif
  {
    double v;
    if (mask.i >= 0) {
      if (linear) {
	double m = interpolate_linear (point, mask, p.x, p.y, p.z);
	if (m < 0.)
	  v = nodata;
	else
	  v = interpolate_linear (point, f, p.x, p.y, p.z);
      }
      else {
	if (mask[] < 0.)
	  v = nodata;
	else
	  v = f[];
      }
    }
    else if (linear)
      v = interpolate_linear (point, f, p.x, p.y, p.z);
    else
      v = f[];

    double pu = (axis == 0) ? p.y : p.x;
    double pv = (axis == 2) ? p.y : p.z;
    int i = (int)((pu - u0)/(u1 - u0)*nu);
    int j = (int)((pv - v0)/(v1 - v0)*nv);
    if (i < 0) i = 0;
    else if (i >= nu) i = nu - 1;
    if (j < 0) j = 0;
    else if (j >= nv) j = nv - 1;
    Color ** alias = ppm; // so that qcc considers ppm a local variable
    alias[nv - 1 - j][i] = colormap_color (cmap, v, min, max);
  }

  if (pid() == 0) {
    if (file)
      fp = open_image (file, opt);

    fprintf (fp, "P6\n%d %d 255\n", nu, nv);
    fwrite (ppm0, sizeof(unsigned char), len, fp);

    if (file)
      close_image (file, fp);
    else
      fflush (fp);

    if (checksum) {
      Adler32Hash hash;
      a32_hash_init (&hash);
      a32_hash_add (&hash, ppm0, sizeof(unsigned char)*len);
      fputs ("# ", checksum);
      if (file)
	fprintf (checksum, "%s: ", file);
      fprintf (checksum, "checksum: %08lx\n", (unsigned long) a32_hash (&hash));
    }
  }

  matrix_free (ppm);
}

#endif // OUTPUT_PPM_SLICE_H
