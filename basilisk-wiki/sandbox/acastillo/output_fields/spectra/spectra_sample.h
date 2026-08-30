/**
# Sampling a plane onto a regular lattice

`sample_scalar_plane()` evaluates every field in `list` on a regular m1 x m2
lattice in the plane $z = h$ and assembles the result, across ranks, into one
contiguous array. 3D only.

This is the piece a horizontal FFT needs: `foreach_region()` returns the cell
containing each query point, so on a tree grid the lattice is a uniform view of
an adaptive field. With `xmin = X0`, `xmax = X0 + L0` and `m1` equal to the
finest grid size, the query points are the finest cells' centres and they cover
$[X_0, X_0+L_0)$ exactly once -- one period, no repeated endpoint.

Follows `output_field()` in [output.h](/src/output.h): the array starts at
`nodata` and is combined with a *min* reduction, so a point no rank owns stays
`nodata` and one located on several ranks resolves to the same value either
way. Values are read as `s[]` rather than interpolated, matching
[profiles_scalar.h](../profiles/profiles_scalar.h); away from the finest level that means
the containing cell's value, so a lattice finer than the local grid low-passes
the field.

`h` on a cell face is ambiguous for a $z$-dependent field -- `locate()` picks
one side. Offset by half a cell to sample cell centres.

The caller allocates `plane` with `m1*m2*len` doubles, laid out as
`plane[(i*m2 + j)*len + k]` for field `k` at lattice point `(i,j)`. Returns the
number of unfilled points, 0 for a complete plane.
*/

#if dimension != 3
  #error sample_scalar_plane() is 3D only
#endif

int sample_scalar_plane (scalar * list, double * plane, double h,
                         double xmin, double xmax,
                         double ymin, double ymax,
                         int m1, int m2)
{
  int len = list_len (list), n = m1*m2*len;
  coord box[2] = {{xmin, ymin, h}, {xmax, ymax, h}};
  coord nsamples = {m1, m2, 1};
  for (int i = 0; i < n; i++)
    plane[i] = nodata;

  coord p;
  foreach_region (p, box, nsamples, reduction(min:plane[:n])){
    double * alias = plane; // so that qcc considers 'plane' a local variable
    int i = (p.x - box[0].x)/(box[1].x - box[0].x)*nsamples.x;
    int j = (p.y - box[0].y)/(box[1].y - box[0].y)*nsamples.y;
    int k = 0;
    for (scalar s in list)
      alias[(i*m2 + j)*len + k++] = s[];
  }

  int holes = 0;
  for (int i = 0; i < n; i++)
    if (plane[i] == nodata)
      holes++;
  return holes;
}

/**
Snap a target height to the nearest cell centre. `foreach_region` returns the
cell containing the point, so a height on a face resolves to one side
arbitrarily -- $z = 0$ is a face whenever the grid size is even. */

double snap_to_cell (double h, int m)
{
  double del = L0/m;
  return Z0 + (floor ((h - Z0)/del) + 0.5)*del;
}
