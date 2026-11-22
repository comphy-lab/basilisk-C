/**
# Periodic `face vector` normal components are not wrapped on GPU

On CPU `grid/multigrid.h`, `periodic (right)` enforces periodicity of the
normal component of a `face vector` (`uf.x[N] == uf.x[0]`) in
`box_boundary_level_normal()`, even though `init_face_vector()` leaves
`uf.x.boundary[left] = uf.x.boundary[right] = NULL`.

On GPU `grid/gpu/grid.h`, `apply_bc()` only updates normal face values when
`s.boundary_right` / `s.boundary_top` / `s.boundary_front` is non-NULL:

~~~c
if (point.i == N*Dimensions.x + GHOSTS)
  for (scalar s in apply_bc_list)
    if (s.face && s.i == s.v.x.i && s.boundary_right)
      foreach_blockf(s)
        s[] = s.boundary_right (neighborp(bc_period_x,0,0), point, s, &data);
~~~

Two issues occur:

1. Because `uf.x.boundary_right == NULL`, `apply_bc()` never copies `uf.x[0]`
   to `uf.x[N]` (nor `uf.y[N]` on `periodic(top)`, nor `uf.z[N]` on 3D
   `periodic(front)`).
2. Even if `s.boundary_right` were set to `periodic_bc`,
   `neighborp(bc_period_x, 0, 0)` uses `+bc_period_x` (`+N*Dimensions.x` when
   periodic) instead of `-bc_period_x`, indexing out of bounds instead of
   wrapping back to `GHOSTS`.

Compare:

~~~bash
# CPU multigrid (passes: left = 1, right = 1, diff = 0)
qcc -O2 -grid=multigrid gpu_periodic_face_bc.c -lm && ./a.out

# GPU multigrid (fails: left = 1, right = 2, diff = 1)
qcc -autolink -O2 gpu_periodic_face_bc.c -lm && ./a.out
~~~ */

#include "grid/gpu/multigrid.h"

int main()
{
  periodic (right);
  init_grid (16);
  face vector uf[];
  foreach_face (x)
    uf.x[] = 1. + x/L0;

  double right_val = 0., left_val = 0.;
  foreach (reduction(+:left_val) reduction(+:right_val)) {
    if (x < Delta)
      left_val += uf.x[];
    if (x > L0 - Delta)
      right_val += uf.x[1];
  }
  left_val /= N;
  right_val /= N;
  double err = fabs (right_val - left_val);
  fprintf (stderr, "left = %.6g, right = %.6g, diff = %.6g (expected 0)\n",
           left_val, right_val, err);
  assert (err < 1e-6);
}
