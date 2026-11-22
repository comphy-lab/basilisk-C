/**
# `apply_bc()` zeroes out `right` / `top` / `front` boundary values of `vertex scalar` fields on GPU

In `src/grid/gpu/grid.h`, `boundary_top()`, `boundary_bottom()`,
`boundary_front()`, `boundary_back()`, and the `left` / `right` centered blocks
in `apply_bc()` do not guard on `b.boundary_* != NULL`:

~~~c
if (point.i == N*Dimensions.x + GHOSTS - 1) { // right
  for (scalar s in apply_bc_list)
    if (!s.face || s.i != s.v.x.i)
      foreach_blockf(s)
        s[-bc_period_x,0,0] = s.boundary_right (point, neighborp(-bc_period_x,0,0), s, &data);
~~~

For a `vertex scalar` (such as `phi[]` used by `fractions (phi, f)`),
`phi.boundary_*` is `NULL`, which compiles to `(0)` in the GPU shader.
Consequently, `apply_bc()` silently overwrites the `right` (`i = N + GHOSTS`),
`top` (`j = N + GHOSTS`), and 3D `front` (`k = N + GHOSTS`) boundary vertices
with `0.0` whenever `foreach_vertex()` finishes.

This corrupts `fractions (phi, f)` on any interface touching the `right`,
`top`, or `front` domain boundary.

Compare:

~~~bash
# CPU multigrid (passes: min(phi) = 1)
qcc -O2 -grid=multigrid gpu_vertex_scalar_bc_zero.c -lm && ./a.out

# GPU multigrid (fails: min(phi) = 0)
qcc -autolink -O2 gpu_vertex_scalar_bc_zero.c -lm && ./a.out
~~~ */

#include "grid/gpu/multigrid.h"

int main()
{
  init_grid (16);
  vertex scalar phi[];
  foreach_vertex()
    phi[] = 1.;

  double pmin = HUGE;
  foreach_vertex (reduction(min:pmin))
    if (phi[] < pmin)
      pmin = phi[];

  fprintf (stderr, "min(phi) = %.6g (expected 1)\n", pmin);
  assert (fabs (pmin - 1.) < 1e-6);
}
