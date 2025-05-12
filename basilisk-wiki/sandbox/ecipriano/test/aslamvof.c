/**
# Aslam Extrapolations

This test case is similar to [aslam.c](aslam.c). However,
instead of starting from the level set field, we start from
a vof field, which is converted into level set, and then
Aslam extrapolations are performed.
*/

#include "grid/multigrid.h"
#include "utils.h"
#include "aslam.h"
#include "redistance.h"
#include "view.h"

/**
We define a function that writes a picture with the
map and isolines of the extended scalar fields.
*/

void write_picture (char* name, scalar u) {
  vertex scalar phi[];
  foreach_vertex()
    phi[] = (u[] + u[-1] + u[0,-1] + u[-1,-1])/4.;
  clear();
  isoline ("levelset", val = 0., lw = 2.);
  isoline ("phi", n = 30);
  squares ("phi", spread = -1);
  box();
  save (name);
}

void write_levelset (void) {
  clear();
  isoline ("levelset", val = 0., lw = 2.);
  squares ("levelset", spread = -1);
  box();
  save ("levelset.png");
}

#define ufunc(x,y)(x*y)
#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

/**
We declare the level set field *levelset*, and
the field to extrapolate *u*. */

scalar levelset[], u[], f[];

int main (void) {

  /**
  We set the domain geometry and we initialize
  the grid. */

  size (2.*pi);
  origin (-pi,-pi);
  init_grid (1 << 8);
  double R0 = 2.;

  /**
  We initialize the volume fraction field. */

  fraction (f, (sq(R0) - sq(x) - sq(y)));

  /**
  We reconstruct the levelset field. */

  vof_to_ls (f, levelset, imax=300);
  write_levelset();

  /**
  We initialize the function *u* to be extrapolated
  and we call the *constant_extrapolation()* function.
  The extrapolations are perfomed using a $\Delta t$
  of 0.01 (it must be small enough to guarantee the
  stability of the explicit in time discretization),
  and using a total number of time steps equal to 300,
  in order to obtain a steady-state solution. */

  foreach()
    u[] = ufunc(x,y)*f[];
  write_picture ("initial.png", u);

  constant_extrapolation (u, levelset, 0.5, 300, c=f);
  write_picture ("constant.png", u);
  fprintf (stderr, "constant = %g\n", statsf(u).sum);

  /**
  We re-initialize the function *u* and we apply the
  *linear_extrapolation()*. */

  foreach()
    u[] = ufunc(x,y)*f[];
  linear_extrapolation (u, levelset, 0.5, 300, c=f);
  write_picture ("linear.png", u);
  fprintf (stderr, "linear = %g\n", statsf(u).sum);
}

/**
## Results

![Level Set](aslamvof/levelset.png)

![Initial field](aslamvof/initial.png)

![Constant extrapolation](aslamvof/constant.png)

![Linear extrapolation](aslamvof/linear.png)

## References

~~~bib
@article{aslam2004partial,
  title={A partial differential equation approach to multidimensional extrapolation},
  author={Aslam, Tariq D},
  journal={Journal of Computational Physics},
  volume={193},
  number={1},
  pages={349--355},
  year={2004},
  publisher={Elsevier}
}
~~~
*/

