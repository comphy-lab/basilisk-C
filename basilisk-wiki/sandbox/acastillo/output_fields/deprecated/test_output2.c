#include "grid/octree.h"
#include "output_vtu_foreach.h"

int main()
{
  L0 = 1;
  X0 = Y0 = Z0 = -L0 / 2;
  N = 1 << 3;
  periodic(left);
  periodic(top);

  init_grid(N);

  scalar s[];
  vector u[];
  foreach ()
  {
    s[] = pid();
    foreach_dimension()
    {
      u.x[] = sq(x) + sq(y) + sq(z);
    }
  }
  boundary({s, u});

  output_vtu((scalar *){s}, (vector *){u}, "test");
}
