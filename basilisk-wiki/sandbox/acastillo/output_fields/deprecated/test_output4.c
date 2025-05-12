#include "grid/octree.h"
#include "output_xmf.h"

#define r2 (sq(x) + sq(y))
int main()
{
  L0 = 1;
  X0 = Y0 = Z0 = -L0 / 2;
  N = 1 << 6;
  init_grid(N);

#if TREE
  refine(((r2 < sq(0.25)) && (r2 > sq(0.1))) && level < 8);
#endif

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

  scalar *slist = (scalar *){s};
  vector *vlist = (vector *){u};

  create_output_h5(H5FILE_NAME);
  write_xml_h5_head("test");
  output_xmf(slist, vlist, "test", "3D000");
  write_xml_h5_tail("test");

#if dimension > 2
  create_output_h5(H5FILE_NAME2);
  write_xml_h5_head("test_x");
  output_slice_xmf(slist, vlist, "test_x", "SX000", (coord){1, 0, 0}, 0);
  write_xml_h5_tail("test_x");

  write_xml_h5_head("test_y");
  output_slice_xmf(slist, vlist, "test_y", "SY000", (coord){0, 1, 0}, 0);
  write_xml_h5_tail("test_y");

  write_xml_h5_head("test_z");
  output_slice_xmf(slist, vlist, "test_z", "SZ000", (coord){0, 0, 1}, 0);
  write_xml_h5_tail("test_z");
#endif
}
