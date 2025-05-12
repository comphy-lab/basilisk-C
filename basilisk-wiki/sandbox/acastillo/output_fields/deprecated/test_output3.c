#include "grid/octree.h"
#include "output_vtu_foreach.h"

int main()
{
	L0 = 1;
	X0 = Y0 = Z0 = -L0 / 2;
	N = 1 << 6;
	periodic(left);
	periodic(top);
	periodic(front);
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

	scalar *slist = (scalar *){s};
	vector *vlist = (vector *){u};

	output_vtu(slist, vlist, "test");

	output_slice_vtu(slist, vlist, "test_x", (coord){1, 0, 0}, 0.0);
	output_slice_vtu(slist, vlist, "test_y", (coord){0, 1, 0}, 0.0);
	output_slice_vtu(slist, vlist, "test_z", (coord){0, 0, 1}, 0.0);
}
