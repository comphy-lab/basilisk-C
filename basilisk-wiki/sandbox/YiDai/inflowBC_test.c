/**
# test of inflow BC
Here we test the function [inflowBC.h](http://www.basilisk.fr/sandbox/YiDai/inflowBC.h).
*/
#include "grid/octree.h"
#include "inflowBC.h"

int main()
{
        L0 = 4.;
        X0 = Y0 = Z0 = 0;
        init_grid(4);
        scalar v[];
        scalar vr[];
        foreach ()
        {
                v[] = x + y;
        }
        v[left] = y;
        boundary({v});
        char nameSliceB[91];
        snprintf(nameSliceB, 90, "Vfield_f");
        write_inflow(v, nameSliceB, 2);

        read_inflow(vr, nameSliceB, 2);
        foreach_boundary(left)
        {
                printf("%g %g %g\n", y, z, vr[-1, 0]);
        }
}
/**
The output is as expected:
*/
/**
~~~bash
0.5 0.5 0.5
2.5 0.5 2.5
2.5 1.5 2.5
2.5 2.5 2.5
2.5 3.5 2.5
3.5 0.5 3.5
3.5 1.5 3.5
3.5 2.5 3.5
3.5 3.5 3.5
0.5 1.5 0.5
1.5 0.5 1.5
1.5 1.5 1.5
0.5 2.5 0.5
1.5 2.5 1.5
0.5 3.5 0.5
1.5 3.5 1.5
~~~
*/

