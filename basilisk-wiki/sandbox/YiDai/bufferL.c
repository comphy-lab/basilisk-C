/**
# Buffer layer
To avoid periodic boundary conditions, here we present a toy example of a gradient damping layer, instead of constant damping. Constant damping can lead to abnormal grid resolving features...
This is a very dull treatment. Is there a smarter way to handle this?

![buffer layer](bufferL/ff.png)
*/

// #include "grid/octree.h"
#include "view.h"

int main()
{
    init_grid(32);
    origin(0, 0, 0);
    scalar f[];
    double bf_L = 0.2;

#if dimension == 2
    foreach ()
    {
        if (y <= -x + L0 && x <= bf_L)
        {
            f[] = bf_L - x;
        }
        if (y <= x && x >= L0 - bf_L)
        {
            f[] = x - (L0 - bf_L);
        }
        if (y >= -x + L0 && y >= x && y >= L0 - bf_L)
        {
            f[] = y - (L0 - bf_L);
        }
    }
    view(fov = 30, width = 600, height = 600, tx = -0.5, ty = -0.5);
    squares("f", min = 0, max = bf_L);
#endif
#if dimension == 3
    foreach ()
    {
        if (y <= -x + L0 && x <= bf_L && z >= x && z <= L0 - x)
        {
            f[] = bf_L - x;
        }
        if (y <= x && x >= L0 - bf_L && z <= x && z >= L0 - x)
        {
            f[] = x - (L0 - bf_L);
        }
        if (y <= -z + L0 && z <= bf_L && x >= z && x <= L0 - z)
        {
            f[] = bf_L - z;
        }
        if (y <= z && z >= L0 - bf_L && x <= z && x >= L0 - z)
        {
            f[] = z - (L0 - bf_L);
        }
        if (y >= -x + L0 && y >= x && y >= -z + L0 && y >= z && y >= L0 - bf_L)
        {
            f[] = y - (L0 - bf_L);
        }
    }
    view(fov = 30, width = 800, height = 800, theta = 0, phi = 1.5, tx = -0.5, ty = 0.5, bg = {0.3, 0.4, 0.6});
    squares("f", n = {0, 1, 0}, alpha = L0 / 2, min = 0, max = bf_L);
#endif

    box();
    save("ff.png");

    FILE *fp1 = fopen("ffile", "w");
    output_field({f}, fp1);
}