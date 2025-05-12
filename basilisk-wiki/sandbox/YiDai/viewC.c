/**
The view in view2D mixed with view in view3D. clear() doesn't make it work.

![t = 0](viewC/view2D0.png)

![t = 1](viewC/view2D1.png)

*/


#include "grid/octree.h"
#include "view.h"
#include "run.h"

int main()
{
    init_grid(1);
    origin(-0.5, -0.5, -0.5);
    run();
}

event view2D(t += 1; t < 2)
{
    view();
    box();
    char name2[90];
    snprintf(name2, 90, "view2D%g.png", t);
    save(name2);
    clear();
}

event view3D(t += 1; t < 2)
{
    view(phi = 0.4);
    box();
    char name3[90];
    snprintf(name3, 90, "view3D%g.png", t);
    save(name3);
    clear();
}
