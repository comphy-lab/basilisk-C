/**
# fiber simulation without flow interaction
Here is a test file to check if a filament will swing under gravity. This is a collaboration work with a master student Julian Sauerbier. 

*/

#include "run.h"

#include "fiber_MG.h"
#include "view.h"
#include "fiber_watch.h"

double TEND = 3.;

int main()
{
    L0 = 6.;
    N = 256;
    X0 = Y0 = Z0 = 0;
    DT = 1E-4;
    run();
}

event movie(t += 0.1)
{
    view(tx = -0.5, ty = -0.5);
    fiberW_node(&FB(0), pc = {0, 128. / 256., 0});
    box();
    save("fiber.mp4");
}

event end(t = TEND)
{
}

/**
![flexible fiber](fiber_gravity/fiber.mp4)
*/