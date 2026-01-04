// dump() is unavailable in cartesian
#include "grid/multigrid.h"
#include "run.h"

void scalar_call()
{
    scalar f[] = 0;

    foreach()
    {
        fprintf(stdout, "%f \t %f \t %f \n", x, y, f[]);
        f[]++;
    }

    delete({f});
}

int main()
{
    X0 = 0.;
    Y0 = 0.;
    L0 = 10;
    N = 10;
    run();
}

event call1 (t = 0)
{
    scalar_call();
    /**
    The code `foreach(){f[] = 0;}` is invalid here because `f` is out of scope, rendering its data inaccessible. However, the data may still reside in memory.
     */
    scalar g[];
    fprintf(stdout, "Call outside the function");
    foreach()
    {
        fprintf(stdout, "%f \t %f \t %f \n", x, y, g[]);
    }
    /**
    The newly defined scalar `g` retains old data. It should be initialized each time to ensure consistent behavior.
     */
}