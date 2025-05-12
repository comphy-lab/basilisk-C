/**
# 2D advection equation

$$\partial_{t}T = c \partial_{x}T + c \partial_{y}T$$

*/
#include "grid/cartesian.h"
#include "advection.h"
// #include "view.h"

scalar T[];
scalar *tracers = {T};
double alpha = 1;

u.n[left] = 0.;
u.n[right] = 0.;
u.n[top] = 0.;
u.n[bottom] = 0.;

int main()
{
    L0 = 4 * pi;
    X0 = -L0 / 2.;
    Y0 = -L0 / 2.;
    
    DT = 1e-4;
    CFL = 0.8;
    for (N = 64; N <= 256; N *= 2)
        run();
}

T[left] = dirichlet(0.);
T[right] = dirichlet(0.);
T[top] = dirichlet(0.);
T[bottom] = dirichlet(0.);

event init(t = 0)
{
    foreach ()
    {
        T[] = 1 * (sqrt(sq(x) + sq(y)) < 1);
    }
}

event velocity(i++){
    trash({u});
    foreach_face(){
        u.x[] = alpha;
        u.y[] = alpha;
    }
}

event movie(t = 0; t <= 0.4; t += 0.01)
{
    if (N == 256)
        output_ppm(T, file = "T1_ap.mp4", min = 0, max = 1, linear = true, map = cool_warm);
}

/**
The basilisk solver shows some error near the sharp edges

![temp video](advec2D_basi/T1_ap.mp4)

*/