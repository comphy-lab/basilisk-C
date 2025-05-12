/**
# Relaxation function in basilisk

The [relaxation function](http://basilisk.fr/src/poisson.h#relax) in poisson solver is used here to solver 2D poisson equation same as [here](http://basilisk.fr/sandbox/YiDai/BASI/relax2D.c)*/

// #include "grid/multigrid.h"
#include "utils.h"
#include "poisson.h"

double solution(double x, double y) { return (sin(pi * x) + cos(pi * y)); }
scalar U[], f[], ee[];
// int nrelax = 10000;
int depth = 8;
// void * data;
struct Poisson pp;
const scalar lambda[] = 0;

U[right] = dirichlet(solution(x, y));
U[left] = dirichlet(solution(x, y));
U[top] = dirichlet(solution(x, y));
U[bottom] = dirichlet(solution(x, y));

int main()
{
    pp.a = U;
    pp.b = f;
    pp.alpha = unityf;
    pp.lambda = lambda;
    L0 = 2 * pi;
    X0 = -L0 / 2.;
    Y0 = -L0 / 2.;

    init_grid(1 << depth);
    foreach ()
    {
        U[] = 0;
        f[] = -sq(pi) * sin(pi * x) - sq(pi) * cos(pi * y);
    }
    double Tolerance = 2;
    double sumerror = 0.;
    foreach ()
    {
        sumerror += fabs(solution(x, y));
    }
    printf("%g\n", sumerror);
    int j = 0;
    while (sumerror > Tolerance)
    {
        double e = 0;
        relax({U}, {f}, depth, &pp);
        foreach ()
        {
            e += fabs(U[] - solution(x, y));
        }
        if (e < sumerror)
            sumerror = e;
        printf("%d %g %g\n", j, sumerror, e);
        j++;
    }

    foreach ()
    {
        ee[] = U[] - solution(x, y);
    }
    FILE *fp1 = fopen("outputfield_U", "w");
    output_field({U}, fp1);
    FILE *fp2 = fopen("outputfield_e", "w");
    output_field({ee}, fp2);
}

/**
~~~gnuplot
reset
plot[-1:7] 'log' using 1:2 with points t "solver", sin(3.*pi*x) with lines t "analytical"
~~~

*/

/**
reset
~~~gnuplot
reset
plot[-1:7][-0.35:0.35] 'log' using 1:3 with l t "error"
~~~
*/

/**
Using brutal force iteration does not work in this case
*/