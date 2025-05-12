/**
# Relaxation function in basilisk

The [relaxation function](http://basilisk.fr/src/poisson.h#relax) in poisson solver is used here to solver 1D poisson equation same as [here](http://basilisk.fr/sandbox/YiDai/BASI/relax1D.c)*/

#include "grid/multigrid1D.h"
#include "utils.h"
#include "poisson.h"


double solution (double x) { return sin(3.*pi*x); }
scalar U[], f[];
int nrelax=100;
int depth = 8;
// void * data;
struct Poisson pp;
const scalar lambda[] = 0;

U[right] = dirichlet (solution(x));
U[left] = dirichlet (solution(x));

int main(){
    pp.a = U; pp.b = f; pp.alpha = unityf; pp.lambda = lambda;
    L0 = 2*pi;
    X0 = 0;

    init_grid(1<<depth);
    foreach(){
        U[]=0;
        f[]= -sq(3*pi)*sin(3*pi*x);
    }
    for (int j = 1; j <nrelax; j+=1){
         relax ({U}, {f}, depth, &pp);
    }
    double max = 0;
    foreach() {
        double e = U[] - solution(x);
        if (fabs(e) > max) max = fabs(e);
            fprintf (stderr, "%g %g %g\n", x, U[], e);
    }
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
Oh what's wrong with basilisk solver, do i need to update the boundary condition?
*/

