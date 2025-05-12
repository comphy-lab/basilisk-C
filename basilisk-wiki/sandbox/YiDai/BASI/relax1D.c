/**
# A novice trial on relaxation method for solving elliptic PDEs.

Take 1D poisson equation as an example

$$\frac{\partial^{2}U}{\partial{x^2}} = f$$

If we set out known solution as $sin(3\pi x)$, then $f = -(3\pi)^2 sin(3\pi x)$.
We use [relaxation method](https://en.wikipedia.org/wiki/Relaxation_(iterative_method)#Model_problem_of_potential_theory) to solve this 1D poisson equation.

by initially assume f = 1, we can iterate until the solution converge.

*/


#include "utils.h"

double solution (double x) { return sin(3.*pi*x); }
scalar U[], f[];
int nrelax=100;  // number of iteration for relaxation

// set BC from solution
U[right] = dirichlet (solution(x));
U[left] = dirichlet (solution(x));

int main(){
    L0 = 2*pi;
    X0 = 0;
    init_grid(1<<8);
    foreach(){
        U[]=0;
        f[]= -sq(3*pi)*sin(3*pi*x);
    }

    for (int j = 1; j <nrelax; j+=1){
        foreach(){
            U[] = (- sq(Delta) * f[] + U[1] + U[-1])/2;
            printf ("%g %d\n", U[], j);
        }
        printf("\n\n");
    }
    double max = 0;
    foreach() {
        double e = U[] - solution(x);
        if (fabs(e) > max) max = fabs(e);
            fprintf (stderr, "%g %g %g\n", x, U[], e);
    }
}

/**
This is actually similar to generic basilisk [relexation function](http://basilisk.fr/src/poisson.h#relax), which use weighted Jacobi relaxation and doesn't store a new field (simple reuse)

This is applied in [here](http://basilisk.fr/sandbox/YiDai/BASI/relax_basi.c)

~~~gnuplot
reset
plot[-1:7] 'log' using 1:2 with points t "solver", sin(3.*pi*x) with lines t "analytical"
~~~
*/

/**
~~~gnuplot
reset
plot[-1:7][-0.05:0.05] 'log' using 1:3 with points t "error"
~~~
*/