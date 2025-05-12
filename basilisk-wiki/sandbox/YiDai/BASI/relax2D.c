/**
# A novice trial on relaxation method for solving elliptic PDEs.

Take 2D poisson equation as an example

$$\frac{\partial^{2}U}{\partial{x^2}} + \frac{\partial^{2}U}{\partial{y^2}} = f$$

If we set out known solution as $sin(\pi x) + cos(\pi y)$, then $f = -(\pi)^2 sin(\pi x) --(\pi)^2 cos(\pi x)$.
We use [relaxation method](https://en.wikipedia.org/wiki/Relaxation_(iterative_method)#Model_problem_of_potential_theory) to solve this 1D poisson equation.

by initially assume f = 1, we can iterate until the solution converge.

*/

#include "utils.h"
#include "view.h"

double solution (double x, double y) { return (sin(pi*x) + cos(pi*y)); }
scalar U[], f[];

// set BC from solution
U[right] = dirichlet(solution(x, y));
U[left] = dirichlet(solution(x, y));
U[top] = dirichlet(solution(x, y));
U[bottom] = dirichlet(solution(x, y));

int main()
{
    L0 = 2*pi;
    X0 = -L0 / 2.;
    Y0 = -L0 / 2.;
    init_grid(1<<8);
    foreach(){
        U[]=0;
        f[]= -sq(pi)*sin(pi*x) -sq(pi)*cos(pi*y);
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
        foreach ()
        {
            
            U[] = (U[1, 0] + U[0, 1] + U[-1, 0] + U[0, -1] - sq(Delta) * f[]) / 4;
            e += fabs(U[] - solution(x, y));
        }
        boundary({U});
        if (e < sumerror)
            sumerror = e;
        printf("%d %g %g\n", j, sumerror, e);
        j++;
    }
    FILE *fp1 = fopen("outputfield_U", "w");
    output_field({U}, fp1);
}



/**
In the end, it took **22972 times** iteration to converge, reaching field sum error of 2. 

~~~gnuplot
reset
set terminal png size 1200,600 enhanced font 'Times-Roman,16'
set key samplen 2 spacing 1.5 font 'Times-Roman,16'

set multiplot layout 1,1
dataU = "outputfield_U"
set size ratio -1
set pm3d map
set palette defined (-3 "blue", 0 "white", 1 "red")
splot[][] dataU u 1:2:3 t "num results"
unset multiplot
~~~
*/