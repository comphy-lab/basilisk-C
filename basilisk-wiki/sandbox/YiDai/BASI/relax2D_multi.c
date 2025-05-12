/**
# Here we use basilisk multigrid [possion solver](http://basilisk.fr/src/poisson.h) to solve simple 2D steady heat equation

*/

#include "utils.h"
#include "poisson.h"

#define dimension 2

scalar U[], f[];

double solution(double x, double y) 
{ 
    return (sin(pi * x) + cos(pi * y)); 
}

U[right] = dirichlet(solution(x, y));
U[left] = dirichlet(solution(x, y));
U[top] = dirichlet(solution(x, y));
U[bottom] = dirichlet(solution(x, y));

int main ()
{
    L0 = 2 * pi;
    X0 = -L0 / 2.;
    Y0 = -L0 / 2.;
    int depth = 7;
    init_grid(1 << depth);
    foreach ()
    {
        U[] = 0;
        f[] = -sq(pi) * sin(pi * x) - sq(pi) * cos(pi * y);
    }
    TOLERANCE = 1e-4;
    poisson(U, f);
    double max = 0;
    foreach ()
    {
        double e = U[] - solution(x, y);
        if (fabs(e) > max)
            max = fabs(e);
        printf("%g %g %g %g\n", x, y, U[], e);
    }
    FILE *fp1 = fopen("outputfield_U", "w");
    output_field({U}, fp1);
}
/**

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
~~~

*/


