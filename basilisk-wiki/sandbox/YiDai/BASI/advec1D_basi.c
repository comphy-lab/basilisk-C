/**
# 1D advection equation

$$\partial_{t}T = c \partial_{x}T$$

use basilisk advection solver to numerically solve this equation
*/
#include "grid/cartesian1D.h"
#include "advection.h"

scalar T[];
scalar T2[];
scalar *tracers = {T, T2};

double alpha = 2;

u.n[left] = 0.;
u.n[right] = 0.;

int main()
{
    L0 = 4 * pi;
    X0 = -L0 / 2.;
    N = 1 << 12;

    DT = 1e-4;
    CFL = 0.8;
    run();
}

T[left] = dirichlet(0.);
T[right] = dirichlet(0.);

event init(t = 0)
{
    foreach ()
    {
        T[] = 0.5 * (fabs(x) < 1);
        T2[] = exp(-sq(x)/4);
    }
    boundary({T, T2});
}

event velocity(i++)
{
    trash({u});
    foreach_face()
    {
        u.x[] = alpha;
    }
}

// event printdata (t = 0; t <= 1000 * DT; t += 100 * DT) {
event printdata(t = 0; t <= 1; t += 0.1)
{
    // event printdata (t = 0.5) {
    static FILE *fp = fopen("AD1D", "w");
    foreach ()
        fprintf(fp, "%g %g %g %g\n", x, T[], T2[], t);
    fprintf(fp, "\n\n");
    fflush(fp);
}

/**
The results after some time steps become unstable a bit

~~~gnuplot
reset
file="AD1D"
set terminal png size 1200,600 enhanced font 'Times-Roman,16'
set key samplen 2 spacing 1.5 font 'Times-Roman,16'

set key bottom right
set multiplot layout 1,2
plot[][-0.05:1.05] file u ($1):($2) t "T1" w l
plot[][-0.05:1.05] file u ($1):($3) t "T2" w l
unset multiplot
~~~

*/