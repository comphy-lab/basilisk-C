/**
# 1D advection equation

$$\partial_{t}T + c \partial_{x}T=0$$

Lax - Friedrichs method is used to avoid the instablity of Forward-Time Central-Space (FTCS) Method.
*/
#include "grid/cartesian1D.h"
#include "run.h"

scalar T[], dT[];
double dt;
double alpha = 1;

int main()
{
    // periodic(left);
    L0 = 4 * pi;
    X0 = -L0 / 2.;
    N = 1 << 9;
    // stable condition of Lax - Friedrichs method dt < dx/c
    DT = L0 / N / alpha;
    run();
}

T[left] = dirichlet(0.);
T[right] = dirichlet(0.);

event init(t = 0)
{
    foreach ()
    {
        T[] = 0.5 * (fabs(x) < 1);
    }
    boundary({T});
}

// event printdata (t = 0; t <= 1000 * DT; t += 100 * DT) {
event printdata(t = 0; t <= 1; t += 0.1)
{
    // event printdata (t = 0.5) {
    static FILE *fp = fopen("AD1D", "w");
    foreach ()
        fprintf(fp, "%g %g %g %g\n", x, T[], dT[], t);
    fprintf(fp, "\n\n");
    fflush(fp);
}

event integration(i++)
{
    double dt = DT;
    dt = dtnext(dt);
    foreach ()
        // Lax - Friedrich method
        dT[] = -alpha * (T[1] - T[-1]) / (2 * Delta);
    foreach ()
        T[] = 0.5 * (T[1] + T[-1]) + dt * dT[];
    boundary({T});
}

/**
The results after some time steps become unstable a bit

~~~gnuplot
reset
file="AD1D"
set terminal png size 1200,600 enhanced font 'Times-Roman,16'
set key samplen 2 spacing 1.5 font 'Times-Roman,16'

set key bottom right
plot[][-0.05:0.7]  file u ($1):($2) t "T" w l
~~~

*/