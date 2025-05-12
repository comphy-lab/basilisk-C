/**
Being inspired by [solve.h](http://www.basilisk.fr/sandbox/Antoonvh/kuramato.c), here we test if multigrid method can solve a bit non-linearity, i.e., preserve the shock-wave features. We compare the Godunov method with [solve.h](http://www.basilisk.fr/src/solve.h). The case of visid-burgers equation

In invisid case, the Godunov method seems to outperform multigrid (which is not designed for this case). Here we add a bit diffusivity. 

*/

#include "grid/multigrid1D.h"
#include "solve.h"
#include "run.h"

double alpha = 0.01;
scalar u[];
int sim;
void solve_GD(scalar T, double dt)
{
    foreach ()
    {
        T[] = T[] - dt * (((T[] >= T[1]) ? \ 
        (((T[] + T[1]) / 2 > 0) ? (sq(T[]) / 2) : (sq(T[1]) / 2)): \ 
        ((T[] > 0) ? (sq(T[]) / 2) : ((T[1] < 0) ? (sq(T[1]) / 2) : 0))) \ 
        - ((T[-1] >= T[]) ? (((T[-1] + T[]) / 2 > 0) ? (sq(T[-1]) / 2) : (sq(T[]) / 2)) : ((T[-1] > 0) ? (sq(T[-1]) / 2) : ((T[] < 0) ? (sq(T[]) / 2) : 0)))) /Delta + alpha * dt * (T[-1] - 2. * T[] + T[1]) / sq(Delta);
    }
}

int main()
{
    init_grid(256);
    L0 = 2. * pi;
    X0 = -L0 / 2;
    periodic(right);
    DT = 1e-3;
    TOLERANCE = 1e-6;
    sim = 0;
    run();
    sim++;
    run();
}

event init(t = 0)
{
    foreach ()
        u[] = exp(-sq(2 * x));
}

event solveB(i++)
{
    dt = dtnext(DT);
    if (sim == 0)
    {
        scalar b[];
        foreach ()
            b[] = u[];
        solve(u,
              u[] + dt / 2 * u[] * (u[1] - u[-1]) / (2. * Delta) - alpha * dt / 2 * (u[-1] - 2. * u[] + u[1]) / sq(Delta),
              b);
        foreach ()
            u[] += (u[] - b[]);
    }
    if (sim == 1)
    {
        solve_GD(u, dt);
    }
}

event printdata(t = 0.; t <= 4.; t += 30 * dt)
{
    char names[30];

    if (sim == 0)
    {
        sprintf(names, "outputfile_SH");
    }
    if (sim == 1)
    {
        sprintf(names, "outputfile_BG");
    }
    static FILE *fp = fopen(names, "w");
    foreach ()
    {
        fprintf(fp, "%g %g %g\n", t, x, u[]);
    }
    fputs("\n", fp);
    fflush(fp);
}

/**
After adding some diffusivity, somehow multigrid method is able to preserve the discontinuity. [Or]()? change viscosity value...

~~~gnuplot

set terminal svg size 1200,600 enhanced font 'Times-Roman,16'
set key samplen 2 spacing 1.5 font 'Times-Roman,16'

set multiplot layout 1,2
set pm3d map
set xlabel 'x'
set ylabel 't'
set xrange [-3.5:3.5]
set yrange [0:4]
set cbrange [0:1]
set title "solve.h"
splot 'outputfile_SH' u 2:1:3
set title "Godunov"
splot 'outputfile_BG' u 2:1:3 
unset multiplot

~~~

*/