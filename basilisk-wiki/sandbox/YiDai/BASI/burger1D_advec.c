/**
Here bcg.h (Bell-Collela-Glaz scheme) is used to solve the burgers equation, we also compare it with manual Godunov.
*/

#include "grid/bitree.h"
#include "run.h"
#include "bcg.h"


scalar u[];
face vector uf[];
int sim;

void solve_GD(scalar T, double dt)
{
    scalar usta[], F[];
    foreach ()
    {
        usta[] = (T[] >= T[1]) ? 
        (((T[] + T[1]) / 2. > 0) ? (T[]) : (T[1])):
        ((T[] > 0) ? (T[]) :((T[1] < 0) ? (T[1]):0));
        F[] = sq(usta[]) / 2.;
        T[] = T[] - dt * (F[] - F[-1])/Delta;
    }
}

int main()
{
    init_grid(512);
    L0 = 2. * pi;
    X0 = -L0 / 2;
    periodic(right);
    CFL = 10;
    DT = 1e-2;
    sim = 0;
    run();
    sim++;
    run();
}

event init(t = 0)
{
    foreach ()
        u[] = exp(-sq(2 * x));
    trash({uf});
    foreach_face()
        uf.x[] = face_value(u, 0);
}

event solveB(i++)
{
    dt = dtnext(DT);
    
    if (sim==0){
    scalar du = new scalar;
    if (u.gradient)
    {
        foreach()
            du[] = u.gradient(u[-1], u[], u[1]) / Delta;
    }
    else
    {
        foreach ()
            du[] = (u[1] - u[-1]) / (2. * Delta);
    }
    trash({uf});
    foreach_face ()
    {
        double un = dt * (u[] + u[-1]) / (2. * Delta), s = sign(un);
        int i = -(s + 1.) / 2.;
        uf.x[] = u[i] + s * (1. - s * un) * du[i] * Delta/2;
    }
    // since we are using face value of uf at Delta/2, the time advancement should be dt/2, check eq3.4 in Bell-Collela-Glaz 1989
    advection({u}, uf, dt/2);
    }
    if (sim==1){
        solve_GD(u, dt);
    }
}

event printdata(t = 0.; t <= 4.; t += 30 * dt)
{
    char fname[30];
    if (sim == 0){
        sprintf(fname, "advec");
    }
    if (sim == 1)
    {
        sprintf(fname, "BG");
    }
    static FILE *fp = fopen(fname, "w");
    foreach ()
    {
        fprintf(fp, "%g %g %g\n", t, x, u[]);
    }
    fputs("\n", fp);
    fflush(fp);
}

/**

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

set title "bcg.h"
splot 'advec' u 2:1:3
set title "manual Godunov"
splot 'BG' u 2:1:3
unset multiplot



~~~



*/