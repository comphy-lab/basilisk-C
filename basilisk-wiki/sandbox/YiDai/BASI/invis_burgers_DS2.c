/**
# solve inviscid burgers equation 
Here presenting a different initial condition, check out the method description [here](http://basilisk.fr/sandbox/YiDai/BASI/invis_burgers_DS.c)

*/


#include "grid/cartesian1D.h"
// #include "grid/bitree.h"
#include "run.h"

scalar T[], Ts[], dT[];
double dt;
#define vis 1.0
int j;

void solve_LF(scalar T, double dt)
{
    foreach ()
        T[] = 0.5 * (T[1] + T[-1]) - dt * 0.5 * (sq(T[1]) - sq(T[-1])) / (2 * Delta);
}

void solve_MC(scalar T, double dt)
{
    scalar Ts[];
    foreach ()
    {
        Ts[] = T[] - dt * 0.5 * (sq(T[1]) - sq(T[0])) / Delta;
        T[] = 0.5 * (T[] + Ts[]) - dt * 0.5 * (sq(Ts[]) - sq(Ts[-1])) / (2 * Delta);
    }
}

void solve_GD(scalar T, double dt)
{
    foreach ()
    {
        T[] = T[] - dt * (((T[] >= T[1]) ? \ 
        (((T[] + T[1]) / 2 > 0) ? (sq(T[]) / 2) : (sq(T[1]) / 2))
                                         : \ 
        ((T[] > 0) ? (sq(T[]) / 2) : ((T[1] < 0) ? (sq(T[1]) / 2) : 0))) \ 
        - ((T[-1] >= T[]) ? (((T[-1] + T[]) / 2 > 0) ? (sq(T[-1]) / 2) : (sq(T[]) / 2)) : ((T[-1] > 0) ? (sq(T[-1]) / 2) : ((T[] < 0) ? (sq(T[]) / 2) : 0)))) /
                        Delta;
    }
}

int main()
{
    L0 = 4;
    X0 = -1.;
    N = 1 << 9;
    // DT = sq(L0/N)/4;
    DT = 2e-4;
    for (j = 1; j < 4; j++)
        run();
}

T[left] = dirichlet(0.);
T[right] = dirichlet(0.);

event init(t = 0)
{
    foreach ()
    {
        T[] = exp(-sq(2 * (x - 1)));
    }
    boundary({T});
}

event integration(i++)
{
    // DT = sq(L0/(1 << grid->maxdepth))/3.;    // smaller time step when the grid is refined
    double dt = DT;
    dt = dtnext(dt);
    // centeral difference
    if (j == 1)
    {
        solve_LF(T, dt);
    }

    if (j == 2)
    {
        solve_MC(T, dt);
    }

    if (j == 3)
    {
        solve_GD(T, dt);
    }
    boundary({T});
}

// event printdata (t = 0; t <= 1000 * DT; t += 100 * DT) {
// event printdata (t = 0; t <= 1; t += 0.01) {

event printdata(t = {0, 0.2, 0.5, 1, 1.5, 2})
{
    char name[80];
    if (j == 1)
        sprintf(name, "BG_LF%g", t * 10);
    if (j == 2)
        sprintf(name, "BG_MC%g", t * 10);
    if (j == 3)
        sprintf(name, "BG_GD%g", t * 10);
    FILE *fp = fopen(name, "w");
    foreach ()
        fprintf(fp, "%g %g %g\n", x, T[], t);
    fprintf(fp, "\n\n");
    fclose(fp);
}


event end(t = 2) {}

/**

~~~gnuplot
reset
file1="BG_LF"
file2="BG_MC"
file3="BG_GD"
set terminal svg size 1200,600 enhanced font 'Times-Roman,12'
set key samplen 2 spacing 1.5 font 'Times-Roman,12'

set key bottom right
set multiplot layout 1,3
set key right top
plot[][-0.05:1.4] file1.'0'  u ($1):($2) t "t=0"   w lp, \
file1.'2'  u ($1):($2) t "t=0.2" w lp, \
file1.'5'  u ($1):($2) t "t=0.5" w lp, \
file1.'10' u ($1):($2) t "t=1.0" w lp, \
file1.'15' u ($1):($2) t "t=1.5" w lp, \
file1.'20' u ($1):($2) t "t=2.0" w lp

set nokey
plot[][-0.05:1.4] file2.'0'  u ($1):($2)  w lp, \
file2.'2'  u ($1):($2) w lp, \
file2.'5'  u ($1):($2) w lp, \
file2.'10' u ($1):($2) w lp, \
file2.'15' u ($1):($2) w lp, \
file2.'20' u ($1):($2) w lp

set nokey
plot[][-0.05:1.4] file3.'0'  u ($1):($2)  w lp, \
file3.'2'  u ($1):($2) w lp, \
file3.'5'  u ($1):($2) w lp, \
file3.'10' u ($1):($2) w lp, \
file3.'15' u ($1):($2) w lp, \
file3.'20' u ($1):($2) w lp
unset multiplot
~~~
*/
