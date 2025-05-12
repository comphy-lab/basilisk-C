/**
# solving inviscid burgers equation
The equation 
$$\partial_t U + U \partial_x U = 0$$

Here we compare three numerical methods: 

* 1. Lax-Friedrichs
$$U_j^{n+1}=\frac{1}{2}\left(U_{j-1}^n+U_{j+1}^n\right)-\frac{k}{2 h}\left[\frac{1}{2}\left(U_{j+1}^n\right)^2-\frac{1}{2}\left(U_{j-1}^n\right)^2\right]$$

* 2. MacCormack

$$\begin{aligned}
U_j^* &=U_j^n-\frac{k}{h}\left[f\left(U_{j+1}^n\right)-f\left(U_j^n\right)\right] \\
U_j^{n+1} &=\frac{1}{2}\left(U_j^n+U_j^*\right)-\frac{k}{2 h}\left[f\left(U_j^*\right)-f\left(U_{j-1}^*\right)\right]
\end{aligned}$$

* 3. Godunov
$$U_j^{n+1}=U_j^n-\frac{k}{h}\left[F\left(U_j^n, U_{j+1}^n\right)-F\left(U_{j-1}^n, U_j^n\right)\right] $$

where $F(U, V)=\frac{\left(u^*\right)^2}{2}$ 

if $U >= V$, 
$$u^*=\left\{\begin{array}{ll}
U, & \text { if } \frac{U+V}{2}>0 \\
V, & \text { in other case. }
\end{array}\right.$$

if $U<V$,
$$u^*=\left\{\begin{array}{ll}
U, & \text { if } U>0 \\
V, & \text { if } V<0 \\
0, & \text { if } U \leq 0 \leq V
\end{array}\right.$$

*/


#include "grid/cartesian1D.h"
// #include "grid/bitree.h"   // check out if it makes a difference
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
        (((T[] + T[1]) / 2 > 0) ? (sq(T[]) / 2) : (sq(T[1]) / 2)) : \ 
        ((T[] > 0) ? (sq(T[]) / 2) : ((T[1] < 0) ? (sq(T[1]) / 2) : 0))) \ 
        - ((T[-1] >= T[]) ? \
        (((T[-1] + T[]) / 2 > 0) ? (sq(T[-1]) / 2) : (sq(T[]) / 2)) : \
        ((T[-1] > 0) ? (sq(T[-1]) / 2) : ((T[] < 0) ? (sq(T[]) / 2) : 0)))) \
        /Delta;
    }
}

int main()
{
    L0 = 1;
    // X0 = -L0/2.;
    N = 1 << 9;
    // DT = sq(L0/N)/4;
    DT =2e-4;
    for (j = 1; j<4; j++)
        run();
}

T[left] = dirichlet(0.);
T[right] = dirichlet(0.);

event init(t = 0){
    foreach(){
        T[] = sin(2*pi*x);
    }
    boundary ({T});
}

event integration(i++)
{
    // DT = sq(L0/(1 << grid->maxdepth))/3.;    // smaller time step when the grid is refined
    double dt = DT;
    dt = dtnext(dt);
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

event end(t = 2){}

/**
oh interesting, Godunov method is able to preserve the shock wave structure. Check out a [differnet initial condition](http://basilisk.fr/sandbox/YiDai/BASI/invis_burgers_DS2.c) here
~~~gnuplot
reset
file1="BG_LF"
file2="BG_MC"
file3="BG_GD"
set terminal png size 1200,600 enhanced font 'Times-Roman,16'
set key samplen 2 spacing 1.5 font 'Times-Roman,16'

set key bottom right
set multiplot layout 1,3
set key right top
plot[0:1][-1.1:1.1] file1.'0'  u ($1):($2) t "t=0"   w lp, \
file1.'2'  u ($1):($2) t "t=0.2" w lp, \
file1.'5'  u ($1):($2) t "t=0.5" w lp, \
file1.'10' u ($1):($2) t "t=1.0" w lp, \
file1.'15' u ($1):($2) t "t=1.5" w lp, \
file1.'20' u ($1):($2) t "t=2.0" w lp

set nokey
plot[][-1.1:1.1] file2.'0'  u ($1):($2)  w lp, \
file2.'2'  u ($1):($2) w lp, \
file2.'5'  u ($1):($2) w lp, \
file2.'10' u ($1):($2) w lp, \
file2.'15' u ($1):($2) w lp, \
file2.'20' u ($1):($2) w lp

set nokey
plot[][-1.1:1.1] file3.'0'  u ($1):($2)  w lp, \
file3.'2'  u ($1):($2) w lp, \
file3.'5'  u ($1):($2) w lp, \
file3.'10' u ($1):($2) w lp, \
file3.'15' u ($1):($2) w lp, \
file3.'20' u ($1):($2) w lp
unset multiplot
~~~
*/


/**
~~~gnuplot
reset
file1="BG_LF"
file2="BG_MC"
file3="BG_GD"
set terminal svg size 1200,600 enhanced font 'Times-Roman,14'
set key samplen 2 spacing 1.5 font 'Times-Roman,14'

set key bottom right
set multiplot layout 1,3 spacing 0.01, 0.01
set key right top
plot[0.4:0.6][-0.5:0.5] file1.'0'  u ($1):($2) t "t=0"   w lp, \
file1.'2'  u ($1):($2) t "t=0.2" w lp, \
file1.'5'  u ($1):($2) t "t=0.5" w lp, \
file1.'10' u ($1):($2) t "t=1.0" w lp, \
file1.'15' u ($1):($2) t "t=1.5" w lp, \
file1.'20' u ($1):($2) t "t=2.0" w lp

set nokey
unset ytics
unset ylabel
plot[0.4:0.6][-0.5:0.5] file2.'0'  u ($1):($2)  w lp, \
file2.'2'  u ($1):($2) w lp, \
file2.'5'  u ($1):($2) w lp, \
file2.'10' u ($1):($2) w lp, \
file2.'15' u ($1):($2) w lp, \
file2.'20' u ($1):($2) w lp

set nokey
unset ytics
unset ylabel
plot[0.4:0.6][-0.5:0.5] file3.'0'  u ($1):($2)  w lp, \
file3.'2'  u ($1):($2) w lp, \
file3.'5'  u ($1):($2) w lp, \
file3.'10' u ($1):($2) w lp, \
file3.'15' u ($1):($2) w lp, \
file3.'20' u ($1):($2) w lp
unset multiplot

~~~
*/