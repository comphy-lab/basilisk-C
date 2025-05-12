/**
The goal of this script is to give an input discharge on the upstream boundary (left) thanks to "discharge.h", with a one-dimensional grid. Unfortunately, it seems that "foreach_boundary()" (called by "discharge.h") does not work together with "grid/cartesian-1D.h", which causes a bug. Other bugs exist in this script, but the are apparently not linked with the bug caused by "discharge.h". I published this bug report because I found out that someone had the same problem on the forum : [https://groups.google.com/g/basilisk-fr/c/AdkEEEVTXno/m/32Dxwc-XBQAJ].
*/

#include <math.h>

#include "grid/cartesian1D.h"
#include "saint-venant.h"

#include "discharge.h"

#define TMAX 12500. // time length of the simulation
#define QBASE 1.65 // base flow

#define MAX(a,b) (a > b ? a :b)

/**
Main constants and call to the solver.
*/

int main()
{
    G  = 9.81;
    L0 = 3000.;
    X0 = 0.;
    N  = 256; 
    run();
}

/**
The initial conditions : base flow ; the boundary conditions : closed upstream, open downstream.
*/

event init(t=0)
{

    double h_base = 1.0291;

    foreach()
    {
        zb[]  = 0.;
        h[]   = h_base;
        u.x[] = QBASE / h_base;
     }

    h[left] = dirichlet(0);
    u.n[left] = dirichlet(0);

    h[right] = neumann(0);
    u.n[right] = neumann(0);

    boundary(all);

}

/**
This paragraph gives the sinusoidal input discharge (left boundary).
It uses "eta_b" from "discharge.h", which itself calls "foreach_boundary()" and causes the main bug.
*/

event flows(i++) // mise a jour des CL
{

    double q_inj, eta_inj;
    q_inj = QBASE + 0.8 * cos(t*6.28/50.);
    eta_inj = eta_b(q_inj, left);
    h[left] = MAX(eta_inj - zb[], 0.);
    u.x[left] = q_inj / h[];

}

/**
The next paragraphs are the outputs.
*/

event outputs_jx(i=1)
{
    static FILE * tab_j_x = fopen( "compScilabDischarge.h.tab-j-x.csv", "write");

    int j = 0; // indice d'espace, j = 1:N

    foreach()
    {
        j += 1;
        fprintf(tab_j_x, "%d; %g;\n", j, x); // tableau avec le temps et les nos. d'iteration
    }


}

event outputs_tab(i++)
{
    int j; // indice d'espace, j = 1:N

    static FILE * tab_i_t = fopen( "compScilabDischarge.h.tab-i-t.csv", "write");

    fprintf(tab_i_t, "%d; %g;\n", i, t); // tableau avec le temps et les nos. d'iteration

    static FILE * tab_h = fopen( "compScilabDischarge.h.tab-h.csv", "write");

    // ecriture des abscisses
    j = 0;
    foreach()
    {
        j += 1;
        fprintf(tab_h, (j < N ? "%g; " : "%g") , h[]);
    }
    fprintf(tab_h, "%s", "\n");

    static FILE * tab_u = fopen( "compScilabDischarge.h.tab-u.csv", "write");

    // ecriture des abscisses
    j = 0;
    foreach()
    {
        j += 1;
        fprintf(tab_u, (j < N ? "%g; " : "%g"), u.x[]);
    }
    fprintf(tab_u, "%s", "\n");


    static FILE * tab_zb = fopen( "compScilabDischarge.h.tab-zb.csv", "write");

    // ecriture des abscisses
    j = 0;
    foreach()
    {
        j += 1;
        fprintf(tab_zb, (j < N ? "%g; " : "%g"), zb[]);
    }
    fprintf(tab_zb, "%s", "\n");
}

/**
End of the program.
*/

event end(t = TMAX)
{
    printf("Programme termine : iteration no. %d ; temps : %g\n", i, t);
}
