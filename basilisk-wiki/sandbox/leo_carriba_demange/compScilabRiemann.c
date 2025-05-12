#include <math.h>

#include "grid/cartesian1D.h"
#include "saint-venant.h"
#include "headers/manning-tilt.h"
#include "headers/Interpol_Txt_2.h"
#include "headers/Riemann_Invariants_Upstream.h"

#define TMAX 10000. // duree de la simulation [s]
#define QBASE 1.65// debit de base unitaire [m^2/s]
#define S0 0. // topographie apparente ("detrendee" par le tilt)
#define TILT 0.001 // rempace la pente pour "manning-tilt.h"
#define NMANNING 0.02 // coefficient de Manning [S.I.]
#define X_INJECTION 1200. // abscisse d'injection d'un hydrogramme recupere au meme endroit dans "projet-GE4-reference"
#define NOM_FICHIER_HYDROGRAMME "data/hydrogramme_q_x1200.txt" // nom du fichier texte contenant l'hydrogramme
#define DELTA_T_MIN 1. // pas de temps par defaut [s] pour l'evaluation des derivees (ex : dans un hydrogramme) si dt = 0.

scalar topo_reelle[]; // topographie reelle, en tenant compte de la pente

double hauteur_normale(double manning, double debit, double pente)
{
    return( manning > 0 && pente > 0 ? pow( manning * manning * debit * debit / pente, 3./10. ) : -1 );
}

#define NB_ITER_OUTPUTS 1 // nombre d'iterations entre chaque output pour les tableaux
/*
/!\ adapter pour ne pas creer des .csv trop grands /!\
*/

#define MAX(a,b) (a > b ? a : b) // function "maximum"



int main()
{
    G  = 9.81; // acceleration de pesanteur [m/s^2]
    L0 = 3000. - X_INJECTION; // longueur du canal/bief [m]
    X0 = X_INJECTION; // abscisse de debut du canal [m]
    N  = 256; // nombre de cellules sur le domaine
    tilt.x = TILT;
    run();
}

event init(t=0)
{
    foreach()
    {
        nmanning[] = NMANNING;
    }

    double h_base = hauteur_normale(NMANNING, QBASE, TILT);

    foreach()
    {
        zb[] = 0.;
        h[]  = h_base;
        u.x[] = QBASE / h_base;
        topo_reelle[] = - TILT * x;
    }

    // Conditions aux limites

    h[left] = dirichlet(0);
    u.n[left] = dirichlet(0);

    h[right] = neumann(0);
    u.n[right] = neumann(0);

    boundary(all);

}


double q0, u_interior, h_interior, Sf_interior, S0_interior;
double u_ghost, h_ghost;
event flows(i++) // mise a jour des CL
{
    // injection des debits
    q0 = interpolation_lineaire_textfile(NOM_FICHIER_HYDROGRAMME, t + dt);

    u_interior = interpolate(u.x, X_INJECTION + L0/N/2.);
    h_interior = interpolate(h,   X_INJECTION + L0/N/2.);

    S0_interior = TILT;
    Sf_interior =  (h_interior > 0. ? u_interior * u_interior * NMANNING * NMANNING / pow(h_interior, 4./3.) : 0. );

    Solve_Riemann_Invariants_Upstream(q0, G, dt, u_interior, h_interior, Sf_interior, S0_interior, &u_ghost, &h_ghost);

    u.n[left] = u_ghost;
    h[left]   = h_ghost;
}

event outputs_jx(i=1)
{
    static FILE * tab_j_x = fopen( "compScilabRiemann.tab-j-x.csv", "write");

    int j = 0; // indice d'espace, j = 1:N

    foreach()
    {
        j += 1;
        fprintf(tab_j_x, "%d; %g;\n", j, x); // tableau avec le temps et les nos. d'iteration
    }


}

event outputs_tab(i+=NB_ITER_OUTPUTS)
{
    int j; // indice d'espace, j = 1:N

    static FILE * tab_i_t = fopen( "compScilabRiemann.tab-i-t.csv", "write");

    fprintf(tab_i_t, "%d; %g;\n", i, t); // tableau avec le temps et les nos. d'iteration

    static FILE * tab_h = fopen( "compScilabRiemann.tab-h.csv", "write");

    // ecriture des abscisses
    j = 0;
    foreach()
    {
        j += 1;
        fprintf(tab_h, (j < N ? "%g; " : "%g") , h[]);
    }
    fprintf(tab_h, "%s", "\n");

    static FILE * tab_u = fopen( "compScilabRiemann.tab-u.csv", "write");

    // ecriture des abscisses
    j = 0;
    foreach()
    {
        j += 1;
        fprintf(tab_u, (j < N ? "%g; " : "%g"), u.x[]);
    }
    fprintf(tab_u, "%s", "\n");

    static FILE * tab_topotilt = fopen( "compScilabRiemann.tab-topo-tilt.csv", "write");

    // ecriture des abscisses
    j = 0;
    foreach()
    {
        j += 1;
        fprintf(tab_topotilt, (j < N ? "%g; " : "%g"), zb[]);
    }
    fprintf(tab_topotilt, "%s", "\n");

    static FILE * tab_toporeelle = fopen( "compScilabRiemann.tab-topo-reelle.csv", "write");

    // ecriture des abscisses
    j = 0;
    foreach()
    {
        j += 1;
        fprintf(tab_toporeelle, (j < N ? "%g; " : "%g"), topo_reelle[]);
    }
    fprintf(tab_toporeelle, "%s", "\n");
}


event end(t = TMAX)
{
    printf("Programme termine : iteration no. %d ; temps : %g\n", i, t);
}



