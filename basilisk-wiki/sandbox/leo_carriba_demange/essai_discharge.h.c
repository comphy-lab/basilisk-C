//13.05.2021
//Stage Injection-Basilisk
//compScilabDischarge.h.c = comparaison des resultats obtenus avec Scilab (projet GE4) et discharge.h (inclus dans le code source de Basilisk)

#include <math.h>

#include "grid/cartesian1D.h"
#include "saint-venant.h"
#include "manning-tilt.h"
//#include "Interpol_Txt_2.h"
#include "discharge.h"

#define TMAX 12500. // duree de la simulation [s] : il ne faut pas la faire durer strictement jusqu'au bout de l'hydrogramme ! (la lecture echoue)
#define QBASE 1.65// debit de base unitaire [m^2/s]
#define S0 0. // topographie apparente ("detrendee" par le tilt)
#define TILT 0.001 // rempace la pente pour "manning-tilt.h"
#define NMANNING 0.02 // coefficient de Manning [S.I.]
#define X_INJECTION 1200 // abscisse d'injection d'un hydrogramme recupere au meme endroit dans "projet-GE4-reference"
//#define NOM_FICHIER_HYDROGRAMME "hydrogramme_q_x1200.txt" // nom du fichier texte contenant l'hydrogramme
#define DELTA_T_MIN 1. // pas de temps par defaut [s] pour l'evaluation des derivees (ex : dans un hydrogramme) si dt = 0.
// sert aussi pour la resolution implicite des invariants de Riemann


scalar Ec[];
scalar q[];
scalar hc[];
scalar hn[];

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

double q_injection(double t) // debit injecte en amont au temps t
{
    double q_max = 10; // debit de pointe unitaire atteint par l'hydrogramme [m^2/s] (Qmax=50 [m3/s] dans le projet de fin de semestre)
    double t_0 = 0.; // temps ou commence la crue [s]
    double t_pic = 20.*60.; // temps ou le debit est maximum [s]
    double t_decrue = 60.*60.; // temps ou la decrue est terminee [s]
    double pente_crue = (q_max - QBASE) / (t_pic - t_0); // gradient de debit pendant la crue
    double pente_decrue = (QBASE - q_max) / (t_decrue - t_pic); // gradient de debit pendant la decrue

    //return ( QBASE + MAX( 0, ( t < t_pic ? pente_crue * (t - t_0) : pente_decrue * (t - t_pic) + q_max - QBASE) ) ); // debit au temps t
    return(QBASE);
}

int main()
{

    //mask (x < X_INJECTION ? left : none); // on masque ce qui est en amont du point d'injection : ne fonctionne pas !

    G  = 9.81; // acceleration de pesanteur [m/s^2]
    L0 = 3000.; // longueur du canal/bief [m]
    X0 = 0.; // abscisse de debut du canal [m]
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
        zb[] = (x < X_INJECTION ? 3.5 : 0.);
        h[]  = ( x < X_INJECTION ? 0. : h_base );
        u.x[] = ( x < X_INJECTION ? 0. : QBASE / h_base );
        topo_reelle[] = (x < X_INJECTION ? 3.5 - TILT * x : 0. - TILT * x);
    }

    // Conditions aux limites

    h[left] = dirichlet(0);
    u.n[left] = dirichlet(0);

    h[right] = neumann(0);
    u.n[right] = neumann(0);

    boundary(all);



}

event disp_iter_and_time(i++)
{
    printf("Simulation en cours... i = %d ; t = %g\n", i, t);
}

event flows(i++) // mise a jour des CL
{
    //injection des debits
    foreach()
    {
        if ( x >= X_INJECTION && x < X_INJECTION + Delta)
        {
            double q_inj, eta_inj;
            //q_inj = interpolation_lineaire_textfile(NOM_FICHIER_HYDROGRAMME, t + dt);
            q_inj = q_injection(t);
            eta_inj = eta_b(q_inj, left);
            h[] = MAX(eta_inj - zb[], 0.);
            u.x[] = q_inj / h[];
        }
    }

    // mise a jour des champs
    foreach()
    {
        Ec[] = u.x[] * u.x[] / (2. * G);
        q[]  = h[] * u.x[];
        hc[] = pow( q[] * q[] / G, (1./3.) );
        hn[] = (NMANNING > 0 && tilt.x > 0 ? hauteur_normale(NMANNING, q[], tilt.x) : -1);
    }
}

event end(t = TMAX)
{
    printf("Programme termine : iteration no. %d ; temps : %g\n", i, t);
    printf("acceleration G = %g\n", G);
    printf("nombre de cellules N = %d\n", N);
    printf("S0 = %g\n", S0);
    printf("tilt = %g\n", TILT);
}



