/** Reprise du cas test "TIC6.c" : canal de pente constante, débit constant, recherche d'un état stationnaire uniforme. "TIC6.c" a montré qu'un tel régime n'était pas atteignable avec de simples conditions de bord "ouvertes". Une première solution est d'intégrer la pente dans le terme source réservé d'habitude aux frottements ("TIC6-tilt.c"). Ceci est un essai pour trouver une autre solution, avec une méthode basée sur les invariants de Riemann.
*/

//26.03.2021
//Stage Injection-Basilisk
//TIC13.c (TIC = Test-Injection-Canal)
//But : reprendre TIC6.c en implementatnt les invariants de Riemann (version "semi-implicite" cad avec q_N+1^k+1 = h_N+1^k+1 * u_N^k) en aval.

#include <math.h>

#include "grid/cartesian1D.h"
#include "saint-venant.h"
#include "manning_n_cst.h"

#define TMAX 13000. // duree de la simulation [s]
#define QBASE 1.65 // debit de base unitaire [m^2/s]
#define S0 0.001 // oppose de la pente du canal [m/m]
#define NMANNING 0.02 // coefficient de Manning [S.I.]

/* Le coefficient de Manning sert à calculer Sf ET a calculer h_base (hauteur normale pour le debit de base).
Si on ne veut pas tenir compte des frottements :
 + commenter " #include "SOURCE/manning.h" " ;
 + commenter " n  = NMANNING; " ;
 + LAISSER " h_base = pow( NMANNING*NMANNING* QBASE*QBASE / S0 , 3./10. ); "
  --> car on a besoin de h_base par la suite, comme hauteur correspondant au debit de base.
*/


#define NB_ITER_OUTPUTS 52 // nombre d'iterations entre chaque output pour les tableaux

#define MAX(a,b) (a > b ? a : b) // function "maximum"

double h_base;

scalar Ec[];
scalar q[];
scalar hc[];
scalar hn[];

double q_injection(double t) // debit injecte en amont au temps t
{
    return(QBASE);
}

double h_CL_right(double h_Nk, double zb_N1k, double zb_Nk, double Delta_Nk, double n_Nk, double u_Nk)
{
    double c_Nk  = sqrt(G*h_Nk);
    double S0_Nk = (zb_N1k-zb_Nk)/Delta_Nk;
    double Sf_Nk = (n_Nk*n_Nk)*(u_Nk*u_Nk)/pow(h_Nk,4./3.);
    double S_Nk  = S0_Nk - Sf_Nk;

    double h_N1k1 = c_Nk*(c_Nk/G+S_Nk*dt);

    return(h_N1k1);
}

double u_CL_right(double u_Nk, double h_Nk, double h_N1k1)
{
    double q_Nk = u_Nk*h_Nk;

    double u_N1k1 = q_Nk/h_N1k1;

    return(u_N1k1);
}

int main()
{
    G  = 9.81; // acceleration de pesanteur [m/s^2]
    L0 = 3000.; // longueur du canal/bief [m]
    X0 = 0.; // abscisse de debut du canal [m]
    N  = 256; // nombre de cellules sur le domaine
    run();
}

h[left] = neumann(0);


event init(t=0)
{

    // Conditions aux limites "invariables"

    zb[left] = neumann(+S0);
    zb[right] = neumann(-S0);

    // conditions initiales
    foreach()
    {
        nmanning[] = NMANNING; // coefficient de Manning [S.I.]
    }

    h_base = pow( NMANNING*NMANNING * QBASE*QBASE / S0 , 3./10. );

    foreach()
    {
        zb[] = -S0 * x;
        h[]  = h_base;
        u.x[] = QBASE/h_base;
    }

    // mise a jour des champs
    foreach()
    {
        Ec[] = u.x[] * u.x[] / (2. * G);
        q[]  = h[] * u.x[];
        hc[] = pow( q[] * q[] / G, (1./3.) );
        hn[] = pow( nmanning[]*nmanning[] * q[]*q[] / S0 , 3./10. );

    }

    boundary(all);

}

double h_end,zb_end,h_ghost,zb_ghost,n_end,u_end;

event flows(i++) // mise a jour des CL variables et des champs donnant des indications sur l'ecoulement
{   
    //  CL aval "invariants de Riemann"
//    Point point = locate(L0-L0/N/2.);
  foreach(){
    h_end = h[];
    zb_end = zb[];
    h_ghost = h[1];
    zb_ghost = zb[1];
    n_end = nmanning[];
    u_end = u.x[];
  }
//    h[right] = h_CL_right(h[], zb[right], zb[], Delta, nmanning[], u.x[]);
//    u[right] = u_CL_right(u.x[], h[], h[right]);
    h[right] = h_CL_right(h_end,zb_ghost,zb_end,L0/N,n_end,u_end);
    u.n[right] = u_CL_right(u_end,h_end,h_ghost);

    //  CL amont "classiques"
    u.n[left] =  q_injection(t)/h[left];
    h[left] = neumann(0);

    // mise a jour des champs
    foreach()
    {
        Ec[] = u.x[] * u.x[] / (2. * G);
        q[]  = h[] * u.x[];
        hc[] = pow( q[] * q[] / G, (1./3.) );
        hn[] = pow( nmanning[]*nmanning[] * q[]*q[] / S0 , 3./10. );

    }


    boundary(all);
}

event outputs_tab(i+=NB_ITER_OUTPUTS)
{
    static FILE * tab_i_t = fopen( "TIC13.tab-i-t.csv", "write");

    fprintf(tab_i_t, "i,t; %d; %g\n", i, t); // tableau avec le temps et les nos. d'iteration

    static FILE * tab_outputs = fopen( "TIC13.tab-outputs.csv", "write");

    // ecriture des abscisses
    fprintf(tab_outputs, "X_i%d", i); // temps
    foreach()
    {
        fprintf(tab_outputs, ";%g", x);
    }
    fprintf(tab_outputs, "%s", "\n");

    // ecriture de la topographie
    fprintf(tab_outputs, "Zb_i%d", i); // temps (le meme qu'a la ligne precedente)
    foreach()
    {
        fprintf(tab_outputs, ";%g", zb[]);
    }
    fprintf(tab_outputs, "%s", "\n");

    // ecriture des surfaces libres
    fprintf(tab_outputs, "Eta_i%d", i); // temps (le meme qu'a la ligne precedente)
    foreach()
    {
        fprintf(tab_outputs, ";%g", eta[]);
    }
    fprintf(tab_outputs, "%s", "\n");

    // ecriture des energies cinetiques
    fprintf(tab_outputs, "Ec_i%d", i); // temps (le meme qu'a la ligne precedente)
    foreach()
    {
        fprintf(tab_outputs, ";%g", Ec[] );
    }
    fprintf(tab_outputs, "%s", "\n");

    // ecriture des hauteurs critiques
    fprintf(tab_outputs, "hc_i%d", i); // temps (le meme qu'a la ligne precedente)
    foreach()
    {
        fprintf(tab_outputs, ";%g", hc[] );
    }
    fprintf(tab_outputs, "%s", "\n");

    // ecriture des hauteurs normales
    fprintf(tab_outputs, "hn_i%d", i); // temps (le meme qu'a la ligne precedente)
    foreach()
    {
        fprintf(tab_outputs, ";%g", hn[] );
    }
    fprintf(tab_outputs, "%s", "\n");

    static FILE * tab_bound_cond = fopen( "TIC13.tab-cond-lim.csv", "write");
    if (i==0)
    {
        //fprintf(tab_bound_cond, "i; t; zb_r; h_r; u.n_r; eta_r\n");
        fprintf(tab_bound_cond, "i; t; zb_L0; zb_r_ghost; h_L0; h_r_ghost; u.n_L0; u.n_r_ghost; eta_L0; eta_r_ghost\n");
    }

    // ecriture des conditions aux limites

    Point point = locate(L0-L0/N/2.);
    //fprintf(tab_bound_cond, "%d; %g; %g; %g; %g; %g\n", i, t, zb[+1], h[+1], u.x[+1], eta[+1]);
    fprintf(tab_bound_cond, "%d; %g; %g; %g; %g; %g; %g; %g; %g; %g\n", i, t, zb[0], zb[right], h[0], h[right], u.x[0], u.x[right], eta[0], eta[right]);

}

/*
Gauge gauging_stations[] =
{

    {"TIC13.ga-Zb-Eta-x0250.txt", 250., 0., "Zb_x0250 Eta_x0250"}, // VIRGULE !
    {"TIC13.ga-Zb-Eta-x1500.txt", 1500., 0., "Zb_x1500 Eta_x1500"},
    {"TIC13.ga-Zb-Eta-x2750.txt", 2750., 0., "Zb_x2750 Eta_x2750"},
    {NULL} // PAS de ";" (pourquoi ? --> car c'est une liste (et pas un evenement) ! )

};

event gauging_event(i+=NB_ITER_OUTPUTS)
{
    output_gauges( gauging_stations, {h} );
}
*/

event end(t = TMAX)
{
    printf("Programme termine : iteration no. %d ; temps : %g\n", i, t);
    printf("acceleration G = %g\n", G);
    printf("nombre de cellules N = %d\n", N);
    printf("hauteur normale correspondant au debit de base h_base = %g\n", h_base);
}



