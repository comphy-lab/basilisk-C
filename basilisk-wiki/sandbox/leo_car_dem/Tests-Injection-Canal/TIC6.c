/**
INTRODUCTION :
La simulation qui suit reproduit un "canal de laboratoire" : il s'agit d'un modèle 1D avec une topographie à pente constante. 
Mon but est de tester plusieurs manières d'implémenter une injection de débit en amont et de vérifier si la réaction à l'intérieur du canal est conforme à des situations de référence.
Comme premier essai, je tente de reproduire un modèle simple vu à l'université, avec un débit imposé en amont et des conditions de bord libres en aval. Je teste deux cas : régime stationnaire uniforme (débit constant et hauteur normale) ou hydrogramme triangulaire en entrée.
*/

/**
Je commence par importer les différents packages nécessaires. Le modèle est en une dimension. Le terme source de frottements est celui de G. Kirstetter (http://basilisk.fr/sandbox/geoffroy/sourceterm/manning.h). Le solveur est celui du système de Saiut-Venant.
*/

#include <math.h>
#include "grid/cartesian1D.h"
#include "saint-venant.h"
#include "SOURCE/manning.h"

/**
Trois paramètres sont facilement modifiables pour tester plusieurs configurations : la durée de simulation, le débit unitaire de base et la pente du canal.
*/

#define TMAX 13000. // duree de la simulation [s]
#define QBASE 1.65 // debit de base unitaire [m^2/s]
#define S0 0.001 // oppose de la pente du canal [m/m]

/** Le coefficient de Manning sert à calculer Sf ET a calculer h_base (hauteur normale pour le debit de base).
Si on ne veut pas tenir compte des frottements : 1°) commenter " #include "SOURCE/manning.h" " ; 2°) commenter " n  = NMANNING; " ; 3°) LAISSER " h_base = pow( NMANNING*NMANNING* QBASE*QBASE / S0 , 3./10. ); " --> Tout ceci car on a besoin de h_base par la suite, comme hauteur correspondant au debit de base.
*/

#define NMANNING 0.02 // coefficient de Manning [S.I.]

/**Les variables globales sont la hauteur normale correspondant au débit de base, et des champs scalaires donnant des informations sur la nature de l'écoulement.
*/

double h_base;

scalar Ec[]; // energie cinetique de l'ecoulement
scalar q[]; // debit unitaire
scalar hc[]; // hauteur critique
scalar hn[]; // hauteur normale

/**
On peut faire varier la fréquence des sorties (graphiques ou textuelles). Une macro simple permet de calculer un maximum.
*/

#define NB_ITER_OUTPUTS 52 // nombre d'iterations entre chaque output pour les tableaux

#define MAX(a,b) (a > b ? a : b) // function "maximum"

/**
La fontion principale regroupe des constantes géométriques et physiques.
*/

int main()
{
    G  = 9.81; // acceleration de pesanteur [m/s^2]
    n  = NMANNING; // coefficient de Manning [S.I.]
    L0 = 3000.; // longueur du canal/bief [m]
    X0 = 0.; // abscisse de debut du canal [m]
    N  = 256; // nombre de cellules sur le domaine
    run();
}

/**
Une fonction calcule le débit à injecter à chaque instant. Deux choix sont possibles : hydrogramme triangulaire ou débit de base constant.
*/

double q_injection(double t) // debit injecte en amont au temps t
{
    double q_max = 10; // debit de pointe unitaire atteint par l'hydrogramme [m^2/s] (Qmax=50 [m3/s] dans le projet de fin de semestre)
    double t_0 = 0.; // temps ou commence la crue [s]
    double t_pic = 20.*60.; // temps ou le debit est maximum [s]
    double t_decrue = 60.*60.; // temps ou la decrue est terminee [s]
    double pente_crue = (q_max - QBASE) / (t_pic - t_0); // gradient de debit pendant la crue
    double pente_decrue = (QBASE - q_max) / (t_decrue - t_pic); // gradient de debit pendant la decrue

    //return ( QBASE + MAX( 0, ( t < t_pic ? pente_crue * (t - t_0) : pente_decrue * (t - t_pic) + q_max - QBASE) ) ); // debit au temps t (hydrogramme triangulaire)
    return(QBASE); // debit de base constant
}

/** 
Les conditions aux limites sont un débit imoposé en amont et une hauteur d'eau et une vitesse libres en aval, quel que soit le régime d'écoulement.
Les conditions initiales correspondent à un régime stationnaire uniforme.
*/

event init(t=0)
{

    // Conditions initiales

    zb[left] = neumann(+S0);
    zb[right] = neumann(-S0);

    h[left] = neumann(0);
    h[right] = neumann(0);

    u.n[left] =  q_injection(t)/h[left];
    u.n[right] = neumann(0);


    h_base = pow( NMANNING*NMANNING * QBASE*QBASE / S0 , 3./10. );
    //h_base = 1.026;

    foreach()
    {
        zb[] = -S0 * x;
        h[]  = h_base;
        u.x[] = QBASE/h_base;
    }

    boundary(all);

}


/** 
A chaque iteration, les conditions aux limites et les champs donnant des renseignements sur la nature de l'écoulement sont mis à jour.
*/

event flows(i++) 
{

    boundary(all); // mise a jour des CL

    // mise a jour des champs
    foreach()
    {
        Ec[] = u.x[] * u.x[] / (2. * G);
        q[]  = h[] * u.x[];
        hc[] = pow( q[] * q[] / G, (1./3.) );
        hn[] = pow( NMANNING*NMANNING * q[]*q[] / S0 , 3./10. );
    }
}

/** 
Des sorties sont possibles sous forme d'images ".ppm".
*/

event outputs_movies(NB_ITER_OUTPUTS)
{
    output_ppm(zb, file = "TIC6.vue-dessus-zb.mp4");
    output_ppm(h, file = "TIC6.vue-dessus-h.mp4");
    output_ppm(eta, file = "TIC6.vue-dessus-eta.mp4");
}

/**
D'autres sorties sont possibles sous forme de fichiers textes (à traiter par exemple sous R).
*/

event outputs_tab(i+=NB_ITER_OUTPUTS)
{
    static FILE * tab_i_t = fopen( "TIC6.tab-i-t.csv", "write");

    fprintf(tab_i_t, "i,t; %d; %g\n", i, t); // tableau avec le temps et les nos. d'iteration

    static FILE * tab_outputs = fopen( "TIC6.tab-outputs.csv", "write");

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

    static FILE * tab_bound_cond = fopen( "TIC6.tab-cond-lim.csv", "write");
    if (i==0)
    {
        //fprintf(tab_bound_cond, "i; t; zb_r; h_r; u.n_r; eta_r\n");
        fprintf(tab_bound_cond, "i; t; zb_L0; zb_r_ghost; h_L0; h_r_ghost; u.n_L0; u.n_r_ghost; eta_L0; eta_r_ghost\n");
    }

    // ecriture des conditions aux limites

    Point point = locate(L0-L0/N/2.);
    //fprintf(tab_bound_cond, "%d; %g; %g; %g; %g; %g\n", i, t, zb[+1], h[+1], u.x[+1], eta[+1]);
    fprintf(tab_bound_cond, "%d; %g; %g; %g; %g; %g; %g; %g; %g; %g\n", i, t, zb[0], zb[+1], h[0], h[+1], u.x[0], u.x[+1], eta[0], eta[+1]);

    /*
    fprintf(tab_bound_cond, "i; t; zb_l; zb_r; h_l; h_r; u.n_l; u.n_r; eta_l; eta_r\n");
    fprintf(tab_bound_cond, "%d; %g; %g; %g; %g; %g; %g; %g; %g; %g\n", i, t, zb[left], zb[right], h[left], h[right], u.n[left], u.n[right], eta[left], eta[right]);
    */

    /*
    fprintf(tab_bound_cond, "i; t; zb_l\n");
    fprintf(tab_bound_cond, "%d; %g; %g\n", i, t, zb[left]);
    */
}

/**
Troisième type de sorties : on peut enregistrer des "hydrogrammes" en un point donné.
*/

/*
Gauge gauging_stations[] =
{

    {"TIC6.ga-Zb-Eta-x0250.txt", 250., 0., "Zb_x0250 Eta_x0250"}, // VIRGULE !
    {"TIC6.ga-Zb-Eta-x1500.txt", 1500., 0., "Zb_x1500 Eta_x1500"},
    {"TIC6.ga-Zb-Eta-x2750.txt", 2750., 0., "Zb_x2750 Eta_x2750"},
    {NULL} // PAS de ";" (pourquoi ? --> car c'est une liste (et pas un evenement) ! )

};

event gauging_event(i+=NB_ITER_OUTPUTS)
{
    output_gauges( gauging_stations, {h} );
}
*/

/**
Le dernier événement met fin aux itérations.
*/

event end(t = TMAX)
{
    printf("Programme termine : iteration no. %d ; temps : %g\n", i, t);
    printf("acceleration G = %g\n", G);
    printf("nombre de cellules N = %d\n", N);
    printf("hauteur normale correspondant au debit de base h_base = %g\n", h_base);
}

/**
CONCLUSION : 
Le débit est injecté correctement. Cependant, les conditions aval ne se comportent pas de la manière attendue. On s'attend à un retour à un régime stationnaire uniforme. Au lieu de cela,  on obient en régime torrentiel une courbe de remous de type S2 et en régime fluvial, une courbe M1. Il est important de trouver une explication car la zone impactée par les conditions aval remonte relatirment loin dans le canal, ce qui perturbe la solution recherchée.
*/

