/**
Introduction :
La simulation qui suit reproduit un "canal de laboratoire" : il s'agit d'un modèle 1D avec une topographie à pejnte constante. 
Mon but est de tester plusieurs manières d'implémenter une injection de débit en amont et de vérifier si la réaction à l'intérieur du canal est conforme à des situations de référence.
Comme premier essai, je tente de reproduire un modèle simple vu à l'université, avec un débit imposé en amont et des conditions de bord libres en aval. Je teste deux cas : régime stationnaire uniforme (débit constant et hauteur normale) ou hydrogramme triangulaire en entrée.
*/

/**
Je commence par importer les différents packages nécessaires. Le modèle est en une dimension. Le terme source de frottements est celui de G. Kirstetter (http://basilisk.fr/sandbox/geoffroy/sourceterm/manning.h). Le solveur est celui du système de Saiut-Venant.
*/

#include "grid/cartesian1D.h"
//#include "grid/bitree.h"
#include "saint-venant.h"
#include "manning_n_cst.h"

/**
Trois paramètres sont facilement modifiables pour tester plusieurs configurations : la durée de simulation, le déboit unitaire de base et lapente du canal.
*/

#define TMAX 13000. // duree de la simulation [s]
#define QBASE 1.65 // debit de base unitaire [m^2/s]
#define S0 0.001 // oppose de la pente du canal [m/m]

/** Le coefficient de Manning sert à calculer Sf ET a calculer h_base (hauteur normale pour le debit de base).
Si on ne veut pas tenir compte des frottements :
 + commenter " #include "SOURCE/manning.h" " ;
 + commenter " n  = NMANNING; " ;
 + LAISSER " h_base = pow( NMANNING*NMANNING* QBASE*QBASE / S0 , 3./10. ); "
  --> car on a besoin de h_base par la suite, comme hauteur correspondant au debit de base.
*/

#define NMANNING 0.02 // coefficient de Manning [S.I.]

/**Les vriables globales sont la hauteur normale correspondant au débit de base, et des champs scalaires donnant des informations sur la nature de l'écoulement.
*/

double h_base;

scalar Ec[]; // energie cinetique de l'ecoulement
scalar q[]; // debit unitaire
scalar hc[]; // hauteur critique
scalar hn[]; // hauteur normale

/**
On peut faire varier la fréquence des sorties (graphiques ou textuelles). Une macro simple permet de calculer maximum.
*/

#define NB_ITER_OUTPUTS 52 // nombre d'iterations entre chaque output pour les tableaux

/**
## Main
La fonction principale regroupe des constantes géométriques et physiques.
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

zb[left] = neumann(+S0);
zb[right] = neumann(-S0);
h[left] = neumann(0);
h[right] = neumann(0);
u.n[right] = neumann(0);

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
##Initialisation 
Les conditions aux limites sont un débit imoposé en amont et une hauteur d'eau et une vitesse libres en aval, quel que soit le régime d'écoulement.
Les conditions initiales correspondent à un régime stationnaire uniforme.
*/

event init(t=0)
{

    // Conditions aux limites

    zb[left] = neumann(+S0);
    zb[right] = neumann(-S0);

    h[left] = neumann(0);
    h[right] = neumann(0);

    u.n[left] =  q_injection(t)/h[left];
    u.n[right] = neumann(0);


    h_base = pow( NMANNING*NMANNING * QBASE*QBASE / S0 , 3./10. );
    //h_base = 1.026;

    // Conditions initiales 
  
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

    h[left] = neumann(0);
    h[right] = neumann(0);

    u.n[left] =  q_injection(t)/h[left];
    u.n[right] = neumann(0);
    boundary(all);
  
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
## Sorties texte
*/

event outputs_tab(t+=60.)
{
  FILE * fp;
  char name[80];

  sprintf(name,"profile-%d.dat",(int)t);
  fp = fopen(name,"w");

  foreach()
    fprintf (fp, "%g %g %g %g %g\n",x,zb[],h[],h[]+zb[],u.x[]);
  
  Point point = locate(L0-L0/N/2.);
  fprintf (fp, "%g %g %g %g %g\n",L0+L0/N/2.,zb[1],h[1],h[1]+zb[1],u.x[1]);
 
  fclose(fp);
}

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
~~~gnuplot
   reset
   set xlabel "x (m)"
   set ylabel "elevation (m)"
   ydim = 800 
   xdim = 480
   plot 'profile-0.dat' u 1:2 w lines lc rgb "black" lw 2 title 'zb', \
        'profile-0.dat' u 1:4 w lines lc rgb "blue" title 'eta(t=0)', \
        'profile-60.dat' u 1:4 w lines lc rgb "red" title 'eta(t=60s)', \
        'profile-120.dat' u 1:4 w lines lc rgb "green" title 'eta(t=120s)', \
        'profile-3600.dat' u 1:4 w lines lc rgb "yellow" title 'eta(t=3600s)'
~~~
~~~gnuplot
   reset
   set xlabel "x (m)"
   set ylabel "velocity (m/s)"
   set key left bottom
   ydim = 800 
   xdim = 480
   plot 'profile-0.dat' u 1:5 w lines lc rgb "blue" title 'u(t=0)', \
        'profile-60.dat' u 1:5 w lines lc rgb "red" title 'u(t=60s)', \
        'profile-120.dat' u 1:5 w lines lc rgb "green" title 'u(t=120s)', \
        'profile-3600.dat' u 1:5 w lines lc rgb "yellow" title 'u(t=3600s)'
~~~
*/



