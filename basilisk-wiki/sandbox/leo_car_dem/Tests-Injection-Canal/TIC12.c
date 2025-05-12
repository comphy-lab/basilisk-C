//25.03.2021
//Stage Injection-Basilisk
//TIC12.c (TIC = Test-Injection-Canal)
//But : reprise de TIC6.c (la version avec manning_n_cst.h). Je veux adapter ce script a hydro.h, pour voir si ça corrige le "bug M1".
//Je veux reproduire un regime stationnaire uniforme.

#include <math.h>

#include "grid/multigrid1D.h"
//#include "grid/cartesian1D.h"
//#include "saint-venant.h" // obsolete
#include "layered/hydro.h"
#include "layered/implicit.h"
#include "manning_n_cst.h"

#define NMANNING 0.02 // coefficient de manning [s/m^(1/3)]
#define TMAX 13000. // duree de la simulation [s]
#define QBASE 1.65 // debit de base unitaire [m^2/s]
#define S0 0.001 // oppose de la pente du canal [m/m]

double h_base;

double f = 0.05; // adapte de "lake flowing into itself"

#define NB_ITER_OUTPUTS 52 // nombre d'iterations entre chaque output pour les tableaux

#define MAX(a,b) (a > b ? a : b) // function "maximum"

int main()
{
    nl = 1; // nombre de couches
    G  = 9.81; // acceleration de pesanteur [m/s^2]
    L0 = 3000.; // longueur du canal/bief [m]
    X0 = 0.; // abscisse de debut du canal [m]
    N  = 256; // nombre de cellules sur le domaine
    run();
}


double q_injection(double t) // debit injecte en amont au temps t
{
    /*
    double q_max = 10; // debit de pointe unitaire atteint par l'hydrogramme [m^2/s] (Qmax=50 [m3/s] dans le projet de fin de semestre)
    double t_0 = 0.; // temps ou commence la crue [s]
    double t_pic = 20.*60.; // temps ou le debit est maximum [s]
    double t_decrue = 60.*60.; // temps ou la decrue est terminee [s]
    double pente_crue = (q_max - QBASE) / (t_pic - t_0); // gradient de debit pendant la crue
    double pente_decrue = (QBASE - q_max) / (t_decrue - t_pic); // gradient de debit pendant la decrue
    */

    //return ( QBASE + MAX( 0, ( t < t_pic ? pente_crue * (t - t_0) : pente_decrue * (t - t_pic) + q_max - QBASE) ) ); // debit au temps t
    return(QBASE);
}


event init(t=0)
{
    theta_H = 0.55; // adapte de "lake flowing into itself"
  
    foreach()
    {
        nmanning[] = NMANNING; // coefficient de Manning [S.I.]
    }

    // Conditions initiales

    zb[left] = neumann(+S0);
    zb[right] = neumann(-S0);
    eta[left] = neumann(+S0);
    eta[right] = neumann(-S0);
  
    h[left] = neumann(0);
    h[right] = neumann(0);

    u.n[left] =  q_injection(t)/h[left];
    u.n[right] = neumann(0);


    h_base = pow( NMANNING*NMANNING * QBASE*QBASE / S0 , 3./10. );

    foreach()
    {
        zb[] = -S0 * x;
        h[]  = h_base;
        u.x[] = QBASE/h_base;
    }

    boundary(all);

}

event source (i++) // adapte de "lake flowing into itself"
{
  foreach()
    u.x[] = (u.x[] + G*S0*dt)/(1. + dt*f/8.*u.x[]/sq(h[]));
  boundary ((scalar *){u});
}



event outputs_movies(NB_ITER_OUTPUTS)
{
    output_ppm(zb, file = "TIC6.vue-dessus-zb.mp4");
    output_ppm(h, file = "TIC6.vue-dessus-h.mp4");
    output_ppm(eta, file = "TIC6.vue-dessus-eta.mp4");
}



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


event end(t = TMAX)
{
    printf("Programme termine : iteration no. %d ; temps : %g\n", i, t);
    printf("acceleration G = %g\n", G);
    printf("nombre de cellules N = %d\n", N);
    printf("hauteur normale correspondant au debit de base h_base = %g\n", h_base);
}


/**
## Sorties texte
*/

event outputs_tab(t = TMAX)
{
  FILE * fp;
  char name[80];

  sprintf(name,"profile-TMAX.dat");
  fp = fopen(name,"w");

  foreach()
    fprintf (fp, "%g %g %g %g %g\n",x,zb[],h[],h[]+zb[],u.x[]);
  
  fclose(fp);
}

/**
~~~gnuplot
   reset
   set xlabel "x (m)"
   set ylabel "elevation (m)"
   ydim = 800 
   xdim = 480
   plot 'profile-TMAX.dat' u 1:2 w lines lc rgb "black" lw 2 title 'zb', \
        'profile-TMAX.dat' u 1:4 w lines lc rgb "blue" title 'eta(t=TMAX)'
~~~
~~~gnuplot
   reset
   set xlabel "x (m)"
   set ylabel "velocity (m/s)"
   set key left bottom
   ydim = 800 
   xdim = 480
   plot 'profile-TMAX.dat' u 1:5 w lines lc rgb "blue" title 'u(t=TMAX)'
~~~
*/






