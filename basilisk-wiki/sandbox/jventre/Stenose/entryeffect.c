/**

# Resolution de Poiseuille dans un cylindre

Equation à résoudre :

$$0 = -\frac{\partial p}{\partial x} + \frac{1}{r} \frac{\partial }{\partial r}  r \frac{\partial u}{\partial r}$$

Ce code résout Poiseuille avec la méthode chinoise (deux phases) pour un Reynolds de 100. 
On impose le gradient de pression qui va pour obtenir un profil de Poiseuille avec un débit 1/2

*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include <sys/stat.h>

int Re = 100;  // Reynolds
double dR = 0; // degrée d'obstruction
double xst = 2.5 ; // position de la sténose
scalar un[];

int main() {
  
  TOLERANCE = 1e-6;
  L0 = 10.;        
  
  // two-phase gives f=1 to 1 and f=0 to phase 2  
   rho1 = 1, rho2 = 1.;
   mu1 = 1./Re ;
   mu2 = 0.;
 
  N = 256;
    run();  
}

/**

## Conditions aux limites
 
 - vitesse nulle dans la phase du haut 
 - neumann sur la vitesse en entrée et sortie du tube
 - dirichlet sur la pression à l'entrée et à la sortie (8*L/Re en entrée, 0 en sortie pour avoir un gradient = 8/Re)
*/
u.t[top] = dirichlet(0);
u.n[top] = dirichlet(0);

u.n[left] = ( y < 1 ? neumann(0.):dirichlet(0.));
u.n[right] =  ( y < 1 ? neumann(0.) : dirichlet(0.) ) ;

p[left] = dirichlet(8*L0/Re);

p[right]   = dirichlet(0.);

pf[left]   = neumann(0.);
pf[right]  = dirichlet(0.);

/**
## Définition de l'interface

*/
event init (t = 0) {
 
  // géométrie de l'interface
  fraction (f, (1. - dR *exp(-(x-xst)*(x-xst)))- y  );
    
  // on sauve l'interface 
  FILE * fp2 = fopen ("interface.csv", "a");
  output_facets (f, fp2);
  fprintf(fp2,"\n\n");
  fclose (fp2);

  // initialisation de un pour la convergence
  foreach()
    un[] = u.x[];
}

/**
## Test de la convergence sur la vitesse
*/
event conv (t += 1) {
    double du = change (u.x, un);
    // fprintf(stdout,"t= %g %g \n",t, du);
    if (i > 0 && du < 1.0e-5)
        return 1; /* stop */
}

event bound (i++,last){

  // on doit redéfinir à chaque itération l'interface
  fraction (f, (1. - dR *exp(-(x-xst)*(x-xst)))- y  );

  // a chaque itération on impose la vitesse nulle en haut de l'interface 
  foreach() {
    if (f[]>1e-6 && f[]<1.-1e-6){
     u.x[]=u.y[]=p[]=0.0;
    } else {
    u.x[] = u.x[]*(f[]);
     u.y[] = u.y[]*(f[]);
     p[] = p[]*(f[]);
   }
  }

}

/**
## Enregistrement des vitesses/pressions finales 
*/
event boucle(t+=1;t<=200){
  fprintf(stdout, " t = %g \n", t);
}


event sauve(t=end)
{
  char filename[80];
  sprintf (filename, "p_final_Re_%d.csv", Re);
  FILE * fp = fopen (filename, "a");
  output_field ({p}, fp, linear = true);
  fprintf(fp,"\n\n");
  fclose (fp);

  char filename2[80];
  sprintf (filename2, "u_final_Re_%d.csv", Re);
  FILE * fp2 = fopen (filename2, "a");
  output_field ({u.x}, fp2, linear = true);
  fprintf(fp2,"\n\n");
  fclose (fp2);

  // pour avoir la pression dans bview 
  scalar press[];
  foreach(){
    press[]=p[];
  }
  dump (file = "snapshot");

}
/**

# Run

pour lancer 

~~~bash

make entreepoiseuille.tst
make entreepoiseuille/plots
make entreepoiseuille.c.html

~~~


# Results
 
## Exemples 

~~~gnuplot  pressure

stats 'p_final_Re_100.csv'  u 1

 reset
 set xlabel "x"
 set ylabel "p"
 plot 'p_final_Re_100.csv' i (STATS_blocks-2) us 1:($2==0?$3:NaN)
 
~~~

~~~gnuplot  velocity

stats 'u_final_Re_100.csv' u 1

plot 'u_final_Re_100.csv' i STATS_blocks-2 us ($2<=1.01?$2:NaN):($1==0?$3:NaN) w linespoints title 'x=0'
replot 'u_final_Re_100.csv' i STATS_blocks-2 us ($2<=1.01?$2:NaN):($1==9.99999?$3:NaN) w linespoints title 'x=L'

set xlabel 'y'
set ylabel 'u'

rep

~~~



# bibliography


  

*/

