/**

# Resolution de Poiseuille dans un cylindre

Equation à résoudre :

$$0 = -\frac{\partial p}{\partial x} + \frac{1}{r} \frac{\partial }{\partial r}  r \frac{\partial u}{\partial r}$$

*/
 
#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include <sys/stat.h>

int Re ;  // Reynolds
double dR = 0.; // degrée d'obstruction
double xst = 2.5 ; // position de la sténose
scalar un[];

double stenose(double x){
    double xb,xl;
    xl=1.;
    xb=xst-.5*xl;
    //return (fabs((x-xb)/xl-.5)<.5?  (1+cos(pi+2*pi*(x-xb)/xl))/2.   : 0 );
    return  exp(-sq((x-xst)/xl));
}

int main() {
  
  TOLERANCE = 1e-6;
  L0=25;
  Re= 100;
  // two-phase gives f=1 to 1 and f=0 to phase 2  
   rho1 = 1, rho2 = 1.;
   mu1 = 1./Re ;
   mu2 = 0.;
 
  N = 512 ; 

  system("rm p_final_Re_100.csv");
  system("rm u_final_Re_100.csv");
  run();  

  // Re = 100;
  // L0 = 18;
  // system("rm p_final_Re_100.csv");
  // system("rm u_final_Re_100.csv");
  // mu1 = 1./Re ;
  // run(); 

  // Re = 200;
  // L0=25;
  // system("rm p_final_Re_200.csv");
  // system("rm u_final_Re_200.csv");
  // mu1 = 1./Re ;
  // run(); 

  // Re = 300;
  // L0=30;
  // N=1024;
  // system("rm p_final_Re_300.csv");
  // system("rm u_final_Re_300.csv");
  // mu1 = 1./Re ;
  // run(); 

  //  Re = 500;
  //  L0=35;
  //  N=1024;
  // system("rm p_final_Re_500.csv");
  // system("rm u_final_Re_500.csv");
  // mu1 = 1./Re ;
  // run(); 

}

/**

# Conditions aux limites
 
 - vitesse nulle dans la phase du haut 
 - dirichlet 1 sur la vitesse en entrée
 - dirichlet sur la pression à l'entrée et à la sortie (8*L/Re en entrée, 0 en sortie pour avoir un gradient = 8/Re)
*/
u.t[top] = dirichlet(0);
u.n[top] = dirichlet(0);

u.n[left] = ( y < 1 ? dirichlet(1.):dirichlet(0.));
p[left] = neumann(0.);
pf[left]   = neumann(0.);

u.n[right] =  ( y < 1 ? neumann(0.) : dirichlet(0.) ) ;
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

/**
# Définition de l'interface

*/
event init (t = 0) {
 
  // géométrie de l'interface
  mask (y >  1 - dR*stenose(x)? top :    none);
    
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
# Test de la convergence sur la vitesse
*/
event conv (t += 1) {
    double du = change (u.x, un);
    fprintf(stdout,"t= %g %g \n",t, du);
    if (i > 0 && du < 1.0e-4)
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
# Enregistrement des vitesses/pressions finales 
*/
event boucle(t+=1;t<=300){
  // printf(" t = %g \n", t);
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

~~~gnuplot velocity

stats 'u_final_Re_100.csv' u 1

plot 'u_final_Re_100.csv' i STATS_blocks-2 us 1:($2==0?$3:NaN) w linespoints title 'Re=100'

set xlabel 'y'
set ylabel 'u'

rep

~~~

~~~gnuplot pressure

stats 'p_final_Re_100.csv' u 1

p'p_final_Re_100.csv' u 1:($2==0?$3:NaN) title 'entry effect'
rep 8.*25/100.-8./100.*x, title "poiseuille pressure"

set xlabel 'y'
set ylabel 'p'

rep

~~~


# bibliography


  

*/

