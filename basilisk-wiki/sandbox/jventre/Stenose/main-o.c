/**
   Ceci est le code de Jeanne */

#include "axi.h"
#include "navier-stokes/centered.h"

#include "two-phase.h"

#define period(a) ((a) - floor(a))  // macro pour prendre que les decimales



int main() {
 
  TOLERANCE = 1e-6;
  L0 = 10.;         // la boite originale est 10x10
  double Re = 50.;  // Reynolds

  
  // two-phase gives f=1 to 1 and f=0 to phase 2  
   rho1 = 1, rho2 = 1.;
   mu1 = 1./Re ;
   mu2 = 0.;
 
  // for (N = 8; N <= 64; N *= 2)

  // Nombre de oints par coté  delta x = L/N
  N = 128;
    run();  
}

u.t[top] = dirichlet(0);
//u.n[left] = ( y < 1 ? (period(t) < 0.4 ? 1.*(sq(1.)-sq(y))*sin(period(t)/0.4*pi) : 0) : 0.);
u.n[left] = ( y < 1 ? dirichlet(1.) : dirichlet(0.));
p[left]    = neumann(0.);
pf[left]   = neumann(0.);


u.n[right] =  ( y < 1 ? neumann(0.) : dirichlet(0.) ) ;
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

event init (t = 0) {
  
  // refine before do fraction!
  refine ( y < 2. & level < 6 );

    fraction (f,  (1. - 0.3*exp(-(x-2.5)*(x-2.5))) -  y );

}

// at eahc iteratio do zero everywhere in phase 2
event bound (i++,last){
  
fraction (f,  (1. - 0.3*exp(-(x-2.5)*(x-2.5))) -  y );
 
  foreach() {
    if (f[] > 1e-6 && f[] < 1.-1e-6){
      u.x[] = 0.0;
      u.y[] = 0.0 ;
      p[] = 0.0 ;
    } else {
    u.x[] = u.x[]*f[];
    u.y[] = u.y[]*f[];
    p[] = p[]*f[];
    }
  }
}

event movies (t += 0.1; t <= 1.) {
  scalar omega[];

  /* vorticity (u, omega); */

  // static FILE * fp = popen ("/Users/jmf/sources/basilisk_0/src/ppm2mp4 > ux.mp4", "w");
  //  static FILE * fp = fopen ("level.ppm", "w");
  //output_ppm (u.x, file = "name.png", box = {{0,0},{10,1.2}});
     
  /*  static FILE * fp = popen ("ppm2mpeg > vort.mpg", "w"); */
  
  /*  output_ppm (omega, file = name, box = {{0,0},{10,1.2}} , linear = true); */


} 

event sauve(t += 1; t <= 100)
{

  // int res;
  // char dir[80];
  // sprintf(dir, "Resultats/Pression/dR=%f/xst=%f/Re=%f", dR,xst,Re)
  // res = mkdir(dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  char filename[80];
  sprintf (filename, "p.csv");
  // remove(filename);
  FILE * fp = fopen (filename, "a");
  output_field ({p}, fp, linear = true);
  fprintf(fp,"\n\n");
  fclose (fp);

  // printf("t = %g\n",t );

  // FILE * fp3 = fopen ("Resultats/fraction.txt", "a");
  // output_field ({f}, fp3, linear = true);
  // fprintf(fp3,"\n\n");
  // fclose (fp3);

}

event print (t=end) {

  dump (file = "snapshot");

}





