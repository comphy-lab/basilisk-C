#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h" 
#include "tension.h" 
#include "reduced.h"
#include "navier-stokes/perfs.h"

double h;
double h_2;
double U0;
double R_d;
double grav;

#define MAXLEVEL 9

FILE * fpmax; 

int main() {

  R_d=0.005; //Rayon du jet initial
  L0=0.4; //Taille de la boite  
  rho2 = 10.;
  rho1=1000.;
  mu1 = 0.1;
  mu2 = 0.01*mu1;
  U0=0.2;
  h=0.12;
  h_2=0.05;
  grav=9.81;

  TOLERANCE = 1e-3 [*];

  u.n[bottom] = dirichlet ((x > -R_d && x <R_d)*U0*(1-(x/R_d)*(x/R_d)));
  u.t[bottom] = dirichlet(0.);
  f[bottom] = !((x > -R_d && x <R_d));

  u.n[top] = u.n[] > 0. ? neumann(0) : dirichlet(0);
  p[top] = dirichlet(0.);
  pf[top] = dirichlet(0.);

  u.n[left] = dirichlet(0.);
  u.n[right] = dirichlet(0.);

  u.n[embed]=dirichlet(0.);
  u.t[embed]=dirichlet(0.);
 
  G.y = -grav;
 
  N=256;
  origin (-L0/2, 0);//set the origin
  init_grid(N);//Maillage, doit être de la forme 2^n

  fpmax =  fopen("log.dat", "w"); //crée un fichier où on va mettre les infos qu'on veut

  f.sigma=0.2;

  run();
}

event init (t = 0) {
  fraction (f, union(union(intersection((y<h+h_2),(y>h_2)),intersection((x<-R_d),(y<h_2))),intersection((x>R_d),(y<h_2))));
  solid(cs,fs, !union(intersection((x<-R_d),(y<h_2)),intersection((x>R_d),(y<h_2))));
}

event logfile (i++) {
  fprintf (stderr, "%d %g \n", i, t);
  fprintf (fpmax, "%d %g \n", i, t);
}


event adapt (i++){
  adapt_wavelet((scalar*){u}, (double[]){0.05, 0.05}, MAXLEVEL,5);
}

/**
We generate videos:
*/
event ppm_output (t = 0; t += 0.05; t <= 20) {
  char name[80];
  sprintf (name, "f.mp4");
  output_ppm (f, file = name, n = 512, min = 0, max = 1, linear = true);
  
  char name1[80];
  sprintf (name1, "uY.mp4");
  output_ppm (u.y, file = name1, n = 512, linear = true);
  
  /**
  char name2[80];
  sprintf (name2, "uX.mp4");
  output_ppm (u.x, file = name2, n = 512, linear = true);
  */
  
  scalar omega[];
  vorticity (u, omega);
  
  char name3[80];
  sprintf (name3, "omega.mp4");
  output_ppm (omega, file = name3, n = 512, linear = true);
}

/**
![Free-surface](bubble_injection/f.mp4)
![Vertical velocity](bubble_injection/uY.mp4)
![Horizontal velocity](bubble_injection/uX.mp4)
![Vorticity](bubble_injection/omega.mp4)
*/