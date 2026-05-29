#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h" 
#include "tension.h" 
#include "reduced.h"
#include "navier-stokes/perfs.h"
#include "tracer.h"
#include "tag.h"

double U0;
double U1;
double W;
double H;
double Do;
double Di;
double Li;
double h;
double Lo;
double Re;

scalar s[];
scalar * tracers = {s};

FILE * fpmax; //

int main() {

  Re=120;
  
  Di=1.;
  H=5*Di;
  W=25*Di;
  Do=Di/2;
  Li=2*Di;
  Lo=2*Do;
  
  L0=W+2*Lo; 
  rho2 = 1.3;
  rho1=1000.;
  U0=7;
  U1=U0*Di/(2*Do);
  mu1 = rho1*U0*Di/Re;
  mu2 = 0.01*mu1;

  TOLERANCE = 1e-3 [*];
  
  u.n[embed] = dirichlet(0.);
  u.t[embed] = dirichlet(0.);

  u.n[bottom] = dirichlet((x > -Di/2 && x <Di/2)*U0*(1-(2*x/Di)*(2*x/Di)));
  u.t[bottom] = dirichlet(0.);
  s[bottom] = dirichlet (U0*(x > -Di/2 && x <Di/2));

  u.n[left] = (y > Li && y < Li+Do) ? dirichlet(-U1*(1-(2*(y-Li-Do/2)/Do)*(2*(y-Li-Do/2)/Do))) : dirichlet(0.);
  u.n[right] = (y > Li && y < Li+Do) ? dirichlet(U1*(1-(2*(y-Li-Do/2)/Do)*(2*(y-Li-Do/2)/Do))) : dirichlet(0.);
  
  u.t[left] = (y > Li && y < Li+Do) ? neumann(0.) : dirichlet(0.);
  u.t[right] = (y > Li && y < Li+Do) ? neumann(0.) : dirichlet(0.);
  
  u.n[top] = u.n[] > 0. ? neumann(0) : dirichlet(0);
  p[top] = dirichlet(0.);
  pf[top] = dirichlet(0.);
 
  G.y = -9.81;

  N=256;
  origin (-L0/2, 0);
  init_grid(N);
  
  fpmax =  fopen("log.dat", "w");
  
  f.sigma=0.5;

  run();
}

event init (t = 0) {
  fraction(f, y<Li+H);
  solid (cs, fs, !union(union(intersection((x<-Di/2),(y<Li)),intersection((x>Di/2),(y<Li))),union(intersection((x<-W/2),(y>(Li+Do))),intersection((x>W/2),(y>(Li+Do))))));
}

event logfile (i++) { 
  fprintf (stderr, "%d %g \n", i, t); 
  fprintf (fpmax, "%d %g \n", i, t);
}


event remove_droplets (i++) {
  remove_droplets (f, threshold=0.01, bubbles=true);
}

/*
int isave1 = 1;
event res_save (t += 0.5; t <= 10) {
  char name[80];
  
  sprintf (name, "interface-%d.txt", isave1);
  FILE * fpfacet = fopen(name, "w");
  output_facets (f, fpfacet);
  fclose(fpfacet);
  
  isave1++;
}
**/

event ppm_output (t = 0; t += 0.05; t <= 70) {
  char name[80];
  sprintf (name, "f.mp4");
  output_ppm (f, file = name, n = 512, min = 0, max = 1, linear = true);
  
  char name0[80];
  sprintf (name0, "uX.mp4");
  output_ppm (u.x, file = name0, n = 512, min = -U0, max = U0, linear = true);
  
  char name1[80];
  sprintf (name1, "uY.mp4");
  output_ppm (u.y, file = name1, n = 512, min = -U0, max = +U0, linear = true);

  char name2[80];
  sprintf (name2, "s.mp4");
  output_ppm (s, file = name2, n = 512, min = 0., max = U0, linear = true);

  char name3[80];
  sprintf (name3, "p.mp4");
  output_ppm (p, file = name3, n = 512, linear = true);
  
  char name4[80];
  sprintf (name4, "cs.mp4");
  output_ppm (cs, file = name4, min=0, max=1, n = 512, linear = true);
}

event profile (t = end) {
  printf ("-----END-----\n");
}

/**
## Movies

![Free surface motion](for_G2L/f.mp4)
![Passive tracer](for_G2L/s.mp4)
![Vertical velocity](for_G2L/uY.mp4)
![Horizontal velocity](for_G2L/uX.mp4)
![Pressure](for_G2L/p.mp4)
![Geometry](for_G2L/cs.mp4)
*/