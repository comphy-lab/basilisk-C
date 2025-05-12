//CC='mpicc -D_MPI=1' make edgeTone.tst 
// CC99='mpicc -std=c99' qcc -Wall -O2 -D_MPI=1 edgeTone.c -o edgeTone -lm
/** 
## Edge tone

 */
#include "embed.h"
#include "navier-stokes/centered.h"
#include "utils.h"

double Reynolds = 375;
int maxlevel = 12;
double b = 1.; // jet diameter
double U0 = 1.; // for top hat
double Tad=1.;
double alpha_b = 5*M_PI/180.;
double W=10.; // edge-distance
double L=2.; // edge faces length
char name[80];

scalar omega[];
face vector muv[];
scalar s[];


double sech2(double y)
{
  return 1./(sq(cosh(y)));
}


int main(int argc, char * argv[])
{
  L0 = 100. [0]; 
  origin (0., -L0/2.);
  N = 4096;
  mu = muv;
  init_grid (N);					
  run();
}

event init (i = 0)
{
  vertex scalar phi[];  
  foreach_vertex() {  
    phi[]=1.;
    if (x > W && x < 15.)
      phi[] = -0.1+fabs(y);
  }
  fractions (phi, cs, fs);
  fractions_cleanup (cs, fs);
  foreach()
    u.x[]=0.;
}


event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]*(b*U0)/Reynolds; // top hat. Remplacé U0 par 1. pour augmentation graduelle de la vitesse via U0. Reference à Re=100
  vorticity (u, omega);  
}

u.n[left]  = dirichlet(U0*sech2(y/(b*0.25)));// sech2 - attention, demi largeur pour Bickley. y=0.5 avec b=1 => y/(0.25*b)=2 et donc U=0
u.t[left]  = dirichlet(0.) ; 

p[left]    = neumann(0.);
pf[left]   = neumann(0.);

u.n[right] = neumann(0.);
u.t[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

u.n[top] = neumann(0.);
u.t[top] = neumann(0.);
p[top]   = dirichlet(0.);
pf[top]  = dirichlet(0.);

u.n[bottom] = neumann(0.);
u.t[bottom] = neumann(0.);
p[bottom]   = dirichlet(0.);
pf[bottom]  = dirichlet(0.);

u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);

event movie(t += 1){
  static FILE * fp3 = fopen ("movie.ppm", "w");
  output_ppm(omega, fp3, n = 256, min = -1, max = 1, box = {{0, -W/2.},{2*W, W/2.}});
}

event adapt (i++) { 
  adapt_wavelet ({cs,u}, (double[]){1e-3,1e-2,1e-2}, 10);
} 

event logfile (i++, t<=2000)
{    
    fprintf (stderr, "%d %g %g \n",
	     i, t, dt);
}
