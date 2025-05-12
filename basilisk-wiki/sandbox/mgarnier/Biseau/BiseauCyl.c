/**
#Biseau Cylindrique
*/

#include "embed.h"
#include "poisson.h"
#include "view.h"

double W = 1.;
double R = 0.8; //Rayon du bout du biseau
double Vitesse = 1.;


static double dirichlet_homogeneous_bc (Point point, Point neighbor,
				 scalar s, void * data) {
  return 0.;
}

static double source (double x,  double y) {
  return Vitesse*4*W*(atan(y/(x-W))-M_PI)/(atan(0/W)-M_PI);
}

int main()
{
  N=128;
  size (24*W); 
  origin (-12*W, -4*(W)-R);
  init_grid (N);

  vertex scalar phi[];// Définition du domaine
    foreach_vertex() {
     double r = sqrt (sq(x-(R+W)) + sq(y));

    if(r-R <0 || (fabs(y)<R && x >R+W))
    {phi[] = -1;
    }
    else
    phi[]=1;

    }
  fractions (phi, cs, fs);  					
  cm = cs;
  fm = fs;

scalar a[], b[]; 

a[bottom] = 0;
a.boundary_homogeneous[bottom] = dirichlet_homogeneous_bc;
a[left] = 0;
a.boundary_homogeneous[left] = dirichlet_homogeneous_bc;
a[top] = 0;
a.boundary_homogeneous[top] = dirichlet_homogeneous_bc;
a[right] = (y + Delta/2< -R) ? Vitesse*(y+(4*W+R))+Delta/2 : (y + Delta/2 > R) ? Vitesse*(y+(4*W+R))+Delta/2 : 0.; 
a[embed] = 4*W*Vitesse; 
//a.boundary_homogeneous[embed]=dirichlet_homogeneous_bc;

a.third = true;


foreach(){
  if(cs[]==0)
  a[] = nodata;
}


  struct Poisson p;
  p.alpha = fs;
  p.lambda = zeroc;
  p.embed_flux = embed_flux;
  scalar res[]; 
  double maxp = residual ({a}, {b}, {res}, &p), maxf = 0.;  
  foreach()  
    if (cs[] == 1. && fabs(res[]) > maxf)  
      maxf = fabs(res[]);  
  fprintf (stderr, "maxres %d %.3g %.3g\n", N, maxf, maxp);  
  timer t = timer_start();
  mgstats s = poisson (a, b, alpha = fs, tolerance = 1e-4, minlevel = 4);
  double dt = timer_elapsed (t);
  printf ("%d %g %d %d\n", N, dt, s.i, s.nrelax);

  FILE * fp = fopen ("streamlines.dat","w");
  foreach() 
    fprintf(fp,"%g %g %g %g %g %g \n", x, y, a[]);
  fflush (fp);

dump ("dump");

}