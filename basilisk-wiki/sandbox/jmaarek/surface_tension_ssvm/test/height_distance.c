#include  "grid/quadtree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "height_distance.h"

#include <vofi.h>
#pragma autolink -L$HOME/local/lib -lvofi

static double xc=0, yc=0, zc=0, Rc=0.7;

static double sphere (creal p[dimension])
{
#if dimension == 3
  return (sq(p[0] - xc) + sq(p[1] - yc) + sq(p[2] - zc)) - sq(Rc);
#else
  return (sq(p[0] - xc) + sq(p[1] - yc)) - sq(Rc);
#endif
}

static void vofi (scalar c, int levelmax)
{
  double fh = Get_fh (sphere, NULL, 1./(1 << levelmax), dimension, 0);
  foreach() {
    creal p[3] = {x - Delta/2., y - Delta/2., z - Delta/2.};
    c[] = Get_cc (sphere, p, Delta, fh, dimension);
  }
}


int LEVEL;
int LEVEL_MAX = (dimension < 3 ? 11 : 8);

int main(){
 size(2);
 origin(-L0/2, -L0/2, -L0/2);

 for (LEVEL = 5; LEVEL <= LEVEL_MAX;  LEVEL++){
 init_grid(1 << LEVEL);
 run();
 }
  return 1;
}

scalar phi[];
vertex scalar phi2[];
scalar f2[];
vector n[];

event init(t = 0){
   //fraction (f, 0.7 - sqrt(sq(x) + sq(y) + sq(z)));
   vofi (f, LEVEL);
   boundary({f});
   height_distance(f, phi, imax = 16);
   #if dimension == 2
   foreach_vertex()
      phi2[] = (phi[] + phi[-1] + phi[0,-1] + phi[-1,-1])/4.;
      fractions(phi2, f2);
   #else
     scalar ff = phi;
     foreach_vertex()
      phi2[] = (ff[] + ff[-1] + ff[0,-1] + ff[-1,-1] +
            ff[0,0,-1] + ff[-1,0,-1] + ff[0,-1,-1] + ff[-1,-1,-1])/8.;

     fractions(phi2, f2);
   #endif

   double e_phi = 0., e_phi_max = 0., vol = 0.;
foreach(reduction(max:e_phi_max) reduction(+:e_phi) reduction(+:vol)) {
  double phi_exact = 0.7 - sqrt(sq(x) + sq(y) + sq(z));       // 2D; add sq(z) for 3D
  double band = 2.*Delta;
  if (fabs(phi_exact) < band) {
    double err = fabs(phi[] - phi_exact);
    e_phi     += sq(err);   // L2 integral
    e_phi_max  = max(e_phi_max, err);   // Linf
    vol       += 1;
  }
}
e_phi = sqrt(e_phi / vol);

double e_n = 0., e_n_max = 0.;
  vol = 0;
  foreach(reduction(max:e_n_max) reduction(+:e_n) reduction(+:vol)) {
  double r = sqrt(sq(x) + sq(y) + sq(z));
  double phi_exact = 0.7-r;
  double band = 2.*Delta;
  if (fabs(phi_exact) < band && r > 1e-10) {
    // Numerical normal from gradient of d
    coord n = (coord){0,0,0};
    double nn = 0;
    foreach_dimension(){
      //n.x = eno2_gradient_x(point, phi);
        n.x = (phi[1] - phi[-1])/2/Delta;
      nn += sq(n.x);
    }
    nn = sqrt(nn);
    foreach_dimension(){
      n.x /= nn;
    }
    // Analytical normal
    coord n_ex = (coord){-x/r, -y/r, -z/r};
    double err = 0;
    foreach_dimension()
        err += sq(n_ex.x - n.x);
    err = sqrt(err);
    e_n     += sq(err);
    e_n_max  = max(e_n_max, err);    vol     += 1;
  }
}
  e_n = sqrt(e_n / vol);

 fprintf(ferr, "%d %g %g %g %g\n", LEVEL, e_phi,  e_phi_max, e_n,  e_n_max);

}
    