/**
# Rayleigh-Plateau instability example of tracer
modify from output_facet
*/
//multigrid
#include "grid/multigrid.h"
#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "vtknew.h"
#include "view.h"

#define epsilon 0.05 //amplitude de la perturbation initial
#define LL 9. //wave length/RR > 2pi can converge
#define Oh 0.01
//#define LEV 6
#define tend 1.
FILE *fp;
FILE *fp2;
int main() {
  fp=fopen("energy","w");
  fp2=fopen("post","w");
  periodic(right);
  rho1 = 1., mu1 = Oh;
  rho2 = 1., mu2 = Oh;
  f.sigma=1.;  
  N = 2 << 6;
  run(); 
}




scalar XX[];
scalar YY[];

trace
  static double surface_ttracer (struct OutputFacets p)
{
  scalar c = p.c;
  face vector s = p.s;
  if (!s.x.i) s.x.i = -1;
  double err =1.;
  double yyy =0.;
  foreach()
    if (c[] > 1e-6 && c[] < 1. - 1e-6) {
    coord n = facet_normal (point, c, s);
    double alpha = plane_alpha (c[], n);
    coord segment[2];

    if (facets (n, alpha, segment) == 2){
      XX[] = x + segment[0].x*Delta;
      YY[] = y + segment[0].y*Delta; 
      if (fabs(XX[] - L0/2.) < err){
        yyy = YY[];
        err = fabs(XX[] - L0/2.);
      }
    }
  }
  return yyy;
}



/**
## init event */


event init (t = 0) {
  fraction (f,  (1.+ epsilon*cos(2.*M_PI*x))/LL -y );
  output_ppm (f, file = "finit.png", linear = true, n=512, min = 0, max = 1);

  view(fov=10.,tx=-0.5,ty=-0.2,width=1280,height=640);
  box();
  draw_vof("f", lw = 1.5); 
  cells();
  save ("fim.png");

  scalar pos[];
  position (f, pos, {0,1});
  double max = statsf(pos).max;
  double hhh = surface_ttracer(f);
  fprintf(stderr,"%g %g %g\n", t, max,hhh);
}

event test (t = 0.; t <= tend; t += tend/100.) {
  scalar pos[];
  position (f, pos, {0,1});
  double max = statsf(pos).max;
  double hhh = surface_ttracer(f);
  fprintf(stderr,"%g %g %g %g\n", t, max,hhh,fabs(-0.15));
}


/**
event outvtk (t = 0.; t <= tend; t += tend/100.) {
        char name[80];
        sprintf (name, "snapshot-%04g.vtk", t*1000);
        FILE * fpvtk = fopen(name, "w");
		    scalar omega[];
        vorticity (u, omega);
        output_vtk ({f, p, u.x, u.y, omega}, fpvtk);
     		fclose (fpvtk);
}
*/