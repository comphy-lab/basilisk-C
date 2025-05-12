#include "spherical.h"
//#include "saint-venant.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
#include "layered/perfs.h"
#include "terrain.h"

#define nlayers 2

#define MGD mgp.i
#define MAXLEVEL 8
#define MINLEVEL 5
#define ETAE     2e-2 // error on free surface elevation (1 cm)
#define ETAMAXE  5e-2 // error on maximum free surface elevation (5 cm)


scalar h;
vector u;
scalar topo[];
scalar etamax[];

int main()
{
  // Earth radius in metres
  Radius = 6371220.;
  // the domain is 10 degrees squared
  size (4.);
  // centered on... longitude,latitude
  origin (-4. - L0/2., 46.5 - L0/2.);
  // acceleration of gravity in m/min^2
  G = 9.81*sq(60.);
  init_grid (1 << MAXLEVEL);
  nl = nlayers;
  //DT = 1e-5;

  CFL_H = 0.2;
  TOLERANCE = 1e-4 [*];
  
  run();
}

int my_adapt() {
  scalar eta[];
  foreach(){
    double H = 0.;
    foreach_layer()
      H += h[];
    eta[] = H > dry ? zb[] + H : 0.;}
  astats s = adapt_wavelet ({eta}, (double[]){ETAE},
			    MAXLEVEL, MINLEVEL);
  fprintf (stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
  return s.nf;
}

/**
## Initial conditions
*/

event init (i = 0)
{
  terrain (zb, "/work/Simu_basilisk/terrain/etopo2", NULL);
  conserve_elevation();
  
  scalar hpert[];
  
  foreach () {
    hpert[] = 20.*exp(-100.*((x+4)*(x+4)+(y-46.5)*(y-46.5))); //initial gaussienne perturbation
    
    foreach_layer()
    h[] = max(0., hpert[] - zb[])/nl;
  }

  u.n[left]   = - radiation(0);
  u.n[right]  = + radiation(0);
  u.n[bottom] = - radiation(0);
  u.n[top]    = + radiation(0);

}

/**
## Outputs */

/*event logfile (i++) {
  scalar H[];
  foreach() {
    H[] = 0.;
    foreach_layer()
      H[] += h[];
  }
  stats s = statsf (H);
  norm n = normf (u.x);
  if (i == 0)
    fprintf (stderr,
	     "t i h.min h.max h.sum u.x.rms u.x.max dt mgD.i speed tn\n");
  fprintf (stderr, "%g %d %g %g %g %g %g %g %d %g %ld\n",
	   t, i, s.min, s.max, s.sum, n.rms, n.max, dt, MGD,
	   perf.speed, grid->tn);
}*/

/**
## Initial conditions */
event pictures1 (t = 0) {
scalar H[];
  foreach() {
    H[] = 0.;
    foreach_layer()
      H[] += h[];
      eta[] = H[] > dry ? zb[] + H[] : 0.;
  }
  output_ppm (H, file = "h_initial.png",min=0, max=100, n = 2048, linear = true);
  output_ppm (eta, file = "eta_initial.png",min=0, max=20, n = 2048, linear = true);
  output_ppm (zb, file = "new_topo.png",min=-500., max=500, n = 2048, linear = true);
}

/**
## Free-surface and level of refinement */

event figures (t = 0; t <= 30; t += 2)
{
  scalar m[], etam[];
  foreach() {
  double H = 0.;
    foreach_layer()
      H += h[];
    m[] = -zb[];
    etam[] = H < dry ? 0. : (zb[] > 0. ? eta[] - zb[] : eta[]);
  }
  //boundary({m, etam})
  char name[80];
  sprintf (name, "eta-%g.png", t);
  output_ppm (etam, mask = m, min = -2, max = 2, file = name, n = 1024,
	      linear = true, box = {{-6,44.5},{-2,48.5}},
	      opt = "-fill white -opaque black");

  /*sprintf (name, "level-%g.png", t);
  scalar l[];
  foreach()
    l[] = level;
  output_ppm (l, min = MINLEVEL, max = MAXLEVEL, file = name, n = 1024,
	      linear = false, box = {{-22,24},{-14,32}});*/


}

// adaptation
//event adapt (i++) my_adapt();