/** We investigate the volume error caused by fraction.h of a sliced sphere (a 'fan').
note: Thanks to Antoon for making it a more clean example.

![This is the shape](fraction_fan2/fan.png)

We compare the intersection via a multiplication of volume fractions against taking the minimum of volume fractions (As in this example of [Alexis Berny](/sandbox/aberny/cgsbool.c)). 
*/
#include "grid/multigrid3D.h"
#include "utils.h"
#include "fractions.h"
#include "view.h"
/** For the volume of a sliced sphere one can visit [wikipedia](https://en.wikipedia.org/wiki/Spherical_cap) */
#define Va (4./3.*M_PI*pow(R, 3.) - 2.*M_PI*pow(R - w/2., 2.)/3.*(3*R - (R - w/2.)))

scalar fan1[], fan2[];

int main() {
  X0 = Y0 = Z0 = -L0/2.;
  init_grid(1 << 7); //a 128 x 128 x 128 grid 
  FILE * fp = fopen("data", "w");
  for(double R = 0.01; R < L0/2.; R*=1.1){
    double VolEst1 = 0., VolEst2 = 0;	// Estimated volume from fractions.h
    double w = R/5.; // R/w = 5
    scalar sph[], planeup[], planedown[];
    fraction(sph, -sq(x) - sq(y) - sq(z) + sq(R));
    fraction(planeup, -x + w/2.);
    fraction(planedown, x + w/2.);
    foreach (){
      fan1[] = sph[]*planeup[]*planedown[];
      fan2[] = min(sph[], min(planeup[], planedown[]));
    }
    foreach (reduction(+:VolEst1) reduction(+:VolEst2)){ 
      VolEst1 += dv()*fan1[];
      VolEst2 += dv()*fan2[];
    }
    fprintf(fp, "%g\t%g\t%g\t%g\n", R, fabs(VolEst1 - Va)/Va, fabs(VolEst2 - Va)/Va, fabs(0.5*(VolEst1+VolEst2) - Va)/Va);
  }
  view(theta = pi/5, phi = pi/8, width = 300, height = 300);
  draw_vof("fan1");
  save("fan.png");
}

/**
## Result

~~~gnuplot It does not really become appearant what formulation to choose.
  set xr [0.005:0.7]
  set logscale xy
  set xlabel 'R/128{/Symbol D} '
  set ylabel 'Rel. Error'
  set size square
  plot 'data' u 1:2 t 'Product formulation',\
        'data' u 1:3 t 'min formulation',\
        'data'u 1:4 t 'linear combination'
~~~
  */