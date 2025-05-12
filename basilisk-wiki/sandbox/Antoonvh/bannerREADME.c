/**
Eventough the [Rayleigh-Taylor instability](rt.c) is more suitable, we use the Kelvin-Helmholtz (KH) instability to generate a banner for the [readme page](README):

![](bannerREADME/banner.png)
![](bannerREADME/banner_cells.png)

![This does not display](combinedbanner.png)

The setup is inspired by the [2D KH study](kh.c).
*/

#include "navier-stokes/centered.h"
#include "view.h"

void black_body (double cmap[NCMAP][3])
{
  /* black body color map from:
   * http://www.kennethmoreland.com/color-advice/
   */
  static double basemap[32][3] = {
    {0.0,0.0,0.0},
    {0.0857913205762,0.0309874526184,0.0173328711915},
    {0.133174636606,0.0588688899571,0.0346802666087},
    {0.180001956037,0.0730689545154,0.0515393237212},
    {0.22981556179,0.0840603593119,0.0647813713857},
    {0.281397607223,0.093912584278,0.075408501413},
    {0.334521638801,0.102639499627,0.0842454688083},
    {0.388957802186,0.110254429637,0.0927990674821},
    {0.444611925648,0.116732501721,0.101402659637},
    {0.501422312285,0.122025816585,0.110058408122},
    {0.559331322331,0.126067584009,0.118767796491},
    {0.618285970576,0.128767919785,0.127531801155},
    {0.678237857955,0.130007052818,0.136351016263},
    {0.712849583079,0.181721849923,0.13081678256},
    {0.743632057947,0.232649759358,0.120991817028},
    {0.774324938583,0.279315911516,0.108089917959},
    {0.804936242903,0.323627020047,0.0907961686083},
    {0.835473266757,0.366524681419,0.0662363460741},
    {0.865942668698,0.408541395043,0.026029485466},
    {0.876634426153,0.46401951695,0.0173065426095},
    {0.883455346031,0.518983528803,0.0149628730405},
    {0.88905246237,0.572164381169,0.013499801006},
    {0.893375939063,0.624108797455,0.0130334871745},
    {0.89637036663,0.675180034619,0.013680092215},
    {0.897973818846,0.725630730259,0.015555776796},
    {0.898116710502,0.775642817733,0.0187767015864},
    {0.896720396485,0.825350944866,0.023459027255},
    {0.927670131094,0.859991226192,0.319086199143},
    {0.956158602738,0.893933112845,0.503316730316},
    {0.97827065392,0.92856476667,0.671307024002},
    {0.993196411712,0.963913323002,0.83560909192},
    {1.0,1.0,1.0},
  };
  for (int i = 0; i < NCMAP; i++) {
    double x = i*(31 - 1e-10)/(NCMAP - 1);
    int j = x; x -= j;
    for (int k = 0; k < 3; k++)
      cmap[i][k] = (1. - x)*basemap[j][k] + x*basemap[j+1][k];
  }
}
int maxlevel = 11;
double Re = 80000.;
const face vector muc[] = {1./Re, 1./Re};

int main(){
  periodic(left);
  X0 = Y0 = -L0/2.;
  mu = muc;
  init_grid(1 << 8);
  run();
}

event init(t = 0){
  refine (fabs(y) < 0.025 && level <= 8);
  refine (fabs(y) < 0.001 && level <= 10);
  foreach(){
    u.x[] = y > 0 ? -0.5 : 0.5;
    u.y[] += 0.0035*noise();
  }
  boundary(all);
}

event adapt(i++)
  adapt_wavelet ({u.x, u.y}, (double[]){0.01, 0.01}, maxlevel);

event snap (t = 0.2){
  view (fov = 0.75, width = 2048, height = 80, samples = 2);
  cells();
  save("banner_cells.png");
  clear();
  scalar omega[];
  vorticity(u, omega);
  boundary({omega});
  while (adapt_wavelet({omega}, (double[]){0.005}, maxlevel + 1).nf){
    boundary({omega});
  }
  foreach()//Invert colors
    omega[] *= -1;
  boundary({omega});
  dump("dump");
  squares("omega", linear = true, map = black_body, min = -185, max = -1);
  save("banner.png");
  system("rm combinedbanner.png");
  system("ffmpeg -i banner.png -i banner_cells.png -filter_complex vstack combinedbanner.png");
  return 1;

}