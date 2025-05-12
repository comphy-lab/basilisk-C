/**
![At a cloud interface evaporation occurs due to mixing of moisture and temperature. Image via [`wideopen spaces`](https://cdn0.wideopenspaces.com/).](https://cdn0.wideopenspaces.com/wp-content/uploads/2014/11/bigstock-Cumulonimbus-71055340.jpg)

# Subsiding Shells at the edge of a cumulus-cloud's updraft

It is well known that mixing at a cumulus cloud's interface can cause a so-called subsiding shell. Nair et al. (2019) propose a scenario to study these physics in an idealized setting. We follow their setup for the flow configuration.


<div class="figure">
<video controls="" preload="metadata" width="900">
<source src="https://surfdrive.surf.nl/files/index.php/s/xsUF3CMV9etQgmA/download
" type="video/mp4">
Your browser does not support the video tag. </video>
<p class="caption">
Volumetric rendering of the cloud (liquid water field) (via surfdrive)
</p>
</div>
*/

#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"
#include "fractions.h" //For the heavyside function;
/**
A sequence of `DEFINE`s is used to compute the cloud-liquid water
content as a function of pressure and temperature. Note that there is
no field dedicated to store the liquid water specific humidity, its prognostic.
*/
double Rv = 461.5;        //in J/kg/K
double Rd = 287.;         //in J/kg/K
double L = 2500000.;      //in J/kg
double cpd = 1004.;       //in J/Kg/Kq
#define P00 100000.       //A constant reference pressure
double P0 = 101780., Pe = 7500; //Surface pressure and 1/e length scale

#ifndef RHO0            
#define RHO0 1.3 //Base density in kg/m^3  
#endif
#define pres  (P0*(exp(-1500./Pe)))  
#define EXNER (pow(pres/P00, Rd/cpd))
#define TK1   (thl[]*EXNER)
#define ES    (610.78*exp(17.27*(TK1 - 273.16)/(TK1 - 35.86)))
#define QSL   (Rd/Rv*(ES/(pres - (1. - Rd/Rv)*ES)))
#define QSAT  (QSL*((1. + ((sq(L)/(Rv*cpd*sq(TK1)))*qt[]))/(1 + ((sq(L)/(Rv*cpd*sq(TK1)))*QSL))))
#define QC    (max(qt[] - QSAT, 0))
/**
The virtual potential temperature in a cell ($\theta_v$,`THV`) can
be computed from $\theta_l$, $q_c$ and $q_t$.
*/
#define THETA (thl[] + L*QC/(cpd*EXNER)) 
#define THV (THETA*(1. - (1. - Rv/Rd)*(qt[] - qt0) - Rv/Rd*QC))
/**
   The virtual potential temperature is related to the buoyancy,
   */
double g_acc = 9.81;
#define B_CENTERED (g_acc*((THV - T_ref)/T_ref))
/**
## The scenario

The case is defined by in-cloud and environment variables:
*/
double T_ref = 294., qt0 = 0.004;
double qtc = 0.009, thvc = 289;
double Lc = 2, wc = 1;
scalar qt[], thl[], * tracers = {qt, thl};
scalar H[]; //Heavyside function

face vector av[];
double tau; 
int main() {
  H.prolongation = H.refine = fraction_refine; //Heavyside fraction
  periodic (bottom);
  L0 = 30;
  a = av;
  tau = Lc/wc;
  const face vector muc[] = {1./1000., 1./1000., 1./1000.};
  mu = muc;
  N = 128;
  run();
}

event init (t = 0) {
  fraction (H, Lc - x); //heavyside fraction 
  foreach() {
    thl[] = thvc*H[] + T_ref*(1 - H[]);
    qt[]  = qtc*H[]  + qt0*(1 - H[]);
    u.y[] = wc*H[];
    u.x[] = 0.001*noise();
  }
  TOLERANCE = 1e-4;
}
/**
The temperature and humidity are subject to diffusion and a in-cloud forcing.
*/
event tracer_diffusion(i++) {
  diffusion (thl, dt, mu);
  diffusion (qt, dt, mu);
  foreach() {
    thl[] += H[]*dt * (thvc - thl[])/tau;
    qt[]  += H[]*dt * (qtc  - qt[] )/tau;
  }
}
/**
The effects of gravity are including via the prognostic buoyancy field. Further, the velocity in the cloud is subject to a forcing.
*/
event acceleration (i++) {
  scalar b[];
  foreach()
    b[] = B_CENTERED;
  boundary ({b, H});
  foreach_face(y)
    av.y[] = ((b[] + b[0,-1])/2. +
	      (wc - (u.y[] + u.y[0,-1])/2.)/tau * (H[] + H[0,-1])/2.);
}
/**
## Grid Adaptation

In order to focus the computational resources, we adaptively refine
and coarsen the grid based on the discrete representation of the
temperature, humidity and velocity-component fields.
*/
#if (dimension == 2)
event adapt (i++) {
  adapt_wavelet ({thl, qt, u}, 
                 (double[]){6./10., 0.005/10., wc/10., wc/10., wc/10.}, 9);
}
#else
event adapt (i++) { 
  adapt_wavelet ({thl, qt, u},
                 (double[]){1., 0.01/10., wc/8., wc/8., wc/8.}, 8);
}
#endif
/**
## Results
   
We output a movie:

![Looks OK](subsiding-shell/sss.mp4)
   
movie generation code:
*/
#if (dimension == 2)
#include "view.h"

void cloud (double cmap[NCMAP][3]) {
  for (int i = 0; i < NCMAP; i++) {
    cmap[i][0]= 0.2 + 0.8 * (double)i/(double)NCMAP;
    cmap[i][1]= 0.2 + 0.8 * (double)i/(double)NCMAP;
    cmap[i][2]= 0.9 + 0.1 * (double)i/(double)NCMAP; 
  }
}

event movie (t += 0.5) {
  scalar b[], ql[];
  foreach() {
    b[] = B_CENTERED;
    ql[] = QC;
  }
  boundary ({ql, b});
  
  view (fov = 40, width = 900, height = 900);
  
  squares ("b", min = -0.05, max = 0.05, linear = true);  
  mirror ({-1,0})
    cells();
  translate (y = -L0) {
    squares ("u.y", min = -wc, max = wc, linear = true);  
    mirror ({-1,0})
      squares ("ql", map = cloud, min = 0, max = 0.002, linear = true);
  }
  draw_string ("Buoyancy", 2);
  draw_string ("Grid Cells", 1);
  draw_string ("Vertical velocity", 3);
  draw_string ("Liquid water", 0);
  save ("sss.mp4");
}

event stop (t = 120);

#else //dimension == 3 ? use bwatch...
#include "bwatch.h"
event bwatcher (t += .25) {
  static FILE * fp = popen ("ppm2mp4 shell.mp4", "w");
  scalar ql[];
  foreach()
    ql[] = QC;
  boundary ({ql});
  watch (fov = 1.2*L0, O = {50, 10, 60}, poi = {-1, 15, 15},
	 nx = 1024, ny = (24*40));
  lattice (width = 0.04, alpha = 0.1);
  volume (ql, sc = 0.0025, shading = 1, cols = true, min = -0.1,
	  max = 0.1, mval = 0.0001);
  store (fp);
  plain();
}

event stop (t = 250);
#endif

/**
## Reference:

Nair, V., T. Heus, and M. van Reeuwijk, *Dynamics of Subsiding Shells In Actively Growing Clouds With Vertical Updrafts.* J. Atmos. Sci.,  [https://doi.org/10.1175/JAS-D-19-0018.1](https://doi.org/10.1175/JAS-D-19-0018.1)
*/
