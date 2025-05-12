/** 
# Viscous capillary waves - Comparison with Prosperetti & Cortelezzi solution.

We test here the ability of the [multilayer](/src/layered.h) solver and the 
surface tension extension to reproduce correctly viscous capillary waves 
for various regimes.

Numerical results are compared to the solution of [Cortelezzi & 
Propseretti](#prosperetti1982). The reference file was computed thanks to the 
Python code [cortelezzi_analytical.py](cortelezzi_analytical.py).
More details are available in an upcoming article. */

/** 
## Includes
We use a constant-resolution grid with the multilayer the Navier--Stokes solver 
and the [surface tension additive](/crober/2_Implicit/hydro-tension.h). */
#include "grid/multigrid1D.h"
#include "crobert/2_Implicit/hydro-tension.h"
#include "layered/nh.h"
#include "crobert/2_Implicit/viscous_surface.h"
#include "layered/remap.h"

/**
##Parameters */
#define Bo 0.001
#define lay 20
double H0 = 1. [1];
double delta=0.05;
double T0 = 1. [0,1];
double k = 0.1;
double sig_rho = 1.; //the sigma field in the solver actually contains sigma/rho since there is no density variable in the solver

/**
## Initialisation */
int nt = 0;
int main() {
  periodic(right);
  G = Bo;
  N = 64;
  nl = lay;
  TOLERANCE = 1e-8 [*];
  linearised = true;
  visc_activate = false;
  
  nl = 3;
  k = 0.1;
  L0 = 2.*pi/k;
  nu = 2.;
  T0 = (3*nu/(pow(H0,3)*(sq(k)*G + sig_rho*pow(k,4))));
  CFL_H = 10.;
  run();
  
  nl=10;
  nu = 0.01;
  T0 = (3*nu/(pow(H0,3)*(sq(k)*G + sig_rho*pow(k,4))));
  CFL_H = 0.5;
  run();
  
  nl = 15;
  k = 1;
  L0 = 2.*pi/k;
  nu = 0.01;
  T0 = (2*pi/sqrt((G + sig_rho*sq(k))*k*tanh(k*1[1])));
  run();
  
  visc_activate = true;
  run();
}

/**
The initial condition is a small amplitude plane wave of wavelength
unity and uniform surface tension and viscosity. */

event init (t = 0) {
  nt = 0;
  foreach () {
    double H = H0*(1. + delta*cos (k*x));
    foreach_layer () {
      h[] = H/nl;
    }
  }
  const scalar sig[]=sig_rho;
  sigma = sig;
}

/**
##Outputs */
/**
We output the evolution of the film (amplitude and vorticity)
 in a file indexed with the simulation parameters. */

event amplitude (t += T0/30.) {
  double max = statsf(eta).max;
  char name[80];
  if (visc_activate)
    sprintf (name, "ampSimu_Oh_%g_k_%g_corr", nu, k);
  else
    sprintf (name, "ampSimu_Oh_%g_k_%g", nu, k);
  static FILE * fpw = fopen (name, "w");
  fprintf (fpw, "%g %g\n", t/T0, (max-1)/delta);
  fflush (fpw);

}

event vorti (t += 0.1*T0; t < T0*1.05) {
  char name[80];
  if (visc_activate)
  sprintf (name, "vortSimu_Oh_%g_k_%g_corr", nu, k);
  else
  sprintf (name, "vortSimu_Oh_%g_k_%g", nu, k);
  
  static FILE * fpv = fopen (name, "w");
  if (nt == 1 || nt == 2 || nt == 5 || nt == 10) {
    foreach_point(L0/4){
      foreach_layer(){
        double vort = (u.x[0,0,1]- u.x[])/h[] - ((w[1] - w[-1])/(2*Delta));
        if (point.l != nl - 1)
          fprintf (fpv, "%g %g %g\n", eta[]*(point.l + 1)/nl, u.x[], vort*T0/delta);
      }
      fprintf (fpv,"\n\n");
    }
  }
  fflush (fpv);
  nt++;
}

/**
At the end of the simulation, we output on standard error the
simulation parameters, the resolution (number of layers and 
grid points) and the error on the final amplitude.*/
event logfile (t = end){
  double ampEnd = (k == 1 ? 0.745291 : (nu = 2 ? 0.37444 : 0.384595));
  double max = statsf(eta).max;
  double err = fabs(((max-1)/delta - ampEnd));
  fprintf (stderr, "%d %d %g %g %g\n", N, nl, k, nu, err);
}


/**
## Results
~~~gnuplot Capillary waves in three different regimes$
set term svg enhanced size 1000,500 font ",10"
set multiplot layout 2,3 scale 1,1 title "{/:Bold=15 Capillary waves in three different regimes : amplitude and vorticity}"

set xlabel 't{/Symbol \057}{/Symbol t}_{relax}'
set format y '%.1f'; 
set ylabel 'Relative amplitude'
set key at screen 1.,0.7
set size 0.3, 0.45
set origin 0.,0.45
set title "Oh = 2 and k = 0.1"
plot [0:1][0:1] exp(-x) w l lc rgb "#aadc32" t "Lubrication", \
                'ampCorte_Oh_2_k_0_1.txt' w l lc rgb "#472c7a" t "Cortelezzi", \
                'ampSimu_Oh_2_k_0.1' w p lc rgb '#0044a5' pt 2 pointsize 0.8 t "Multilayer"

unset key
set format y ''; unset ylabel
set size 0.3, 0.45
set origin 0.3,0.45
set title "Oh = 0.01 and k = 0.1"
plot [0:1][0:1] exp(-x) w l lc rgb "#aadc32" t "Lubrication", \
                'ampCorte_Oh_0_01_k_0_1.txt' w l lc rgb "#472c7a" t "Cortelezzi", \
                'ampSimu_Oh_0.01_k_0.1' w p lc rgb '#0044a5' pt 2 pointsize 0.8 t "Multilayer"
                

set xlabel 't{/Symbol \057}{/Symbol t}_{osc}'
set size 0.3, 0.45
set origin 0.6, 0.45
set title "Oh = 0.01 and k = 1"
plot [0:1][0:1] 'ampCorte_Oh_0_01_k_1.txt' w l lc rgb "#472c7a" t "Cortelezzi", \
                'ampSimu_Oh_0.01_k_1_corr' w p lc rgb '#0044a5' pt 2 pointsize 0.8 t "Multilayer"


set key at screen 1.,0.3
set xlabel '{/Symbol w} * {/Symbol t}_{relax}'
set size 0.3, 0.45
set origin 0.,0.
set format y '%.1f'; set ylabel 'Height'
set ytics
unset title
plot 'vortCorte_Oh_2_k_0_1.txt' u 2:1 w l lc rgb "#aadc32" t "0.1*T0", \
    'vortCorte_Oh_2_k_0_1.txt' u 3:1 w l lc rgb "#27ad81" t "0.2*T0", \
    'vortCorte_Oh_2_k_0_1.txt' u 4:1 w l lc rgb "#2c718e" t "0.5*T0", \
    'vortCorte_Oh_2_k_0_1.txt' u 5:1 w l lc rgb "#472c7a" t "T0", \
    'vortSimu_Oh_2_k_0.1' u 3:1 w p lc rgb '#0044a5' pt 2 pointsize 0.8 t "Multilayer"

unset key
set format y ''; unset ylabel  
set size 0.3, 0.45
set origin 0.3,0.
plot 'vortCorte_Oh_0_01_k_0_1.txt' u 2:1 w l lc rgb "#aadc32" t "0.1*T0", \
    'vortCorte_Oh_0_01_k_0_1.txt' u 3:1 w l lc rgb "#27ad81" t "0.2*T0", \
    'vortCorte_Oh_0_01_k_0_1.txt' u 4:1 w l lc rgb "#2c718e" t "0.5*T0", \
    'vortCorte_Oh_0_01_k_0_1.txt' u 5:1 w l lc rgb "#472c7a" t "T0", \
    'vortSimu_Oh_0.01_k_0.1' u 3:1 w p lc rgb '#0044a5' pt 2 pointsize 0.8 t "Multilayer"

set xlabel '{/Symbol w} * {/Symbol t}_{osc}'
set size 0.3, 0.45
set origin 0.6,0.
plot 'vortSimu_Oh_0.01_k_1_corr' u 3:1 w p lc rgb '#0044a5' pt 2 pointsize 0.8 t "Multilayer", \
    'vortCorte_Oh_0_01_k_1.txt' u 2:1 w l lc rgb "#aadc32" t "0.1*T0", \
    'vortCorte_Oh_0_01_k_1.txt' u 3:1 w l lc rgb "#27ad81" t "0.2*T0", \
    'vortCorte_Oh_0_01_k_1.txt' u 4:1 w l lc rgb "#2c718e" t "0.5*T0", \
    'vortCorte_Oh_0_01_k_1.txt' u 5:1 w l lc rgb "#472c7a" t "T0"
    
unset multiplot
~~~
The results are presented for three different regimes : relaxation regime (in 
 the lubrication domain), aperiodic regime (low k, low viscosity) and 
 oscillating regime (here, medium k, low viscosity). In all cases, multilayer
 model reproduce very well the results of Cortelezzi and Prosperetti.

~~~gnuplot Capillary waves in three different regimes$
set multiplot layout 2,2 scale 1,1 title "{/:Bold=15 Comparison with or with viscous correction terms}"

set origin 0.,0.45
set xlabel 't{/Symbol \057}{/Symbol t}_{relax}'
set format y '%.1f'; 
set ylabel 'Relative amplitude'
set key at screen 1.,0.7
set size 0.45, 0.45
set title "With correction"
plot [0:1][0:1] 'ampCorte_Oh_0_01_k_1.txt' w l lc rgb "#472c7a" t "Cortelezzi", \
                'ampSimu_Oh_0.01_k_1_corr' w p lc rgb '#0044a5' pt 2 pointsize 0.8 t "Multilayer"
                
unset key
set format y ''; unset ylabel
set xlabel 't{/Symbol \057}{/Symbol t}_{osc}'
set size 0.45, 0.45
set origin 0.45,0.45
set title "Without correction"
plot [0:1][0:1] 'ampCorte_Oh_0_01_k_1.txt' w l lc rgb "#472c7a" t "Cortelezzi", \
                'ampSimu_Oh_0.01_k_1' w p lc rgb '#0044a5' pt 2 pointsize 0.8 t "Multilayer"


set key at screen 1.,0.3
set xlabel '{/Symbol w} * {/Symbol t}_{relax}'
set size 0.45, 0.45
set origin 0.,0.
set format y '%.1f'; set ylabel 'Height'
set ytics
unset title
plot 'vortSimu_Oh_0.01_k_1_corr' u 3:1 w p lc rgb '#0044a5' pt 2 pointsize 0.8 t "Multilayer", \
    'vortCorte_Oh_0_01_k_1.txt' u 2:1 w l lc rgb "#aadc32" t "0.1*T0", \
    'vortCorte_Oh_0_01_k_1.txt' u 3:1 w l lc rgb "#27ad81" t "0.2*T0", \
    'vortCorte_Oh_0_01_k_1.txt' u 4:1 w l lc rgb "#2c718e" t "0.5*T0", \
    'vortCorte_Oh_0_01_k_1.txt' u 5:1 w l lc rgb "#472c7a" t "T0"

unset key
set format y ''; unset ylabel  
set xlabel '{/Symbol w} * {/Symbol t}_{osc}'
set size 0.45, 0.45
set origin 0.45,0.
plot 'vortSimu_Oh_0.01_k_1' u 3:1 w p lc rgb '#0044a5' pt 2 pointsize 0.8 t "Multilayer", \
    'vortCorte_Oh_0_01_k_1.txt' u 2:1 w l lc rgb "#aadc32" t "0.1*T0", \
    'vortCorte_Oh_0_01_k_1.txt' u 3:1 w l lc rgb "#27ad81" t "0.2*T0", \
    'vortCorte_Oh_0_01_k_1.txt' u 4:1 w l lc rgb "#2c718e" t "0.5*T0", \
    'vortCorte_Oh_0_01_k_1.txt' u 5:1 w l lc rgb "#472c7a" t "T0"
    
unset multiplot
~~~
In the case of high k, the terms of the continuity of surface stress must be
 added to ensure a correct reproduction of the oscillations.
 These are added in [viscous_surface.h](../viscous_surface.h)
 The normal projection of the surface stress continuity has an influence on the
 wave attenuation and the tangential projection impacts vorticity. Besides,
 stratification is very important in this cases (as vorticity diffusion is not
 immediate).
 
## References
~~~bib
@article{prosperetti1982,
  title={Small-amplitude waves produced by a submerged vorticity distribution on the surface of a viscous liquid},
  author={Prosperetti, Andrea and Cortelezzi, L},
  journal={The Physics of Fluids},
  volume={25},
  number={12},
  pages={2188--2192},
  year={1982},
  publisher={American Institute of Physics}
}
~~~
*/
