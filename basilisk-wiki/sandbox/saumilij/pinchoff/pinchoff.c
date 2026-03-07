/** 
# Viscous Plateau-Rayleigh Pinchoff

This case contains the code for axisymmetric simulation of the Rayleigh-Plateau instability with finite viscosity.
For very small viscosities the universal similarity scalings for pinch-off dynamics are reproduced:
minimum radius $r_{\min} \sim (t_0 - t)^{2/3}$ and maximum velocity
$u_{\max} \sim (t_0 - t)^{-1/3}$, where $t_0$ is the pinch-off time.
The setup is similar to the [Basilisk test case](https://basilisk.fr/src/test/plateau.c), but with finite viscosity and a different approach to refine the mesh during the pinchoff. In the inviscid limit, the results are identical to those of the [test case](https://basilisk.fr/src/test/plateau.c).

## Physical Setup
A liquid cylinder of initial radius $R_0 = 1$ is perturbed sinusoidally with a wave number $k = 0.5$ and a small amplitude of peturbation $\epsilon = 0.1$. Surface tension drives the Rayleigh-Plateau instability resulting in necking and pinch-off of the cylinerical filament to form satellite drops. 

The initial perturbations is given as:

$$R(x, t=0) = R_0(1 + \epsilon \sin(kx))$$. 

The non-dimensional parameter Ohnesorge number,
 $$Oh = \frac{\eta_1}{\rho_1 \gamma R_0}$$

 denotes the ratio of the inertio-capillary and the visco-capillary timescales. Here $\gamma$, $\rho_1$ and $\eta_1$ are the coefficient if surface tension, density and dynamic viscosity of the liquid.

The air to liquid density ratio is kept fixed $\rho_2/\rho_1 = 10^{-3}$ and the surrounding viscosity is kept negligible. Subscript 1 and 2 are used to indicate the liquid and the gas phase respectively.

## Numerical code
We solve the axisymmetric, incompressible Navier-Stokes equations with two phases and surface tension.
*/

#include "axi.h" //sets the case to axisymmetric
#define FILTERED 1 // Smear density and viscosity jumps
#include "navier-stokes/centered.h"
#include "two-phase.h" //f: 1 is liq, 0 is gas phase
#include "navier-stokes/conserving.h"
#include "tension.h" //include surface tension
#include "view.h"

//error tolerances (for mesh refinement)
#define fErr (1e-3)
#define VelErr (1e-3)
#define KErr (1e-4)

//time-intervals for saving
#define tsnap (0.01)  //snapshot saving interval

/**
 * Id1 :indicates the liquid media formimg the jetted filament
* Id2: indicates the surrounding gas/fluid(Newtonian)
*/

#define Rho21 (1e-3)

//declarations
int MAXlevel, MINlevel;
int maxlimit = 18;//maximum limit for grid MAXlevel at pinchoff 
double tmax, Oh1, Ldomain, epsilon;


//Boundary conditions
//x-axis axisymmetric//top - outflow, others symmetric
u.n[top] = neumann(0);
p[top] = dirichlet(0);


int main(int argc, char const *argv[]){
  //assignments
  MAXlevel = 10; //for grid refinement
  MINlevel = max(6, (MAXlevel-4));//for grid refinement
  tmax = 10.0;
  Ldomain = 2*pi; //domain length for simulation
  Oh1 = 1e-3; //liq film Ohnesorg no
  epsilon = 0.1; //amplitude of perturbation

  fprintf(ferr, "Level %d, tmax %g, Oh1 %3.2e, epsilon %g, Lo %g\n", MAXlevel, tmax, Oh1, epsilon, Ldomain);

  //domain and grid
  L0=Ldomain;
  X0=0.; Y0=0.;
  init_grid (1 << (MINlevel));

  //make a directory "intermediate" to store the snapshots
  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);

  //properties
  rho1 = 1.0; rho2 = Rho21;//phase densities
  f.sigma = 1; //coefficient of surface tension
  mu1 = Oh1; mu2 = 1e-4;//phase viscosities

  //run the simulation
  run();
}

/**
### Initial Condition
We start with a liquid cylinder perturbed with 10% amplitude sinusoidal deformation. Perturbation parameters are wavenumber $k = 0.5$, $\epsilon = 0.1$ and a base radius $R_0 = 1$.
 */
event init (t = 0){
    if(!restore (file = "dump")){
        refine((y < 1+epsilon) && (level < MAXlevel));
        
        double k = 0.5;
        
        vertex scalar phi[];
        foreach_vertex() {
            phi[] = (1+epsilon*cos(k*x))-y;
        }
        fractions (phi, f);
  }
}

/** 
### Adaptive Mesh refinement 
We adapt the mesh according to the errors of the volume fraction field, velocity and the interface curvature. 
We check if the filament breaks or not. If the filament is not broken, we increase MAXlevel as the filament thins down upto a maximum level of 18, until breakup, and after breakup continue using the minimum value of MAXlevel.
*/ 

event adapt (i++){
  //Allocate a scalar field `Y' to store the radial position of the interface relative to the axis of symmetry.
  scalar Y[];
  position (f, Y, {0,1});
  boundary({Y});
  //check if the filament is broken
  static bool broken = false;
  double y_min = statsf(Y).min;
  if (!broken)
    broken = y_min < 1./(1 << MAXlevel);

  if(broken){ MAXlevel = 10;}
  else{
    //If filament is not broken, and the height of the interface is less than ~5*grid cells, we increment the MAXlevel by 1 (upto a maximum limit)
    while(((statsf(Y).min)<(5*L0/(1<<MAXlevel))) && (MAXlevel<maxlimit)){
        MAXlevel = MAXlevel + 1;
    }
  }
  
  scalar KAPPA[]; //field for curvature
  curvature(f, KAPPA);
  boundary({KAPPA});

  adapt_wavelet ((scalar *){f, u.x, u.y, KAPPA},
    (double[]){fErr, VelErr, VelErr, KErr},
    MAXlevel, MINlevel);
}


/** 
### Outputs
We save snapshots of the simulation at regular intervals to restart the simulation or to post-process with bview.
*/
//static
event writingFiles (t = 0, t += tsnap; t <= tmax) {
  dump (file = "dump");
  char nameOut[80];
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump (file = nameOut);
}
/** We log the minimum position of the interface relative to the axis of symmetry, as well as the maximum axial velocity with time.
 */

event logWriting (i++) {
  
  scalar Y[];
  position (f, Y, {0,1});
  boundary({Y});

  double h_min, ux_max;
  h_min = statsf(Y).min;
  ux_max = normf(u.x).max;

  static FILE * fp;
  if (pid() == 0){
    if (i == 0) {
      fprintf (ferr, "i dt t y_min1 x_min1 y_min2 x_min2 ux_max\n");
      fp = fopen ("log", "w");
      fprintf (fp, "Level %d, tmax %g, Oh %3.2e, epsilon %3.2e, Lo %g\n", MAXlevel, tmax, Oh1, epsilon, Ldomain);
      fprintf (fp, "i dt t h_min ux_max\n");
      fprintf (fp, "%d %g %.12f %.12f %.12f\n", i, dt, t, h_min, ux_max);
      fclose(fp);
    } else {
      fp = fopen ("log", "a");
      fprintf (fp, "%d %g %.12f %.12f %.12f\n", i, dt, t, h_min, ux_max);
      fclose(fp);
    }
    fprintf (ferr, "%d %g %.12f %.12f %.12f\n", i, dt, t, h_min, ux_max);
  }
}
/**
Animation Generate movie showing mesh refinement and interface evolution.
*/
event movie (t = 9.6; t += 0.0002; t <= 9.71)
{
  view (quat = {0.000, 0.000, 0.000, 1.000},
      fov = 30, near = 0.01, far = 1000,
      tx = -0.438, ty = -0.171, tz = -0.636,
      width = 1428, height = 962);
  squares (color = "level", min = MINlevel, max = maxlimit, spread = -1);
  draw_vof (c = "f");
   mirror ({0,-1}) {
    squares (color = "level", min = MINlevel, max = maxlimit, spread = -1);
    draw_vof (c = "f");
  }
  save ("movie.mp4");
}
/**
## Results
The animation shows how adaptive mesh refinement tracks high curvatures and short
timescales near pinch-off. Up to 18 levels of refinement capture roughly four
orders of magnitude in spatial scales.

![Mesh and interface evolution](pinchoff/movie.mp4)(width="50%")

The scaling plots below show the theoretical fits. The fit is excellent over
at least four orders of magnitude. Departures from power laws near pinch-off
result from saturation of spatial resolution (grid size $L_{domain}/2^{18}$).

~~~gnuplot Evolution of the minimum radius
reset
# We define the breakup time as the time when the axial velocity is maximum
t0="`awk '{ if (NR>2 && NF == 5 && $5 > max) { max = $5; t0 = $3; } } END{ print t0 }' < log`"
set key spacing 1.5
set grid
set xlabel 't_0 - t' font ', 12'
set ylabel 'R_{min}' font ', 12'
set logscale
set format x '10^{%L}'
set format y '10^{%L}'
fit [1e-5:1e-1] a*x**(2./3.) 'log' u (t0 - $3):4 via a
plot [1e-7:][6.2831/2**18:]'log' u (t0 - $3):4 ps 0.5 t '', a*x**(2./3.) t 'x^{2/3}'
~~~
 */

/**
## See Also
- [Basilisk test case - plateau.c](https://basilisk.fr/src/test/plateau.c)
- [Same test with Gerris](http://gerris.dalembert.upmc.fr/gerris/tests/tests/plateau.html)
- [3D Plateau example with Gerris](http://gerris.dalembert.upmc.fr/gerris/examples/examples/plateau.html)
*/