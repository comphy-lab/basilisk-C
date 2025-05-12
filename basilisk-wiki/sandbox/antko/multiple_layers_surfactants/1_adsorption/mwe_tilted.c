/**
# Diffusion of surfactants
We investigate the diffusion of surfactants in water, including the adsorption at the air/water adsorption.
Surfactants diffusion answers to a Fick diffusion equation of a specie $c$.
$$
\frac{\partial c}{\partial t} =\nabla^2 c 
$$
which can be solved with the reaction--diffusion solver.

The interface with air is taken into account, the adsorption of surfactant is modelled using a Langmuir adsorption.*/

/**
## Includes*/
#include "grid/multigrid.h"
#include "fractions.h"
#include "diffusion.h"
#include "../mixtures_cr.h"
#include "../../qmagdelaine/my_functions.h"
#include "run.h"

/**
## Geometry and resolution
L is the size of the box, alpha the parameters of the line tilt (between -1 and 1)*/
#define L 10.

#define LEVEL 6
#define alpha 0.4
#define MY_TOLERANCE 1e-11

#define F_ERR 1e-10

#define T_END 2
#define DELTA_T (1e-1)
/**
## Physical parameters
In the bulk:
* Diffusion coefficient of the solute in the liquid
* Initial concentration of the solute*/
#define Diff 1.
#define solute0 1.

/**
At the interface:
* Initial surface excess concentration
* Equilibrium constant for surfactant adsorption at the interface*/
#define gamma0 0.

#define omega 30.
#define sigma 10.
#define DT_MAX 1./(omega*50.)
#define security 30.

/**
## Fields and functions*/
/**
Several scalars and vector fields are allocated to describe the interface, the concentration and he surface concentration fields*/
scalar solute[];
scalar surface[];
scalar c_solute[];
scalar c_solute_diagonal[];
scalar f[], L_I[];

/**
Attributes are required for diffusion*/
attribute {
  double D;
  double inverse;
}

/**
During the step of adsorption at the interface, the solute concentration is modified. If we do not take care, at equal adsorption rate, the less liquid fraction is in the cell, the more the solute will be concentrated. This is not physical: a cell should be concentrated more quickly than its neighbour just because its liquid fraction is smaller. To solve this problem we try to define two functions to compute a concentration at the interface which average the concentration over a reconstructed cell, which have always the same volume (or area in 2D). This reconstructed cell share its volume between the cell at the interface and its neighbors accondingly to the normal to the interface.*/
double interfacial_concentration (Point point, scalar f, scalar solute) {
  double local_solute = solute[]*f[];
  coord n = interface_normal (point, f);
  foreach_dimension() {
    int i = (n.x > 0. ? -1 : 1);
    local_solute += fabs(n.x)*solute[i]*(1. - f[]);
  }
  return local_solute;
}

double interfacial_concentration_diagonal (Point point, scalar f, scalar solute) {
  double local_solute = solute[]*f[];
  coord n = interface_normal (point, f);
  double factor = pow((pow(n.y,2) - pow(n.x,2))/(pow(n.y,2) + pow(n.x,2)),2);
  
  foreach_dimension() {
    int i = (n.x > 0. ? -1 : 1);
    local_solute += fabs(n.x)*solute[i]*(1. - f[])*factor;
  }
  int i = (n.x > 0. ? -1 : 1);
  int j = (n.y > 0. ? -1 : 1);
  local_solute += solute[i,j]*(1. - f[])*(1-factor);
  return local_solute;
}

/**
We will store the statistics on the diffusion solvers in `mgd`. */
mgstats mgd;

/**
## Boundary conditions */
/**
* On the top and bottom boundaries the species concentration is imposed. */
solute[bottom] = neumann(0);
solute[top] = neumann(0);

/**
* And there is no flux on the right and left boundaries. */
solute[right] = neumann(0);
solute[left]  = neumann(0);

/**
## Parameters
The domain is the square box. */

int main() {
  L0 = L;
  X0 = -L0/2;
  N = 1 << LEVEL;
  DT = DT_MAX;
  TOLERANCE = MY_TOLERANCE;
  
  run();
}

/**
## Initialisation

Before the first step, we initialise the concentration fields to $c_0$
in the bulk and \gamma_0 at the interface. */

double initial_solute;

event init (i = 0) {
/**
* Volume fraction:*/
  fraction(f, L0/2 - y + alpha * x);
  boundary ({f});

/**
* Bulk and surface concentration: */ 
  foreach() {
    solute[] =  f[]*solute0;
    L_I[] = interface_length(point, f)/Delta;
    surface[] = (interfacial(point, f) ? gamma0*L_I[] : 0.);
    c_solute[] = 0;
    c_solute_diagonal[] = 0;
  }  
  boundary({solute, L_I, surface});
  initial_solute = statsf(solute).sum + statsf(surface).sum;
}

/**
## Integration*/
event integration (i++) {

/**
##Adsorption-Desorption and bulk diffusion of surfactants */
  
  /**
  The quantity fields are turned into concentration fields. */
  foreach() {
    solute[] = (f[] > F_ERR ? solute[]/f[] : 0.);
    L_I[] = interface_length(point, f)/Delta;
    surface[] = (L_I[] > F_ERR ? surface[]/L_I[] : 0.);
  }
  boundary({solute, L_I, surface});
    
  /**
  The adsorption velocity for each cell is computed. The dtmax is corrected if to high. */
  static FILE * fpW = fopen("dataW.txt", "w");
  double ads = 0.;
  double N = 0.;
  scalar adsorbed[], source[];
  double dtt = DT_MAX;
  vector n[];
    
  foreach(){
    if (interfacial(point, f)){  
      c_solute[] = interfacial_concentration(point, f, solute);       
      c_solute_diagonal[] = interfacial_concentration_diagonal(point, f, solute);
      adsorbed[] = omega * (sigma * c_solute[] * (1 - surface[]) - surface[]);

      if (fabs(adsorbed[]) > 0)
        dtt = min(dtt, max(Delta/(sigma*omega), fabs(c_solute[]/(security*adsorbed[]*L_I[]))));
      ads += adsorbed[];    
      N += 1;  
    }
    else
      adsorbed[] = 0.;
  }
  
  boundary({adsorbed});
   
  /**
  The adsorption velocity is transformed into a source term that will be added in the diffusion solver. 
  
  !!!!!! Here, the averaged adsorption flux is taken into account, contrary to the final program !!!!!!!!!
  */ 
  foreach(){
      source[] = - ads/N*L_I[];
  }
  boundary({source});

  /**
  The timestep is set according to the timing of upcoming events and to avoid negative concentration in interfacial cells. */
  dt = dtnext (dtt);
  
  /**
  We use the diffusion solver to advance the system from $t$
  to $t+dt$ : the solute is diffused in the bulk. 
  Adsorbed surfactants are removed from the bulk */
  solute.D = Diff;
  solute.inverse = false;
  mgd = with_flux_diffusion (solute, source, f, dt);

  /**
  Adsorbed surfactants are added to the interface. */
  foreach()
    if (interfacial(point, f)) {  
      surface[] += ads/N*dt;
    }
  boundary({surface});  

  /**
  The concentration field is turned back into a quantity field. */  
  foreach() {
    solute[] *= f[];
    surface[] *= L_I[];
  }
  boundary({solute, surface});
}

/**
## Results*/
event viewing (t += DELTA_T; t <= T_END) {
  static FILE * fpl;
  scalar l[];
  output_ppm (f = solute, fp = fpl, file = "timeline.mp4", min = 0, max = 1, linear = false, map=cool_warm);         


  static FILE * fpI = fopen("profil1.txt", "w");
    foreach()
      if (interfacial(point, f))
        fprintf (fpI, "%g %g %g %g %g %g %g %g %g\n", x, y, (L_I[] > F_ERR ? surface[]/L_I[] : 0.),  (f[] > F_ERR ? solute[]/f[] : 0.), c_solute[], c_solute_diagonal[], L_I[], f[],t);
    fprintf (fpI, "\n");
}   

event outputs(t += DELTA_T; t = DELTA_T; t <= T_END) {
  double total_solute = statsf(solute).sum ;
  double total_surface = statsf(surface).sum ;
  static FILE * fp = fopen("surf_mass.txt", "w");
  fprintf (fp, "%g %g %g %g %g\n", t, total_solute, total_surface, total_solute + total_surface, total_solute + total_surface - initial_solute);
}

/**
The surfactants are adsorbed from bulk to surface.
~~~gnuplot Evolution of mass of surfactant
set output 'plot1.png'
set xlabel "t"
plot \
  './surf_mass.txt' u 1:2 t 'bulk' w l, \
  './surf_mass.txt' u 1:3 t 'surf' w l, \
  './surf_mass.txt' u 1:4 t 'tot' w l
~~~

The error on the total mass of surfactant is checked and should stay under 1e-12
~~~gnuplot Error on the total mass of surfactants.
set terminal @PNG enhanced size 640,640 font ",8"
set output 'plot0.png'
set xlabel "t"
set key left
plot './surf_mass.txt' u 1:5 t 'erreur' w l
~~~

As the adsorption is averaged over the line, surface concentration is homogeneous at a give time.
~~~gnuplot Surface concentration
set output 'plot2.png'
L0=5
set xlabel "x"
set ylabel "gamma"
set yrange [0:1]
plot \
  './profil1.txt' u 1:3 t 'c' w l
~~~

Adsorption of surfactants is equivalent over all the surface but cell size is not a constant. As a consequence, the subsurface concentration is not homogeneous at a given time t.
~~~gnuplot Subsurface concentration
set output 'plot3.png'
L0=5
set xlabel "x"
set ylabel "c_s"
set yrange [0:1]
plot \
  './profil1.txt' u 1:4 t 'c' w l
~~~

~~~gnuplot Subsurface concentration with correction other direct neighbours.
set output 'plot4.png'
L0=5
set xlabel "x"
set ylabel "c_s"
set yrange [0:1]
plot \
  './profil1.txt' u 1:5 t 'c' w l
~~~

~~~gnuplot Subsurface concentration with correction other direct and diagonal neighbours.
set output 'plot5.png'
L0=5
set xlabel "x"
set ylabel "c_s_diag"
set yrange [0:1]
plot \
  './profil1.txt' u 1:5 t 'c' w l
~~~

## Ameliorations
Two problems appear at two different scales. 
* First an overall variation, due to diffusion and boundary conditions. This should be quite reduced in the real problems.
* Then a quick variation due to the size of the cells. The average over the neighbooring cells helps but do not resolve everything...
*/