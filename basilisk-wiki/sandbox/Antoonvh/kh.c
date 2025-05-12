/**
![Clouds can serve as tracers to reveal the occurrence of a shear
 instability in the atmosphere. Image courtesy of
 [EarthSky](http://earthsky.org/earth/kelvin-helmholzt-clouds).](http://en.es-static.us/upl/2016/05/kelvin-hemhotltz-wave-cloud-mashka.jpg)

# The Kelvin-Helmholtz instability in two dimensions

On this page a classical shear instability is simulated. According to
the Kelvin-Helmholtz instability description, a thin shearing layer
may be unstable and can then roll-up on itself creating vortices, or
so-called Kelvin-Helmholtz billows. On this page the effect of the
fluid's viscosity on the observed dynamics is investigated.

## Set-up

A fluid with an infinetly thin shear layer is initialized. The
difference in velocity between the induvidual layers *U* and the
length scale associated with the periodicity of the solution
($\mathcal{L}$) in our set-up can be used to normalize the value of
the fluid's viscosity. Resulting in a Reynolds number,

$$Re=\frac{U\mathcal{L}}{\nu}.$$  

In this study, the Reynolds number is varied to be $Re=\{20000, 40000,
80000\}$. In order to analyze the flow's evolution, the emergence of
coherent vortex structures is monitored. Therefore, the number of
vortices is counted using the marvalous `tag()` function. It is
remarkable how complex seamingly simple things can be.
*/
#include "navier-stokes/centered.h"
#include "tag.h"

int maxlevel;
double vis, Re;
FILE * fpn;
char fnamev[99];
char fnameg[99];

/**
  As mentioned earlier periodic boundaries are used in the stream-wise
  direction. Furthermore, the boundaries in the span-wise direction of
  the square domain are of the free-slip type.
*/
int main(){
  periodic(left);
  X0 = Y0 = -L0/2;
  /**
     A loop is used to run the simulation for the three different Reynolds numbers, 
     increasing the grid resolution each iteration to a maximum of 11 levels for $Re = 80000$.
  */
  maxlevel = 9;
  for (Re = 20000; Re <= 80001; Re = Re*2.){
    init_grid (1 << 4);
    run();
    maxlevel++;
  }
}

event init(t = 0){
  vis = 1./Re;
  const face vector muc[] = {vis, vis};
  mu = muc;
  /**
     The shear layer is initialized and a small random perturbation is
     added to the span-wise-velocity-componenent field to kick-off the
     growth of the instability.
  */
  astats adapting;
  do{
    foreach() 
      u.x[] = 0.5 - (y < 0);
    boundary({u.x});
    adapting = adapt_wavelet ({u.x}, (double[]){0.01}, maxlevel);
  }while (adapting.nf);
  foreach()
    u.y[] = 0.004*noise();
 
  /**
     Output file names are defined and their names carry a
     Reynolds-number identifier.
  */
  sprintf (fnamev, "KH%g.mp4", Re);
  sprintf (fnameg, "KHlev%g.mp4", Re);
  char name1[99];
  sprintf (name1, "nrofvortices%g", Re);
  fpn = fopen (name1, "w");
}
/**
## Output
   
   The usual suspects of animation types are generated. Meaning that
   the evolution of the vorticity ($\omega$) field and the grid
   resolution are visualized. The maximum value of the vorticity field
   ($\|\omega\|_{max}$) is monitored. consequently, a vortex can be
   defined as a connected region with a vorticity value
   $\|\omega\|>\|\omega\|_{max}/3$. Based on this definition, we count
   and log the number of vortices.
*/

event output (t = 0.005; t += 0.01){
  scalar omega[], lev[], f[];
  int n = 0;
  double m = 0;
  foreach(reduction(+:n) reduction(max:m)){
    n++;
    lev[] = level;
    omega[] = (u.x[0,1]-u.x[0,-1] - (u.y[1,0]-u.y[-1,0]))/(2*Delta);
    if (fabs(omega[]) > m)
      m = fabs(omega[]);
  }
  boundary({omega});
  output_ppm (omega, file = fnamev, n = 512, min = -10, max = 10, linear = true);
  output_ppm (lev, file = fnameg, n = 512, min = 2, max = 11);
  foreach()
    f[] = (omega[] > m/3.); //1 or 0
  int nrv = tag(f);
  fprintf (fpn,"%g\t%d\t%g\t%d\t%d\t%d\n",t, nrv, m,
           n,((1 << (maxlevel*dimension)))/(n), i);
}

/**
## The adaptation and the end-of-run events

The grid resolution is adapted based on the wavelet-estimated
discretization error in the representation of the velocity vector
field. In the end-of-run event, the files associated with the run's
output file pointer is closed.
*/

event adapt (i++)
  adapt_wavelet((scalar*){u}, (double[]){0.01, 0.01}, maxlevel);


event end (t = 5.)
  fclose(fpn);

/**
## Results

One may visually inspect the evolution of the various flow set-ups
using the visualizations of the solutions and their grids.

For $Re=20000$,  

![Evolution of the vorticity field](kh/KH20000.mp4)
![Evolution of the grid structure](kh/KHlev20000.mp4)

$Re=40000$,  

![Evolution of the vorticity field](kh/KH40000.mp4)
![Evolution of the grid structure](kh/KHlev40000.mp4)

and $Re=80000$,  

![Evolution of the vorticity field](kh/KH80000.mp4)
![Evolution of the grid structure](kh/KHlev80000.mp4)

It seems as if the fastest growing mode is a function of the Reynolds
number, this is consistent with the linear perturbation analysis that
is presented in many text books on this topic. The ratio of the
wavelength associated with this mode and the periodicity length scale
($\mathcal{L}$) is related to the number of vortices that emerge. We
can plot the evolution of the number of vortices for the the different
runs.

~~~gnuplot The evolution of the number of vortices
set yr [1:50]
set xlabel 't U/L [-]'
set ylabel 'Vortex Count [-]'
set logscale y
set size square
plot 'nrofvortices20000' u 1:2 w l lw 3 t 'Re=20000' ,\
     'nrofvortices40000' u 1:2 w l lw 3 t 'Re=40000' ,\
     'nrofvortices80000' u 1:2 w l lw 3 t 'Re=80000'
~~~   

Apart from the fact that the used vortex detection algorithm
intermittently identifies decaying shear layers as vortices, these
numbers seem to correspond with what we saw in the
visualizations. Well done `tag()` function!

Finally, we check if it was sensible to use an adaptive grid by
plotting the grid-cell-compression ratio ($\Pi$) over time. This ratio
is defined as the number of grid cells required to fill the domain
with an equidistant and static grid with a cell size $\Delta_{min}$
divided by the number of grid cells employed by the adaptation
algoirithm (i.e. $\Pi \geq 1$).

~~~gnuplot The results speak for them selves (note the logaritmic axis).
set yr [1:5000]
set xlabel 'Solver Iteration'
set ylabel 'Pi'
set key left
plot 'nrofvortices20000' u 6:5 w l lw 3 t 'Re=20000' ,\
     'nrofvortices40000' u 6:5 w l lw 3 t 'Re=40000' ,\
     'nrofvortices80000' u 6:5 w l lw 3 t 'Re=80000'
~~~   
*/

   
