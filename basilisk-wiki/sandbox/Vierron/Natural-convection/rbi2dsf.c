/**
# Rayleight Benard instability with stratified fluid

On this page we use the same code rbi2d.c but with stratified fluid.

## General set-up

We include as before :

*/
#include "grid/multigrid.h"
/**

Now we include :

*/
#include "tag.h"
#include "view.h"
/**

Tag.h is use for counting the number of cell which is smaller than an initial tag value (more information [here](http://basilisk.dalembert.upmc.fr/src/tag.h)). Then [view.h](http://basilisk.dalembert.upmc.fr/src/bview#interactive-basilisk-view) is a Interactive basilisk view interface.


## Model equations

*/
#include "convection_boussinesq_buoyancy.h"
#include "global_nusselt.h"
#include "profil5.h"
#include "navier-stokes/perfs.h"
#include "curvature.h"
/**

## Dimensionless parameters

Level of refinement of the adaptive mesh size 2**MINLEVEL, here the adaptive mesh is disabled because grid/multigrid.h is used.

EndTime = Max time step.

A = Initial height of the fluid layer.

B =  buoyancy coefficient.

The domain depends on the number of processors affected in order to be able to run the code in parallel.

We have an aspect ratio of one processor at y and 4 at x. 
DT variable is the maximum time step to help run the code. 

TOLERANCE corresponds to the minimum value to be reached for residues ([poisson.h](http://basilisk.dalembert.upmc.fr/src/poisson.h)).

We set Ra=1e5, Pr=1. and B=2.

*/
#define MINLEVEL 8
#define MAXLEVEL 8

double A = 0.2; 
double EndTime= 300.;

int main() {
  size (npe()); 
  origin (-0.5, -0.5);
  dimensions (ny = 1);
  DT = 0.1;
  TOLERANCE = 1e-6;
  N = 1 << MINLEVEL;
  Ra = 1e5; Pr = 1.; B = 0.75;
  run();
}
/**

## Boundary conditions 

no-slip walls and fixed temperature at the top and the bottom. In addition Neumman condition for the temperature on the side of the domain.

*/
T[top] = dirichlet(-0.5);
T[left] = neumann(0.);
T[right] = neumann(0.);
T[bottom] = dirichlet(0.5);

u.n[top] = dirichlet(0.);
u.n[bottom] = dirichlet(0.);
u.n[right] = dirichlet(0.);
u.n[left] = dirichlet(0.);
u.t[top] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);
u.t[right] = dirichlet(0.);
u.t[left] = dirichlet(0.);
/**

## Initialization

Initialization of the fluid layer and linear temperature profile with no motion.
*/
event init (t=0) {
  fraction (f, - (y - Y0 - A));
  foreach(){
    T[] = - y;
    foreach_dimension()
      u.x[] = 0.;
  }
  boundary ({T,u});
}
/**

## Outputs

Scalar div corresponds to the field of divergence. Stats s0 allows to evaluate statistical quantities on this field ([utils.h](http://basilisk.dalembert.upmc.fr/src/utils.h)).

s0.sum/s0.volume is the average value of the divergence, s0.max is the maximum value of this field.

It is important to check these quantities because the method chosen by Basilisk roughly forces the null divergence.
statsf(div) function returns the minimum, maximum, volume sum, standard deviation and volume for field div.
The mgT is the statistics on the [Poisson.h](http://basilisk.dalembert.upmc.fr/src/poisson.h) solver for the diffusion step which is in the [convection_boussinesq_buoyancy.h](http://basilisk.dalembert.upmc.fr/sandbox/Vierron/convection_boussinesq_buoyancy.h).

Creation of video file mp4 of the temperature field.

*/
event logfile (t += 1.0; t <= EndTime) {

  scalar div[];
  foreach() {
    div[] = 0.;
    foreach_dimension()
      div[] += u.y[1] - u.y[];
    div[] /= Delta;
  }

  stats s0 = statsf (div);
  static FILE * fs = fopen("stat","w");
  fprintf (fs, "%f %.9g %.9g %d %d\n",
	   t, s0.sum/s0.volume, s0.max, mgT.i, mgT.nrelax);

  output_ppm (T, file="temperature.mp4", n = 1024, box = {{-0.5,-0.5},{-0.5 + L0, 0.5}});
/**

Variable initialization nusselt and nusselt calculation ([global_nusselt.h](http://basilisk.dalembert.upmc.fr/sandbox/Vierron/global_nusselt.h)).

Creation of the data file in order to store my variables at each time step (t += 1.0; t <= EndTime).

statsf(velx and vely) returns the minimum, maximum, volume sum, standard deviation and volume for field velx and vely.

Creation of video file mp4 of the  fluid stratification.

*/
  
  double nu_vol=0. , nu_t=0. ,nu_b=0.;
  nu_vol = nusselt_vol(T,u); 
  nu_t = nusselt_top(T); 
  nu_b = nusselt_bot(T);

  static FILE * fy = fopen ("data", "w");
  if (t<1.){
    fprintf (fy, "[1]Ra [2]Pr [3]N [4]nutop [5]nubot [6]nuvol [7]umin [8]umax [9]vmin [10]vmax [11]t\n");
  }
  stats velx = statsf (u.x), vely = statsf (u.y);
  fprintf (fy, "%.9g %.9g %d %.9g %.9g %.9g %.9g %.9g %.9g %.9g %f \n",
	   Ra,Pr,N,nu_t,nu_b,nu_vol, velx.min, velx.max, vely.min, vely.max, t);

  output_ppm (f, file="strat.mp4", n = 1024, box = {{-0.5,-0.5},{-0.5 + L0, 0.5}});
/**

Creation of video file mp4 of the stratification and the temperature field.
*/
#if 1
  view(tx=+2.);
  clear();
  box ();
  draw_vof("f");
  squares ("T", spread=1., linear = true);
  save("temp-strat.mp4");
#endif
}

/**

Calculation of average temperature profile ([Antoonvh Sandbox](http://basilisk.dalembert.upmc.fr/sandbox/Antoonvh/))
saving numerical data temperature field in Temperature file.
Saving of the simulation data set at t = EndTime (dump).
*/
#if 1
event tempfile(t=EndTime){
  profile({T},-0.5,0.5,t,1);
  static FILE * fT ;
  fT= fopen("Temperature","w");
  output_field({T}, fT, N, linear = true);
  dump ();
}
#endif
/**

Counting droplets [atomisation.c](http://basilisk.dalembert.upmc.fr/src/examples/atomisation.c#counting-droplets).
*/
event droplets(t+=1.)
{
  scalar m[];
  foreach()
    m[] = f[] > 0.0001;
  int n = tag (m);

  /**
  Once each cell is tagged with a unique droplet index, we can easily
  compute the volume *v* and position *b* of each droplet. Note that
  we use *foreach_leaf()* rather than *foreach()* to avoid doing a
  parallel traversal when using OpenMP. This is because we don't have
  reduction operations for the *v* and *b* arrays (yet). */

  double v[n];
  coord b[n];
  for (int j = 0; j < n; j++)
    v[j] = b[j].x = b[j].y = 0.;
  foreach_leaf()
    if (m[] > 0) {
      int j = m[] - 1;
      v[j] += dv()*f[];
      coord p = {x,y};
      foreach_dimension()
	b[j].x += dv()*f[]*p.x;
    }

 /**
 When using MPI we need to perform a global reduction to get the
 volumes and positions of droplets which span multiple processes. */

#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, b, 3*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

  /**
  Finally we output the volume and position of each droplet to
  standard output. */
  static FILE * fd = fopen("dist","w");
  if (t<1.){
    fprintf (fd, "[1]i [2]t [3]j [4]v[j] [5]b[j].x/v[j] [6]b[j].y/v[j]\n");
  }
  for (int j = 0; j < n; j++)
    fprintf (fd, "%d %g %d %g %g %g\n", i, t,
	     j, v[j], b[j].x/v[j], b[j].y/v[j]);
  fflush (fd);
  
  /**
  Calculation maximum surface  */
  int Time = EndTime;
  int tt = t;
  double smax[Time];
  static FILE * fmax = fopen("max surface","w");
  if (t<1.){
    fprintf (fmax, "[1]t [2]smax\n");
  }
  for (int tt=0; tt< Time; tt++)
    smax[tt] = 0.;
  for (int j=0; j < n; j++)
    if (smax[tt]<v[j]){
       smax[tt]=v[j]*L0;
    }
  fprintf (fmax, "%g %g \n", t, smax[tt]);
  
 /**
First, we tag the connected bubbles;
  */
  scalar temp[];
  foreach()
    temp[] = f[] > 0.0001;
  int nb = tag (temp);
  /**
The `tag()` function identified regions:
Now we see what regions appear at the border of interest.
   */
  int * bttm = calloc (nb + 1, sizeof(int));;
  foreach_boundary(bottom)
    bttm[(int)temp[]] = 1;
  /**
The cells with the relevant tag values obtain the original volume
fraction.
   */
# if _MPI
  MPI_Allreduce (MPI_IN_PLACE, bttm, nb + 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD); 
#endif
  foreach() 
    temp[] = bttm[(int)temp[]] == 1 ? f[] : 0;
  free (bttm);
  boundary ({temp});

  output_ppm (temp, file="bottom-layer.mp4", n = 1024, box = {{-0.5,-0.5},{-0.5 + L0, 0.5}}, map = gray, min = 0., max= 1.); 
  double area =0.;
  scalar bot[];
  foreach(reduction(+:area))
      if (temp[]==f[]){
        area += temp[]*dv();
      }
  static FILE * fma = fopen("bottom","w");
  if (t<1.){
    fprintf (fma, "[1]t [2]area [3]mixture\n");
  }
  fprintf (fma, "%g %g %g\n", t, area, A*L0-area);
  fflush (fma); 

  //remove_droplets(temp,5);
#if 1
  /** 
  Calculation maximum peak to peak amplitude. Doesn't calculate exactly the peak-to-peak amplitude because of the first destabilization. Small drops of low density are trapped at the bottom of the domain */
  static FILE * fa = fopen ("amplitude", "w");
  if (t<1.){
    fprintf (fa, "[1]t [2]Amplitude [3]max [4]min [5]Ra [6]Pr [7]B [8]A\n");
  }
  scalar pos[];
  position (temp, pos, {0.,1.});
  stats s = statsf(pos);
  fprintf (fa, "%g %g %.9g %.9g %1.1e %.1f %.2f %.1f\n",
	   t, s.max - s.min, s.max, s.min, Ra, Pr, B, A);
  fflush (fa); 
#endif
}

/**

## Results


~~~gnuplot
set output "Nu=f(t).svg"
set xlabel 'times t'
set ylabel 'Nusselts'
Ra=system("awk 'NR==2 {print $1}' data")
Ra=Ra+0.
Pr=system("awk 'NR==2 {print $2}' data")
Pr=Pr+0.
B=system("awk 'NR==2 {print $7}' amplitude")
B=B+0.
A=system("awk 'NR==2 {print $8}' amplitude")
A=A+0.
set title sprintf("Ra=%1.1e Pr=%.1f B=%.2f A=%.2f",Ra, Pr, B, A)
plot 'data' u 11:4 w l title 'Nu-top',\
      '' u 11:5 w l title 'Nu-bot'
~~~

~~~gnuplot
set output "amplitude.svg"
reset
set xlabel 'times t'
set ylabel 'maximum peak to peak amplitude'
Ra=system("awk 'NR==2 {print $1}' data")
Ra=Ra+0.
Pr=system("awk 'NR==2 {print $2}' data")
Pr=Pr+0.
B=system("awk 'NR==2 {print $7}' amplitude")
B=B+0.
A=system("awk 'NR==2 {print $8}' amplitude")
A=A+0.
set title sprintf("Ra=%1.1e Pr=%.1f B=%.2f A=%.2f",Ra, Pr, B, A)
plot 'amplitude' u 1:2 w l title 'bottom-layer'
~~~

~~~gnuplot
set output "area.svg"
reset
set xlabel 'times t'
set ylabel 'area'
Ra=system("awk 'NR==2 {print $1}' data")
Ra=Ra+0.
Pr=system("awk 'NR==2 {print $2}' data")
Pr=Pr+0.
B=system("awk 'NR==2 {print $7}' amplitude")
B=B+0.
A=system("awk 'NR==2 {print $8}' amplitude")
A=A+0.
set title sprintf("Ra=%1.1e Pr=%.1f B=%.2f A=%.2f",Ra, Pr, B, A)
plot 'max_surface' u 1:2 w l title 'Max-Surface', \
     'bottom' u 1:2 w l title 'bottom-layer',\
     'bottom' u 1:3 w l title 'mixture'
~~~
## video

![Temperature field.](rbi2dsf/temperature.mp4)

![Stratification.](rbi2dsf/strat.mp4)

![Bottom layer.](rbi2dsf/bottom-layer.mp4)

## References

*/
