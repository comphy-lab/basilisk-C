/**
# Lamb-Chaplygin dipole collision with a wall using second order accurate interpolation methods
We can test the effect of applying second order accurate interpolation strategies using a simulation case nearly identical to an [earlier version](lamb-dipole.c). You can view these methods [here](additional-attributes.h).

## Numerical Set-up
Only minor things need to be changed compared to the version using the default settings. First we need to include the file that contains the definitions for the new attributes.  
*/

#include "navier-stokes/centered.h"
#include "additional_attributes.h"

scalar omg[],r[],s[],psi[];
int maxlevel = 10;
double xo= 1, yo= 0.00312;
double k=3.83170597;
double temp= 12;
double vis = 1./500.;
const face vector muc[] = {vis,vis}; 
u.t[right]=dirichlet(0);

int main()
{
  L0=10;
  X0=Y0=-L0/2;
  init_grid(1<<(8));
  run();
}
/**
## Initialization
We can set the restriction, coarsen, prolongation and refine attributes of the fields we solve for to the second order accurate formulations.  
*/
event init(t=0)
{
  # if 1
  u.x.restriction=quadratic_coarsening;
  u.x.coarsen=quadratic_coarsening;
  u.x.prolongation=refine_quadratic;
  u.x.refine=refine_quadratic;
  
  u.y.restriction=quadratic_coarsening;
  u.y.coarsen=quadratic_coarsening;
  u.y.prolongation=refine_quadratic;
  u.y.refine=refine_quadratic;
  
  p.restriction=quadratic_coarsening;
  p.coarsen=quadratic_coarsening;
  p.prolongation=refine_quadratic;
  p.refine=refine_quadratic;
  # endif
  mu=muc;
  refine(pow(pow((x-xo),2)+(pow((y-yo),2)),0.5) < 1.5 && level<11);
  foreach() {
    r[] = pow(pow((x-xo),2)+(pow((y-yo),2)),0.5);
    s[] = (y-yo)/r[];
    psi[] = ((r[]>1)*((1/r[]))*s[]) + ((r[]<1)*((-2*j1(k*r[])*s[]/(k*j0(k)))+(s[]*r[])));
  }
  boundary({psi});
  foreach() {
    u.x[] = ((psi[0,1]-psi[0,-1])/(2*Delta));
    u.y[] = -(psi[1,0]-psi[-1,0])/(2*Delta);
  }
  boundary(all);
}
/**
## Adaptation
We do not wish to use the second order accurate formulations for the "adapt_wavelet" function as that this would result in a different grid structure. Therefore, we define two dummy scalar fields uu[] and vv[] that use the default attributes such that the wavelet error estimation is not altered. Note that potentially the second order accurate attributes can be usefull for error estimation as well.        
*/
event adapt(i++)
{
  scalar uu[],vv[];
  foreach()
    {
      uu[]=u.x[];
      vv[]=u.y[];
    }
  boundary({uu,vv});
  adapt_wavelet((scalar *){uu,vv},(double []){0.01,0.01},maxlevel); 
}
double aa=0.1;
event gfsview(t+=aa;t<=temp)
{
  double pz=0;
  double en=0;
  int n=0;
  foreach() 
    omg[]=(u.x[0,1]-u.x[0,-1] - (u.y[1,0]-u.y[-1,0]))/(2*Delta);
  
  foreach(reduction(+:pz) reduction(+:en) reduction(+:n)) {
    pz+=pow(omg[],2)*Delta*Delta;
    n++;
    foreach_dimension()
      en+=sq(u.x[])*Delta*Delta;
  }
  fprintf(ferr,"%d\t%g\t%g\t%g\t%d\n",i,t,pz,en,n);
  
  #if 0
  static FILE * fp =
    popen ("gfsview-batch2D lamb.gfv| ppm2gif > lamb4000.gif", "w");
  output_gfs (fp);
  fprintf (fp, "Save stdout { format = PPM width = 640 height = 480}\n");

  static FILE * fp2 =
    popen ("gfsview-batch2D lambinzoom.gfv | ppm2gif > lamb4000inz.gif", "w");
  output_gfs (fp2);
  fprintf (fp2, "Save stdout { format = PPM width = 1000 height = 500}\n");
  #endif
    aa=0.1/(pz/100);
}

/**
## Results and Performance
It appears that the results become more accurate using the new second order accurate formulated attributes. We can for example compare the total resolved enstrophy ($\Omega$).

![Evolution of enstrophy and number of grid points](http://www.basilisk.fr/sandbox/Antoonvh/2ndordertestlamb.jpg)

The performance appears to be affected as this is what I retrieved from compiling with the "-DTRACE=2" tag

Using the default attributes:

~~~bash
# Quadtree, 3193 steps, 154.649 CPU, 154.7 real, 3.99e+05 points.step/s, 20 var
# 4 procs, MPI: min 38 (25%) avg 41 (27%) max 45 (29%)
   calls    total     self   % total   function
   83053    39.79    28.41     18.4%   boundary():/home/antoon/basilisk/src/grid/cartesian-common.h:274
    6388    48.14    27.48     17.8%   project():/home/antoon/basilisk/src/poisson.h:418
    3194    35.90    22.76     14.7%   viscous_term():/home/antoon/basilisk/src/navier-stokes/centered.h:295
 4062230    16.04    16.04     10.4%   mpi_boundary_level():/home/antoon/basilisk/src/grid/tree-mpi.h:421
 4033794    12.94    12.94      8.4%   mpi_boundary_restriction():/home/antoon/basilisk/src/grid/tree-mpi.h:428
    3202    19.77    12.41      8.0%   balance():/home/antoon/basilisk/src/grid/balance.h:390
    6416     7.87     7.87      5.1%   mpi_boundary_update_buffers():/home/antoon/basilisk/src/grid/tree-mpi.h:929
   31662     8.82     7.56      4.9%   mpi_boundary_coarsen():/home/antoon/basilisk/src/grid/tree-mpi.h:1107
    3194    38.90     6.94      4.5%   advection_term():/home/antoon/basilisk/src/navier-stokes/centered.h:252
   77728     3.95     3.95      2.6%   mpi_all_reduce0():/home/antoon/basilisk/src/common.h:580
    3194    43.54     2.96      1.9%   adapt_wavelet():/home/antoon/basilisk/src/grid/tree-common.h:293
    3194    28.46     1.66      1.1%   projection():/home/antoon/basilisk/src/navier-stokes/centered.h:352
    3202     1.72     1.56      1.0%   z_indexing():/home/antoon/basilisk/src/grid/tree-mpi.h:1389
    .....
    And some more lines with less than 1% total time
~~~
And With the second order accurate attributes:

~~~bash
# Quadtree, 3266 steps, 183.536 CPU, 183.6 real, 3.61e+05 points.step/s, 20 var
# 4 procs, MPI: min 41 (22%) avg 44 (24%) max 46 (25%)
   calls    total     self   % total   function
    6534    64.33    39.21     21.4%   project():/home/antoon/basilisk/src/poisson.h:418
   87943    46.46    34.36     18.7%   boundary():/home/antoon/basilisk/src/grid/cartesian-common.h:274
    3267    45.19    29.45     16.0%   viscous_term():/home/antoon/basilisk/src/navier-stokes/centered.h:295
 4681147    17.31    17.31      9.4%   mpi_boundary_level():/home/antoon/basilisk/src/grid/tree-mpi.h:421
 4652085    15.27    15.27      8.3%   mpi_boundary_restriction():/home/antoon/basilisk/src/grid/tree-mpi.h:428
    3276    19.98    12.83      7.0%   balance():/home/antoon/basilisk/src/grid/balance.h:390
   32372     8.95     7.96      4.3%   mpi_boundary_coarsen():/home/antoon/basilisk/src/grid/tree-mpi.h:1107
    6563     7.83     7.83      4.3%   mpi_boundary_update_buffers():/home/antoon/basilisk/src/grid/tree-mpi.h:929
    3267    40.08     7.23      3.9%   advection_term():/home/antoon/basilisk/src/navier-stokes/centered.h:252
   81029     3.50     3.50      1.9%   mpi_all_reduce0():/home/antoon/basilisk/src/common.h:580
    3267    44.86     3.05      1.7%   adapt_wavelet():/home/antoon/basilisk/src/grid/tree-common.h:293
    3267    45.33     1.76      1.0%   projection():/home/antoon/basilisk/src/navier-stokes/centered.h:352
    3276     1.77     1.59      0.9%   z_indexing():/home/antoon/basilisk/src/grid/tree-mpi.h:1389
    .....
    And some more lines with a low impact on performance
~~~
These results were obtained using my laptop. 

##Further readings
view this topic on the user forum:
[https://groups.google.com/forum/#!topic/basilisk-fr/16E-uYkeF4s](https://groups.google.com/forum/#!topic/basilisk-fr/16E-uYkeF4s) and the links therein.
*/

