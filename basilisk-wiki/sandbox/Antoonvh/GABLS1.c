/**
![The first GABLS intercomparison case was inspired by the simulations
 of a stable atmospheric boundary layer over the Arctic ocean by
 Kosovic and Curry (2000). Image couretesy of
 [pinterest](https://www.pinterest.com/pin/276267758371956563/).](https://i.pinimg.com/originals/5b/75/39/5b7539d06d1b9d22f73c5d278a93e3c1.jpg)

# The GABLS1 case with an adaptive-grid Single Column model

The details of the GABLS1-case set-up are presented.

## General set-up  

The so-called `binary-tree-grid' structure is used to solve the
reaction-diffusion equation. A generic timeloop function is included.
 */
#include "grid/bitree.h"
#include "diffusion.h"
#include "run.h"

/**
## Closures for turbulent transport.  

We use the closures from Louis et al. (1982), Holtslag and Boville
(1992) and England and Mcnider (1995) to parameterize the turbulent
transport. Their forms are defined below:
 */

#define fris(Ri) (sq((1-(Ri/0.20)))*(Ri<0.20))                          //Critical Ri, Short-tail mixing, England and Mcnider (1995)
//#define fris(x) (1/(1+(10*x*(1+8*x))))                                // We do not use out-dated Long tail mixing
#define friu(Ri) (sqrt(1-(18*Ri)))                                      // Holtslag en Boville 1992
#define friubm(Ri,y) ((1-((10*Ri)/(1+75*y*sqrt((x+zo/zo)*fabs(Ri))))))  // Louis 1982
#define friubh(Ri,y) ((1-((15*Ri)/(1+75*y*sqrt((x+zo/zo)*fabs(Ri))))))  // louis 1982
/**
## GABLS1 set-up  

Buoyancy is used as the thermodinamic variable. The cooling surface
buoyancy is prescribed by the case according to $0.25
\mathrm{K/hour}$. The surface roughness is set to 10 cm and the
maximum level of refinement is set to correspond to a $2^6=64$-cell
grid. Some variables are declared that will be useful later.
*/
#define bbottom (-0.25/26.5*(t/3600))
double zo = 0.1;
int maxlevel = 6;
mgstats mgb;
int nn;
double Up[100], uu[100], vv[100], bb[100];
double Cm,Ch;
int m = 0;
/**
   Among which are the fieds for the velocity components $(u,v)$ and
   the buoyancy ($b$);
*/
scalar u[], v[], b[];
/**
   We initialize a grid with $N=64$ cells and set a domain height of
   400 metres. After this the simulation starts to run.
*/

int main(){
  init_grid (1 << maxlevel);
  L0 = 400;
  X0 = 0;
  run();
}	
/**
### Boundary Conditions  

Here we set so-called no-slip the boundary conditions for the velocity
components at the bottom boundary (labelled *left*), and use the
default (stress-free) conditions for the top boundary. For the
buoyancy, we set the bottom boundary condition according to the case
description and set a Neumann boundary condition consistent with the
initialized profile (see next section) at the top boundary (labelled
*right*). The surface transport is entirely described by a
Monon-Obukhov-type closure that (in our implementation) does not
require the definition of ghostcell values (see below).
*/
u[left]  = dirichlet (0.);
v[left]  = dirichlet (0.);
b[left]  = dirichlet (bbottom);
b[right] = neumann (0.01/26.5);
/**
### Initialization  

The solution is initialized according to the prescribed initial
profiles (Cuxart et al. 2006). 
*/
event init (t = 0) {
  DT = 1.;
  foreach() {
    u[] = 8;
    v[] = 0;
    b[] = (x>100)*(0.01/26.5)*(x - 100);
  }
  boundary (all);
  while( adapt_wavelet({u, v, b},(double[]){0.25, 0.25, 0.5/26.5}, maxlevel, 3, {u,v,b}).nc){
    foreach() {
      u[] = 8;
      v[] = 0;
      b[] = (x>100)*(0.01/26.5)*(x - 100);
    }
    boundary (all);
  }
}
/**
## Time integration

   In this event, the time integration is carried out.
*/
event Diffusion (i++) {
  /**
     The tendencies due to the Coriolis force and surface transport
     (`rx, ry`) for the velocity components $(u,v)$, respectively, and
     the the tendency field for the buoyancy (`rb`) is declared.
  */
  scalar rx[], ry[], rb[], bf[];
  /**
     Face vector fields are declared for the turbulent diffusivity
     (`kh`), the gradient Richardson number(`Ri`), the stability
     correction function (`fRi`) and the surface drag coefficient
     (`CN`) is declared.
  */  
  face vector kh[], sqd[], Ri[], fRi[];
  double CN;
  /**
     On adaptive meshes, face-vectors need to be redefined explicitly
     before they appear in any computation when a cell is refined.
  */
  kh.x.refine  = no_restriction;
  sqd.x.refine = no_restriction;
  Ri.x.refine  = no_restriction;
  fRi.x.refine = no_restriction;
  /**
     We first calculate the tendency term due to the Coriolis force
     and the geostrophic forcing.
  */  
  foreach() {
    rx[] = 0.000139*v[];
    rb[] = 0;
    ry[] = 0.000139*(8 - u[]);
    /** 
	For the lowest grid cell, an additional tendency term due to
	the surface transport is calculated.
    */
    if (x < Delta) {
      /**
	 Therefore, we calculate the value of the stability correction
	 function based on the Bulk Richardson number. We distinguish
	 between stable and unstable conditions. Noting that the
	 GABLS1 case is exclusively stable.
      */
      if (b[] > bbottom){
	Cm = sq(0.4/log((x)/zo))*fris(((x - zo)*(b[] - (bbottom))/(sq(u[]) + sq(v[]))));
	Ch = Cm;
      }
      else {
	CN = sq(0.4/log((x)/zo));
	Cm = CN*friubm((x - zo)*(b[] - (bbottom))/(sq(u[]) + sq(v[])),CN);
	Ch = CN*friubh((x - zo)*(b[] - (bbottom))/(sq(u[]) + sq(v[])),CN);
      }
      /**
	 We `add' the surface fluxes to the respective tendency terms:
      */
      rx[] -= (u[]            *Cm*sqrt(sq(u[]) + sq(v[])))/Delta;
      ry[] -= (v[]            *Cm*sqrt(sq(u[]) + sq(v[])))/Delta;
      rb[] -= ((b[] - bbottom)*Cm*sqrt(sq(u[]) + sq(v[])))/Delta;
    }
  }
  /**
     A call the boundary function is needed so that all cells are in
     a locally equidistant neighbourhood (see Van Hooft et al. 2018).
  */
  boundary (all);
  /**
     The turbulent diffusivities (on cell faces) are computed below.
  */
  foreach_face() {
    /**
       We start out with evaluating the gradient Richardson number on
       each face.
    */
    sqd.x[] = (sq((u[] - u[-1])/(Delta)) + sq((v[] - v[-1])/(Delta)));
    Ri.x[] = ((b[] - b[-1])/(Delta))/(sqd.x[] + 0.00001);
    /**
       That is used to calculate the stability-correction function,
       using a different formulation for stable and unstable
       conditions.
    */
    if (Ri.x[] < 0)
      fRi.x[] = friu(Ri.x[]);
    else
      fRi.x[] = fris(Ri.x[]);
    /**
       Finally, *kh* can be evaluated.
    */
    kh.x[] = sq(min(0.4*x, 70))*(sqrt(sqd.x[]))*fRi.x[];
  }
  /**
   *kh* should be defined consistently at resolution boundaries.
   */
  boundary({kh.x});
  /**
     The time is advanced by `dt`. We log the convergence properties of
     the interative multigrid scheme that will be used later.
  */
  dt  = dtnext(DT);
  mgb = diffusion(u,dt,kh,rx);
  nn += mgb.i;
  mgb = diffusion(v,dt,kh,ry);
  nn += mgb.i;
  mgb = diffusion(b,dt,kh,rb);
  nn += mgb.i;
}
/**
## Output  

Every ten minutes we output statistics of our simulation. Most
notably, the used number of cells and profiles of the used resolution.
*/
event output (t += 360) {
  static FILE * fp2 = fopen ("GABLScells.dat","w");
  int nnn = 0;
  foreach()
    nnn++;
  fprintf (fp2, "%g\t%g\t%d\t%d\n", t, dt, i, nnn);
  fflush (fp2);
  double yp = 0;
  static FILE * fp1 = fopen ("prfileGABLS10m.dat","w");
  while (yp<400) {
    Point point = locate(yp);
    yp = x;
    fprintf (fp1,"%g\t%g\t%g\t%g\t%g\t%d\n", yp, u[], v[], b[], sqrt(sq(u[]) + sq(v[])), level);
    yp += Delta/1.5;
  }
  fflush(fp1);
  
  static FILE * fp5 = fopen("gabls1grid.dat","w");
  for (double mm=0.; mm<=400; mm+=3.125) {
    Point point = locate((double)mm);
    fprintf (fp5,"%d\t",level);
  }
  fprintf (fp5,"\n");
  fflush (fp5);
}
/**
   Furthermore, in the last hour of simulation averaged profiles are
   calculated, we do this by evaluating the solution on an equidistant
   grid using interpolation to 67 points within the domain.
*/
event avgprof(t = 8*3600; i += 20) {
  scalar U[];
  U[left] = dirichlet(0);
  int ng = 0;
  foreach() {
    U[] = sqrt(sq(u[])+sq(v[]));
    ng++;
  }
  boundary (all);
  
  static FILE * fp = fopen ("prfileGABLS.dat","w");
  double yp = 0.;
  int j = 0;
  m++;
  while (yp < 400) {
    Up[j] += interpolate (U,yp);
    uu[j] += interpolate (u,yp);
    vv[j] += interpolate (v,yp);
    bb[j] += interpolate (b,yp);
    if (t == 8*3600)
      fprintf (fp,"%g\t%g\t%g\t%g\t%g\n", yp, uu[j]/m, vv[j]/m, bb[j]/m, (Up[j]/(m)));
    j++;
    yp = yp+400./67.;
  }
  fflush (fp);
}
/**
## Adaptation
Each timestep the grid is adaptated. Furthermore, the timestep is adapted based on a *Vertrouwen-komt-te-voet-en-gaat-the-paard* strategy. The simulation is stopped when $t = 9 \mathrm{hours}$    
*/
event adapt (i++; t<=9*3600) {
  adapt_wavelet ({u,v,b}, (double[]){0.25,0.25,0.5/26.5}, maxlevel, 2, {u,v,b});
  if (nn > 14)//Quickly reduce the timestep if things get rough
    DT = max(DT/(1+((double)nn/10.)), 1.);
  if (nn < 8)//Slowly increase the timestep when time integration is easy.
    DT = min(DT*(1 + ((double)nn/100.)),15.);
}

/**
## Results

Here is a visualization of the output,:

![](http://www.basilisk.fr/sandbox/Antoonvh/gabls1.outputgabls.mp4)(width=800
 height=300)

Looks good, we can plot the used number of grid cells over time:

 ~~~gnuplot
set xr [0:9*3600]
set yr [10:25]
 
set size square
set xlabel 't [s]'
set ylabel 'cells'
plot 'GABLScells.dat' u 1:4 w lines
~~~

And for the results, we get sensible profiles for wind ($u,v$) and
buoyancy ($b$), that is now expressed as potential temperature. A
low-level jet at approx. 170 m above the surface and with a 9 m/s
magnitude is observed. If you wish to lower the low-level jet a bit,
one could decrease the value of the critical Richardson number.

![Intercomparison of the obtained average profiles over the eigtht
 hour. The grey shaded region corresponds to the $\pm \sigma$ of the
 LES results](http://www.basilisk.fr/sandbox/Antoonvh/GABLS1SCM.jpg)

The obtained results correspond very well to the LES results presented
by Beare et al. (2006) who performed a LES intercomparison study of
the same case. Cuxart et al. (2006) suggested to use the results from
these models as a benchmark.
 
## References  

Kosović, B., & Curry, J. A. "A large eddy simulation study of a
quasi-steady, stably stratified atmospheric boundary layer". Journal
of the atmospheric sciences, 57(8) (2000): 1052-1068.

Beare, Robert J., et al. "An intercomparison of large-eddy simulations
of the stable boundary layer." Boundary-Layer Meteorology 118.2
(2006): 247-272.

Cuxart, Joan, et al. "Single-column model intercomparison for a stably
stratified atmospheric boundary layer." Boundary-Layer Meteorology
118.2 (2006): 273-303.

Holtslag, A. A. M., and B. A. Boville. "Local versus nonlocal
boundary-layer diffusion in a global climate model." Journal of
Climate 6.10 (1993): 1825-1842.

Louis, J. "A short history of PBL parameterization at ECMWF." paper
presented at the Workshop on Planetary Boundary Layer
Parameterization, Eur. Cent. For Medium-Range Weather Forecasts,
Reading, England, 1982. 1982.

England, D. E. and McNider, R. T.: Stability functions based upon
shear functions, Boundary-Layer Meteorology, 74, 113–130, 1995.
*/ 
