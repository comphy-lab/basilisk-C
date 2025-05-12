/**
# Large Eddy Simulation of isotropic Turbulence. 
On this page we aim to validate the Large Eddy Simulation (LES) formulation as defined in the [SGS.h](SGS.h) header file. Fortunately, benchmark results from a Direct Numerical Simulation (DNS) are provided under /src/examples [here](http://www.basilisk.fr/src/examples/isotropic.c). We only need to slightly alter the formulation presented there to turn the DNS into an LES. This page will highlight these changes and present the results.

##Setup
First we include the Sub-Grid-Scale formulation. Furthermore, we include a function to help visualize the $\lambda_2$ contours. 
*/

#include "grid/multigrid3D.h"
#include "navier-stokes/centered.h"
#include "SGS.h"
#include "lambda2.h"
#include "structurefunction.h"
#include "view.h"
#define BVIEW 1
#include "particles.h"
#define MU 0.01

face vector av[];
FILE * fp1;
/**
Second, Facilitated by our *filtered* approach, we reduce the resolution from $128^3$ to $32^3$.
*/
int maxlevel = 5;
int lin = 0;
int main (int argc, char * argv[]) {
  if (argc > 1)
    maxlevel = atoi(argv[1]);
  L0 = 2.*pi;
  X0 = Y0 = Z0 = -L0/2.;
  foreach_dimension()
    periodic (right);
  a = av;
  run();
}

/**
Third, a relevant value for the SGS-tuning constant is set and we tell the simulation that the molecular viscosity is the *MU* value. Additionally, we open a file to we can display the results on this page.  
*/

event init (i = 0) {
  init_grid(1 << (maxlevel - 1));
  init_particles_in_cells(); //16x16x16 particles
  init_grid(1 << maxlevel);
  N = (1 << maxlevel);
  Csmag=0.2;
  molvis=MU;
  fp1=fopen("LES.dat","w");
  foreach() {
    u.x[] = cos(y) + sin(z);
    u.y[] = sin(x) + cos(z);
    u.z[] = cos(x) + sin(y);
  }
  boundary ((scalar *){u});
}

event acceleration (i++) {
  coord ubar;
  foreach_dimension() {
    stats s = statsf(u.x);
    ubar.x = s.sum/s.volume;
  }
  foreach_face()
    av.x[] += 0.1*((u.x[] + u.x[-1])/2. - ubar.x);
}

/**
Fourth, The dissipation should be calculated with the eddy viscosity. 
*/

event logfile (i+=10; t <= 300) {
  coord ubar;
  foreach_dimension() {
    stats s = statsf(u.x);
    ubar.x = s.sum/s.volume;
  }
  double vd;
  double ke = 0., vdt = 0., vdm = 0.,vde = 0., vol = 0.;
  foreach(reduction(+:ke) reduction(+:vdt) reduction(+:vde) reduction(+:vdm) reduction(+:vol)) {
    vol += dv();
    foreach_dimension() {
      // mean fluctuating kinetic energy
      ke += dv()*sq(u.x[] - ubar.x);
      // viscous dissipation
      vd = dv()*(sq(u.x[1] - u.x[-1]) +
		  sq(u.x[0,1] - u.x[0,-1]) +
		  sq(u.x[0,0,1] - u.x[0,0,-1]))/sq(2.*Delta);
      vde+=(Evis[]-MU)*vd;
      vdm+=MU*vd;
      vdt+=Evis[]*vd;
    }
  }
  ke /= 2.*vol;
  vdm /= vol;
  vde /= vol;
  vdt /= vol;

  if (i == 0)
    fprintf (fp1,
	     "t dissipationt energy disse dissm perf.t perf.speed\n");
  fprintf (fp1, "%g %g %g %g %g %g %g\n",
	   t, vdt, ke, vde, vdm, perf.t, perf.speed);
}
/**
As a concsistency check, we evaluate the second-order structure function at fixed time intervals. Or is it really the structure funtion we aim for? 
*/
event structurefun(t+=50){
  char name2[100];
  lin++;
  sprintf(name2,"structure%d.dat",lin);
  FILE * fp2 = fopen(name2,"w");
  long2structure(fp2,u,1000,L0/2.,100);
  fflush(fp2);
  fclose(fp2);
}

/**
For our visual statisfaction, we also generate an mp4 movie of particle advection.
*/

event viewer (t = 20; t <= 100; t += 0.1) {
  view(phi = 0.5, theta = 0.25);
  clear();
  scatter(loc);
  save("LESiso.mp4");
}

/**
## Results
Movie before plots: 


![Particles!](isotropicLES/LESiso.mp4)

We also look at some statistics. The resolved kinetic energy is shown below:

~~~gnuplot Evolution of kinetic energy
set xlabel 'Time'
set ylabel 'Kinetic energy'
set logscale  y
plot 'LES.dat' u 1:3 w l lw 6 t 'Present LES' ,\
     'isotropicLES.DNSbasilisk.dat' u 1:3 w l lw 4 t 'DNS Basilisk' ,\
     'isotropicLES.hit3d.dat' u 1:($3*3./2.) w l lw 4 t 'Spectral DNS hit3d'
~~~

Compared to the spectral results, the LES 'suffers' from the early onset of turbulence, similar as the DNS run with Basilisk. The fluctuating kinetic energy after this initial transition (i.e. $t>50$) appears to be quite well captured. So atleast the implementation of the SGS formulation did not totally ruin the results. However, one may argue that energy is a rather simple quantity to get correct as it may very well be a result of the construct of the case definition. Therefore, we also look at the resolved dissipation of the simulations.       

~~~gnuplot Evolution of dissipation
set ylabel 'Dissipation'
plot  'LES.dat' u 1:2 w l lw 3 t 'Present LES' ,\
 'isotropicLES.DNSbasilisk.dat' u 1:2 w l lw 3 t 'DNS Baslisk',\
 'isotropicLES.hit3d.dat' u 1:2 w l lw 3 t 'Spectral DNS hit3d'
~~~

Seems OK. However, the initial transition shows a bit of a deveation, even between both Basilsik-based runs. We look at it in a bit more detail:

~~~gnuplot Evolution of the dissipation during the initial instability and its partitioning
set xr [10:25]
set ylabel 'Dissipation'
unset logscale y
plot  'LES.dat' u 1:2 w l lw 3 t 'LES total' ,\
'isotropicLES.DNSbasilisk.dat' u 1:2 w l lw 3 t 'DNS Baslisk' ,\
  'LES.dat' u 1:6 w l lw 3 t 'LES molecular' ,\
  'LES.dat' u 1:4 w l lw 3 t 'LES SGS'

~~~

The second-order longitudional structure function was also computed and it is plotted below.


~~~gnuplot We do observe a range with a scaling characteristic. But is is not the 2/3 law I was expexting.
reset
set  xr [0.1:5]
set  yr [0.005:0.5]
set ylabel 'dv?'
set xlabel 'Distance'
set logscale y
set logscale x  
set key bottom right
plot  'structure5.dat' u 1:2 t'data t = 200',\
'structure6.dat' u 1:2 t'data t = 250',\
'structure7.dat' u 1:2 t'data t = 300',\
  (x**(1.33))/5 w l lw 3 t 'x^{4/3} scaling'
   
~~~

We conclude that the SGS-model is not able to properly capture the dissipation dynamics within the transition. This is not a show stopper perse, more that this serves as a cautionary note;  for an accurate representation of certain flow statistics an LES is not always a suitable approach. The fact that the SGS model was able to preserve large-eddy quantities, like the kinetic energy, is motivating for its future usage. Furthermore, after the transition the turbulent strucures become larger and the spectrum more resolved. The SGS controbutions also vanish in that case. That is a nice feature.      

Finally, note that the reduction from a $128^3$ grid to a $32^3$ grid results is approximately $16^2=256$ times less effort. Meaning that instead of using 512 cores for the [original example](http://www.basilisk.fr/src/examples/isotropic.c), this simulation is reasonably fast on a single processor. So in that respect, well done LES-techniques! Also I checked the overhead of the LES formulation on my own system, using a single core. In DNS mode, I obtained an averaged performance figure of $4.6 \times 10^5$ points/sec during the first $20$ time units, and with the LES formulation this number was reduced to $3.9\times 10^5$ points/sec. Thus a reduction of 15%. These numbers were obtained with no movie output and less frequent calls to the diagnostic event.   



*/
