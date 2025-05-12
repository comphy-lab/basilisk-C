/**
# Internalwaves Using an Adaptive Grid
In a stratified fluid so-called internal waves can exist (also reffered to as gravity waves). An interesting feature of these waves is the so-called dispersion relation between the angle of wave propagation ($\theta$), stratification strength ($N^2$) and the freqency of the wave ($\omega$), according to, 

$$ \omega = N^2 \cos(\theta).$$

## Set-up
The Navier-Stokes equantions under Boussinesq approximation are solved on an adaptive octree grid. In the centre of the domain an oscillating force exites the internal waves with a freqency corresponding to $\theta = 45^o$.
*/

#include "navier-stokes/centered.h"
#include "tracer.h"

scalar b[];
face vector av[];
double sqN = 1,omega=pow(2,0.5)/2;
b[top]=neumann(sqN);
b[bottom]=neumann(-sqN);
scalar * tracers = {b};

int main(){
  L0=30;
  X0=Y0=-L0/2;
  init_grid(256);
  run();
}
/**
## Initialization
We initialize the simulation with a small tolerance for the Poisson problems and a very short timestepping. This is chosen so that the pressure field ($p$) can be 'found' by the solver before the rest of the simulation is done. Note that $p=\frac{N^2}{2}y^2 + c$ with $c$ an arbitrarry constant. In the acceleration event during the 100-th iteration the timestepping and tolerance is altered to more sensible values. 
*/

event init(i=0){
  p.prolongation = refine_linear;
  p.refine = refine_linear;
  b.refine = refine_linear;
  TOLERANCE=1e-10;
  DT=0.000000001; 
  a=av;
  foreach()
    b[]=sqN*y;
}

event acceleration(i++){
  coord del = {0,1};
  foreach_face(){
    av.x[]= del.x*((b[]+b[-1])/2 + 0.1*(sin(omega*t)*((sq(x)+sq(y))<1)));
  }
}

event output(t+=0.5;t<=75){
  fprintf(ferr,"i = %d, t = %g\n",i,t);
  scalar grb[],ddpddy[];
  foreach(){
    ddpddy[]=(p[0,1]+p[0,-1]-2*p[])/(sq(Delta));
    grb[]=0;
    foreach_dimension()
      grb[]+=sq((b[1]-b[-1])/(2*Delta));
    grb[]=pow(grb[],0.5);
  }
  
  static FILE * fp =
    popen ("gfsview-batch2D internalwavesAMR.pres.gfv | ppm2mp4 internalwavesAMRpres2nd.mp4", "w");
  output_gfs(fp);
  fprintf(fp, " Save stdout { format = PPM width = 600 height = 600}\n");
  
  static FILE * fp2 =
    popen ("gfsview-batch2D internalwavesAMR.bf.gfv | ppm2mp4  internalwavesbf2nd.mp4", "w");
  output_gfs (fp2, list = {u,b,p,grb});
  fprintf(fp2, "Save stdout {format = PPM width = 600 height =600}\n");
}

event adapt(i++)
{
  if (t>1)
    adapt_wavelet((scalar *){u,b},(double[]){0.01,0.01,0.005},8);
  if (i==25){
    DT=0.05;
    TOLERANCE=1e-4;
  }
}
/**
## Results
The results appears fine and the grid seems to refine consistenly. One may compare the results against those obtained with a fixed equidistant grid at the maximum resolution via [this](http://basilisk.fr/sandbox/Antoonvh/internalwacesMG.c) link.  

![Magnitude of the gradient  of the buoyancy field ($\|\nabla b\|$)](internalwavesAMR/internalwavesbf2nd.mp4)

![Double vertical derrivative of the pressure field](internalwavesAMR/internalwavesAMRpres2nd.mp4)

*/