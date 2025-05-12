/**
# Internal Waves 
In a stratified fluid so-called internal waves can exist (also reffered to as gravity waves). An interesting feature of these waves is the so-called dispersion relation between the angle of wave propagation ($\theta$), stratification strength ($N^2$) and the freqency of the wave ($\omega$), according to, 

$$ \omega = N^2 \cos(\theta).$$

## Numerical set-up

The Navier-Stokes equantions under Boussinesq approximation are solved on a $256 \times 256$ miltigrid. In the centre of the domain an oscillating force exites the internal waves with a freqency corresponding to $\theta = 45^o$.
*/


#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "tracer.h"

scalar b[];
scalar * tracers = {b};
face vector av[];
double sqN = 1,omega=pow(2,0.5)/2;
b[top]=neumann(sqN);
b[bottom]=neumann(-sqN);

int main(){
  L0=30;
  X0=Y0=-L0/2;
  init_grid(256);
  run();
}
/**
# Initialization
We initialize the simulation with a small tolerance for the Poisson problems and a very short timestepping. This is chosen so that the pressure field ($p$) can be 'found' by the solver before the rest of the simulation is done. Note that $p=\frac{N^2}{2}y^2 + c$ with $c$ an arbitrarry constant. In the acceleration event during the 100-th iteration the timestepping and tolerance is altered to more sensible values. 
*/
event init(i=0){
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
  if (i==100){
    DT=0.05;
    TOLERANCE=1e-5;
  }
}
/**
## Output 
We output a .gif file showing the evolution of the magnitude of the gradient of the buoyancy field (|$\nabla b$|).
*/
event output(t+=0.5;t<=75){
  fprintf(ferr,"i = %d, t = %g\n",i,t);
  scalar grb[];
  foreach(){
    grb[]=0;
    foreach_dimension()
      grb[]+=sq((b[1]-b[-1])/(2*Delta));
    grb[] = pow(grb[],0.5);
  }
  static FILE * fp = popen ("ppm2gif > grbMG2.gif", "w");
  output_ppm (grb, fp, min = 0.8, max = 1.2);
}
/**
## Results
The dispersion relation appears to be statisfied. So thats good.

![Visualization of the internal waves](internalwacesMG/grbMG2.gif)

The next step is to perform this simulation using adaptive grids. [See here](internalwacesAMR.c). 
*/