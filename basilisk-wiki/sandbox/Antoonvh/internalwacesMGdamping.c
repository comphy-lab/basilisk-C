/**
# Damping Internal Waves In a stratified fluid  
so-called internal waves can exist (also reffered to as gravity waves). An interesting
feature of these waves is the so-called dispersion relation between
the angle of wave propagation ($\theta$), stratification strength
($N^2$) and the freqency of the wave ($\omega$), according to,

$$ \omega = N^2 \cos(\theta).$$

This page aims to damp these waves somewhat in the upper part of the
domain so that they do not reflect at the top boundary. The set-up is enherited from this [page](interwacesMG.c) on internal waves.

## Numerical set-up
The Navier-Stokes equantions under Boussinesq approximation are solved
on a $256 \times 256$ miltigrid. In the centre of the domain an
oscillating force exites the internal waves with a freqency
corresponding to $\theta = 45^o$.
*/
#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"
/**
A parabolicly increasing damping stength with height in the top part
of the domain is defined in the line below. Notice that $-15<y<15$ in
this set-up.
*/
#define Damp (-(y>10)*((y-10)/2.5)*((y-10)/2.5))

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
## Initialization

We initialize the simulation with a small tolerance (`TOLERANCE`)
for the Poisson problems and a very small timestep(`DT`). This is chosen
so that the pressure field ($p$) can be "found" by the solver before
the rest of the simulation is done. Note that $p=\frac{N^2}{2}y^2 + c$
with $c$ an arbitrarry constant. In the acceleration event during the
100-th iteration the timestepping and tolerance is altered to more
sensible values.
*/
event init(i=0){
  TOLERANCE=1e-9;
  DT=0.00000001; 
  a=av;
  foreach()
    b[]=sqN*y;
}

event acceleration(i++){
  coord del = {0,1};
  foreach_face(){
    av.x[]= del.x*((b[]+b[-1])/2. + 0.1*(sin(omega*t)*((sq(x)+sq(y))<1)) + (Damp*(u.x[]+u.x[-1])/2.)); // <- The last term concerns the damping of velocity fluctuations
  }
  if (i==100){
    DT=0.05;
    TOLERANCE=1e-5;
  }
}
/**
## Buoyancy fluctuation damping
In the samping layer, the buoyancyfield is nudged towards the initialized stratification.
*/
event tracer_diffusion(i++){
  scalar damping[];
  foreach()
    damping[]=Damp*(b[]-sqN*y);
  boundary({damping});
  diffusion(b,dt,zerof,damping);
}

/**
## Output 
We output two movie files (.mp4) with 200 frames each, that show the evolution of the magnitude of the
gradient of the buoyancy field (i.e. |$\nabla b$|) and the vertical velocity. 
*/
event output(t+=0.5;t<=100){
  fprintf(ferr,"i = %d, t = %g\n",i,t);
  scalar grb[];
  foreach(){
    grb[]=0;
    foreach_dimension()
      grb[]+=sq((b[1]-b[-1])/(2*Delta));
    grb[] = pow(grb[],0.5);
  }
  static FILE * fp = popen ("ppm2mp4 grbMG.mp4", "w");
  output_ppm (grb, fp,n=512, min = 0.8, max = 1.2);
  static FILE * fp2 = popen ("ppm2mp4 uyMG.mp4", "w");
  output_ppm (u.y, fp2,n=512, min = -0.1, max = 0.1);
}
/**
## Results  
We can view the results and compare the behaviour of the solution
near the damped (upper) and undamped (lower) boundary.

![Visualization of the internal waves (buoyancy structure)](internalwacesMGdamping/grbMG.mp4)

![Visualization of the internal waves (vertical velocity)](internalwacesMGdamping/uyMG.mp4)

The results seem statisfactory as the waves appear to be damped before
reflection occurs. Notice however that the exact shape of the `Damp`
function defined above is rather ad hoc and may not translate well to
other problems. Furthermore, it is not clear how the damping itself
influences the solution in the inner part of the domain.
*/