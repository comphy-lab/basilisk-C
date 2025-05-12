/**
# Diffusion equation and Self Similar Solution

## Theory
### Diffusion equation:

$$\frac{\partial c}{\partial t} = D \Delta c$$

with:

* $c$ the concentration field
* $D$ the diffusivity of the chemical
* $x$ & $t$ time and spatial variables


we want to solve it for this initial condition

~~~gnuplot
set output "init.png"
set xlabel 'x'
set ylabel 'c'
p[-0.5:0.5][-0.1:1.1] 'start.txt' u 2:3 t 'initial condition' w l
~~~

### Self Similar Solution:

The self similar solution is:
$$c = \frac{1}{2} \text{erfc}(\frac{\eta}{2})$$
With $\eta = xt^{-1/2}$ the self-similarity variable

This computation gives:
~~~gnuplot
reset
set output "sss.png"
set xlabel 'eta'
set ylabel 'c'
p 'bulk.txt' u ($2/(sqrt($1))):3 t 'evolution',\
   'end.txt'  u ($2/(sqrt($1))):3 t 'end',\
    0.5*erfc(x/2) t 'analytical'
~~~

## Code:
### Setup:
For a simple test case we use the diffusion-reaction solver:
*/

#include "run.h"
#include "diffusion.h"

#define PICTURES 0

//Outputs
FILE * bulk, * start, *end;

//physical Parameters
double dt;
double DT;
double tmax;
scalar A[];
const face vector da[] = {1. , 1.};
//solver stats
mgstats mgd1;

#if PICTURES
char name[30];
#endif //PICTURES

/**
Box of size 1, origin in the center, and initialization of files:
**/
int main(){
  L0= 1.;
  origin(-0.5,-0.5);
  tmax = .1;
  DT = tmax*1e-2;
  N = 256;

  bulk = fopen("bulk.txt","w");
  fclose(bulk);
  start = fopen("start.txt","w");
  fclose(bulk);
  end = fopen("end.txt","w");
  fclose(end);

  run();
}

/**
### Initial condition
**/
event init (i=0){
  // Step function for concentration field
  foreach(){
    A[] = 1*(x<0); // a gauche
  }

#if PICTURES
  sprintf(name,"A-init.png");
  output_ppm(A, file = name, min = 0., max = 1.,linear = true);
#endif //PICTURES
  // Save the initial concentration field
  start = fopen("start.txt","a");
  foreach(){
    if (y>-0.01 && y < 0.01)
      fprintf(start, "%g %g %g\n", t , x , A[]);
  }
  fclose(start);
}

/**
### Time integration
*/

event integration (i++){
  dt = dtnext(DT);

  mgd1 = diffusion(A,dt, da);
}
/**
### Outputs 
*/

#if PICTURES
event pictures (t += 5*DT, t <= tmax){
  sprintf(name,"A-%d.png",i);
  output_ppm(A, file = name, min = 0., max = 1.,linear = true);
}
#endif // PICTURES

event profiles (t += 5*DT, t<= tmax){
  bulk = fopen("bulk.txt","a");
  foreach(){
    if (y>-0.01 && y < 0.01)
      fprintf(bulk, "%g %g %g\n", t , x , A[]);
  }
  fclose(bulk);
}

event end(t = tmax){
  end = fopen("end.txt","a");
  foreach(){
    if (y>-0.01 && y < 0.01)
      fprintf(end, "%g %g %g\n", t , x , A[]);
  }
  fclose(end);
}
/** # Biblio:

* [cours PYL](http://www.lmm.jussieu.fr/~lagree/COURS/M2MHP/SSS.pdf)

*/
