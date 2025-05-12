/**
# Turbulent Neutral Ekman layer

Here we test some closures in a Single-Column model model of the Neutral turbulent Ekman layer. 

We follow the LES set-up of [Qingfang Jiang et
al. (2018)](https://journals.ametsoc.org/doi/pdf/10.1175/JAS-D-17-0153.1)
somewhat.
 */

#include "grid/multigrid1D.h"
#include "diffusion.h"
#include "run.h"

#define Y_0 (0.2)
double k = 0.41;
double lut[20];

scalar u[], v[];
int maxlevel = 9;
double Ugeo =  1;
double f = 0.0001;
face vector muc[];
int main(){
  init_grid( 1 << maxlevel);
  L0 = 2000;
  run();
}

event init(t=0){
  foreach()
    u[] = Ugeo;
  lut[0] = 0.; //level > 0
  for (int m = 1; m <= 19; m++)
    lut[m] = sq(k/(log(L0/((double)(1<<m)*Y_0))-1.));
  DT = 0.1;
}

#define dvdt(u) (sign(u)*lut[level]*sq(u)/Delta)

event step(i++){
  dt = dtnext(DT);
  double S;
  scalar rx[],ry[];
  foreach_face(x){
    S = sqrt(sq((u[]-u[-1])/(Delta))+sq((v[]-v[-1])/(Delta)));
    muc.x[] = sq(min(k*x, 70))*S;
  }
  foreach(){
    rx[] = f*v[];
    ry[] = f*(Ugeo-u[]);
  }
  double ui; //scratch for the mid-point-value estimate
  foreach_boundary(left){
    ui = u[] + (dt*dvdt(u[])/2.); 
    rx[] -= dvdt(ui);
    ui = v[] + (dt*dvdt(v[])/2.);
    ry[] -= dvdt(ui);
  }
  diffusion(u, dt, muc, r = rx);
  diffusion(v, dt, muc, r = ry);
  DT *= 1.0004;
}
/**
## Results

We validate if we retrieve a log profile...
 */
event stop(t = 3600*24){
  FILE * fp = fopen("prof", "w");
  foreach()
    fprintf(fp, "%g\t%g\t%g\t%g\n", x, u[], v[], sqrt(sq(u[]+sq(v[]))));
  fclose(fp);	    
  return 1;
}

/**
Here is a plot of the profile of the absolute velocity ($U = \sqrt{u^2 + v^2}$):

~~~gnuplot
set yr [0.2: 4000]
set xr [0:1.2]
set ylabel 'Height'
set xlabel 'U'
set logscale y; set ytics 0.2, 10, 2000
set key off
set size ratio 1
plot 'prof' u 4:1 w l lw 3
~~~

The results looks OK, most notably, we have a log layer with a vanishing velocity at height$=y_0$. 

We also plot a Hodographic projection:

~~~gnuplot
unset logscale y
set ytics 0, 0.5, 1
set xtics 0, 0.1, 0.3
set yr [-0.01: 1.1]
set xr [-0.01: 0.3]
set size ratio -1
set ylabel 'u.y'
set xlabel 'u.x'

plot 'prof' u 3:2 
~~~

that is quite a funky nearsurface-velocity angle?
 */  