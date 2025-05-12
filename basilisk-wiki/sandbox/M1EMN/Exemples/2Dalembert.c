/**
# Two dimensional Wave Equation


This example is a simple extension of the tutorial `bump` 
and aims to compute a simple 2D cylindrical wave flow with an initial small elevation.
We use Saint-Venant, 
whose linearisation will give the wave equation:
$$\frac{\partial^2 h}{\partial t^2} - c_0^2 (\frac{\partial^2 h}{\partial x^2}+\frac{\partial^2 h}{\partial y^2})$$

with suppose at $t=0$ a Gaussian bump $h =h_0 (1+ a g(x,y))$ with 
$g(x,y)=e^{-(x^2+y^2)}$ and no initial flow. 

 ![animation of the propagation](2Dalembert/eta.mp4)
 
 
## Code

Code is almost the tutorial `bump` ([http://basilisk.fr/Tutorial]())

*/ 

#include "saint-venant.h"
#define LEVEL 9
double tmax;

/**
A small Gaussian bump, a larger one gives slope steepening (a=0.05) and even shock (a=0.5)
*/
event init (t = 0) {
  double a = .025, b = 1.;
  foreach()
      h[] = 1 + a * exp(- b*(x*x + y*y));
}
/**
Domain is large enough in space and time (changed compared to tutorial) 
*/
int main() {
  L0 = 30;
  tmax = 10;
  DT = HUGE;
  origin (-L0/2, -L0/2);
  init_grid (1 << LEVEL);
  run();
}

/**
Output with a gnuplot flag to pipe in, or animations for web interface 
*/
#ifdef gnuX
event output (t += .5; t < tmax) {
  fprintf (stdout, "p[-5:5][:]  '-' u 1:2 ,'' u 1:4 ,'Q.txt' u (-$1):($2*.025/2/sqrt(2)) not  w   l \n");  
   double dx= L0/pow(2,LEVEL);	   
    for (double x = dx ; x < L0/2; x += dx)
    fprintf (stdout, "%g %g %g %g %g\n", x-t, 
    	  (interpolate (h, x, 0)-1)*sqrt(x), sqrt(x)*interpolate (u.x, x, 0.), sqrt(x)*(interpolate (h, x/sqrt(2), x/sqrt(2), t)-1), t);
    fprintf (stdout, "e\n\n");
}
#else 
event graphs (i++) {
  stats s = statsf (h);
  fprintf (stderr, "%g %g %g\n", t, s.min, s.max);
}

event images (t += 4./300.) {
#if 0
  output_ppm (h, linear = true);
  scalar l[];
  foreach()
    l[] = level;
  static FILE * fp = fopen ("grid.ppm", "w");
  output_ppm (l, fp, min = 0, max = LEVEL);
#endif
  output_ppm (h, file = "eta.mp4", n = 512,  linear = true);
}

#endif 

/**
Mesh adaptaion
*/ 
event adapt (i++) {
  adapt_wavelet ({h}, (double []){4.e-6}, maxlevel = LEVEL);
}

/**

A file to monitor during computation 
*/ 
event printdata (t +=1 ) {
    FILE *  fpx = fopen("cut.txt", "w");
    scalar l[];
    foreach()
     l[] = level;
    double dx = L0/256;
    for (double x = dx; x < L0/2; x += dx){
        fprintf (fpx, "%g %g %g %g %g %g \n",
                 x, (interpolate (h, x, 0)), (interpolate (h, x/sqrt(2), x/sqrt(2))), interpolate (u.x, x, 0.),(interpolate (l, x, 0)),t);
        }
    fclose(fpx);
}

/** 
## Results



The final signal $h(x,t_{max})$ at end of computation

~~~gnuplot
set key right 
set xlabel "x"
set ylabel "h(x,t)" 
p'cut.txt' u 1:2 t '0',''u 1:3 t '\pi/4'
~~~



The asymptotic solution is $h(x,t)=\dfrac{1}{\sqrt{x}}q(x-c_0t)$
where $$q(t)=\int_{-\infty}^{t}\frac{g'(t)}{2 \sqrt{2\pi}\sqrt{t-T}}dT$$
(half derivative) 

We plot on to the ray $(x,0)$ and the ray $\pi/4$, we plot the velocity.
They are all superposed with the analitical solution 

~~~gnuplot
set key left
set xlabel "x-t"
set ylabel "(h(x-t)-h_0)(r)^{1/2}" 
a=0.025
p[-5:5][:]'cut.txt' u ($1-$6):(($2-1)*sqrt($1)/a) t'0','' u ($1-$6):(($3-1)*sqrt($1)/a) t'pi/4' ,'' u ($1-$6):(($4)*sqrt($1)/a) t'u','Q.txt' u (-$1):($2/2/sqrt(2)) t'Lighthill 78'  w   l
~~~



## Links

* [http://basilisk.fr/Tutorial]()


## References

* James Lighthill, Waves in Fluids, Cambridge Univ Press 1978, page 21
* G. B. Whitham, Linear and Nonlinear Waves, Wiley-Interscience 1999, page 222



## Annex
   
   Values of 
   $$\int_{-\infty}^{t}\frac{-2 T e^{-T^2}}{\sqrt{\pi}\sqrt{t-T}}dT$$
*/


event end (t = tmax) {
     FILE *  fpx = fopen("Q.txt", "w");
     fprintf (fpx,"-2. 0.0371514 \n");
     fprintf (fpx,"-1.8 0.075587 \n");
     fprintf (fpx,"-1.4 0.242043 \n");
     fprintf (fpx,"-1.  0.54     \n");
     fprintf (fpx,"-0.8 0.71     \n");
     fprintf (fpx,"-0.6 0.84     \n");
     fprintf (fpx,"-0.4 0.898    \n");
     fprintf (fpx,"-0.2 0.848    \n");
     fprintf (fpx,"-0.2 0.848695 \n");
     fprintf (fpx," 0. 0.691367  \n");
     fprintf (fpx," 0.2 0.452983 \n");
     fprintf (fpx," 0.4 0.181756 \n");
     fprintf (fpx," 0.6 -0.06931 \n");
     fprintf (fpx," 0.8 -0.25990 \n"); 
     fprintf (fpx," 1. -0.372496 \n"); 
     fprintf (fpx,"1.2 -0.412126 \n"); 
     fprintf (fpx,"1.6 -0.354543 \n"); 
     fprintf (fpx,"1.8 -0.30061  \n");
     fprintf (fpx,"2. -0.249051  \n");  
     fprintf (fpx,"3. -0.109932  \n");  
     fclose(fpx);    
}
