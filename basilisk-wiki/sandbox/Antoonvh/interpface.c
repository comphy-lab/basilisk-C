/**
# Solution estimates at faces 

We use the finite volume discretization technique for a centered
scalar field $s$. The aim is to *estimate* the corresponding solution
values at the location of a cell's face. The discretized signal and
its derivatives are:

$$s = \mathrm{cos}(x),$$
$$\frac{\mathrm{d}s}{\mathrm{d}x}=-\mathrm{sin}(x),$$
$$\frac{\mathrm{d}^2s}{\mathrm{d}x^2}=-\mathrm{cos}(x).$$

*/
#include "grid/multigrid1D.h"

scalar s[];
face vector fs2[], fs4[], fs[];
/**
We use an exact definition of the cell-averaged value of $s$ to
define `s[]`.
 */
#define sol ((sin(x + Delta/2) - sin(x - Delta/2))/Delta)

int main (){
  FILE * fp = fopen("cs" , "w"); // This will be handy later.
  FILE * fp4 = fopen("cs4", "w"); 
  FILE * fpxx = fopen("cxx", "w"); 
  L0 = 2*pi;
  periodic (left);
  for (N = 4; N <= 2048; N *= 2){
    init_grid (N);
    foreach()
      s[] = sol;
    boundary({s});
    double err20 = 0, err40 = 0, err21 = 0, err41 = 0, err22 = 0, err42 = 0;
    foreach_face(){
      // Interpolated value
      fs2.x[] = (s[] + s[-1])/2.;
      fs4.x[] = (7.*(s[] + s[-1]) - (s[-2] + s[1]))/12.;
      // The errors are weighted with the grid-cell size
      err20 += fabs(fs2.x[] - cos(x))*Delta;
      err40 += fabs(fs4.x[] - cos(x))*Delta;
      // The first derivative
      fs2.x[] = (s[] - s[-1])/Delta;
      fs4.x[] = (15.*(s[] - s[-1]) - (s[1] - s[-2]))/(12.*Delta);
      err21 += fabs(fs2.x[] + sin(x))*Delta;
      err41 += fabs(fs4.x[] + sin(x))*Delta;
    }   
    /** 
    And finally, the second derivative. Note that both methods use a four-point stencil. Can you spot the difference?
    */
    foreach_face(){ 
      fs.x[] = (fs2.x[1] - fs2.x[-1])/(2*Delta);
      fs4.x[] = (s[-2] - s[-1] - s[] + s[1])/(2.*sq(Delta));
      err22 += fabs(fs.x[] + cos(x))*Delta;
      err42 += fabs(fs4.x[] + cos(x))*Delta;
    }
    printf("%d\t%g\t%g\t%g\t%g\t%g\t%g\n", N, err20,
	   err40, err21, err41, err22, err42);
    /**
## Result

We plot the convergence data for the interpolation;

~~~gnuplot Atleast something went OK
set xr[2 : 4096]
set logscale xy 2
set xlabel 'N'
set ylabel 'L1 Error'
set size square
set key box outside 
plot 'out' u 1:2 t 'Two Cell Stencil',\
     'out' u 1:3 t 'Four Cell Stencil',\
     10*x**(-2) w l lw 3 t 'Second order',\
     100*x**(-4) w l lw 3 t 'Fourth order'
~~~

Also we can plot the convergence of the error for the first derivative;

~~~gnuplot Atleast two things went OK
set xr[2 : 4096]
set logscale xy 2
set xlabel 'N'
set ylabel 'L1 Error'
set size square
set key box outside 
plot 'out' u 1:4 t 'Two Cell Stencil',\
     'out' u 1:5 t 'Four Cell Stencil',\
     10*x**(-2) w l lw 3 t 'Second order',\
     100*x**(-4) w l lw 3 t 'Fourth order'
~~~


Finally, we check the error-convergence properties for the second
derivative;

~~~gnuplot We obtain 1st and 2nd order accuracy
set xr[2 : 4096]
set logscale xy 2
set xlabel 'N'
set ylabel 'L1 Error'
set size square
set key box outside 
plot 'out' u 1:6 t 'Naive Weights',\
     'out' u 1:7 t 'Proper Weights',\
     10*x**(-1) w l lw 3 t 'First order',\
     100*x**(-2) w l lw 3 t 'Second order'
~~~

It is tempting to conclude that the implementations are correct, but we should really be more thorough. Since we know, from dimensional analysis, an estimate of $s_i$ ($[s_i]$) with Zth order accuracy,

$$[s_i]_{Z} - s_i = c \frac{\mathrm{d}^Zs}{\mathrm{d}x^Z}\Delta^Z + \mathcal{O}\left(\Delta^{Z+1}\right)$$

We are curions to check if $c$ does indeed exist. Therefore we recall the $Z=2$ estimation of $s$ and compute $c$ at each face location, 
*/
 
    foreach_face()
      fprintf(fp, "%g\t%g\n", x, ((s[] + s[-1])/2. - cos(x)) / (-cos(x)*sq(Delta))); //Singed!
    /**
    
    ~~~gnuplot
       reset
       set xr [0:6.29]
       set yr [-0.33:0.33]
       set key off
       set grid
       set ylabel 'c'
       set xlabel 'x'
       set size square
       plot 'cs' u 1:2
    ~~~
    
    It appears that for all grids, and at virtually all locations, $c\approx 1/6$. This value can also be derived analytically from the Taylor expansion. In fact, the Taylor expansion was used to derive the heigher order estimates. There also exist some outlayers when we devide two small numbers. 
    
    We can carry out the same excersize for the 4th order approximations of $s$:
*/
    foreach_face()
      fprintf(fp4, "%g\t%g\n", x, ((7.*(s[] + s[-1]) - (s[-2] + s[1]))/12. - cos(x)) / (cos(x)*sq(sq(Delta)))); 
   /**
   
    ~~~gnuplot The theoretical value is -1/30 = -0.0333..
       reset
       set xr [0:6.29]
       set yr [-0.1:0.1]
       set key off
       set grid
       set ylabel 'c'
       set xlabel 'x'
       set size square
       plot 'cs4' u 1:2
    ~~~
    
    Finally we compute $c$ for the second order accurate estimate of the second derivative:
    */
    foreach_face()
      fprintf(fpxx, "%g\t%g\n", x, ((s[-2] - s[-1] - s[] + s[1])/(2.*sq(Delta)) + cos(x)) / (cos(x)*sq(Delta))); 
    /**
    ~~~gnuplot The theoretical value is 1/4 = 0.25
       reset
       set xr [0:6.29]
       set yr [-0.5:0.5]
       set key off
       set grid
       set ylabel 'c'
       set xlabel 'x'
       set size square
       plot 'cxx' u 1:2
    ~~~
    */
  }
}