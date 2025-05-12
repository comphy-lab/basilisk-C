/**
$$\int \|\nabla s\| ^2$$

We compute the the integral $A = \int_{\mathrm{domain}} \|\nabla s\| ^2 \mathrm{d}V$ of the following setups:

![5 situaltions](total_gradient/s.mp4)

Here it the data:

~~~gnuplot
set xr [0.5:5.5]
set ylabel 'A'
set xlabel 'n'
set size square
set grid
set key top left
plot 'out' t 'data' , 20*x**2 + 20 t 'parabola'
~~~
*/
#include "utils.h"

int main() {
  
  
  init_grid(256);
  for (int n = 1; n <= 5; n++) {
    scalar s[];
    foreach() {
      s[] = sin(x*2*pi*n)*cos(y*2*pi*n);
    }
    output_ppm (s, file = "s.mp4", opt = "-r 3", n = 300);
    
    double ds = 0;
    foreach(reduction(+:ds)) {
      foreach_dimension()
        ds += dv()*sq((s[1] - s[-1])/(2*Delta));
    }
    printf ("%d %g\n", n , ds);
  }
}
  
