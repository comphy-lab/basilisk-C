/**
# A test for the axisymmetric hydrostratic balance

When $a_r = r$, then in a hydrostatic balance, $\frac{\mathrm{d}p}{\mathrm{d}r}= r$. So:

$$p(r) = c+\frac{1}{2}r^2$$


*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "profile5c.h"

face vector av[];
p[top] = neumann(a.y[]); //Boundary conditions are important

int main(){
  a = av;
  N = 32;
  DT = 0.01;
  run();
}

event acceleration(i++){
  foreach_face(y)
    av.y[] = y;
}

event pres(t = 5)
  profile({p}, fname = "prof");

/**
~~~gnuplot Close Enough
ftitle(b,c,d) = sprintf("Fit: %f^2+%4.4f*x+%4.2f", b, c, d)
f(x) = b*x**2+c*x+d
fit f(x) 'prof' u 1:2 via b,c,d
set xlabel 'r'
set ylabel 'Pressure'
set key top left box
set size square
plot 'prof' u 1:2 t 'Data', f(x) t ftitle(b,c,d)
~~~

The fit convergences to the analytical solution for larger N. 
 */
