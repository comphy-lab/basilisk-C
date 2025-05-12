/**
# A Gaussian vortex

Say, 

$$ v_\theta = r e^{-r^2},$$

Then the centrifugal acceleration is

$$ a_r = \frac{v_\theta^2}{r} = r e^{-2r^2}.$$

This is balance by the pressure gradient

$$-\nabla p = a.$$

Hence, 

$$ p = -\frac{1}{4}e^{-2r^2} + C.$$

Is it true?

~~~gnuplot Pressure profile
set xr [0:2]
set yr [-0.3:0.05]
set grid
set xlabel 'r'
set ylabel 'Pressure'
plot 'pressure', -exp(-2*x**2)/4
~~~

Yes...
*/
#include "navier-stokes/centered.h"
#include "profile6.h"

int main() {
  L0 = 20;
  X0 = Y0 = -L0/2;
  N = 256;
  run();
}

event init (t = 0) {
  TOLERANCE = 1e-6;
  foreach() {
    u.y[] =  x*exp(-(sq(x) + sq(y)));
    u.x[] = -y*exp(-(sq(x) + sq(y)));
  }
}

event stop (i = 1) {
  profile ({p}, sqrt(sq(x) + sq(y)), "pressure");
}
