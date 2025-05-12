/**
# An axisymmetric free surface

The free surface of a fluid in solid body rotation is affected by
gravity. The surface height is a function of the radial distance from
the centre ($r$);

$$z_i = \frac{\omega^2}{2g}r^2 + c,$$

where $\omega$ is the angular frequency (in $\mathrm{Rad./s}$) and $g$
the acceleration due to gravity. 
 */

#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
//p[top] = neumann(a.y[]); // This does not ensure no penetration

uf.n[top] = 0;
#include "tracer.h"
#include "diffusion.h"
scalar ut[];
scalar * tracers = {ut}; 
face vector av[];
event acceleration(i++){
  foreach_face(y)//Centrifugal
    av.y[] += y*sq((ut[] + ut[0,-1]) / 2.);
  foreach_face(x)//Gravity
    av.x[] -= 1.;
}

event tracer_diffusion(i++)
  diffusion(ut, dt, mu);

double omega = 0.5;
ut[left] = dirichlet(omega);
ut[top] = dirichlet(omega);

int main(){
  a = av;
  N = 32;
  TOLERANCE = 10E-4;
  DT = 0.01;
  mu1 = 0.1;
  mu2 = 0.1;
  rho1 = 2;
  rho2 = 1;
  run();
}
/**
We initialize a fluid in solid body rotation with the corresponding
free surface.
 */
event init(t = 0){
  fraction(f, 0.25 + 0.5*sq(omega)*sq(y)  - x);
  foreach()
    ut[] = omega;
  boundary(all);
}

event free_surface(t = 25){
  FILE * fp = fopen("facets", "w");
  output_facets(f, fp = fp);
}
  
/**
~~~gnuplot The initialized interface is maintained-ish
ftitle(b,c,d) = sprintf("%4.3fr^2+%4.3fr+%4.3f", b, c, d)
f(x) = b*x**2+c*x+d
fit f(x) 'facets' u 2:1 via b,c,d
set xlabel 'r'
set ylabel 'z'
set key top left
set size ratio -1
plot 'facets' u 2:1 t 'Interface facets', f(x) t ftitle(b,c,d)
~~~

Noting that $0.132 \neq 0.125$

## Another test

* [Axisymmetric hydrostatic pressure](axihsb.c)
 */
