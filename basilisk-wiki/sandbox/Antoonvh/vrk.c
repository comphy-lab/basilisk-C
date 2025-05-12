/**
# A Kauffman-vortex ring

We test a naive vortex-ring-initialization method and characterize the
resulting structure a bit. It can be compared against a ring strucutre
that emergeces from [injection](vri.c)
 */

#include "axi.h"
#include "navier-stokes/centered.h"

double major_rt = 1, Gamma, R1, compact = 6; 
double xo = 3;

#define KAUFMANN (d*Gamma/((sq(R1) + sq(d))*2*pi)) 

int maxlevel = 9;
int main() {
  L0 = 15;
  R1 = major_rt/compact;
  Gamma = 0.5*major_rt*2.*pi; // V = 0.5
  foreach_dimension()
    u.x.refine = refine_linear;
  run();
}

coord cross (coord a, coord b) {
  coord product;
  foreach_dimension(3) 
    product.x = a.y*b.z - a.z*b.y;
  return product;
}

#define radi (sqrt(sq(x - xo) + sq(y - major_rt)))
event init (t = 0) {
  refine (radi < 3*R1 && level < maxlevel - 1);
  refine (radi < 1.5*R1 && level < maxlevel);
  // A non-solenoidal flow field
  foreach() {
    coord c = {x, y, 0}, pr = {xo, major_rt, 0}; // vortex center
    double d = 0;
    foreach_dimension()
      d += sq(c.x - pr.x);
    d = sqrt(d);
    double vtheta = KAUFMANN;
    coord rc_to_pc = {x - pr.x, y - pr.y, 0}; 
    coord n = cross ((coord){0, 0, 1}, rc_to_pc);
    normalize (&n);
    foreach_dimension()
      u.x[] += vtheta*n.x;

    pr = (coord){xo, -major_rt, 0}; // vortex center
    d = 0;
    foreach_dimension()
      d += sq(c.x - pr.x);
    d = sqrt(d);
    vtheta = KAUFMANN;
    rc_to_pc = (coord){x - pr.x, y - pr.y, 0}; 
    n = cross ((coord){0,0,-1}, rc_to_pc);
    normalize (&n);
    foreach_dimension()
      u.x[] += vtheta*n.x;
  }
  boundary ((scalar*){u});
}

event adapt (i++)
  adapt_wavelet ((scalar*){u}, (double[]){0.005, 0.005}, maxlevel);

event stop (t = 15);

/**
## Output and analysis

A movie is generated, showing the vorticity field, stream lines,
vortex center location and grid cells.

![Coherent and steady](vrk/mov.mp4)
 */

#include "view.h"
int n_part = 1;
#include "scatter.h"
#include "axistream.h"
#define vertex_value(s) ((s[] + s[-1] + s[0,-1] + s[-1,-1])/4.)

coord pos1;
coord position (void) {
  double ax = 0, ay = 0, tot = 0;
  foreach(reduction(+:ax) reduction (+:ay)  reduction (+:tot)) {
    double omg = (u.y[1] - u.y[-1] - u.x[0,1] + u.x[0,-1])/(2.*Delta);
    if (omg > 1) {
      ax += omg*x*dv();
      ay += omg*y*dv();
      tot += omg*dv();
    }
  }
  ax /= tot > 0 ? tot : HUGE;
  ay /= tot > 0 ? tot : HUGE;
  return (coord){ax, ay};
}

event movie (t += 0.25) {
  scalar omg[], psi[];
  psi[left] = dirichlet (0.);   //not true..
  psi[right] = dirichlet (0.);
  psi[bottom] = dirichlet (0.); //?
  foreach()
    omg[] = (u.y[1] - u.y[-1] - u.x[0,1] + u.x[0,-1])/(2.*Delta);
  boundary ({omg});
  axistream (u, psi);
  vertex scalar phi[];
  foreach_vertex()
    phi[] = vertex_value(psi);
  boundary ({omg, phi});
  coord a = position();
  
  view (fov = 15, tx = -0.35, ty = 0.1);
  for (double val = -0.5; val <= 0.5; val += 0.05) 
    isoline ("phi", val);
  squares ("omg", map = cool_warm, min = -1, max = 1, linear = true);
  mirror ({0,1}) {
    scatter (&a);
    squares ("-omg", map = cool_warm, min = -1, max = 1); // Neat!
    translate (y = 4)
      cells();
  }
  save ("mov.mp4");
}

/**
I am also curious about the shape of the separatrix. 

~~~gnuplot Fit by an ellipse by eye.

set xr [-2:2]
set yr [0:3]
set xlabel 'z'
set ylabel 'r'
set size ratio -1
set grid
set obj ellipse size 2.6,3.6 front fc lt 2 lw 3
plot 'out' t 'Separatrix cloud data',\
       NaN t 'Ellipse Rmin = 1.3, Rmax = 1.8' with ellipses  lt 2 lw 3
~~~

The exentricity is $\approx$ 3.6/2.6 = 1.38. Less than the injected
vortex ring.

The Radius of the vortex centre line and translation velocity

~~~gnuplot Speed and Size
reset
set xlabel 't'
set ylabel 'R_c and V'
set size square
set grid
set key top left
plot 'log' w l lw 2 t 'R_c', 'log' u 1:3 w l lw 2 t 'V' 
~~~
 */

event analyze (t = 5; t += 0.5) {
  coord pos2 = position();
  if (pos1.x == 0) 
    pos1 = position();
  else {
    double Vz = pos2.x - pos1.x;
    fprintf (stderr, "%g %g %g\n", t, pos2.y, Vz);
    scalar psi[];
    psi[left] = dirichlet (0.);  
    psi[right] = dirichlet (0.);
    psi[bottom] = dirichlet (0.);
    axistream (u, psi);
    foreach() 
      psi[] += 0.5*sq(y)*Vz;
    boundary ({psi});
    vertex scalar phi[];
    foreach_vertex()
      phi[] = vertex_value(psi);
    boundary ({phi});
    scalar c[];
    fractions (phi, c);
    c.prolongation = fraction_refine;
    boundary ({c});
    face vector s; s.x.i = -1;
    foreach()
      if (c[] > 1e-6 && c[] < 1. - 1e-6) {
	coord n = facet_normal (point, c, s);
	double alpha = plane_alpha (c[], n);
	coord segment[2];
	if (facets (n, alpha, segment) == 2)
	  fprintf (stdout, "%g %g\n%g %g\n\n", 
		   x + segment[0].x*Delta - pos2.x, y + segment[0].y*Delta, 
		   x + segment[1].x*Delta - pos2.x, y + segment[1].y*Delta);
      }
    /**
       For steady vortex structures, $\omega = f(\psi)$ in the
       co-moving frame. we check if this is indeed the case.

       ~~~gnuplot The omega-psi relation appears
       reset
       set size square
       set xlabel 'psi'
       set ylabel 'omega'
       set grid
       f(x) = a*(-x)**b
       b = 3.5
       a = 80
       fit f(x) 'omgpsi' u ($1<0?$1:NaN):($1<0?$2:NaN) via a,b
       plot 'omgpsi' u ($1<0?$1:NaN):2 t 'data',\
	f(x) lw 3 t sprintf("omg = %.2g*psi^{%.2g}", a ,b)
       ~~~

The exponent is larger than the injected vortex ring.
    */
    if (t > 7.99 && t < 8.01) {
      FILE * fp = fopen ("omgpsi", "w");
      foreach() {
	double omg = (u.y[1] - u.y[-1] - u.x[0,1] + u.x[0,-1])/(2.*Delta);
	if (fabs(omg) > 0.01) 
	  fprintf (fp, "%g %g\n", psi[], omg);
      }
	/**
A transact of $u_z(r)$ accros the vortex center is also plotted. 

~~~gnuplot
reset
set xr [-1.:2.1]
set yr [0:4]
set parametric
set tr [0:4]
c = 0
set xlabel 'u.z'
set ylabel 'r'
set size square
set grid
set key off
plot 'prof' w l lw 2, c,t lw 2
~~~
	 */
      FILE * fp2 = fopen ("prof", "w");
      for (double rp = 0; rp < 4; rp += 0.01) 
	fprintf (fp2, "%g %g\n", interpolate (u.x, pos2.x, rp), rp);
      fclose (fp);
      fclose (fp2);
    }
    pos1 = pos2;
}
}
    

