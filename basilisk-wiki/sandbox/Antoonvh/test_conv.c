/**
# Do the adaptive discretizations of `nsf4t.h` work?

We test convergence by tuning the refinement criterion. The test case
is the flow field of a "Gaussian vortex";

$$v_\theta = r e^{-r^2},$$
$$v_r = 0.$$

The effective resolution $N_{\mathrm{eff}}$ is defined as:

$$N_{\mathrm{eff}} = \sqrt{\# Cells}.$$

~~~gnuplot Face-averaged to vertex values
set xr [15:512]
set logscale x 2
set logscale y 10
set xlabel 'N_{effective}'
set ylabel 'error'
set grid
set size square
fit [1:] a*x + b 'log' u (log($1)):(log($2)) via a, b
fit [1:] c*x + d 'log' u (log($1)):(log($3)) via c, d
plot 'log' u 1:2 t 'L_1', '' u 1:3 t 'L_{inf}',\
 exp(b)*x**(a) t sprintf("%.0fN^{%4.2f}", exp(b),a),\
 exp(d)*x**(c) t sprintf("%.0fN^{%4.2f}", exp(d),c)
~~~

~~~gnuplot 1st derivative (Compact vertex scheme)
fit [1:] a*x + b 'log' u (log($1)):(log($4)) via a, b
fit [1:] c*x + d 'log' u (log($1)):(log($5)) via c, d
plot 'log' u 1:4 t 'L_1', '' u 1:5 t 'L_{inf}',\
 exp(b)*x**(a) t sprintf("%.0fN^{%4.2f}", exp(b),a),\
 exp(d)*x**(c) t sprintf("%.0fN^{%4.2f}", exp(d),c)
~~~

~~~gnuplot 2nd derivate (Central vertex scheme)
fit [1:] a*x + b 'log' u (log($1)):(log($6)) via a, b
fit [1:] c*x + d 'log' u (log($1)):(log($7)) via c, d
plot 'log' u 1:6 t 'L_1', '' u 1:7 t 'L_{inf}',\
 exp(b)*x**(a) t sprintf("%.0fN^{%4.2f}", exp(b),a),\
 exp(d)*x**(c) t sprintf("%.0fN^{%4.2f}", exp(d),c)
~~~

~~~gnuplot Vertex to face-averaged values
fit [1:] a*x + b 'log' u (log($1)):(log($8)) via a, b
fit [1:] c*x + d 'log' u (log($1)):(log($9)) via c, d
plot 'log' u 1:8 t 'L_1', '' u 1:9 t 'L_{inf}',\
 exp(b)*x**(a) t sprintf("%.0fN^{%4.2f}", exp(b),a),\
 exp(d)*x**(c) t sprintf("%.0fN^{%4.2f}", exp(d),c)
~~~

~~~gnuplot Derivative on faces (projection)
fit [1:] a*x + b 'log' u (log($1)):(log($10)) via a, b
fit [1:] c*x + d 'log' u (log($1)):(log($11)) via c, d
plot 'log' u 1:10 t 'L_1', '' u 1:11 t 'L_{inf}',\
 exp(b)*x**(a) t sprintf("%.0fN^{%4.2f}", exp(b),a),\
 exp(d)*x**(c) t sprintf("%.0fN^{%4.2f}", exp(d),c)
~~~

![level](test_conv/lev.png)

The convergence is OK when the PDE has a small prefactor for the term
with second derivatives (i.e. the viscosity for the Navier-Stokes
eqs.)
 


*/
#include "nsf4t.h"
scalar * tracers = NULL;
#define funcux (-y*exp(-sq(x) - sq(y)))
#define funcuy ( x*exp(-sq(x) - sq(y)))

double funu_x (double x, double y) {
  return funcux;
}

double funu_y (double x, double y) {
  return funcuy;
}

double dfunux_x (double x, double y) {
  return 2.*x*y*exp(-sq(x) - sq(y));
}
double dfunux_y (double x, double y) {
  return (2.*sq(y) - 1)*exp(-sq(x) - sq(y));
}

double dfunuy_x (double x, double y) {
  return (1 - 2.*sq(x))*exp(-sq(x) - sq(y));
}

double dfunuy_y (double x, double y) {
  return -2*x*y*exp(-sq(x) - sq(y));
}

double d2funux_x (double x, double y) {
  return (2*y - 4.*sq(x)*y)*exp(-sq(x) - sq(y));
}

double d2funux_y (double x, double y) {
  return (6*y - 4.*cube(y))*exp(-sq(x) - sq(y));
}

double d2funuy_x (double x, double y) {
  return (4.*cube(x) - 6*x)*exp(-sq(x) - sq(y));
}

double d2funuy_y (double x, double y) {
  return (4.*sq(y)*x - 2*x)*exp(-sq(x) - sq(y));
}

/** 
The Corresponding proportional pressure field,
[see.](gaussianvortex.c)

*/
double pressure (double x, double y) {
  return exp(-2*(sq(x) + sq(y)))/4.;
}

double dp_x (double x, double y) {
  return -x*exp(-2*(sq(x) + sq(y)));
}

double dp_y (double x, double y) {
  return -y*exp(-2*(sq(x) + sq(y)));
}

double ue, expo = 1.;
int main() {
  L0 = 15.;
  X0 = -L0/2. + pi/25.;
  Y0 = -L0/2. - exp(1.)/25.;
  N = 8;
  run();
}

event check_errors (i = 0) {
  for (ue = 1e-1; ue > 1e-9; ue /= 4) {
    do {
      foreach_face() 
	u.x[] = Gauss6_x (x, y, Delta, funu_x);
      boundary ((scalar*){u});
    } while (adapt_flow (ue, 99, expo).nf > 1);
    double e1 = 0, einf = -1;
    vector v[];
    foreach_dimension() {
      v.x.prolongation = refine_vert5;
      v.x.restriction = restriction_vert;
    }
    // Face to vertex
    foreach() {
      foreach_dimension() {
	v.x[] = FACE_TO_VERTEX_4(u.x);
	double el = fabs(funu_x(x - Delta/2, y - Delta/2) - v.x[]);
	e1 += el*sq(Delta);
	if (einf < el)
	  einf = el;
	v.x[] = funu_x(x - Delta/2, y - Delta/2); //...
      }
    }
    boundary ((scalar*){v});
    vector * grads = NULL;
    for (scalar s in {v}) {
    vector dsd = new_vector ("gradient");
    grads = vectors_add (grads, dsd);
    }
    for (vector grad in grads) {
      foreach_dimension() {
	grad.x.prolongation = refine_vert5;
	grad.x.restriction = restriction_vert;
      }
    }
    compact_upwind ((scalar*){v}, grads, v);
    double u1d = 0, uinfd = -1;
    // Gradients
    foreach() {
      vector gradux = grads[0]; 
      vector graduy = grads[1]; 
      foreach_dimension() {
	double el = fabs (gradux.y[] -
			  dfunux_y (x - Delta/2., y - Delta/2));
	
	u1d += sq(Delta)*el;
	if (uinfd < el) {
	  uinfd = el;
	}
	el = fabs (graduy.y[] -
		   dfunuy_y (x - Delta/2., y - Delta/2));
	u1d += sq(Delta)*el;
	if (uinfd < el) {
	  uinfd = el;
	}
      }
    }
    double u1d2 = 0, uinfd2 = -1;
    // double derivatives
    foreach() {
      scalar s = v.x;
      foreach_dimension() {
	double el = fabs (D2SDX2/sq(Delta) -
			  d2funux_x (x - Delta/2., y - Delta/2.));
	u1d2 += sq(Delta)*el;
	if (uinfd2 < el)
	  uinfd2 = el;
      }
      s = v.y;
      foreach_dimension() {
	double el = fabs (D2SDX2/sq(Delta) -
			  d2funuy_x (x - Delta/2., y - Delta/2.));
	u1d2 += sq(Delta)*el;
	if (uinfd2 < el)
	  uinfd2 = el;
      }
    }
    double u1r = 0, uinfr = -1;
    // Vertex to face
    foreach_face() {
      double el =  fabs (Gauss6_x(x, y, Delta, funu_x)
			 - VERTEX_TO_FACE_4(v.x));
      u1r += sq(Delta)*el;
      if (uinfr < el)
	uinfr = el;
    }
    // p.prolongation = refine_5th; for 4th order L_{inf}
    foreach()
      p[] = Gauss6 (x, y , Delta, pressure);
    boundary ({p});
    double uf1 = 0, ufinf = -1;
    // Face gradient
    foreach_face() {
      double el = fabs ((p[-2,0] - 15.*p[-1,0] +
			 15.*p[0,0] - p[1,0])/(12*Delta) -
			Gauss6_x (x, y, Delta, dp_x));
      uf1 += sq(Delta)*el;
      if (ufinf < el)
	ufinf = el;
    }
    free (grads); grads = NULL;
    fprintf (stderr,"%g %g %g %g %g %g %g %g %g %g %g\n",sqrt(grid->tn),
	     e1, einf, u1d, uinfd, u1d2, uinfd2, u1r, uinfr, uf1, ufinf);
  }
  scalar lev[];
  foreach()
    lev[] = level;
  output_ppm (lev, file = "lev.png", n = 350, max = depth() + 0.5);
  return 1;
}
