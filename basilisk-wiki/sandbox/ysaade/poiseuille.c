/**
# Starting Poiseuille flow in a cylindrical pipe
![Axial velocity.](poiseuille/mov.mp4)(width="1200" height="150")

At $t = 0$, a gradient of pressure $G$ is applied between both ends of
a cylindrical pipe, the fluid in which is initially at rest. The
converged solution has the well-known Poiseuille velocity profile. The
aim of this test case, however, is to check the transient development
of the flow, and to validate it against the transient analytical
solution.
*/

#include "grid/multigrid.h"
#include "axi.h"
#include "navier-stokes/centered.h"
#include "bessel-roots.h"

/**
We will be using the multigrid in order to construct a rectangular
domain. The aspect ratio `AR` is used in the `dimensions()`
function. The total size of the domain is also set to `AR`. This fixes
the radial dimension to $R = 1$, thus simplifying the analytical
expressions. Space and time are dimensionless. This is necessary to be
able to use the `mu = fm` trick. The latter means that the dynamic
viscosity is $\mu = 1$ (with the metric incorporated under the hood),
and therefore, G is chosen accordingly so as to yield a maximum
velocity of $u_z(r = 0, t\rightarrow\infty) = 1$.

When ran in parallel, the domain is divided into square subdomains,
each assigned to a different processor. The `LEVEL` of refinement
herein pertains to the subdomains. The total grid is then initialized
as a multiple of the number of cells in the smaller domains.
*/

#define LEVEL 5
#define AR 8
#define G 4.

int main() {
  dimensions (nx = AR, ny = 1);
  size (AR [0]);

  DT = 1e-3 [0];
  
  TOLERANCE = 1e-6;

  init_grid (AR*(1 << LEVEL));
  
  run();
}

/**
The top boundary is a no-slip rigid wall. The bottom boundary is the
axis of symmetry. A constant pressure boundary condition is applied at
both ends so that $\Delta p/L_0 = G$.
*/

u.n[top] = dirichlet(0.);
u.t[top] = dirichlet(0.);

u.n[left] = neumann(0.);
p[left] = dirichlet(G*L0);
pf[left] = dirichlet(G*L0);

u.n[right] = neumann(0.);
p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);

scalar un[];

/**
We define global variables that store the first `nr` roots $\lambda_n$
of the zeroth-order Bessel function of the first kind $J_0(\lambda) =
0$. These will be used in the analytical solution of the transient
velocity profile defined as $$\displaystyle u(r,t) =
\frac{G}{4\mu}(R^2 - r^2) - \frac{2GR^2}{\mu}\sum_{n =
1}^{\infty}\frac{1}{\lambda_n^3}\frac{J_0(\lambda_n
r/R)}{J_1(\lambda_n)}\exp{(-\lambda_n^2\nu t/R^2)},$$ where $J_0$ and
$J_1$ are the zeroth- and first-order Bessel functions of the first
kind, and $\nu$ is the kinematic viscosity. In our case, the upper
bound of the summation is `nr`. We start from reasonbale approximation
of the roots and use a Newton-Raphson method to refine them to a great
degree of accuracy. This is done in
[`bessel-roots.h`](bessel-roots.h).
*/

int nr;
double * roots;

event init (t = 0) {
  mu = fm;

  nr = 50;
  roots = calloc (nr, sizeof(double));
  zeros (nr, roots);
  
  /**
  We initialize the reference velocity. */
  
  foreach()
    un[] = u.x[];
}

/**
We save the veloctiy profiles $u(r = L_0/2,t)$ as well as the
respective analytical solution calculated using the equation above.
*/

event sampling (t += 0.1) {
  int s = 3*N/AR;
  int k = 0;
  double * v = calloc (s, sizeof(double));
  
  foreach_face(x, reduction(max:k) reduction(max:v[:s])) {
    if (x == L0/2.) {
      v[k] = y;
      v[k + N/AR] = (u.x[] + u.x[-1])/2.;
      double temp = 0;
      for (int i = 0; i < nr; i++)
	temp += j0(y*roots[i])*exp(-t*sq(roots[i]))/
	  j1(roots[i])/cube(roots[i]);
      temp *= -2.*G;
      temp += G*(1 - sq(y))/4.;
      v[k + 2*N/AR] = temp;
      k++;
    }
  }

  if (pid() == 0) {
    char name[80];
    sprintf (name, "ux-%g.dat", t);
    FILE * fp = fopen (name, "a");
    for (int i = 0; i < N/AR; i++)
      fprintf (fp, "%.6f\t%.6f\t%.6f\n",
	       v[i], v[i + N/AR], v[i + 2*N/AR]);
    fclose(fp);
  }  
  free(v);
}

/**
We generate a movie of the flow development.
*/

event movie (t += 0.05) {
  output_ppm (u.x, file = "mov.mp4",
	      min = 0., max = 1., linear = true);
}

/**
We check for a stationary solution.
*/

event convergence (t += 0.1; i <= 10000) {
  double du = change (u.x, un);
  if (i > 0 && du < 1e-12)
    return 1; 			/* stop */
}

/**
We also save the converged velocity profile to be checked against the
Poiseuille velocity profile $$\displaystyle u(r) = \frac{G}{4\mu}(R^2
- r^2).$$
*/

event error (t = end) {
  int s = 2*N/AR;
  int k = 0;
  double * v = calloc (s, sizeof(double));
  
  foreach_face(x, reduction(max:k) reduction(max:v[:s])) {
    if (x == L0/2.) {
      v[k] = y;
      v[k + N/AR] = (u.x[] + u.x[-1])/2.;
      k++;
    }
  }
  
  if (pid() == 0) {
    FILE * fp = fopen ("final.dat","a");
    for (int i = 0; i < N/AR; i++)
      fprintf (fp, "%.6f\t%.6f\n", v[i], v[i + N/AR]);
    fclose(fp);
  }
  free(v); free(roots);
}

/**
~~~gnuplot Converged velocity profile
set xlabel 'r'
set ylabel 'u(r)'
set size ratio -1
set xrange[0:1]
set yrange[0:1]
f(x) = 1 - x*x
plot 'final.dat' u 1:2 w p title 'Basilisk' ,\
f(x) w l title 'Analytical'
~~~
~~~gnuplot Transient velocity profiles at t = [0:0.1:1]
set xlabel 'r'
set ylabel 'u(r)'
set size ratio -1
set xrange[0:1]
set yrange[0:1]
l = 10
array styles[2]
styles[1] = "p pt 7"
styles[2] = "l"
cmd = "plot "
cmd = cmd . "keyentry w p pt 7 lc 1 title \"Basilisk\","
cmd = cmd . "keyentry w l lc 1 title \"Analytical\","
do for [i = 0:l] {
do for [c = 1:2] {
cmd = cmd . sprintf("'ux-%g.dat' u 1:%d w %s lc %d notitle",\
0.1*i, c + 1, styles[c], i)
if (i*c != 2*l) {cmd = cmd . ", "}
}
}
eval cmd
~~~
*/
