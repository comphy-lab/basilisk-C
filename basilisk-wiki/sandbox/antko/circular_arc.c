/**
# Circular arc

We solve the equilibrium equation for a capillary surface of constant curvature $\kappa$:
$$
\partial_x\left(\frac{\partial_x \eta}{\sqrt{1+\partial_x \eta^2}}\right) = -\kappa
$$

To solve this nonlinear equation we use the multigrid solver and Newton iterations.
*/

/**
# Circles are not parabola

For constant $\kappa$ the solution of the (full) equation is a circular arc, when the linearised version of the equation yields parabola (of constant second derivatives, rather than constant curvature). For steep slopes the different between the two can be perceptible (illustration here with a "contact angle" of 45$^\circ$).
~~~gnuplot
set term @SVG size 960,240 font ',10'
set size ratio -1
set style line 1 \
    linecolor rgb '#774F38' \
    linetype 1 linewidth 2 \
    pointtype 7 pointsize 1.
set style line 2 \
    linecolor rgb '#E08E79' \
    linetype 1 linewidth 2 \
    pointtype 7 pointsize 1.
set style line 5 \
    linecolor rgb '#3498DB' \
    linetype 1 dashtype 2 linewidth 3 \
    pointtype 7 pointsize 1.

set ylabel 'Surface elevation'
set xlabel "x"
set key left
plot -cos(pi/4)/sqrt(2) + sqrt(1./2-(x-1./2.)**2) t 'exact shape' linestyle 5, './shape' u 1:2 t 'nonlinear solution' with p ls 2, './shape' u 1:3 t 'linear solution' with p ls 1
~~~
*/

#include "grid/multigrid1D.h"
#include "poisson.h"

/**
# Newton iterations
We build a residual and a relaxation functions to make use of the Basilisk multigrid solver.
*/

static double residual_circle (scalar * al, scalar * bl,
                        scalar * resl, void * data)
{
  scalar deta = al[0], b = bl[0], res = resl[0];
  scalar eta = *((scalar *)data);
  double maxres = 0.;
  double lengthl, lengthr;
  foreach (reduction(max:maxres)) {
    lengthl = sqrt (1. + sq ((eta[]  - eta[-1])/Delta));
    lengthr = sqrt (1. + sq ((eta[1] - eta[]  )/Delta));
    res[] = b[] - (deta[1] - deta[])/(sq(Delta)*sq(lengthr)*lengthr)
      + (deta[] - deta[-1])/(sq(Delta)*sq(lengthl)*lengthl);
    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);
  }
  boundary (resl);
  return maxres;
}

/**
The relaxation iterations are performed by computing a local jacobian
*/

static void relax_circle (scalar * al, scalar * bl, int l, void * data)
{
  scalar deta = al[0], b = bl[0];
  scalar eta = *((scalar *)data);
  double lengthl, lengthr;
  foreach_level_or_leaf (l) {
    lengthl = sqrt (1. + sq ((eta[]  - eta[-1])/Delta));
    lengthr = sqrt (1. + sq ((eta[1] - eta[]  )/Delta));
    deta[] = (b[]
              - deta[1] / (sq(Delta)*sq(lengthr)*lengthr)
              - deta[-1]/ (sq(Delta)*sq(lengthl)*lengthl))/
      (-1./(sq(Delta)*sq(lengthr)*lengthr) -
       1./(sq(Delta)*sq(lengthl)*lengthl));
  }   
}

mgstats solve (scalar deta, scalar eta, double kappa, scalar curvature)
{
  scalar b[];
  foreach()
    b[] = -kappa - curvature[];
  boundary({b});
  return mg_solve ({deta}, {b}, residual_circle, relax_circle, &eta);
}

/**
# Linear version
For comparison purposes, we also compute the classic linearised version of the relaxation operator 
*/
static double residual_lin (scalar * al, scalar * bl,
                            scalar * resl, void * data)
{
  scalar eta_lin = al[0], b = bl[0], res = resl[0];
  double maxres = 0.;
  foreach (reduction(max:maxres)) {
    res[] = b[] - (eta_lin[1] - eta_lin[])/sq(Delta)
      + (eta_lin[] - eta_lin[-1])/sq(Delta);
    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);
  }
  boundary (resl);
  return maxres;
}

static void relax_lin (scalar * al, scalar * bl, int l, void * data)
{
  scalar eta_lin = al[0], b = bl[0];
  foreach_level_or_leaf (l) {
    eta_lin[] = (b[]
              - eta_lin[1] / sq(Delta)
              - eta_lin[-1]/ sq(Delta)) /
      (-1./sq(Delta) -
       1./sq(Delta));
  }   
}

mgstats solve_lin (scalar eta_lin, double kappa)
{
  scalar b[];
  foreach()
    b[] = -kappa;
  boundary({b});
  return mg_solve ({eta_lin}, {b}, residual_lin, relax_lin, NULL);
}

/**
# Main function
We compute the shape of a constant curvature surface (hence, a circle, if you followed so far) displaying a $45^\circ$ contact angle on the boundaries.
*/

int main()
{
  init_grid(16);
  L0 = 1.;
  scalar eta[], deta[], eta_lin[], curvature[];
  double kappa = sqrt(2.), maxerror = HUGE;//, contact_angle = asin(L0*kappa/2.);
  double lengthl, lengthr;
  int niter = 0;
  eta[left] = dirichlet(0.);
  eta[right] = dirichlet(0.);
  eta_lin[left] = dirichlet(0.);
  eta_lin[right] = dirichlet(0.);
  deta[left] = dirichlet(0.);
  deta[right] = dirichlet(0.);
  foreach()
    eta[] = eta_lin[] = deta[] = curvature[] = 0.;
  boundary({eta,eta_lin,deta});
  TOLERANCE = 1e-8;
  while (maxerror > 1e-12) {
    niter++;
    maxerror = 0.;
    solve (deta, eta, kappa, curvature);
    foreach() {
      eta[] += deta[];
    }
    boundary ({eta});
    foreach() {
      lengthl = sqrt (1. + sq ((eta[]  - eta[-1])/Delta));
      lengthr = sqrt (1. + sq ((eta[1] - eta[]  )/Delta));
      curvature[] = (eta[1] - eta[])/(sq(Delta)*lengthr)
        - (eta[] - eta[-1])/(sq(Delta)*lengthl); // actually this is minus the curvature
      maxerror = max(maxerror, fabs(-kappa-curvature[]));
    }
  }
  fprintf (stderr, "%i Newton iterations needed. Final true residual : %g\n", niter, maxerror);
  /**
  We also compute the linearised solution for reference.
  */
  solve_lin (eta_lin, kappa);
  char name[80];
  sprintf (name, "shape");
  FILE * fpout = fopen (name, "w");
  foreach() {
    fprintf (fpout, "%g %g %g\n", x, eta[], eta_lin[]);
  }
  fclose (fpout);
}

/** 
# Epilogue 
For this particular example, several tens of Newton iterations were required (note that the first Newton iteration captures here exactly the linearised solution).

`20 Newton iterations needed. Final true residual : 8.08464e-13`

Of course, at each iteration, multigrid iterations are performed to get the solution of the linear system, thereby weighing seriously the overall computational cost.
*/
