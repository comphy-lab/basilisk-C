/**
# Inviscid steady state

We solve the inviscid steady state of the 1D blood flow equations in
an artery with variable properties. A steady state verifies the
following system of equations:
$$
\left\{
\begin{aligned}
& Q = cst \\
& \frac{1}{2}\rho U^2 + p = cst\\
\end{aligned}
\right.
$$
From a numerical standpoint, these steady states are not easily
obtained in an artery with variable properties. Indeed, a numerical
error between the treatment of the flux and the topography source term
can lead to the apparition of spurious velocities. This test case is
therefore designed to test the ability of different the low-Shapiro
hydrostatic reconstruction technique to capture inviscid steady
states. */

/**
We choose here to use an adaptive cartesian grid. */

#include "grid/bitree.h"
#include "../bloodflow-hr.h"

#define lmin (4) // N = 16
int lmax;
#define cmax (1.e-3)

/**
We define the artery's geometrical and mechanical properties. */

#define SH (1.e-3)

#define R0 (1.)
#define DR (-1.e-1)
#define XS (3./10.*(L0))
#define XE (7./10.*(L0))
#define shape(x,delta) ((XS) <= (x) && (x) <= (XE) ? 1. + (delta)/2.*(1. + cos(pi + 2.*pi*((x) - (XS))/((XE) - (XS)))) : 1.)

#define K0 (1.e4)
#define DK (1.e-1)

int main()
{
  /**
  The domain is 10.*/

  L0 = 10.;
  size (L0);
  origin (0.);

  /**
  We run the computation for maximum level of refinements *lmax*. */
  
  for (lmax = 5; lmax <= 8; lmax++) {
    N = (1 << (lmin));
    init_grid (N);
    run();
  }
}

/**
## Boundary conditions

We impose the flow rate at the inlet and use homogeneous Neumann
boundary conditions on all other variables. */

#define celerity(k,a) (sqrt (0.5*(k)*sqrt ((a))))
#define ai(r,s) (pi*sq ((r)*(1 + (s))))
#define qi(k,a,s) ((s)*(a)*(celerity ((k),(a))))

a[left] = neumann (0.);
q[left] = dirichlet (qi ((K0),(ai ((R0),(SH))), (SH)));

a[right] = neumann (0.);
q[right] = neumann (0.);

/**
## Defaults conditions
*/

event defaults (i = 0)
{
  gradient = zero;
}

/**
## Refinement of topography

On a *TREE* we need to make sure that refined cells are initialised
with the correct topography. */

#if TREE
void refine_k (Point point, scalar k)
{
  foreach_child()
    k[] = (K0)*(shape (x, (DK)));
}

void refine_zb (Point point, scalar zb)
{
  foreach_child()
    zb[] = (K0)*(shape (x, (DK)))*sqrt (pi)*(R0)*(shape(x, (DR)));
}
#endif

/**
## Initial conditions

We define *am1* and *qm1* to store the *a* and *q* from the previous
time step, in order to compute the temporal convergence error. */

scalar am1[], qm1[];

event init (i = 0)
{
  /**
  We ensure that $\eta$ is preserved when reconstructing $a$ on a
  *TREE* and that the topography is properly initialized. */

#if TREE
  k.refine = refine_k;
  zb.refine = refine_zb;
  conserve_elevation();
#endif // TREE
  
  /**
  We initialize the variables *k*, *zb*, *a* and *q*. */
  
  foreach() {   
    k[] = (K0)*(shape (x, (DK)));
    zb[] = k[]*sqrt (pi)*(R0)*(shape(x, (DR)));
    a[] = sq (zb[]/k[]);
    q[] = (qi ((K0),(ai ((R0),(SH))), (SH)));

    am1[] = a[];
    qm1[] = q[];
  }
}

/**
## Post-processing

We first compute the temporal convergence error. */

scalar t_err_a[], t_err_q[];

event t_error (i++)
{
  if (lmax == 7) {
    
    foreach() {
      t_err_a[] = (a[] - am1[]);
      t_err_q[] = (q[] - qm1[]);

      am1[] = a[];
      qm1[] = q[];

    }
    boundary ((scalar *) {am1, qm1, t_err_a, t_err_q});
    
    norm na = normf (t_err_a);
    norm nq = normf (t_err_q);

    char name[80];
    sprintf (name, "t_err.dat");
    static FILE * ft = fopen (name, "w");
    fprintf (ft, "%g %.14f %.14f\n",
	     t, na.rms, nq.rms);
  }
}
/**
Next, we compute the spatial error for the flow rate. */

event error (t = end)
{
  scalar err_q[];
  foreach() {
    err_q[] = fabs (q[] - (qi ((K0),(ai ((R0),(SH))), (SH))));
  }
  boundary ((scalar *) {err_q});
  
  norm nq = normf (err_q);

  fprintf (ferr, "%d %g %g %g\n",
	   lmax,
	   nq.avg, nq.rms, nq.max);
}

/**
Finally, we plot the arterial properties as well as the
cross-sectional area and flow rate at the final time. */

event artery_properties (t = end)
{  
  if (lmax == 7) {
    
    char name[80];
    sprintf (name, "properties.dat");
    static FILE * fp = fopen (name, "w");
    foreach() {
      fprintf (fp, "%g %g %g %g %g %g %d\n",
	       x,
	       sq(zb[]/k[])/(pi*(R0)*(R0)), // a0/A0
	       k[]/K0, 
	       fabs (a[] - sq(zb[]/k[]))/(ai ((R0),(SH))), // (a - a0)/ai
	       (qi ((K0),(ai ((R0),(SH))), (SH))), q[],  // qi, q
	       level
	       );
    }
  }
}

/**
## Mesh adaptation */

event adapt (i++)
{
  adapt_wavelet ({zb,eta,q}, (double[]){(cmax),(cmax),(cmax)},
		 maxlevel = (lmax), minlevel = (lmin));
}

/**
## End of simulation */

event stop_run (t = 1.5)
{
  return 0;
}

/**
## Results for first order

#### Mesh adaptation

We first the distribution of the mesh levels throughout the artery.

~~~gnuplot levels for $lmax=7$. 
reset
set xlabel 'x'
set ylabel 'l'
set yrange[0:15]
plot 'properties.dat' u 1:7 w l lw 2 lc rgb 'blue' t 'level'
~~~

#### Arterial properties

~~~gnuplot $a_0$ and $k$ for $lmax=7$
reset
set key top right
set xlabel 'x'
set ylabel 'a_0,k'
plot 'properties.dat' u 1:2 w l lw 2 lc rgb 'blue' t 'a_0/a_0(0)', \
     'properties.dat' u 1:3 w l lw 2 lc rgb 'red' t 'k/k(0)'
~~~

#### Flow rate

We now plot the flow rate, we should be conserved by the well-balanced
scheme. We We compare the results with those obtained with a [first-order
scheme](steady-state.c) for $N=128$.

~~~gnuplot $q$ for $lmax=7$ 
set xlabel 'x'
set ylabel 'q'
plot 'properties.dat' u 1:5 w l lw 3 lc rgb 'black' t 'analytic', \
     '../steady-state/properties.dat' u 1:6 w l lw 2 lc rgb 'blue' t 'uniform', \
     'properties.dat' u 1:6 w l lw 2 lc rgb 'red' t 'adaptive'
~~~

#### Convergence

We now plot the time evolution of the relative $L_2$ error between two
consecutive time steps for the cross-sectional area $a$ and the flow
rate $q$.

~~~gnuplot Temporal convergence for $a$ and $q$
set xlabel 't'
set ylabel 'L_2(a),L_2(q)'
set format y '%.1e'
set logscale y
plot 't_err.dat' u 1:2 w l lw 2 lc rgb 'blue' t 'a', \
     't_err.dat' u 1:3 w l lw 2 lc rgb 'red' t 'q'   
~~~

Finally, we plot the evolution of the error for the flow rate $q$ with
the number of cells $N_{max}$. We compare the results with those obtained
with a [uniform grid](steady-state.c).

~~~gnuplot Spatial convergence for $q$
reset
set xlabel 'N_{max}'
set ylabel 'L_1(q),L_2(q),L_{max}(q)'
set format y '%.1e'
set logscale

ftitle(a,b) = sprintf('order %4.2f', -b)

f1(x) = a1 + b1*x
f2(x) = a2 + b2*x
f3(x) = a3 + b3*x
fit f1(x) '../steady-state/log' u (log($1)):(log($2)) via a1, b1
fit f2(x) '../steady-state/log' u (log($1)):(log($3)) via a2, b2
fit f3(x) '../steady-state/log' u (log($1)):(log($4)) via a3, b3

f11(x) = a11 + b11*x
f22(x) = a22 + b22*x
f33(x) = a33 + b33*x
fit f11(x) 'log' u (log(2**$1)):(log($2)) via a11, b11
fit f22(x) 'log' u (log(2**$1)):(log($3)) via a22, b22
fit f33(x) 'log' u (log(2**$1)):(log($4)) via a33, b33

plot '../steady-state/log' u 1:2 w p pt 6 ps 1.5 lc rgb "blue" t '|q|_1, '.ftitle(a1, b1), \
     exp (f1(log(x))) ls 1 lc rgb "red" notitle, \
     'log' u (2**$1):2 w p pt 6 ps 1.5 lc rgb "sea-green" t '|q|_1, '.ftitle(a11, b11), \
     exp (f11(log(x))) ls 1 lc rgb "coral" notitle, \
     '../steady-state/log' u 1:3 w p pt 7 ps 1.5 lc rgb "navy" t '|q|_2, '.ftitle(a2, b2), \
     exp (f2(log(x))) ls 1 lc rgb "red" notitle, \
     'log' u (2**$1):3 w p pt 7 ps 1.5 lc rgb "dark-green" t '|q|_2, '.ftitle(a22, b22), \
     exp (f22(log(x))) ls 1 lc rgb "coral" notitle, \
     '../steady-state/log' u 1:4 w p pt 5 ps 1.5 lc rgb "skyblue" t '|q|_{max}, '.ftitle(a3, b3), \
     exp (f3(log(x))) ls 1 lc rgb "red" notitle, \
     'log' u (2**$1):4 w p pt 5 ps 1.5 lc rgb "forest-green" t '|q|_{max}, '.ftitle(a33, b33), \
     exp (f33(log(x))) ls 1 lc rgb "coral" notitle
~~~
*/


