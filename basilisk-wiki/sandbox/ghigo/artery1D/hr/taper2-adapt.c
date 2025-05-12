/**
# Inviscid wave propagation in a tapered artery

We solve the propgation of an inviscid wave in a tapered artery using
the inviscid 1D blood flow equations.

## Analytic solution

To capture the perturbations induced by vessel tapering, we introduce
the following non-dimensional variables:
$$
t = T \bar{t},\, x = X \bar{x},\,R_0 = R_0 \bar{f},\,R = R_0 \left[
\bar{f} + \Delta_R \bar{R} \right],\,K = K_0 \bar{g},\,Q = Q
\bar{Q},\,p = p_0,\,\Pi \tilde{p}.
$$
Injecting these non-dimensional variables in the 1D blood flow
equations which we then linearize, we obtain the following simplified
equation:

$$
\frac{1}{\bar{c}^2} \frac{\partial^2 \bar{Q}}{\partial \bar{t}^2} -
\frac{\partial^2 \bar{Q} }{\partial \bar{x}^2 } = \left[
\frac{1}{\bar{g}}\frac{\partial \bar{g} }{\partial \bar{x} } -
\frac{1}{\bar{f}}\frac{\partial \bar{f} }{\partial \bar{x} } \right]
\frac{\partial \bar{Q} }{\partial \bar{x} },
$$
where $\bar{c} = \sqrt{\bar{f}\bar{g}}$ is the dimensionless speed. We
then search for a solution of the form:
$$
\bar{Q} = \tilde{Q}\left( \bar{x} \right)\exp\left( i \omega \bar{t}
\right), \quad \omega \in \mathbb{R},
$$
and therefore rewrite the previous equation as:
$$
\frac{\mathrm{d}^2 \tilde{Q} }{\mathrm{d} \bar{x}^2 } +
\frac{\omega^2}{\bar{c}^2} \tilde{Q} = - \left[
\frac{1}{\bar{g}}\frac{\mathrm{d} \bar{g} }{\mathrm{d} \bar{x} } -
\frac{1}{\bar{f}}\frac{\mathrm{d} \bar{f} }{\mathrm{d} \bar{x} }
\right] \frac{\mathrm{d} \tilde{Q} }{\mathrm{d} \bar{x} }.
$$
To keep track of the slowly varying neutral radius $R_0$ and arterial
wall rigidity $K$, we use the following change of variables, to place
ourselves at long $x$ while keeping track of local variations of the
wave speed:
$$
\frac{d \xi}{d \bar{x}} = \Phi^{\prime}\left( X \right) , \quad X =
\epsilon \bar{x},
$$
where $\epsilon$ is the small parameter characterizing the slow
variations of the neutral radius $R_0$ and the arterial wall rigidity
$K$. The function $\Phi^\prime$ represents the wave distortion. Using
this change of variables, we have obtain the following solution at
first order:
$$
\tilde{Q_0} = \frac{B}{\sqrt{\omega}} \frac{
\bar{f}^{\frac{3}{4}}}{\bar{g}^{\frac{1}{4}}} \exp \left( i\omega
\left[ \bar{t} - \frac{1}{\epsilon} \int_{0}^{X} \frac{1}{c}
\mathrm{d}X \right] \right) + C.C., \qquad B = \mathrm{cst}.
$$

We choose here to use an adaptive cartesian grid. */

#include "grid/bitree.h"
#include "../bloodflow-hr.h"

#define BGHOSTS 2

#define lmin (4) // N = 16
#define lmax (7) // N = 128
#define cmax (1.e-3)

/**
We define the artery's geometrical and mechanical properties. */

#define SH (1.e-2)

#define R0 (1.)
#define DR0 (-0.1)

#define K0 (1.e4)
#define DK0 (0.1)

#define shape(x,delta) (1. + (delta)*(x))

#define T0 (1.)
#define celerity(k,a) (sqrt (0.5*(k)*sqrt ((a))))
#define DR ((DR0)/(T0)/celerity((K0), pi*sq (R0)))
#define DK ((DK0)/(T0)/celerity((K0), pi*sq (R0)))

int main() {

  /**
  The domain is 200.*/

  L0 = 200.;
  size (L0);
  origin (0.);
  
  DT = 1.e-5;

  /**
  We run the computation for a maximum level of refinement *lmax*. */
  
  N = (1 << (lmin));
  init_grid (N);
  run();
}

/**
## Boundary conditions

We impose the the flow rate at the inlet and Neumann conditions on all
other variables. */

#define AI (pi*sq ((R0)*(1. + (SH))))
#define QI ((SH)*(AI)*(celerity ((K0),(AI))))

k[left] = (K0)*(shape (x, (DK)));
zb[left] = (K0)*(shape (x, (DK)))*sqrt (pi)*(R0)*(shape (x, (DR)));
q[left] = (QI)*(t < (T0) ? max (0., 0.5*(1. + cos (pi + 2.*pi/(T0)*t))) : 0.);

k[right] = (K0)*(shape (x, (DK)));
zb[right] = (K0)*(shape (x, (DK)))*sqrt (pi)*(R0)*(shape (x, (DR)));

/**
## Defaults conditions
*/

event defaults (i = 0)
{
  gradient = minmod;
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
## Initial conditions */
 
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
    zb[] = k[]*sqrt (pi)*(R0)*(shape (x, (DR)));
    a[] = sq (zb[]/k[]);
    q[] = 0.;
  }
}

/**
## Post-processing

We output the computed fields. */

event field (t = {0., 0.5, 1., 1.5, 2., 2.5})
{
  char name[80];
  sprintf (name, "fields-%.1f-pid-%d.dat", t, pid());
  FILE * ff = fopen (name, "w");

  foreach()
    fprintf (ff, "%g %g %g %g %g %d\n",
	     x, k[]/(K0), sq (zb[]/k[])/(pi*sq ((R0))),
	     (a[] - sq (zb[]/k[]))/(pi*sq ((R0))),
	     q[],
	     level
	     );
}

/**
## Mesh adaptation */

event adapt (i++)
{
  adapt_wavelet ({q}, (double[]){(cmax)},
		 maxlevel = (lmax), minlevel = (lmin));
}

/**
## End of simulation */

event stop_run (t = 2.5)
{
  return 0;
}

/**
# Results for second order

#### Mesh adaptation

We first plot the distribution of the mesh levels throughout the
artery at $t={0, 0.5, 1., 1.5, 2., 2.5}$ for $N_{max}=128$.

~~~gnuplot levels for $lmax=14$. 
reset
set xlabel 'x'
set ylabel 'l'
set yrange[:15]
plot '< cat fields-0.0-pid-*' u 1:6 w l lw 2 lc rgb "blue" t 't=0', \
     '< cat fields-0.5-pid-*' u 1:6 w l lw 2 lc rgb "red" t 't=0.5', \
     '< cat fields-1.0-pid-*' u 1:6 w l lw 2 lc rgb "sea-green" t 't=1', \
     '< cat fields-1.5-pid-*' u 1:6 w l lw 2 lc rgb "coral" t 't=1.5', \
     '< cat fields-2.0-pid-*' u 1:6 w l lw 2 lc rgb "dark-violet" t 't=2', \
     '< cat fields-2.5-pid-*' u 1:6 w l lw 2 lc rgb "skyblue" t 't=2.5'
~~~

#### Arterial properties

~~~gnuplot $a_0$ and $k$ for $Nmax=128$
reset
set key top right
set xlabel 'x'
set ylabel 'a_0,k'
plot '< cat fields-0.0-pid-*' u 1:3 w l lw 2 lc rgb 'blue' t 'a_0/a_0(0)', \
     '< cat fields-0.0-pid-*' u 1:2 w l lw 2 lc rgb 'red' t 'k/k(0)'
~~~

#### Cross-sectional area and flow rate

Next, we plot the spatial evolution of the cross-sectional area $a$
and flow rate $q$ at $t={0, 0.5, 1., 1.5, 2., 2.5}$ for
$N_{max}=128$. We compare the results with those obtained with a
uniform grid with $N=128$.

~~~gnuplot $a/a_0$ for $Nmax=128$. 
reset
set xlabel 'x'
set ylabel '(a - a_0)/a_0'
plot '< cat ../taper2/fields-0.0-pid-*' u 1:4 w l lw 3 lc rgb "black" t 'uniform', \
     '< cat ../taper2/fields-0.5-pid-*' u 1:4 w l lw 3 lc rgb "black" notitle, \
     '< cat ../taper2/fields-1.0-pid-*' u 1:4 w l lw 3 lc rgb "black" notitle, \
     '< cat ../taper2/fields-1.5-pid-*' u 1:4 w l lw 3 lc rgb "black" notitle, \
     '< cat ../taper2/fields-2.0-pid-*' u 1:4 w l lw 3 lc rgb "black" notitle, \
     '< cat ../taper2/fields-2.5-pid-*' u 1:4 w l lw 3 lc rgb "black" notitle, \
     '< cat fields-0.0-pid-*' u 1:4 w l lw 2 lc rgb "blue" t 't=0',	\
     '< cat fields-0.5-pid-*' u 1:4 w l lw 2 lc rgb "red" t 't=0.5', \
     '< cat fields-1.0-pid-*' u 1:4 w l lw 2 lc rgb "sea-green" t 't=1', \
     '< cat fields-1.5-pid-*' u 1:4 w l lw 2 lc rgb "coral" t 't=1.5', \
     '< cat fields-2.0-pid-*' u 1:4 w l lw 2 lc rgb "dark-violet" t 't=2', \
     '< cat fields-2.5-pid-*' u 1:4 w l lw 2 lc rgb "skyblue" t 't=2.5'
~~~

~~~gnuplot $q$ for $Nmax=128$. 
reset
set xlabel 'x'
set ylabel 'q'
plot '< cat ../taper2/fields-0.0-pid-*' u 1:5 w l lw 3 lc rgb "black" t 'uniform', \
     '< cat ../taper2/fields-0.5-pid-*' u 1:5 w l lw 3 lc rgb "black" notitle, \
     '< cat ../taper2/fields-1.0-pid-*' u 1:5 w l lw 3 lc rgb "black" notitle, \
     '< cat ../taper2/fields-1.5-pid-*' u 1:5 w l lw 3 lc rgb "black" notitle, \
     '< cat ../taper2/fields-2.0-pid-*' u 1:5 w l lw 3 lc rgb "black" notitle, \
     '< cat ../taper2/fields-2.5-pid-*' u 1:5 w l lw 3 lc rgb "black" notitle, \
     '< cat fields-0.0-pid-*' u 1:5 w l lw 2 lc rgb "blue" t 't=0',	\
     '< cat fields-0.5-pid-*' u 1:5 w l lw 2 lc rgb "red" t 't=0.5', \
     '< cat fields-1.0-pid-*' u 1:5 w l lw 2 lc rgb "sea-green" t 't=1', \
     '< cat fields-1.5-pid-*' u 1:5 w l lw 2 lc rgb "coral" t 't=1.5', \
     '< cat fields-2.0-pid-*' u 1:5 w l lw 2 lc rgb "dark-violet" t 't=2', \
     '< cat fields-2.5-pid-*' u 1:5 w l lw 2 lc rgb "skyblue" t 't=2.5'
~~~
*/
