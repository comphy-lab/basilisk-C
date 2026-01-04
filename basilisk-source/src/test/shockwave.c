/**
# Shock tube problem for a single ideal gas (strong shock wave)

This test verifies that the method captures the correct propagation
speed of shock waves propagating in an adiabatic perfect gas (with
known polytropic coefficient $\gamma$) given the post-shocked (L) and
pre-shocked conditions (R).
$$
\left(\begin{array}{c}
    p_R\\
    \rho_R\\
    u_R
  \end{array}\right)
  =
\left(\begin{array}{c}
    0.1\\
    1\\
    0
  \end{array}\right)  
$$
$$
\left(\begin{array}{c}
    p_L\\
    \rho_L\\
    u_L
  \end{array}\right)
  =
\left(\begin{array}{c}
    10\\
    \rho_L\\
    u_L
  \end{array}\right)  
$$

See also section 4.1.2 in [Fuster & Popinet,
2018](/src/compressible/two-phase.h#fuster2018).

We use the two-phase flow formulation. */

#include "compressible/two-phase.h"
#include "compressible/Mie-Gruneisen.h"

/**
Parameters of the problem. */

double rhoR = 1.;
double pL = 10., pR = 0.1;
double tend = 1.;
double rhoL, uL, gr, ushock;

/**
Boundary conditions: */

p[left]     = dirichlet (pL);
frho1[left] = dirichlet (rhoL);
q.n[left]   = dirichlet (uL*rhoL);

/**
We obtain the post-shocked state using the Rankine-Hugoniot relations
$$
\rho_L = \rho_R \frac{\gamma_R \frac{p_L}{p_R} + 1}{\gamma_R + \frac{p_L}{p_R}}
$$
$$
u_L = \frac{\sqrt{\gamma p_R/\rho_R}}{\gamma} \frac{\frac{p_L}{p_R} -
1}{\sqrt{\frac{\gamma + 1}{2 \gamma} \left(\frac{p_L}{p_R} - 1\right)
+ 1}}
$$
where $\gamma_R = \frac{\gamma + 1}{\gamma - 1}$. */

int main()
{

  /**
  The EOS for an adiabatic perfect gas is defined by its polytropic
  coefficient $\Gamma = \gamma = 1.4$. */
  
  gamma1 = 1.4;

  gr = (gamma1 + 1.)/(gamma1 - 1.);
  rhoL = rhoR*(gr*pL/pR + 1.)/(gr + pL/pR);
  uL = sqrt(gamma1*pR/rhoR)/gamma1*(pL/pR - 1.)/sqrt((gamma1 + 1.)/2./gamma1*
						     (pL/pR - 1.) + 1.);
  ushock = sqrt(gamma1*pR/rhoR)*sqrt((gamma1 + 1.)/2./gamma1*(pL/pR - 1.) + 1.);
  
  /**
  Size of the domain: */
  
  size (10. [1]);
  origin (-L0/2.);
  
  /**
  We use an upwind method for the tracer advection associated to the
  VOF tracer f. */
  
  f.gradient = zero; 

  /**
  We perform a convergence study. */
  
  for (N = 256; N >= 32; N /= 2)
    run();
}

/**
Variable initialization of the conservative variables: density,
momentum and energy The shock is initially placed at $x = 0$. */

event init (i = 0)
{   
  foreach() {
    double m = (x < 0.);
    p[] = m*pL + (1. - m)*pR;
    frho1[] = rhoR*(gr*p[]/pR + 1.)/(gr + p[]/pR);
    q.x[] = m*frho1[]*uL;
    fE1[] = p[]/(gamma1 - 1.) + 0.5*sq(q.x[])/frho1[];
  }
}

/**
Grid adaptation. */

#if TREE
event adapt (i++) {
  adapt_wavelet ((scalar *){p}, (double[]){0.01}, maxlevel = log(N)/log(2.));
}
#endif

/**
At the end of each simulation we output the relative position of the
shock with respect to the exact theoretical position together with the
conservative variables and pressure to verify that the wave structure
is not distorted by the numerical method The theoretical shockwave
speed is
$$
u_{shock} = \sqrt{\gamma \frac{p_R}{\rho_R} \left( \frac{\gamma + 1}{2
\gamma} \left(\frac{p_L}{p_R} - 1 \right) + 1 \right)}
$$
*/

event endsim (t = tend)
{
  foreach()
    printf ("%i %g %g %g %g %g \n",
	    N, (x - ushock*t)/(ushock*t), p[], frho1[], q.x[]/frho1[], fE1[]);
  printf ("\n");

  /**
  We also compute the L_1 error norm to check convergence. */
  
  double perr = 0., rhoerr = 0., uerr = 0., vol = 0.;
  foreach (reduction(+:vol) reduction(+:perr)
	   reduction(+:rhoerr) reduction(+:uerr)) {
    vol += Delta;
    if (x - ushock*tend <= 0.) {
      perr += fabs(p[] - pL)*Delta;
      rhoerr += fabs(frho1[] - rhoL)*Delta;
      uerr += fabs(q.x[]/frho1[] - uL)*Delta;
    }
    else {
      perr += fabs(p[] - pR)*Delta;
      rhoerr += fabs(frho1[] - rhoR)*Delta;
      uerr += fabs(q.x[]/frho1[])*Delta;
    }
  }
  fprintf (stderr, "%i %g %g %g\n",
	   N, perr/vol/pL, rhoerr/vol/rhoL, uerr/vol/uL);
}

/**
~~~gnuplot Error convergence
set xlabel 'N'
set ylabel 'Average error'
set log xy
set xtics 16,2,512
set xrange [16:512]
set grid
set key bottom left
plot "log" u 1:2 t 'p (adaptive)' w p, "" u 1:3 t 'rho (adaptive)' w p, \
     "" u 1:4 t 'u (adaptive)' w p, \
     "clog" u 1:2 t 'p (multigrid)' w p, "" u 1:3 t 'rho (multigrid)' w p,   \
     "" u 1:4 t 'u (multigrid)' w p, x**(-1.) t '1/x' w l lc 0
~~~ 

~~~gnuplot Velocity, density and pressure profiles
reset
set palette rgb 33,13,10;
set xrange[-0.5:0.5]
set xlabel '(x - s*t)/(s*t)'
set ylabel 'p'
set multiplot
set size 1.,0.33
set origin 0.,0.
unset colorbox
set xtics -2,0.5
set pointsize 0.5
unset ytics
p "out" u 2:3:(log($1)) not w lp palette pt 7 
set origin 0,0.33
set ylabel '{/Symbol r}'
p "out" u 2:4:(log($1)) not w lp palette pt 7
set origin 0, 0.66
set ylabel 'u'
p "out" u 2:5:(log($1)) not w lp palette pt 7
set origin 0,0.75
unset multiplot
~~~ 
*/
