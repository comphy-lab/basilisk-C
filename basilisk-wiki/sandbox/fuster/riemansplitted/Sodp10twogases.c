/**
# Shock tube problem in a system containing two ideal gases (strong shock wave) using the split Riemann solver

*/

#include "grid/multigrid.h"
#include "vofsplitriemannsolver.h"
#include "Sodtheory.h"

double CFLacous = 0.5;

double rhoL = 10., rhoR = 0.125;
double pL = 10., pR = 0.1;
double gammaL = 1.4, gammaR = 1.6;
double tend = 1.;
double dtmax;

FILE * fp;

void flux (const double * s, double * f, double e[2])
{

  /**
  We first recover each value ($\rho$, $E$, $w_x$ and $w_y$) and then
  compute the corresponding fluxes (`f[0]`, `f[1]`, `f[2]` and
  `f[3]`). */

  double rho = s[0], E = s[1], wn = s[2 + sweep_dir], w2 = 0.;
  for (int i = 2; i < 2 + dimension; i++)
    w2 += sq(s[i]);
  double un = wn/rho;
  double gammao = 1./(cf/(gammaL - 1.) + (1. - cf)/(gammaR - 1.)) + 1.;

  double p = (E - 0.5*w2/rho)*(gammao - 1.);
  
  f[0] = wn;
  f[1] = un*(E + p);
  for (int i = 2; i <= 2 + dimension; i++)
    f[i] = un*s[i];
  f[2+sweep_dir] += p;

  /**
  The minimum and maximum eigenvalues for the Euler system are the
  characteristic speeds $u \pm \sqrt(\gamma p / \rho)$. */

  double c = sqrt(1./(rho*cf/(gammaL*p) + (1 - cf)*rho/(gammaR*p)));
  e[0] = un - c; // min
  e[1] = un + c; // max
}

int main() {
  L0 = 10;
  X0 = Y0 = -L0/2.;
  CFL = 0.1;
  
  for (N = 32; N <= 256; N *= 2) {
    dtmax = CFLacous*L0/N/sqrt(fmax(gammaL*pL/rhoL,gammaR*pR/rhoR));
    run();
  }
  
}

event init (i = 0) {
  foreach() {
    c[] = (x < 0.);
    double p= c[]*pL + (1. - c[])*pR;
    rho[] = c[]*rhoL + (1. - c[])*rhoR;
    Etot[] = c[]*p/(gammaL - 1.) + (1. - c[])*p/(gammaR - 1.);
    foreach_dimension ()
      q.x[] = 0.;
  }
  boundary ({rho,Etot});
}

event outputdata ( t = tend )
{

  FILE * fp    = fopen ("numerical.dat", "w");
  FILE * fp1    = fopen ("theory.dat", "w");
  
  scalar perr[], rhoerr[], uerr[];
  
  foreach () {

      double w2 = 0.;
      foreach_dimension ()
        w2 += sq(q.x[]);
      double p = (Etot[] - 0.5*w2/rho[])/(c[]/(gammaL - 1.) + (1. - c[])/(gammaR - 1.));

      fprintf(fp, "%g %g %g %g %g \n", x, t, p, rho[], q.x[]/(rho[]));
  
      struct SodSol sol = Sod_theory (x/t);
      fprintf(fp1, "%g %g %g %g \n", x/t, sol.pe, sol.rhoe, sol.ue);
       
      perr[] = fabs(p - sol.pe);
      rhoerr[] = fabs(rho[] - sol.rhoe);
      uerr[] = fabs(q.x[]/rho[]  - sol.ue);
    
  }

  printf ("%i %g %g %g \n", N, normf(perr).avg, normf(rhoerr).avg, normf(uerr).avg);
  
  fclose(fp);
  fclose(fp1);

}

/**
 *
~~~gnuplot Pressure profile
set output 'p.png'
set log y
set xrange[-5:5]
set xlabel 'x/t'
set ylabel 'p'
set cblabel 't'
p "theory.dat" u 1:2 t 'Theory' w l lc 0 lw 3,  "./numerical.dat" u ($1/$2):3:2 every 500 t 'Sim' w p pt 7 lc 1
~~~ 
 
~~~gnuplot Density profile
set output 'r.png'
set log y
set xrange[-5:5]
set xlabel 'x/t'
set ylabel '{/Symbol r}'
set cblabel 't'
p "theory.dat" u 1:3 t 'Theory' w l lc 0 lw 3,  "./numerical.dat" u ($1/$2):4:2 every 500 t 'Sim' w p pt 7 lc 1
~~~ 

~~~gnuplot Velocity profile
set output 'u.png'
unset log y
set xrange[-5:5]
set xlabel 'x/t'
set ylabel 'u'
set cblabel 't'
p "theory.dat" u 1:4 t 'Theory' w l lc 0 lw 3,  "./numerical.dat" u ($1/$2):5:2 every 500 t 'Sim' w p pt 7 lc 1
~~~ 

~~~gnuplot Convergence study
set output 'conv.png'
set log xy
set au xy
set xlabel 'N'
set ylabel 'err'
p "out" u 1:2 t 'Err(p)' w p,  "out" u 1:3 t 'Err(rho)' w p, "out" u 1:4 t 'Err(u)' w p, 10*x**(-1) t '10/x' w l
~~~ 
 */  