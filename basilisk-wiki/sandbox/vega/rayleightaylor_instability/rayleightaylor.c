/**
# Rayleigh--Taylor instability 

## Acknowledgment

This code was inspired and can be considered as an upgrade of this Basilisk [test case](https://basilisk.fr/src/test/rt.c) , and this [code in the sandbox of haouche](https://basilisk.fr/sandbox/haouche/Hydrodynamic_instabilities/rayleigh-taylor.c)

## How to run ?

DTENSION can be set to zero (or undefined) if you want to run a simulation without surface tension.

~~~
rm -r rt/
mkdir rt
qcc -autolink -O2 -DMTRACE=3 -DTREE=1 -DREDUCED=1 -DTENSION=1 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -o rt/rt rt.c -lm

cd rt
./rt > out 2> log
~~~

## The Code

*/


#include "navier-stokes/centered.h"
#include "two-phase.h"
#if TENSION
  #include "tension.h"
#endif
#include "navier-stokes/conserving.h"
#if REDUCED
  # include "reduced.h"
#endif

#define LEVEL 7

#if TENSION
  double Bo = 50.;
#else
  double Bo = HUGE;
#endif
double r = 0.5;

double dt_output = 0.1;

double Tmax = 5.;
double eps = 0.01;

bool do_fraction = false;
bool do_levels = false;
bool do_velocity = false;

bool do_plots = false;
bool do_movie = false;
bool do_snap = true;



void gnu_title (FILE * fp, char * title) {
  fprintf(fp,
          "set title \" %s \\nBo = %g, r = %g , {/Symbol e} = %g \\n LEVEL = %d, DT = %g \"\n",
          title, Bo, r, eps, LEVEL, DT);
//  fprintf(fp,
//          "set title \" %s - %s \\nRe = %g, {/Symbol a} = %g, Fr = %g , {/Symbol e} = %g \\n n_res = %d, DT = %g \"\n",
//          title, simulation_type, Re, ratio, Fr, epsilon, N_res, DT);
}

int main() {

  

  if (true) { //WARNING : the use of periodic BC breaks the resolution of the Poisson equation for the pressure
    uf.n[left]   = 0; //uf meaning u face
    p[left] = neumann(0.);
    uf.n[right]  = 0;
    p[right] = neumann(0.);
  }
  else {
    periodic (right);
  }

  size (4);
  origin (-2, -2);
  init_grid (1 << LEVEL);
  mu1 = mu2 = 0.;//00313;
  rho1 = 1. + r, rho2 = 1 - r; //1.694
#if TENSION
  f.sigma = 1./Bo;
#endif
#if REDUCED
  G.x = 0.;
  G.y = -1;
#else
    a[] = {0.,-1.};
#endif
  DT = 5e-3;
  TOLERANCE = 1e-10;
  run();
}

event init (t = 0) {
#if TREE  
  mask (x < -0.5 ? left : x > 0.5 ? right : none);
#endif
  fraction (f, eps*cos (2.*pi*x) + y); // defined by y >< 0 - eps cos (k x)
  //initialize with a droplet
  //double delta = 0.1;
  //fraction (f, 2*exp(-(x*x + y*y)/(delta*delta)) - 1);// < 0 : 0 (fluid 1) / > 0 : 1 (fluid 2)

  foreach() {
    u.x[] = 0.;
    u.y[] = 0.;
  }
}

event logfile (i++) {
  if (i==0){
    fprintf (stderr, "t, i, dt, s.sum - 2., s.min, s.max - 1., mgp.i, mgpf.i, mgu.i\n");
  }
  stats s = statsf (f);
  fprintf (stderr, "%g %d %g %g %g %g %d %d %d\n", 
	  t, i, dt, s.sum - 2., s.min, s.max - 1., mgp.i, mgpf.i, mgu.i);
  //assert (s.min >= -1e-10 && s.max <= 1. + 1e-10);
  //assert (fabs (s.sum - 2.) < 9e-7);
}


event output_ux (t=0.;t+=dt_output) {
  if (do_velocity) {
    char file_name[32];
    snprintf(file_name, sizeof(file_name), "data_ux_%06d", i);
    FILE * ff = fopen(file_name, "w");
    output_field ({u.x}, ff);
    fclose (ff);
  }
}

event output_levels (t=0.;t+=dt_output) {
  if (do_levels) {
    char file_name[32];
    snprintf(file_name, sizeof(file_name), "data_levels_%06d", i);
    FILE * ff = fopen(file_name, "w");
    scalar l[];
    foreach()
      l[] = level;
    output_field ({l}, ff);
    fclose (ff);

    if (false) {
      static FILE * fp = popen("gnuplot 2> /dev/null", "w");
      fprintf(fp,
                  "set pm3d map\n"
                  "set palette defined(0 \"red\", 1 \"blue\")\n"
                  "set cbrange [0:1]\n"
                  "sp 'data_levels_%06d' u (abs($1) < 0.5 ? $1 : NaN):2:3\n" //$1 x-axis
                  "set term pngcairo size 1024,720\n"
                  "set output 'plot_levels_%06d.png'\n"
                  "replot\n",
                  i, i); 
    }
  }
}

event output_fraction (t=0.;t+=dt_output) {
  if (do_fraction) {
    char file_name[32];
    snprintf(file_name, sizeof(file_name), "data_fraction_%06d", i);
    FILE * ff = fopen(file_name, "w");
    output_field ({f}, ff);
    fclose (ff);

    static FILE * fp = popen("gnuplot 2> /dev/null", "w");

    char title[64];
    snprintf(title, sizeof(title), "Volume Fraction (t = %g)", t);
    gnu_title (fp, title);
    fprintf(fp,
            "set xlabel 'x'\n"
            "set ylabel 'y'\n"
            );
    fprintf(fp,
                  "set pm3d map\n"
                  "set palette defined(0 \"red\", 1 \"blue\")\n"
                  "set cbrange [0:1]\n"
                  "sp 'data_fraction_%06d' u (abs($1) < 0.5 ? $1 : NaN):2:3 notitle\n" //$1 x-axis
                  "set term pngcairo size 1024,720\n"
                  "set output 'plot_fraction_%06d.png'\n"
                  "replot\n",
                  i, i);
    //fflush(fp); // ?
  }
}

event output_interface (t=0.;t+=dt_output) {
  // outputing the interface in a file
  char file_name[32];
  snprintf(file_name, sizeof(file_name), "data_interface_%06d", i);
  FILE * ff = fopen(file_name, "w");
  output_facets (f, ff);
  fclose (ff);
  // using this file to plot
  static FILE * fp = popen("gnuplot 2> /dev/null", "w");
  fprintf(fp,
              "reset\n"
              "set fit errorvariables\n"
              "FIT_CONVERGED = 0\n"
              "a = %g\n"
              "f(x) = -a*cos(2*pi*x)\n"
              "fit f(x) 'data_interface_%06d' u 1:2 via a\n"
              "set print 'amplitude' append\n" // give some error values if not converged
              "if (!FIT_CONVERGED) {a = NaN; a_err = NaN; FIT_STDFIT = NaN}\n"
              "print %g, a, a_err, FIT_STDFIT\n"
              "unset print\n",
              eps, i, t);
  if (false) {
    char title[64];
    snprintf(title, sizeof(title), "Interface {/Symbol n}(t = %g)", t);
    gnu_title (fp, title);
    fprintf(fp,
            "set xlabel 'x'\n"
            "set ylabel 'y'\n"
            );

    fprintf(fp,
                "plot 'data_interface_%06d' u 1:2 w p t 'data', 'data_interface_%06d' u 1:(f($1)) w l t 'fit a*cos'\n"
                "set term pngcairo size 1024,720\n"
                "set output 'snap_%06d.png'\n"
                "replot\n",
                i, i, i);
  }
}


#if TREE
event adapt (i++) {
  coord momb = {0};
  foreach()
    foreach_dimension()
      momb.x += dv()*rho(f[])*u.x[];
  
  adapt_wavelet ({f}, {5e-3}, LEVEL);

  coord moma = {0};
  foreach()
    foreach_dimension()
      moma.x += dv()*rho(f[])*u.x[];

  /**
  We check that each component of the momentum is conserved. */
  
//  foreach_dimension()
//    assert (fabs(momb.x - moma.x) < 1e-10);
}
#endif

double sigma_theory (double r, double Bo) {
  return 2*pi*r - 8*pow(pi,3)/(2*Bo);
}

event fit_amplitude(t = end) {
  static FILE * fp = popen("gnuplot 2> /dev/null", "w");
  fprintf(fp,
              //"plot 'data_%06d' u 1:2\n"
              "set fit errorvariables\n"
              "a = 0.05\n" //theory ?
              "b = log(%g)\n" // ln(eps) ?
              "f(x) = a*x+b\n"
              "fit f(x) 'amplitude' u ($1 > 2 ? $1 : NaN):($3 < 0.01 ? log($2) : NaN) via a,b\n" // via a ?
              "set print 'sigma' append\n" // give some error values if not converged
              "print a, FIT_STDFIT, %g, %g, %d\n" //add Bo, r and LEVEL ?
              "unset print\n",
              eps , Bo, r, LEVEL);

  fprintf(fp,
              //"plot 'data_%06d' u 1:2\n"
              "set fit errorvariables\n"
              "w = 1.\n" //theory ?
              "A = %g\n" 
              "f(x) = A*cos(w*x)\n"
              "fit f(x) 'amplitude' u 1:2 via w\n"
              "set print 'omega' append\n" // give some error values if not converged
              "print w, FIT_STDFIT, %g, %g, %d\n" //add Bo, r and LEVEL ?
              "unset print\n",
              eps , Bo, r, LEVEL);
}

event stop (t = Tmax) {
  return 1.;
}