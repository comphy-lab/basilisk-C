/**
# Marangoni flow vs Laplace Pressure 

We test here the ability of the multilayer solver to reproduce
 correctly the Benard-Marangoni instability.

For a periodic wave with a height profil and a periodic surface
 tension (phisically, such a gradient could be obtained for a 
 thin film with a temperature gradient) : 
$$
h = h_0 (1 + \alpha \cos(kx)) \textrm{ and } \sigma = \sigma_0 (1 + \beta \cos(kx)) 
$$

Marangoni flows tend to accumulate liquid on the low surface tension 
areas whereas Laplace pressure tends to flatten the surface. 
In lubrication theory, we can demonstrate that this case is stable if :
$$
\beta = \alpha \frac{2}{3} h^2 k^2 g
$$
*/ 

/**
## Include and parameters
This test uses the multilayer solver and the surface tension additive. */
#include "grid/multigrid1D.h"
#include "crobert/2_Implicit/hydro-tension.h"
#include "layered/implicit.h"
#include "layered/remap.h"

const double g = 100. [1,-2];
const double k =1;
const double h_0 = 0.1;                   //epsilon
#define alpha 0.1                         //relative height variation
#define beta (2./3.*(h_0*h_0)*alpha*g)    //Surface tension variation

#define T_END 200.

/**
## Main function
The boundary conditions are set to periodic. */
int main()
{
  L0 = 2*pi;
  periodic(right);
  
  G = g;
  nu = 1.;
  
  nl = 8;
  TOLERANCE = 1e-6 [*];
  CFL_H = 0.5;
  
  for (N = 4; N <= 128; N*=2)
    run();

  N = 32;  
  for (nl = 1; nl <= 16; nl*=2)
    run();
}

/**
## Initialisation 
Surface tension is initialised as a sinusoidal in x. The initial shape of the
 interface is a sinusoidal wave, corresponding to the theoretical equilibrium 
 shape. */
scalar Hi[], sig[];
double maxerror;
event init (i = 0)
{
  foreach() {
    Hi[] = h_0*(1 + 0.5*cos(2.*x*k));
    foreach_layer() {
      h[] = Hi[]/nl;
      sig[] = 1. + beta * cos(x*k);
    }
  }
  sigma = sig;
  boundary({h, sigma, Hi});
  maxerror = 0;
}

/**
## Outputs */
event output (t += T_END/100.; t <= T_END)
{
  double dH = change (eta, Hi);
  maxerror = max(maxerror, dH);
  
  char name[80];
  sprintf (name, "height-nl%d-N%d", nl, N);  
  static FILE * fph = fopen (name, "w");
  if (t > T_END/100.)
    fprintf(fph, "%g %g %g\n", t, statsf(eta).sum, dH);
}

event logfile (t = end) {
  /**
  We recompute the reference solution. */
  scalar Href[];
  foreach()
    Href[] =  h_0 *(1 + alpha * cos(x*k));
  boundary({Href});
  
  /**
  We compute the shape error *ec* and output these in the *log* file. */
  scalar ef[];
  foreach() {
    ef[] = (eta[] - Href[])/h_0;
  }
  norm ne = normf (ef);
  fprintf(stderr, "%d %d %g %g %g %g\n", N, nl, maxerror, ne.avg, ne.rms, ne.max);

  /**
  We output the velocity profile. */  
  char name[80];
  sprintf (name, "center-nl%d-N%d", nl, N);
  FILE * fpc = NULL;
  fpc = fopen (name, "w");
  int J = 0;
  foreach() {
    double z = zb[];
    if (J == N/4)
      foreach_layer() {
        z += h[]/2;
	      fprintf (fpc, "%g %g %g %g\n", z, u.x[], - h_0*alpha*sin(x*k) * g/nu*(sq(z)/(2.) - z*h_0/3.), x);
        z += h[]/2;
      }
    J += 1;
  }
  fflush(fpc);
  fclose(fpc);
  
  /**
  We save the final shape of the interface :*/ 
  sprintf (name, "final_shape-nl%d-N%d", nl, N);
  FILE * fpf = fopen (name, "w");
  foreach() {
    fprintf(fpf, "%g %g %g %g\n", x, eta[], Href[], (Href[] - eta[])/Href[]);
  }
  fflush(fpf);
  fclose(fpf);
}

/**# Movie
Plot the evolution in time with the max speed*/
#if 1
void setup (FILE * fp)
{
  fprintf (fp,
	   "set pm3d map interpolate 2,2\n"
	   "# jet colormap\n"
	   "set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25"
	   " 0 0.5647 1, 0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392, 0.625"
	   " 1 0.9333 0, 0.75 1 0.4392 0, 0.875 0.9333 0 0, 1 0.498 0 0 )\n"
	   "unset key\n"
	   "set cbrange [-2e-3:2e-3]\n"
	   "set xlabel 'x'\n"
	   "set ylabel 'height'\n"
	   "set xrange [0:%g]\n"
	   "set yrange [%g:%g]\n"
     "set lmargin at screen 0.1\n"
     "set rmargin at screen 0.9\n", L0, h_0 - beta, h_0 + beta
	   );
}

void plot (FILE * fp)
{
  fprintf (fp,
	   "set title 't = %.2f'\n"
	   "sp '-' u 1:2:3\n", t);
  foreach (serial) {
    double z = zb[];
    fprintf (fp, "%g %g %g\n", x, max(z,h_0 - beta), u.x[]);
    foreach_layer() {
      z += h[];
      fprintf (fp, "%g %g %g\n", x, min(max(z,h_0 - beta), h_0 + beta), u.x[]);
    }
    fprintf (fp, "\n");
  }
  fprintf (fp, "e\n\n");
  //  fprintf (fp, "pause 1\n");
  fflush (fp);  
}

event gnuplot (t += T_END/100.; t <= T_END)
{
  if (N == 64 && nl == 8) {
    static FILE * fp = popen ("gnuplot 2> /dev/null", "w");
    fprintf (fp,
  	   "set term pngcairo font \",10\" size 900,300\n"
  	   "set output 'plot-%04d.png'\n", i/50);
    if (i == 0)
      setup (fp);
    plot (fp);
  }
}

event moviemaker (t = end)
{

  if (N == 64 && nl == 8)
    system ("for f in plot-*.png; do convert $f ppm:- && rm -f $f; done | "
	    "ppm2mp4 movie.mp4");
}
#endif

/**
# Results
~~~gnuplot Final shape using multilayer
set style line 1 pt 7 ps 0.7

set terminal @SVG enhanced size 640,640 font ",14"
set xlabel "x"
set ylabel "z"

plot \
  'final_shape-nl8-N64' u 1:2 w l t 'multilayer', \
  'final_shape-nl8-N64' u 1:3 w l t 'model'
~~~

~~~gnuplot Relative error on final shape
set ylabel "relative error"

plot 'final_shape-nl8-N64' u 1:4 w l t 'relative error'
~~~

~~~gnuplot Numerical and analytical velocity profiles at the center depending of the number of layers.
set xlabel "velocity"
set ylabel "height"
plot  \
  'center-nl1-N32' u ($2):1 w p t 'nl = 1' ps 3 , \
  'center-nl2-N32' u ($2):1 w p t 'nl = 2' ps 2, \
  'center-nl4-N32' u ($2):1 w p t 'nl = 4' ps 1, \
  'center-nl16-N32' u ($2):1 w p t 'nl = 16' ps 1, \
  'center-nl16-N32' u ($3):1 w l t 'model' lw 2
~~~
We get the coorect interface and velocity profile.


~~~gnuplot multi changes
set xlabel "t"
set ylabel "dH"
set logscale y
plot 'height-nl8-N64' u 1:3 w p
~~~
There is a slight error in each cases. The interface exponentially 
relaxes to a new equilibrium position

~~~gnuplot Final error with horizontal resolution
set xlabel "N"
set ylabel "error"
set logscale x
set logscale y
plot [4:128] 'log' every ::0::5 u 1:5 w l ls 1
~~~
The error is strongly dependant on the horizontal resolution 
and very few on the vertical resolution.
*/