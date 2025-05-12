/**
# Marangoni flow vs Laplace Pressure 

We test here the ability of the multilayer solver to reproduce
 correctly the Benard-Marangoni instability.

For a periodic wave with a height profil and a periodic surface
 tension (phisically, such a gradient could be obtained for a 
 thin film with a temperature gradient) : 
$$
h = h_0 (1 + \beta * cos(kx)) \textrm{ and } \sigma = \sigma_0 (1 + \alpha * cos(kx)) 
$$

Marangoni flows tend to accumulate liquid on the low surface tension 
areas whereas Laplace pressure tends to flatten the surface. 
In lubrication theory, we can demonstrate that this case is stable if :
$$
\alpha = \beta * \frac{2}{3} h^2 k^2
$$
*/ 

/**
## Include and parameters
This test uses the multilayer solver and the surface tension additive. */
#include "grid/multigrid1D.h"
#include "../hydro.h"
#include "layered/remap.h"
#include "../nh_PhiS.h"
#include "../tension.h"


#define Re 10.
#define sigma0 1. //Adimensionned
#define _g 100. //To be superieur to sigma
#define h_0 0.1 //adimensionned (epsilon)
#define alpha 0.01 //relative height amplitude
#define beta (2./3. * (h_0*h_0) * alpha * _g) //Surface tension amplitude

#define T_END 10.
#define DELTA_T (T_END/100.)
#define H_ERR 1e-10

/**
## Main function
The boundary conditions are set to periodic. The test is done with zero-gravity. */
int main()
{
  corr_dux = true;
  corr_dwx = true;
  L0 = 2*pi;
  periodic(right);
  G = _g;
  nu = 1/Re;
  nl = 16;
  for (N = 4; N <= 128; N*=2)  
    run();

  N = 64;  
  for (nl = 1; nl <= 32; nl*=2)
    run();
}

/**
## Initialisation 
The initial shape of the interface is a sinusoidal wave, corresponding 
to the theoretical equilibrium shape. 
Surface tension is initialised as a sin as well. */
scalar Hi[];
double maxerror;
event init (i = 0)
{
  foreach() {
    double H = h_0 * (1 + alpha * cos(x));
    Hi[] = H;
    foreach_layer() {
      h[] = H/nl;
      sigma[] = sigma0 *(1 + beta * cos(x));
    }
  }
  boundary({h, sigma, Hi});
  maxerror = 0;
}

/**
## Outputs 
*/

event output (t += DELTA_T; t = 0; t <= T_END)
{
  /**
  At every timestep, we check whether the height field has
  converged. */
  double dH = change (eta, Hi)*N/L0;
  if (i > 1 && dH < H_ERR)
    return 1; /* stop */
  maxerror = max(maxerror, dH);
  
  char name[80];
  sprintf (name, "height-nl%d-N%d", nl, N);  
  static FILE * fph = fopen (name, "w");
  if (t > DELTA_T)
    fprintf(fph, "%g %g %g\n", t, statsf(eta).sum, dH);
}

event logfile (t = end) {
  /**
  We recompute the reference solution. */
  scalar Href[];
  foreach()
    Href[] =  h_0 *(1 + alpha * cos(x));
  boundary({Href});
  
  /**
  We compute the shape error *ec* and output these in the *log* file. */
  scalar ef[];
  foreach() {
    ef[] = (eta[] - Href[]);
  }
  norm ne = normf (ef);
  fprintf(stderr, "%d %d %g %g %g %g %g %g %g\n", N, nl, h_0, G, nu, maxerror, ne.avg, ne.rms, ne.max);
  
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
	      fprintf (fpc, "%g %g %g %g\n", z, u.x[], - h_0*alpha*sin(x) * _g/nu*(sq(z)/(2.) - z*h_0/3.), x);
        z += h[]/2;
      }
    J += 1;
  }
  fflush(fpc);
  fclose(fpc);
  
  /**
  At the end, we save the final shape of the interface :*/ 
  sprintf (name, "final_shape-nl%d-N%d", nl, N);
  FILE * fpf = fopen (name, "w");
  foreach() {
    fprintf(fpf, "%g %g %g\n", x, eta[], Href[]);
  }
  fflush(fpf);
  fclose(fpf);
}

/**# Movie
Plot the evolution in time with the max speed*/
#if 0
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
	   "set yrange [0:%g]\n"
     "set lmargin at screen 0.1\n"
     "set rmargin at screen 0.9\n", L0, 2*h_0
	   );
}

void plot (FILE * fp)
{
  fprintf (fp,
	   "set title 't = %.2f'\n"
	   "sp '-' u 1:2:3\n", t);
  foreach (serial) {
    double z = zb[];
    fprintf (fp, "%g %g %g\n", x, z, u.x[]);
    foreach_layer() {
      z += h[];
      fprintf (fp, "%g %g %g\n", x, z, u.x[]);
    }
    fprintf (fp, "\n");
  }
  fprintf (fp, "e\n\n");
  //  fprintf (fp, "pause 1\n");
  fflush (fp);  
}

event gnuplot (i+=30)
{
  if (N == 64 && nl == 32) {
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

  if (N == 64 && nl == 32)
    system ("for f in plot-*.png; do convert $f ppm:- && rm -f $f; done | "
	    "ppm2mp4 movie.mp4");
}
#endif

/**
# Results
~~~gnuplot Final shape using multilayers
set style line 1 pt 7 ps 0.7

set terminal @SVG enhanced size 640,640 font ",14"
set output '_final_shape.svg'
set key 
set border
set tics
set xlabel "x"
set ylabel "y"
set size ratio 0.2

plot \
  'final_shape-nl16-N64' u 1:2 w l t 'nl = 16'
~~~

~~~gnuplot Numerical and analytical velocity profiles at the center depending of the number of layers.
set output '_velocity.svg'
set xlabel "velocity"
set ylabel "height"
set key left top
set size ratio 1
plot  \
  'center-nl1-N64' u ($2):1 w p t 'nl = 1' ps 3 , \
  'center-nl2-N64' u ($2):1 w p t 'nl = 2' ps 2, \
  'center-nl4-N64' u ($2):1 w p t 'nl = 4' ps 1, \
  'center-nl16-N64' u ($2):1 w p t 'nl = 16' ps 1, \
  'center-nl32-N64' u ($2):1 w l t 'model' lw 2
~~~
We get the right interface and the correct velocity profile

~~~gnuplot multi changes
set output '_multi_changes.svg'
set xlabel "t"
set ylabel "dH"
plot 'height-nl1-N64' u 1:3 w p t 'h_0 = 0.5' 
~~~
There is a slight error in each cases. The interface exponentially 
relaxes to a new equilibrium position

~~~gnuplot Final error with horizontal resolution
set output '_errorN.svg'
set xlabel "N"
set ylabel "error"
set logscale x
set logscale y
set xrange [4:128]
unset key
plot 'log' every ::0::5 u 1:9 w l ls 1 t 'nl = 16'
~~~
The error is strongly dependant on the horizontal resolution 
and very few on the vertical resolution.
*/