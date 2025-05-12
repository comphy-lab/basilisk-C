/**
# Kapitza waves
This examples tests the ability of the multilayer solver and the surface 
tension additive to reproduce Kapitza waves. The parameters are chosen to
 reproduce the experiments of Liu and Gollub (1993) and the simulations 
 of Malamataris(2002).*/

#include "grid/multigrid1D.h"
#include "../hydro.h"
#include "layered/remap.h"
#include "../nh_PhiS.h"
#include "../tension.h"

/** A flow $Q$ of liquid of viscosity $\nu$ and surface tension $\sigma$ is 
imposed on the top of an inclined plan, which makes an angle $\beta$ 
with the horizontal. 

Without any perturbation, this problem has a analytical solution, known as 
the Nusselt film. The film thickness is related to the Reynolds number 
$Re = Q/\nu$ and is given by:
$$
H_{Nusselt} = \left(\frac{3\,\nu^2\,Re}{g\,\sin{\beta}}\right)^{1/3}
$$

The mean Nusselt velocity $U$ verifies $Q = UH$ and is thus:
$$
U_{Nusselt} = \frac{g\,\sin{\beta}\,H_{Nusselt}^2}{3\,\nu}
$$*/

#define beta (6.4*pi/180.)
#define Re 19.33
#define nu0 0.0000063
#define H0 (pow((3.*sq(nu0)*Re)/(9.81*sin(beta)), 1./3.))
#define U0 ((9.81*sin(beta)*sq(H0))/(3.*nu0))


/**
In this simulation, we choose to non-dimesionalise the parameters with the 
density $\rho$, the Nusselt height $H_{Nusselt}$ and velocity $U_{Nusselt}$.
We obtain the Froude and Weber numbers:
$$
We = \frac{\sigma}{\rho U^2 H}
$$
and
$$
Fr^2 = \frac{U}{\sqrt{gH}} = \frac{Re}{3} \sin{\beta}
$$*/
#define Fr2 (Re*sin(beta)/3.)
#define We 5.43
#define T0 (H0/U0)

/**A periodic perturbation of frequency $F$ is imposed. This leads to non linear
 instability called Kapitza waves.
![Surface Profil](Kapitza_7/movie.mp4)
![Surface Profil of Malamataris](Kapitza/Malamataris7.png)
*/

#define F 7. //dimensionned frequency (Hz)
#define k (2.*pi*F*T0)

#define T_END (0.6/H0)
#define DELTA_T (T_END/500.)

/**
## Main function */
int main()
{
  corr_dux = true;
  corr_dwx = true;
  horvisco = true;
  L0 = 1./H0;
  N = 1024;
  nl = 50;
  
  G = cos(beta)/Fr2;
  nu = 1./Re;
  
  system ("rm -f plot-*.png");
  run();
}

/**
## Initialisation 
The inflow is a parabolic profile with a total flow rate Q. The
function below computes the height zc of the middle of the layer and
returns the corresponding velocity. For the multilayer solver, we need
some trickery to define the inflow velocity profile as a function of
zc. We need to sum the layer thicknesses to get the total depth H and
the height zc of the current layer (of index point.l). */



double Qparab (Point point)
{
  double H = 0.;
  double zc = 0.;
  for (int l = 0; l < nl; l++) {
    H += h[0,0,l];
    if (l < point.l)
      zc += h[0,0,l];
  }
  zc += (h[]/2.);
  return 3./2.*(2.*(zc/H) - sq(zc/H));
}

event init (i = 0)
{
  foreach() {
    sigma[] = We;
    foreach_layer()
      h[] = 1./nl;
    foreach_layer()
      u.x[] = Qparab (point);
  }
  
  u.n[left] = Qparab (point);
  eta[left] = dirichlet(1 + 0.01*sin(k*t)); 
    
  u.n[right] = neumann(0.);
  eta[right] = dirichlet(1.);
  phi[right] = dirichlet(0.);
  boundary(all);
}

/**
## Slope

This imposes a slope without changing the topography. */

event source (i++)
{
  double slope = sin(beta);
  foreach()
    foreach_layer()
      u.x[] += (1./Fr2)*slope*dt;
  boundary ((scalar *){u});
}

/**
## Outputs*/
event logfile(t = end) {
  foreach() {
    double H = 0;
    double Q = 0;
    foreach_layer() {
      H += h[];
      Q += h[]*u.x[];
    }
    fprintf (stderr, "%g %d %g %g\n", x, point.l, H, Q);  
  }
}

/**
# Movie*/
void setup (FILE * fp)
{
  fprintf (fp,
	   "set pm3d map interpolate 2,2\n"
	   "# jet colormap\n"
	   "set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25"
	   " 0 0.5647 1, 0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392, 0.625"
	   " 1 0.9333 0, 0.75 1 0.4392 0, 0.875 0.9333 0 0, 1 0.498 0 0 )\n"
	   "unset key\n"
	   "set cbrange [0:3]\n"
	   "set xlabel 'x'\n"
	   "set ylabel 'height'\n"
	   "set xrange [%g:%g]\n"
	   "set yrange [0.6:1.2]\n"
     "set lmargin at screen 0.1\n"
     "set rmargin at screen 0.9\n", (X0*H0*1000), (X0 + L0)*H0*1000
	   );
}

void plot (FILE * fp)
{
  fprintf (fp,
	   "set title 't = %.2f'\n"
	   "sp '-' u 1:2:3\n", t*T0);
  foreach (serial) {
    double z = 0;
    fprintf (fp, "%g %g %g\n", x*H0*1000, max(0.6, z), u.x[]);
    foreach_layer() {
      z += h[];
      fprintf (fp, "%g %g %g\n", x*H0*1000, min(max(0.6, z), 1.2), u.x[]);
    }
    fprintf (fp, "\n");
  }
  fprintf (fp, "e\n\n");
  //  fprintf (fp, "pause 1\n");
  fflush (fp);  
}

event gnuplot (t += DELTA_T; t <= T_END)
{
  static FILE * fp = popen ("gnuplot 2> /dev/null", "w");
  fprintf (fp,
	   "set term pngcairo font \",10\" size 1024,320\n"
	   "set output 'plot-%04d.png'\n", (int) (t/DELTA_T));
  if (i == 0)
    setup (fp);
  plot (fp);
}

event moviemaker (t = end) {
    system ("for f in plot-*.png; do convert $f ppm:- && rm -f $f; done | "
	    "ppm2mp4 movie.mp4");
}
