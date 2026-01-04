/**
# Tidally-induced internal lee waves

This example illustrates the importance of non-hydrostatic effects for
the generation of internal waves by barotropic tides flowing over a
sill, as first studied by [Xing & Davies, 2007](#xing2007).

A close up of the full domain close to the sill is illustrated in the
movie below. The sill geometry (hidden in [Berntsen et al,
2009](#berntsen2009)) is given by
$$
\begin{aligned}
H(x) & = -50 + \frac{35}{1+(x/500)^4} & \text{ if }x < 0 \\
H(x) & = -100 + \frac{85}{1+(x/500)^4} & \text{ if }x > 0 \\
\end{aligned}
$$
The barotropic tide is imposed as inflow on the left
boundary. Breaking, non-hydrostatic internal waves are generated on
the lee side.

![Evolution of the temperature field](lee/movie.mp4)

The results at $t=2/8 T$ and $t = 3/8 T$ below, with $T$ the M2 tidal
period, can be compared, for example, with the corresponding results
of [Klingbeil & Burchard, 2013](#klingbeil2013), Figure 9, a1 and b1.

![Temperature field at $t=2/8 T$](lee/T-2.png)

![Temperature field at $t=3/8 T$](lee/T-3.png)

## References

~~~bib
@article{davies2007,
  title={On the influence of stratification and tidal forcing upon 
         mixing in sill regions},
  author={Davies, Alan M and Xing, Jiuxing},
  journal={Ocean Dynamics},
  volume={57},
  number={4-5},
  pages={431--451},
  year={2007},
  publisher={Springer},
  doi = {10.1007/s10236-007-0114-5}
}

@article{xing2007,
  title = {On the importance of non-hydrostatic processes in determining 
           tidally induced mixing in sill regions},
  journal = {Continental Shelf Research},
  volume = {27},
  number = {16},
  pages = {2162 - 2185},
  year = {2007},
  issn = {0278-4343},
  doi = {https://doi.org/10.1016/j.csr.2007.05.012},
  url = {http://www.sciencedirect.com/science/article/pii/S0278434307001379},
  author = {Jiuxing Xing and Alan M. Davies}
}

@article{berntsen2009,
  title={Numerical studies of flow over a sill: sensitivity of the 
         non-hydrostatic effects to the grid size},
  author={Berntsen, Jarle and Xing, Jiuxing and Davies, Alan M},
  journal={Ocean Dynamics},
  volume={59},
  number={6},
  pages={1043},
  year={2009},
  publisher={Springer},
  doi={10.1007/s10236-009-0227-0}
}

@article{klingbeil2013,
  title = {Implementation of a direct nonhydrostatic pressure 
           gradient discretisation into a layered ocean model},
  journal = {Ocean Modelling},
  volume = {65},
  pages = {64 - 77},
  year = {2013},
  issn = {1463-5003},
  doi = {10.1016/j.ocemod.2013.02.002},
  url = {http://www.sciencedirect.com/science/article/pii/S1463500313000280},
  author = {Knut Klingbeil and Hans Burchard}
}
~~~

## Code

We use the 1D (horizontal) non-hydrostatic multilayer solver. */

#include "grid/multigrid1D.h"
#include "layered/hydro.h"
#include "layered/nh.h"

/**
The Boussinesq density perturbation is given as a function of the
"temperature" field $T$. */

double drho0 = 1e-3, T0 = 13.25, T1 = 8., H0 = 100.;
#define drho(T) (drho0*(T0 - T)/(T0 - T1))
#define T0(z) (T1 + (T0 - T1)*(z + H0)/H0)
#include "layered/dr.h"

/**
The layer positions are remapped to $\sigma$ levels and performance
statistics are displayed. */

#include "layered/remap.h"
#include "layered/perfs.h"
// #include "profiling.h"

/**
The horizontal viscosity is set to 0.1 m/s^2^. */

double nu_H = 0.1;

int main()
{
  L0 = 21500 [1];
  X0 = - 6500;
  G = 9.81;

  /**
  100 layers are used and the horizontal resolution of ~10 metres
  matches that of [Klingbeil & Burchard,
  2013](#klingbeil2013). Monotonic limiting is used for vertical
  remapping. */

  nl = 100;
  N = 2048;
  cell_lim = mono_limit;

  /**
  The maximum timestep is set to 100 seconds. The actual timestep is
  limited to about 5 seconds due to the CFL condition based on the
  maximum horizontal velocity and spatial resolution. Note that this
  is still much larger than the timestep used by Klingbeil & Burchard
  (0.56 seconds) and Bernsten et al., 2009 (0.3 seconds). An
  explanation for such a small timestep could be that the CFL
  restriction due to vertical motions can be quite restrictive for a
  vertically-Eulerian discretisation. For this (non-hydrostatic)
  example, the vertical velocities are of the same order as the
  horizontal velocities and since the vertical resolution is
  approx. 10 times larger than the horizontal resolution, the vertical
  CFL criterion is correspondingly smaller. This restriction is
  avoided with the vertically-Lagrangian solver. */
  
  DT = 100. [0,1];
  
  nu = 1e-3;
  
  system ("rm -f plot-*.png"); // fixme: should not be necessary
  run();
}

/**
The M2 tidal period (in seconds). */

const double M2 = 12.*3600. + 25.2*60.;

/**
The temperature profile $T0(z)$ is imposed at inflow. The slightly
complicated function below computes the vertical coordinate of a layer
and returns the corresponding temperature. */

double Tleft (Point point)
{
  double zc = zb[];
  for (int l = - point.l; l < nl - point.l; l++)
    if (l < 0)
      zc += h[0,0,l];
  zc += h[]/2.;
  return T0(zc);
}

/**
The M2 tide with an amplitude of 0.3 m/s is imposed at inflow (left
boundary) as well as the temperature profile. The outflow (right
boundary) is free. */

double U2 = 0.3;

event init (i = 0)
{
  u.n[left]  = dirichlet (U2*sin(2.*pi*(t/M2)));
  u.n[right] = neumann(0.);
  h[right] = dirichlet(H0/nl);
  T[left] = Tleft(point);
  T[right] = Tleft(point);

  /**
  The sill geometry, initial layer depths and initial temperature
  profile. */

  double a = 35., b = 85., c = 500. [1];
  foreach() {
    zb[] = x < 0. ?
      -H0/2. + a/(1. + pow(x/c,4)) :
      -H0 + b/(1. + pow(x/c,4));
    double z = zb[];
    foreach_layer() {      
      h[] = - zb[]/nl;
      z += h[]/2.;
      T[] = T0(z);
      z += h[]/2.;
    }
  }
}

/**
A naive discretisation of the horizontal viscosity. */

event viscous_term (i++)
{
  if (nu_H > 0.) {
    scalar d2u[];
    foreach_layer() {
      foreach()
	d2u[] = (u.x[1] + u.x[-1] - 2.*u.x[])/sq(Delta);
      foreach()
	u.x[] += dt*nu_H*d2u[];
#if NH
      foreach()
	d2u[] = (w[1] + w[-1] - 2.*w[])/sq(Delta);
      foreach()
	w[] += dt*nu_H*d2u[];
#endif // NH
    }
  }
}

/**
## Outputs */

// fixme: plotting is (almost) the same as overflow.c

void setup (FILE * fp)
{
  fprintf (fp,
#if ISOPYCNAL
	   "set pm3d map corners2color c2\n"
#else
	   "set pm3d map interpolate 2,2\n"
#endif
	   "# jet colormap\n"
	   "set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25"
	   " 0 0.5647 1, 0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392, 0.625"
	   " 1 0.9333 0, 0.75 1 0.4392 0, 0.875 0.9333 0 0, 1 0.498 0 0 )\n"
	   "unset key\n"
	   "set cbrange [8:13.5]\n"
	   "set xlabel 'x (m)'\n"
	   "set ylabel 'depth (m)'\n"
	   "set xrange [-1500:2000]\n"
	   "set yrange [-100:1]\n"
	   );
}

void plot (FILE * fp)
{
  fprintf (fp,
	   "set title 't = %.2f M2/8'\n"
	   "sp '-' u 1:2:4\n",
	   t/(M2/8.));
  foreach (serial) {
    double z = zb[];
    fprintf (fp, "%g %g %g %g\n", x, z, u.x[], T[]);
    foreach_layer() {
      z += h[];
      fprintf (fp, "%g %g %g %g\n", x, z, u.x[], T[]);
    }
    fprintf (fp, "\n");
  }
  fprintf (fp, "e\n\n");
  //  fprintf (fp, "pause 1\n");
  fflush (fp);  
}

event gnuplot (t += M2/1024.)
{
  static FILE * fp = popen ("gnuplot 2> /dev/null", "w");
  fprintf (fp, "set term x11 size 1024,300\n");
  if (i == 0)
    setup (fp);
  plot (fp);
  fprintf (fp,
	   "set term pngcairo font \",10\" size 1024,300\n"
	   "set output 'plot-%04d.png'\n"
	   "replot\n", i);
}

event figures (t <= M2/2.; t += M2/8.)
{
  FILE * fp = popen ("gnuplot 2> /dev/null", "w");  
  fprintf (fp,
	   "set term pngcairo font \",10\" size 1024,300\n"
	   "set output 'T-%g.png'\n", t/(M2/8.));
  setup (fp);
  plot (fp);
}

event moviemaker (t = end)
{
  system ("for f in plot-*.png; do convert $f ppm:- && rm -f $f; done | "
	  "ppm2mp4 movie.mp4");
}
