/**
# A bubble shrinking due to thermal effects

This reproduces Figure 3.a of [Saade et al, 2023](#saade2023) where
more detailed explanations can be found (section "4.1 Epstein-Plesset
like problem for temperature").

~~~gnuplot Bubble radius as a function of time
set xlabel 't/tau'
set ylabel 'R/R_0'
plot "log" u 1:2 w l t 'Spherical', "../shrinking-axi/log" u 1:2 w l t 'Axisymmetric'
~~~

![Animation of the pressure and temperature
 fields.](shrinking-axi/T2.mp4)( width = "50%" )

~~~bib
@hal{saade2023, hal-03950917}
~~~

We run both the axisymmetric and the spherically-symmetric versions. */

#if AXIS
# include "axi.h"
# include "compressible/thermal.h"
# include "view.h"
# define LEVEL 8
#else
# include "spherisym.h"
# include "compressible/thermal.h"
# define LEVEL 9
#endif
#include "compressible/NASG.h"

/**
The initial density of the gas is chosen such that the initial
temperature inside the bubble is twice the far-field temperature
$T_\infty$. */

double rhoL = 1., rhoR = 0.02011771644;
double p0L = 1.;
double p0 = 1.;
double tend = 1.;
double R0 = 1.;
double tau;

/**
The problem is rendered dimensionless using the ambient pressure, the
liquid density, the far-field temperature and the bubble initial
radius. The values employed for this simulation are respectively
listed. */

double pdim = 5e6;
double rhodim = 975.91;
double Tdim = 350;
double Rdim = 1e-4;

p[right]   = dirichlet(p0L);
q.n[right] = neumann(0.);

#if AXIS
p[left]    = dirichlet(p0L);
q.n[left]  = neumann(0.);

p[top]    = dirichlet(p0L);
q.n[top]  = neumann(0.);
#endif

/**
Although the thermal solver is implicit and unconditionally stable, a
diffusive CFL condition is employed for better accuracy. */

event stability (i++) {
  dtmax = rhoR*cp2*sq(L0/pow(2,LEVEL))/kappa2/2.;
}

int main()
{
  L0 = 8.;
#if AXIS  
  X0 = -L0/2.;
#endif

  /**
  Liquid water parameters in the Noble-Abel Stiffened Gas (NASG)
  equation of state. */

  gamma1 = 1.187;
  PI1 = 7028e5/pdim;
  b1 = 6.61e-4*rhodim;
  q1 = -1177788*rhodim/pdim;

  /**
  Specific heats and thermal conductivity of the fluids. */
  
  cv1 = 3610*rhodim*Tdim/pdim; cv2 = 729.1*rhodim*Tdim/pdim;
  cp1 = 4285*rhodim*Tdim/pdim; cp2 = 1063*rhodim*Tdim/pdim;

  kappa1 = 0.6705/(Rdim/Tdim*sqrt(cube(pdim)/rhodim));
  kappa2 = 0.03153/(Rdim/Tdim*sqrt(cube(pdim)/rhodim));

  mu1 = 3.7e-4/(Rdim*sqrt(rhodim*pdim));
  mu2 = mu1*1e-2;

  /**
  The diffusive time scale $\tau$ based on the gas properties. */

  tau = rhoR*cp2/kappa2;
  
  tend *= tau;

#if TREE
  N = 1 << 4;
#else
  N = 1 << LEVEL;
#endif
  
  run();
}

event init (t = 0)
{
  if (!restore (file = "restart")) {
    
    /**
    The static mesh refinement. */

#if TREE
  for (int l = 4; l <= LEVEL; l++)
    refine (level < l && sqrt(sq(x) + sq(y)) < (2.5*R0 + 4.*sqrt(2.)*L0/(1 << (l - 1))));
#endif
  
    /**
    Initialization of a bubble with initial radius `R0`. */
    
    fraction (f, - (sq(R0) - sq(x) - sq(y)));
    
    foreach() {
      frho1[]  = f[]*rhoL;
      frho2[]  = (1. - f[])*rhoR;
    
      double pL = p0L;    
      p[] = pL*f[] + p0*(1. - f[]);
      T[] = average_temperature (point, p[]);
    
      fE1[] = (pL + gamma1*PI1)/(gamma1 - 1.)*(f[] - frho1[]*b1) + frho1[]*q1;
      fE2[] = (1. - f[])*(p0/(gamma2 - 1.));
    }
  }
}

/**
We log the evolution of the bubble radius. */

event centroid (i += 20)
{
  double volume = 0.;
  foreach(reduction(+:volume))
    volume += dv()*(1. - f[]);
#if AXIS
  volume /= 2.;
#endif
  fprintf (stderr ,"%g %g\n", t/tau, pow(3.*volume,1./3.));
}

/**
Output of some statistics about the fields. */

event logfile (i++)
{
  stats sp = statsf (p), su = statsf (q.x), sT = statsf (T);
  if (i == 0)
    fprintf (stdout, "t dt max(p) max(T) max(u)\n");
  fprintf (stdout, "%g %g %g %g %g\n",
	   t/tau, dt/tau, sp.max, sT.max, su.max);
}

/**
On the fly movie generation. */

#if AXIS
event movie (t += 0.01*4737.81)
{
  view (fov = 12.5, quat = {0,0,-cos(M_PI/4.),cos(M_PI/4.)}, width = 640, height = 990);
  draw_vof ("f");
  squares ("p", min = 1., map = cool_warm);
  char s[80];
  sprintf (s, "t = %.2f", t/4737.81);
  mirror({0,1}) {
    draw_vof ("f");
    squares ("T", min = 1.0896, max = 2.1792, map = cool_warm);
    draw_string (s, pos = 2, size = 16, lc = {255,255,255}, lw = 4);
    draw_string ("Temperature", pos = 3, size = 25, lc = {255,255,255}, lw = 4);
    draw_string ("Pressure", size = 25, lc = {255,255,255}, lw = 4);
  }
  save ("T2.mp4");

  /**
  Saving dump files for post-processing. (Uncomment) */
  
#if 0
  char name[80];
  sprintf (name,"dump-%g",t/4737.81);
  dump (name, list = (scalar *){f,p,T});
#endif
}
#endif // AXIS

event ending (t = tend);
