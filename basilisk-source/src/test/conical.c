/**
# Runup of a solitary wave on a conical island

In this test case, we compare experimental data with solutions of the
[Saint-Venant](/src/saint-venant.h),
[Green-Naghdi](/src/green-naghdi.h) and
[layered](/src/layered/hydro.h) equations. For details on the test
case and experimental data many references are available (e.g. [Hou el
al, 2013](/src/references.bib#hou2013), [Nikolos and Delis,
2009](/src/references.bib#nikolos2009), [Lannes and Marche,
2014](/src/references.bib#lannes2014), [Kazolea et al,
2012](/src/references.bib#kazolea2012) etc...).

We also test both the explicit and the implicit versions of the
Saint-Venant solver. */

#if SAINT_VENANT
# if IMPLICIT
#  include "saint-venant-implicit.h"
# else
#  include "saint-venant.h"
# endif
# define MGD 0
#elif ML
# include "layered/hydro.h"
# include "layered/nh.h"
 // fixme: remap does not work
 // # include "layered/remap.h"
# define MGD mgp.i
#else
# include "green-naghdi.h"
# define MGD mgD.i
#endif

/**
The domain is 27.6 metres squared, resolved with a maximum of 256^2^
grid points. */

#define MAXLEVEL 8
#define MINLEVEL 0

int main()
{
  G = 9.81;
  size (27.6);

  /**
  We set an "acoustic" CFL for the implicit solver, so that the
  timestep is comparable to that of the explicit Saint-Venant
  solver. With larger timesteps the scheme is too dissipative. */
  
#if IMPLICIT
  CFLa = 0.25;
#endif

#if ML
  breaking = 0.07;
  // nl = 1;
  CFL_H = 0.5;
  TOLERANCE = 1e-4;
#endif
  
  init_grid (1 << MAXLEVEL);
  run();
}

/**
We reproduce case B i.e. a relative amplitude of the solitary wave of
0.09. */

double H0 = 0.32, H = 0.32*0.09, T = 2.45;
#define C sqrt(G*(H0 + H))

double sech2 (double x)
{
  double e = exp (-x);
  double s = 2.*e/(1. + e*e);
  return s*s;
}

double eta0 (double t)
{
  return H*sech2 (sqrt(3.*H/(4.*cube(H0)))*C*(t - T));
}

double uleft (double t)
{
  double e = eta0 (t);
  return C*e/(H0 + e);
}

/**
This is the definition of the topography of the conical island. */

double island (double x, double y)
{
  double r0 = 3.6, r1 = 1.1, H0 = 0.625;
  x -= 25.92/2.; y -= L0/2.;
  double r = sqrt(x*x + y*y);
  return r > r0 ? 0. : r < r1 ? H0 : (r - r0)/(r1 - r0)*H0;
}

/**
On the tree we need to make sure that refined cells are
initialised with the correct topography. */

#if TREE
void refine_zb (Point point, scalar zb)
{
  foreach_child()
    zb[] = island (x, y);
}
#endif

event init (i = 0)
{
  
  /**
  We use the definition of the solitary wave to impose conditions on the
  left side (i.e. the "wave paddle"). */

#if ML
  eta[left] = dirichlet (H0 + eta0(t));
  h[left] = dirichlet((H0 + eta0(t))/nl);
  u.n[left] = uleft (t)/nl;
  u.t[left] = 0.;
#else
  h[left] = H0 + eta0(t);
#if IMPLICIT
  q.n[left] = (H0 + eta0(t))*uleft (t);
  q.t[left] = 0.;
#else
  u.n[left] = uleft (t);
  u.t[left] = 0.;
#endif
#endif
  
#if TREE
  /**
  We setup the refinement function and choose to conserve surface
  elevation rather than volume. */

  zb.refine = refine_zb;
  default_sea_level = H0;
  conserve_elevation();
#endif

  /**
  These are the initial conditions i.e. conical island and
  lake-at-rest. */

  foreach() {
    zb[] = island (x, y);
#if ML
    foreach_layer()
      h[] = max(0., H0 - zb[])/nl;
#else
    h[] = max(0., H0 - zb[]);
#endif
  }
}

/**
We define the wave gauges corresponding with the experimental
locations. */

Gauge gauges[] = {
  {"WG3",   6.82, 13.05},
  {"WG6",   9.36, 13.80},
  {"WG9",  10.36, 13.80},
  {"WG16", 12.96, 11.22},
  {"WG22", 15.56, 13.80},
  {NULL}
};

/**
We output snapshots of fields at various times. */

event outputfile (t = {9, 12, 13, 14, 20})
{
  static int nf = 0;
  printf ("file: conical-%d\n", nf);
#if ML
  scalar H[];
  foreach() {
    H[] = 0.;
    foreach_layer()
      H[] += h[];
  }
#else
  scalar H = h;
#endif
  output_field ({H,zb}, stdout, N, linear = true);

  /**
  Here we output the level of refinement of the grid. */

  printf ("file: level-%d\n", nf);
  scalar l[];
  foreach()
    l[] = level;
  output_field ({l}, stdout, N);
  nf++;
}

/**
Which gives the following wave (left column) and mesh (right column)
evolutions. 

~~~gnuplot Evolution of free surface (left) and mesh (right)
set term @PNG enhanced size 640,1120 font ",8"
set output 'free.png'

! awk '{ if ($1 == "file:") file = $2; else print $0 > file; }' < out

unset key
unset xtics
unset ytics
  
dry=1e-4
  
set size ratio -1
set pm3d
set pm3d map
# set contour base
set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25 0 0.5647 1, \
0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392, 0.625 1 0.9333 0,	\
0.75 1 0.4392 0, 0.875 0.9333 0 0, 1 0.498 0 0 )
  
set multiplot layout 4,2 scale 1.35,1.35
set cbrange [0.3:0.4]
splot './conical-0' u 1:2:($3>dry?$3+$4:1e1000)
set cbrange [4:8]
splot './level-0'
set cbrange [0.3:0.4]
splot './conical-1' u 1:2:($3>dry?$3+$4:1e1000)
set cbrange [4:8]
splot './level-1'
set cbrange [0.3:0.4]
splot './conical-2' u 1:2:($3>dry?$3+$4:1e1000)
set cbrange [4:8]
splot './level-2'
set cbrange [0.3:0.4]
splot './conical-3' u 1:2:($3>dry?$3+$4:1e1000)
set cbrange [4:8]
splot './level-3'
unset multiplot
~~~
*/

event logfile (i++) {

  /**
  We output various diagnostics. */

#if ML
  scalar H[];
  foreach() {
    H[] = 0.;
    foreach_layer()
      H[] += h[];
  }
#else
  scalar H = h;
#endif
  stats s = statsf (H);

  #if ML
  norm n = normf (u.x);
  #elif IMPLICIT
  scalar u[];
  foreach()
    u[] = h[] > dry ? q.x[]/h[] : 0.;
  norm n = normf (u);
  #else
  norm n = normf (u.x);
  #endif
  
  if (i == 0)
    fprintf (stderr, "t i h.min h.max h.sum u.x.rms u.x.max dt mgD.i\n");
  fprintf (stderr, "%g %d %.8f %g %.8f %g %.4g %g %d\n", 
	   t, i, s.min, s.max, s.sum, n.rms, n.max, dt, MGD);

  /**
  Here we output surface elevation at the various gauges... */

  #if IMPLICIT
  scalar eta[];
  foreach()
    eta[] = h[] > dry ? h[] + zb[] : nodata;
  #endif
  
  output_gauges (gauges, {eta});

  /**
  We also use a simple implicit scheme to implement quadratic bottom
  friction i.e.
  $$
  \frac{d\mathbf{u}}{dt} = - C_f|\mathbf{u}|\frac{\mathbf{u}}{h}
  $$
  with $C_f=10^{-4}$. */
  
  foreach() {
#if IMPLICIT
    double a = h[] < dry ? HUGE : 1. + 1e-4*dt*norm(q)/sq(h[]);
    foreach_dimension()
      q.x[] /= a;
#else
    double a = h[] < dry ? HUGE : 1. + 1e-4*dt*norm(u)/h[];
    foreach_dimension()
      u.x[] /= a;
#endif
  }
}

/**
... and compare it with experimental data for the Green-Naghdi solver
(blue), the Saint-Venant solver (magenta) and layered solver
(black). The results for the Green-Naghdi solver are very similar to
those of [Kazolea et al, 2012](/src/references.bib#kazolea2012),
Figure 14 and [Lannes and Marche,
2014](/src/references.bib#lannes2014), Figure 20 and significantly
better than the results of the Saint-Venant solver (see also
e.g. Figure 18 of [Hou et al, 2013](/src/references.bib#hou2013) and
Figure 19 of [Nikolos and Delis,
2009](/src/references.bib#nikolos2009)) which have too sharp fronts.

~~~gnuplot Timeseries of surface elevation. Numerical results and experimental data (symbols).
reset
set term svg enhanced size 640,800 font ",10"
set multiplot layout 5,1 scale 1.,1.
set tmargin 0.5
set xrange [3:20]
t0=22.
h0=0.32
plot '../ts2b.txt' u ($1-t0):($4-0.001) pt 7 ps 0.25 t 'WG3',	\
'../conicalsv/WG3' u 1:($2-h0) w l lt 4 t 'Saint-Venant explicit',	\
'../conical-implicit/WG3' u 1:($2-h0) w l lt 5 t 'Saint-Venant implicit', \
'./WG3' u 1:($2-h0) w l lt 3 t 'Green-Naghdi', \
'../conical-ml/WG3' u 1:($2-h0) w l lw 2 lt -1 t 'layered'
plot '../ts2b.txt' u ($1-t0):($6-0.0004) pt 7 ps 0.25 t 'WG6',	\
'../conicalsv/WG6' u 1:($2-h0) w l lt 4 t '',			\
'../conical-implicit/WG6' u 1:($2-h0) w l lt 5 t '',		\
'./WG6' u 1:($2-h0) w l lt 3 t '', \
'../conical-ml/WG6' u 1:($2-h0) w l lw 2 lt -1 t ''
plot '../ts2b.txt' u ($1-t0):($7) pt 7 ps 0.25 t 'WG9',	\
'../conicalsv/WG9' u 1:($2-h0) w l lt 4 t '',		\
'../conical-implicit/WG9' u 1:($2-h0) w l lt 5 t '',	\
'./WG9' u 1:($2-h0) w l lt 3 t '', \
'../conical-ml/WG9' u 1:($2-h0) w l lw 2 lt -1 t ''
plot '../ts2b.txt' u ($1-t0):($8-0.0017) pt 7 ps 0.25 t 'WG16',	\
'../conicalsv/WG16' u 1:($2-h0) w l lt 4 t '',		\
'../conical-implicit/WG16' u 1:($2-h0) w l lt 5 t '',		\
'./WG16' u 1:($2-h0) w l lt 3 t '', \
'../conical-ml/WG16' u 1:($2-h0) w l lw 2 lt -1 t ''
plot '../ts2b.txt' u ($1-t0):($9+0.0015) pt 7 ps 0.25 t '',	\
'../conicalsv/WG22' u 1:($2-h0) w l lt 4 t 'WG22',		\
'../conical-implicit/WG22' u 1:($2-h0) w l lt 5 t '',	\
'./WG22' u 1:($2-h0) w l lt 3 t '', \
'../conical-ml/WG22' u 1:($2-h0) w l lw 2 lt -1 t ''
unset multiplot
  
! rm -f conical-?
~~~

We adapt the mesh based on the wavelet-estimated error on the free
surface elevation. */

#if TREE
event adapt (i++) {
  scalar eta1[];
  foreach() {
#if ML
    double H = 0.;
    foreach_layer()
      H += h[];
    eta1[] = H > dry ? zb[] + H : 0.;
#else
    eta1[] = h[] > dry ? zb[] + h[] : 0.;
#endif
  }

  astats s = adapt_wavelet ({eta1}, (double[]){3e-4}, MAXLEVEL, MINLEVEL);
  fprintf (stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
}
#endif
