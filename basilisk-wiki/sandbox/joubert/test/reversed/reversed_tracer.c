/**
# Time-reversed VOF advection with tracer in a vortex

This is a modified version of the [reversed](http://www.basilisk.fr/src/test/reversed.c)
test case. We add a VOF tracer and a Godunov tracer (standard tracer) center on the right
side of the circle interface and test its advection.
There is a similar test case using Gerris
([see](http://gerris.dalembert.upmc.fr/gerris/tests/tests/shear.html)).
*/

#define ADAPT 0
#include "advection.h"
#include "vof.h"
#include "view.h"

#define VOFTHR 1e-10
/**
  The volume fraction is stored in scalar field `f` which is listed as
  an *interface* for the VOF solver. We advect tracer $C$ attached to the *interface* with
  the VOF scheme and $G$ with the default (diffusive) advection scheme. */

scalar f[], C[], G[];
scalar * interfaces = {f}, * tracers  = {G};

int MAXLEVEL;
double initial_C,initial_Cconc,initial_G;
FILE * fp = NULL;

/**
  We center the unit box on the origin and set a maximum timestep of 0.1 */

int main() {
  origin (-0.5, -0.5);
  fp = fopen ("tracer_stats.dat", "w");
  /**
    We then run the simulation for different levels of refinement. */

  for (MAXLEVEL = 5; MAXLEVEL <= 7; MAXLEVEL++) {
    init_grid (1 << MAXLEVEL);
    run();
    fprintf (fp,"\n\n");
  }
  fclose (fp);
}

/**
  The initial interface is a circle of radius 0.2 centered on
  (-0.2,-0.236338) (for historical reasons). We use the levelset
  function `circle()` to define this interface. The initial tracer distribution is
  a gaussian centered on the right side of the circle interface. We use the levelset
  function `gaussian()` to define the tracer distribution. 

  The period of the stretching cycle is set to 10, which will lead to
  not too strong stretching. Milder conditions can be obtained by decreasing it. */
#define circle(x,y) (sq(0.2) - (sq(x + 0.2) + sq(y + .236338)))
#define T 10.
#define gaussian(x,y) (exp(-100.*(sq(x) + sq(y + .236338))))

/**
  We define the levelset function $\phi$ on each vertex of the grid and
  compute the corresponding volume fraction field. */

event init (i = 0) {

  fraction (f, -circle(x,y));
#if ADAPT
  adapt_wavelet ({f}, (double[]){5e-3}, MAXLEVEL, list = {f});
#endif
  foreach() {
    C[]=gaussian(x,y)*f[];
    G[]=gaussian(x,y);
  }
  boundary ((scalar *){C,G});
  scalar Cconcentration[];
  foreach()
    Cconcentration[] = (f[] > VOFTHR ? C[]/f[] : C[]);
  boundary({Cconcentration});
  initial_C = statsf(C).sum;
  initial_G = statsf(G).sum;
  initial_Cconc = statsf(Cconcentration).sum;
}

/**
  This is used to enable the VOF scheme for the tracer advection*/
static scalar * interfaces1 = NULL;

event vof (i++) {

  C.gradient = zero; //minmod2;

  /**
    We associate the transport of $C$ with $f$ and transport
    all fields consistenTfy using the VOF scheme. */

  f.tracers = (scalar *){C};
  vof_advection ({f}, i);

  /**
    We set the list of interfaces to NULL so that the default *vof()*
    event does nothing (otherwise we would transport $f$ twice). */

  interfaces1 = interfaces, interfaces = NULL;
}

/**
  We set the list of interfaces back to its default value. */

event tracer_advection (i++) {
  interfaces = interfaces1;
}

event velocity (i++) {

  /**
    This event defines the velocity field.

    On trees we first adapt the grid so that the estimated error on
    the volume fraction is smaller than $5\times 10^{-3}$. We limit the
    resolution at `MAXLEVEL` and we refine the volume fraction field
    `f` and the tracers `C` and `G`. */

#if ADAPT
  adapt_wavelet ({f,C,G}, (double[]){5e-3, 1e-2, 1e-2}, MAXLEVEL, list = {f,C,G});
#endif

  /**
    The velocity field is defined through a streamfunction $\psi$, defined
    on the vertices of the grid. */

  vertex scalar psi[];
  foreach_vertex()
    psi[] = - 1.5*sin(2.*pi*t/T)*sin((x + 0.5)*pi)*sin((y + 0.5)*pi)/pi;

  /**
    We can then differentiate the streamfunction to get the velocity
    components. This guarantees that the velocity field is exactly
    non-divergent. */

  trash ({u});
  struct { double x, y; } f = {-1.,1.};
  foreach_face()
    u.x[] = f.x*(psi[0,1] - psi[])/Delta;
  boundary ((scalar *){u});
}

/**
  To compute the error, we reinitialise field `e` at the end of the
  simulation with the initial shape and compute the difference with the
  final shape. We output the norms as functions of the maximum
  resolution `N`. */

event field (t = T) {
  scalar e[], eC[], eG[];
  fraction (e, -circle(x,y));
  foreach() {
    eC[]=gaussian(x,y)*e[];
    eG[]=gaussian(x,y);
  }
  boundary ((scalar *){eC,eG});
  foreach() {
    eC[] -= C[];
    eG[] -= G[];
  }
  boundary ((scalar *){eC,eG});
  norm n = normf (eC);
  norm o = normf (eG);
  fprintf (stderr, "%d %g %g %g %g %g %g\n", N, n.avg, n.rms, n.max, o.avg, o.rms, o.max);
}

/**
  We save some snapshot of the simulation (but only on the finest grid considered). */

event shape (t += T/2.) {
  if (N == 128) {
    scalar Cconcentration[];
    foreach()
      Cconcentration[] = (f[] > VOFTHR ? C[]/f[] : C[]);
    boundary({Cconcentration});
    char name[80];
    sprintf (name, "snapshot-%g", t);
    dump (name);
  }
}

/**
  At $t+=T/4$ we check some stats on the tracer. the sum, min and max
  values of the volume fraction field. The sum must be constant to
  within machine precision and the volume fraction should be bounded by
  zero and one. */

event statistic (t +=T/10.) {
  scalar Cconcentration[];
  foreach()
    Cconcentration[] = (f[] > VOFTHR ? C[]/f[] : C[]);
  boundary({Cconcentration});
  stats u = statsf (C);
  stats v = statsf (G);
  stats s = statsf (Cconcentration);
  double variation_C = (initial_C-u.sum)/initial_C;
  double variation_Cconc = (initial_Cconc-s.sum)/initial_Cconc;
  double variation_G = (initial_G-v.sum)/initial_G;
  double Cliq = 0., Cgas = 0., Cconcliq = 0., Cconcgas = 0.,Gliq = 0., Ggas = 0., volgas = 0., volliq = 0.;
  foreach(reduction(+:Cliq) reduction(+:Cgas) reduction(+:Cconcliq) reduction(+:Cconcgas) reduction(+:Gliq) reduction(+:Ggas) reduction(+:volgas) reduction(+:volliq)) {
    if (f[] >= (1.-VOFTHR)) {
      Cliq += C[]*dv();
      Cconcliq += Cconcentration[]*dv();
      Gliq += G[]*dv();
      volliq += dv();
    }
    else if (f[] <= VOFTHR) {
      Cgas += C[]*dv();
      Cconcgas += Cconcentration[]*dv();
      Ggas += G[]*dv();
      volgas += dv();
    }
  }
  fprintf (fp, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n", 
      t, variation_C, variation_Cconc, variation_G, volgas, volliq, Cliq/volliq, Cgas/volgas, Cconcliq/volliq, Cconcgas/volgas, interface_area(f), Gliq/volliq, Ggas/volgas, dt, s.min, s.max, s.sum/s.volume);
}

/**
  Some movies are generated. */

event mov (t += T/200.) {
  if (N == 128) {
    scalar Cconcentration[];
    foreach()
      Cconcentration[] = (f[] > VOFTHR ? C[]/f[] : C[]);
    boundary({Cconcentration});
    clear();
    view (fov = 24, quat = {0,0,0,1}, tx = 0, ty = 0, bg = {1,1,1}, width = 600, height = 600);
    box ();
    draw_vof ("f", lw=2);
    squares ("C", min=0, max=1);
    save ("C.mp4");
    clear();
    box ();
    draw_vof ("f", lw=2);
    squares ("Cconcentration", min=0, max=1);
    save ("Cconc.mp4");
    clear();
    box ();
    draw_vof ("f", lw=2);
    squares ("G", min=0, max=1);
    save ("G.mp4");
    clear();
    box ();
    cells();
    save ("cells.mp4");
  }
}

/**
## Results

VOF tracer C field | Tracer G field
:-------------------------:|:-------------------------:
![VOF tracer C](reversed_tracer/C.mp4)(width="800" height="600") | ![Tracer G](reversed_tracer/G.mp4)(width="800" height="600") 

~~~gnuplot maximimum, average and minimum deviation of the concentration value in function of time
T=10
reset
set xlabel 'time in T unit'
set ylabel 'concentration value'
plot 'tracer_stats.dat' i 2:2 u ($1/T):16 w l t 'maximum',\
'tracer_stats.dat' i 2:2 u ($1/T):17 w l t 'average',\
'tracer_stats.dat' i 2:2 u ($1/T):15 w l t 'minimum'
~~~

~~~gnuplot Comparison of convergence of error between initial value and final value of $C$ and $G$ in function of cells number
reset
ftitle(a,b) = sprintf("%.0f/x^{%4.2f}", exp(a), -b)
f(x)=a+b*x
fit f(x) 'log' u (log($1)):(log($4)) via a,b
f2(x)=a2+b2*x
fit f2(x) 'log' u (log($1)):(log($7)) via a2,b2
C(x)=ac+bc*x
fit C(x) 'log' u (log($1)):(log($2)) via ac,bc
C2(x)=ac2+bc2*x
fit C2(x) 'log' u (log($1)):(log($5)) via ac2,bc2
set xlabel 'Maximum resolution'
set ylabel 'Maximum error'
set logscale 
set xrange [16:256]
set xtics 16,2,256
set grid ytics
set cbrange [1:1]
plot 'log' u 1:4 t 'max (C)', exp(f(log(x))) t ftitle(a,b), \
'log' u 1:7 t 'max (G)', exp(f2(log(x))) t ftitle(a2,b2), \
'log' u 1:2 t 'norm1 (C)', exp(C(log(x))) t ftitle(ac,bc), \
'log' u 1:5 t 'norm1 (G)', exp(C2(log(x))) t ftitle(ac2,bc2)
~~~

~~~gnuplot Mass conservation (volume sum) of C in function of time
reset
T=10
set xlabel 'time in T unit'
set ylabel 'C relative difference between t=0 and t=t'
plot 'tracer_stats.dat' i 0:0 u ($1/T):2 w lp t 'lvl 5',\
'tracer_stats.dat' i 1:1 u ($1/T):2 w lp t 'lvl 6',\
'tracer_stats.dat' i 2:2 u ($1/T):2 w lp t 'lvl 7'
~~~

~~~gnuplot Mass conservation (volume sum) of C in function of time
reset
T=10
set xlabel 'time in T unit'
set ylabel 'Concentration relative difference between t=0 and t=t'
plot 'tracer_stats.dat' i 0:0 u ($1/T):3 w lp t 'lvl 5',\
'tracer_stats.dat' i 1:1 u ($1/T):3 w lp t 'lvl 6',\
'tracer_stats.dat' i 2:2 u ($1/T):3 w lp t 'lvl 7'
~~~

~~~gnuplot Comparison of total mass conservation (volume sum) of C in function of time
reset
T=10
set xlabel 'time in T unit'
set ylabel 'absolute relative difference between t=0 and t=t'
set logscale y
plot 'tracer_stats.dat' i 2:2 u ($1/T):(abs($2)) w lp t 'C lvl 7',\
'tracer_stats.dat' i 2:2 u ($1/T):(abs($3)) w lp t 'Cconc lvl 7',\
'tracer_stats.dat' i 2:2 u ($1/T):(abs($4)) w lp t 'G lvl 7'
~~~
*/
