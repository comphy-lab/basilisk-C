/**
# Advection-diffusion of a VOF tracer

In this test case we compare the advection and diffusion of passive scalars $C$
transported by VOF scheme and passive scalar $G$ transported with the default
scheme. This test is inspired by Quentin's
[test](http://www.basilisk.fr/sandbox/qmagdelaine/phase_change/2_diffusion_within_a_domain/translation_mwe.c).*/

#define ADAPT 0
#define ADVECTION 1
#define DIFFUSION 1
#define NODIFF 0
#define VOLCORRECTION 0

#if ADVECTION
#include "grid/multigrid.h"
#include "../../qmagdelaine/phase_change/advection_Q.h"
#include "vof.h"
#else
#include "grid/multigrid.h"
#include "run.h"
#include "fractions.h"
double dt;
#endif

#if DIFFUSION
#include "diffusion.h"
#if NODIFF
#define D1 5e-3
#define D2 0.
#ifndef D
# define D(f)  (clamp(f,0.,1.)*(D1 - D2) + D2)
#endif
#else
#define D1 5e-3
#endif
mgstats mgd1, mgd2;
face vector Diff[];
#endif
#include "view.h"

/**
  The volume fraction is stored in scalar field `f` which is listed as
  an *interface* for the VOF solver. We advect tracer $C$ with either
  the default (diffusive) advection scheme of the advection solver or
  the VOF scheme. */

#define VOFTHR 1e-10
#define VELOCITY 0.2
scalar f[], C[], G[];

int MAXLEVEL;
double initial_C, initial_Cconc, initial_G;
FILE * fp = NULL;

scalar * interfaces = {f};
#if ADVECTION
scalar * tracers  = {G};
uf.n[right] = dirichlet(fm.n[ghost]*VELOCITY);
uf.n[left] = dirichlet(fm.n[ghost]*VELOCITY);
#endif

/**
  We center the unit box on the origin and set a maximum timestep of 5e-3 */

int main() {
  origin (-0.5, -0.5);
#if DIFFUSION
  DT=1e-3;
  TOLERANCE=1e-4;
#endif
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
  in the middle of the box. We use the levelset
  function `circle()` to define this interface. 
  The maxtime is set to 1. */
#define circle(x,y) (sq(0.2) - (sq(x) + sq(y)))
#define T 1.

/**
  We define the levelset function $\phi$ on each vertex of the grid and
  compute the corresponding volume fraction field. */

event init (i = 0) {
  fraction (f, -circle(x,y));
  fraction (G, -circle(x,y));
  foreach()
      C[]=f[];
  boundary ((scalar *){C});
scalar Cconcentration[];
  foreach()
    Cconcentration[] = (f[] > VOFTHR ? C[]/f[] : 0.);
  boundary({Cconcentration});
#if ADVECTION
  foreach_face(x)
    uf.x[] = fm.x[]*VELOCITY;
  boundary((scalar *){uf});
#endif
  initial_C = statsf(C).sum;
  initial_Cconc = statsf(Cconcentration).sum;
  initial_G = statsf(G).sum;
}

#if ADVECTION
event stability (i++) {
  foreach_face(x)
    uf.x[] = fm.x[]*VELOCITY;
  boundary((scalar *){uf});
}

/**
  This is used to enable the VOF scheme for the tracer advection*/
static scalar * interfaces1 = NULL;
event vof (i++) {
  /**
    We associate the transport of $C$ with $f$ and transport
    all fields consistenTfy using the VOF scheme. */

  C.gradient = zero; //minmod2;
#if ADAPT
  C.refine = refine_linear;
  C.restriction = restriction_volume_average;
#endif
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
#endif

#if DIFFUSION 
#if ADVECTION 
event tracer_diffusion (i++) {
  /**
    Here we compute  the Diff value for two phase case with average and
    conditions on $f$ the face vector value. */
#if NODIFF
  foreach_face() {
    if (f[] > 1.- VOFTHR && f[-1] >1.-VOFTHR)
      Diff.x[]=D1*fm.x[];
    else 
      Diff.x[]= D2*fm.x[];
  }
#else
  foreach_face()
    Diff.x[]=D1*fm.x[];
#endif
  boundary ((scalar *){Diff});
  foreach()
    C[] = (f[] > VOFTHR ? C[]/f[] : 0.);
  boundary({C});
#if VOLCORRECTION
  scalar volume_correction[];
#if TREE
  volume_correction.prolongation = volume_correction.refine = fraction_refine;
#endif
  foreach()
    volume_correction[] = cm[]*max(f[], VOFTHR);
  boundary ({volume_correction});
  mgd1 = diffusion (C, dt, Diff, theta = volume_correction);
#else
  mgd1 = diffusion (C, dt, Diff);
#endif
  foreach()
    C[] *= f[];
  boundary({C});
  mgd2 = diffusion (G, dt, Diff); 
}
#else
event time_integration (i++) {
  dt = dtnext (DT);
  /**
    Here we compute  the Diff value for two phase case with average and
    conditions on $f$ the face vector value. */
  const face vector Diff[] = {D1, D1};
  foreach()
    C[] = (f[] > VOFTHR ? C[]/f[] : 0.);
  boundary({C});
  mgd1 = diffusion (C, dt, Diff); 
  foreach()
    C[] *= f[];
  boundary({C});
  mgd2 = diffusion (G, dt, Diff); 
}
#endif
#endif

#if ADAPT
event adapt(i++) {
  adapt_wavelet ({f,C,G}, (double[]){5e-3, 1e-2, 1e-2}, MAXLEVEL, list = {f,C,G});
}
#endif

/**
  At $t+=T/4$ we check some stats on the tracer. the sum, min and max
  values of the volume fraction field. The sum must be constant to
  within machine precision and the volume fraction should be bounded by
  zero and one. */

event statistic (t +=T/50.; t <=T) {
  scalar Cconcentration[];
  foreach()
    Cconcentration[] = (f[] > VOFTHR ? C[]/f[] : 0.);
  boundary({Cconcentration});
  stats s = statsf (Cconcentration);
  stats u = statsf (C);
  stats v = statsf (G);
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
#if DIFFUSION 
  fprintf (fp, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %d %d\n", 
      t, variation_C, variation_Cconc, variation_G, volgas, volliq, Cliq/volliq, Cgas/volgas, Cconcliq/volliq, Cconcgas/volgas, interface_area(f), Gliq/volliq, Ggas/volgas, dt, mgd1.i, mgd2.i);
#else
  fprintf (fp, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g\n", 
      t, variation_C, variation_Cconc, variation_G, volgas, volliq, Cliq/volliq, Cgas/volgas, Cconcliq/volliq, Cconcgas/volgas, interface_area(f), Gliq/volliq, Ggas/volgas, dt);
#endif
}

/**
  We also output the shape of the reconstructed interface at regular
  intervals and some snapshot of the simulation (but only on the finest grid considered). */

event shape (t += T/2.) {
  if (N == 128) {
scalar Cconcentration[];
  foreach()
    Cconcentration[] = (f[] > VOFTHR ? C[]/f[] : 0.);
  boundary({Cconcentration});
    char name[80];
    sprintf (name, "snapshot-%g", t);
    dump (name);
  }
}

/**
  Some movies are generated. */

event mov (t += T/200.) {
  if (N == 128) {
scalar Cconcentration[];
  foreach()
    Cconcentration[] = (f[] > VOFTHR ? C[]/f[] : 0.);
  boundary({Cconcentration});
    clear();
    translate (z = 0.01)
      box ();
    draw_vof ("f", lw=2);
    squares ("C", min=0, max=1);
    save ("C.mp4");
    clear();
      box ();
    draw_vof ("f", lw=2);
    squares ("Cconcentration", min=0, max=1);
    save ("Cconcentration.mp4");
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

![G field advection + diffusion](tracer_advection_diffusion_mwe/G.mp4)(width="800" height="600")

![Concentration field advection + diffusion](tracer_advection_diffusion_mwe/Cconcentration.mp4)(width="800" height="600")

~~~gnuplot Comparison mass conservation (volume sum) of C, outside the drop.
reset
T=1
set key bottom left
set xlabel 'time in T unit'
set ylabel 'Tracer quantity outside drop [-]'
plot 'dtracer_stats.dat' i 2:2 u ($1/T):7 w lp t 'C diffusion',\
'tracer_stats.dat' i 2:2 u ($1/T):7 w lp t 'C advection-diffusion',\
'tracer_stats.dat' i 2:2 u ($1/T):12 w lp t 'G advection-diffusion'
~~~

~~~gnuplot Comparison mass conservation (volume sum) of C, inside the drop
reset
T=1
set xlabel 'time in T unit'
set ylabel 'Tracer quantity inside drop [-]'
set key top left
plot 'dtracer_stats.dat' i 2:2 u ($1/T):8 w lp t 'C diffusion',\
'tracer_stats.dat' i 2:2 u ($1/T):8 w lp t 'C advection-diffusion',\
'tracer_stats.dat' i 2:2 u ($1/T):13 w lp t 'G advection-diffusion'
~~~

~~~gnuplot Comparison mass conservation (volume sum) of Cconc, outside the drop
reset
T=1
set xlabel 'time in T unit'
set ylabel 'concentration quantity outside drop [-]'
plot 'dtracer_stats.dat' i 2:2 u ($1/T):9 w lp t 'Cconc diffusion',\
'tracer_stats.dat' i 2:2 u ($1/T):9 w lp t 'Cconc advection-diffusion'
~~~

~~~gnuplot Comparison mass conservation (volume sum) of Cconc, inside the drop
reset
T=1
set xlabel 'time in T unit'
set ylabel 'Concentration quantity inside drop [-]'
plot 'dtracer_stats.dat' i 2:2 u ($1/T):10 w lp t 'Cconc diffusion',\
'tracer_stats.dat' i 2:2 u ($1/T):10 w lp t 'Cconc advection-diffusion'
~~~

~~~gnuplot Comparison of total mass conservation (volume sum) of C in function of time
reset
T=1
set xlabel 'time in T unit'
set ylabel 'C relative difference between t=0 and t=t'
set key top left
plot 'dtracer_stats.dat' i 2:2 u ($1/T):2 w lp t 'C diffusion',\
'tracer_stats.dat' i 2:2 u ($1/T):2 w lp t 'C advection-diffusion',\
'tracer_stats.dat' i 2:2 u ($1/T):4 w lp t 'G advection-diffusion'
~~~

~~~gnuplot Comparison of total mass conservation (volume sum) of Cconcentration in function of time
reset
T=1
set xlabel 'time in T unit'
set ylabel 'Cconcentration relative difference between t=0 and t=t'
set key bottom right
plot 'dtracer_stats.dat' i 2:2 u ($1/T):3 w lp t 'Cconc diffusion',\
'tracer_stats.dat' i 2:2 u ($1/T):3 w lp t 'Cconc advection-diffusion'
~~~

~~~gnuplot Comparison mass conservation (volume sum) of C, outside the circle
reset
T=1
set key bottom left
set xlabel 'time in T unit'
set ylabel 'Tracer total quantity outside circle [-]'
plot 'atracer_stats.dat' i 2:2 u ($1/T):7 w lp t 'C advection',\
'dtracer_stats.dat' i 2:2 u ($1/T):7 w lp t 'C diffusion',\
'nodifftracer_stats.dat' i 2:2 u ($1/T):7 w lp t 'C advection-diffusion NODIFF',\
'tracer_stats.dat' i 2:2 u ($1/T):12 w lp t 'G advection-diffusion'
~~~

~~~gnuplot Comparison mass conservation (volume sum) of C, inside the circle
reset
T=1
set key top left
set xlabel 'time in T unit'
set ylabel 'Tracer total quantity outside circle [-]'
plot 'atracer_stats.dat' i 2:2 u ($1/T):8 w lp t 'C advection',\
'dtracer_stats.dat' i 2:2 u ($1/T):8 w lp t 'C diffusion',\
'nodifftracer_stats.dat' i 2:2 u ($1/T):8 w lp t 'C advection-diffusion NODIFF',\
'tracer_stats.dat' i 2:2 u ($1/T):13 w lp t 'G advection-diffusion'
~~~

*/
