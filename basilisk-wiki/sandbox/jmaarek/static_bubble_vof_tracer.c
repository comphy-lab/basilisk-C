#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "reduced.h"

#include "diffusion.h"

#include "navier-stokes/perfs.h"

#define MAXLEVEL 9

//#include "view.h"

#define VOFTHR 0.5


double D1 =0.1;
double D2 = 1.;
face vector Diff[];
scalar C1[];
double initial_C1;

double c_alpha = 0.001;

/**
The boundary conditions are slip lateral walls (the default) and
no-slip on the right and left walls. */

u.t[right] = dirichlet(0);
u.t[left]  = dirichlet(0);

/**
We make sure there is no flow through the top and bottom boundary,
otherwise the compatibility condition for the Poisson equation can be
violated. */

uf.n[bottom] = 0.;
uf.n[top] = 0.;

int main() {

  /**
  The domain will span $[0:2]\times[0:0.5]$ and will be resolved with
  $256\times 64$ grid points. */

  size (10.);
  init_grid (1 << (MAXLEVEL));

  rho1 = 1.0;
  rho2 = 0.01;

  TOLERANCE = 1e-4;

  run();
}

event init (t = 0) {
  //mask (y > 0.5 ? top : none);
  fraction (f, sq(x - 5.0) + sq(y) + sq(z) - sq(1.0));
  foreach()
        C1[] = c_alpha*f[];
  boundary ((scalar *){C1});
}

static scalar * interfaces1 = NULL;
event vof (i++) {
  C1.gradient = zero; //minmod2;
  C1.refine = refine_linear;
  C1.restriction = restriction_volume_average;
  f.tracers = (scalar *){C1};
  vof_advection ({f}, i);
  interfaces1 = interfaces, interfaces = NULL;
}

event tracer_advection (i++) {
  interfaces = interfaces1;
}

event tracer_diffusion (i++) {

  foreach_face() {
  double ff = (f[] + f[-1])/2.;
  Diff.x[] = fm.x[]*D1*D2/(D1*(1. - ff) + ff*D2);
  }
  boundary ((scalar *){Diff});
  foreach()
    C1[] = (f[] > VOFTHR ? C1[]/f[] : 0.0);
  boundary ({C1});
  scalar volume_correction[];
  scalar sharp_interface[];
  volume_correction.prolongation = volume_correction.refine = fraction_refine;
  sharp_interface.prolongation = sharp_interface.refine = fraction_refine;

  foreach(){
    volume_correction[] = cm[]*max(f[], VOFTHR);
    sharp_interface[] = (f[] < VOFTHR ? -1e12 : 0.0);
  }
  boundary ({volume_correction, sharp_interface});
  diffusion (C1, dt, Diff, theta = volume_correction, beta = sharp_interface);
  foreach()
    C1[] *= f[];
  boundary({C1});

}

event extract (t = 0; t += 1.; t <= 2.)
{
  double co = 0.0;
  foreach (reduction(+:co)){
  if (f[] >= VOFTHR)
    	co += C1[]*dv();
  }

  fprintf(fout,"%e %.12e %.12e %.12e\n",t, c_alpha - interpolate(C1,5.0,1.5), co, dt);
}