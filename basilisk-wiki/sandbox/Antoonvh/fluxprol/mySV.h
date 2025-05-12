/**
# A solver for the Saint-Venant equations

This is an adaptation of [`saint-venant.h`](/src/saint-venant.h).

1. Centered gradients are computed without halo-ghost values.

 */
#include "mytree-utils.h"
scalar zb[], h[], eta[];
vector u[];

double G = 1.;
double dry = 1e-10;

#if !LAYERS
int nl = 1;
#endif

#include "multilayer.h"
#include "predictor-corrector.h"

scalar * evolving = NULL;

trace
static void advance_saint_venant (scalar * output, scalar * input, 
				  scalar * updates, double dt) {
  scalar hi = input[0], ho = output[0], dh = updates[0];
  vector * uol = (vector *) &output[1];
  foreach() {
    double hold = hi[];
    ho[] = hold + dt*dh[];
    eta[] = zb[] + ho[];
    if (ho[] > dry) {
      for (int l = 0; l < nl; l++) {
        vector uo = vector(output[1 + dimension*l]);
      	vector ui = vector(input[1 + dimension*l]),
	  dhu = vector(updates[1 + dimension*l]);
	foreach_dimension()
	  uo.x[] = (hold*ui.x[] + dt*dhu.x[])/ho[];
      }
      if (nl > 1)
	vertical_viscosity (point, ho[], uol, dt);
    }
    else // dry
      for (int l = 0; l < nl; l++) {
        vector uo = vector(output[1 + dimension*l]);
        foreach_dimension()
	  uo.x[] = 0.;
      }
  }
  scalar * list = list_concat ({ho, eta}, (scalar *) uol);
  multigrid_restriction (list);
  free (list);
}

#if TREE
static void refine_eta (Point point, scalar eta) {
  foreach_child()
    eta[] = zb[] + h[];
}

static void restriction_eta (Point point, scalar eta) {
  eta[] = zb[] + h[];
}
#endif

#include "riemann.h"

trace
double update_saint_venant (scalar * evolving, scalar * updates, double dtmax)
{
  scalar h = evolving[0], dh = updates[0];
  vector u = vector(evolving[1]);
  face vector Fh[], S[];
  tensor Fq[];
  vector gh[], geta[];
  tensor gu[];
  for (scalar s in {gh, geta, gu}) {
    s.gradient = zero;
    #if TREE
      s.prolongation = refine_linear;
    #endif
  }
  /**
Two lines are modified here `gradients()` $\rightarrow$ `my_gradients()`$
   */
  my_gradients ({h, eta, u}, {gh, geta, gu});
  for (int l = 0; l < nl; l++) {
    vector u = vector (evolving[1 + dimension*l]);
    if (l > 0)
      my_gradients ((scalar *) {u}, (vector *) {gu});
    foreach_face (reduction (min:dtmax)) {
      double hi = h[], hn = h[-1];
      if (hi > dry || hn > dry) {
	double dx = Delta/2.;
	double zi = eta[] - hi;
	double zl = zi - dx*(geta.x[] - gh.x[]);
	double zn = eta[-1] - hn;
	double zr = zn + dx*(geta.x[-1] - gh.x[-1]);
	double zlr = max(zl, zr);
	
	double hl = hi - dx*gh.x[];
	double up = u.x[] - dx*gu.x.x[];
	double hp = max(0., hl + zl - zlr);
	
	double hr = hn + dx*gh.x[-1];
	double um = u.x[-1] + dx*gu.x.x[-1];
	double hm = max(0., hr + zr - zlr);
	
	/**
	#### Riemann solver
	
	We can now call one of the approximate Riemann solvers to get
	the fluxes. */
	
	double fh, fu, fv;
	kurganov (hm, hp, um, up, Delta*cm[]/fm.x[], &fh, &fu, &dtmax);
	fv = (fh > 0. ? u.y[-1] + dx*gu.y.x[-1] : u.y[] - dx*gu.y.x[])*fh;
	
	/**
	#### Topographic source term
      
	In the case of adaptive refinement, care must be taken to ensure
	well-balancing at coarse/fine faces (see [notes/balanced.tm]()). */
	
        #if TREE
	if (is_prolongation(cell)) {
	  hi = coarse(h);
	  zi = coarse(zb);
	}
	if (is_prolongation(neighbor(-1))) {
	  hn = coarse(h,-1);
	  zn = coarse(zb,-1);
	}
        #endif
	
	double sl = G/2.*(sq(hp) - sq(hl) + (hl + hi)*(zi - zl));
	double sr = G/2.*(sq(hm) - sq(hr) + (hr + hn)*(zn - zr));
	
	/**
	#### Flux update */
      
	Fh.x[]   = fm.x[]*fh;
	Fq.x.x[] = fm.x[]*(fu - sl);
	S.x[]    = fm.x[]*(fu - sr);
	Fq.y.x[] = fm.x[]*fv;
      }
      else // dry
	Fh.x[] = Fq.x.x[] = S.x[] = Fq.y.x[] = 0.;
    }

    boundary_flux ({Fh, S, Fq});

  
    vector dhu = vector(updates[1 + dimension*l]);
    foreach() {
      double dhl =
	layer[l]*(Fh.x[1,0] - Fh.x[] + Fh.y[0,1] - Fh.y[])/(cm[]*Delta);
      dh[] = - dhl + (l > 0 ? dh[] : 0.);
      foreach_dimension()
	dhu.x[] = (Fq.x.x[] + Fq.x.y[] - S.x[1,0] - Fq.x.y[0,1])/(cm[]*Delta);
      
      if (l < nl - 1) {
	scalar div = wl[l];
	div[] = dhl;
      }

      double dmdl = (fm.x[1,0] - fm.x[])/(cm[]*Delta);
      double dmdt = (fm.y[0,1] - fm.y[])/(cm[]*Delta);
      double fG = u.y[]*dmdl - u.x[]*dmdt;
      dhu.x[] += h[]*(G*h[]/2.*dmdl + fG*u.y[]);
      dhu.y[] += h[]*(G*h[]/2.*dmdt - fG*u.x[]);
    }
  }

   if (nl > 1)
    vertical_fluxes ((vector *) &evolving[1], (vector *) &updates[1], wl, dh);
    
  return dtmax;
}

event defaults (i = 0)
{
  assert (ul == NULL && wl == NULL);
  assert (nl > 0);
  ul = vectors_append (ul, u);
  for (int l = 1; l < nl; l++) {
    scalar w = new scalar;
    vector u = new vector;
    foreach_dimension()
      u.x.l = l;
    w.l = l;
    ul = vectors_append (ul, u);
    wl = list_append (wl, w);
  }

  evolving = list_concat ({h}, (scalar *) ul);
  foreach()
    for (scalar s in evolving)
      s[] = 0.;
  multigrid_restriction (evolving);
  
 
  layer = qmalloc (nl, double);
  for (int l = 0; l < nl; l++)
    layer[l] = 1./nl;
  
  advance = advance_saint_venant;
  update = update_saint_venant;
#if TREE
  for (scalar s in {h,zb,u,eta}) {
    s.refine = s.prolongation = refine_linear;
    s.restriction = restriction_volume_average;
  }
  eta.refine  = refine_eta;
  eta.restriction = restriction_eta;
#endif
}

event init (i = 0) {
  foreach()
    eta[] = zb[] + h[];
  multigrid_restriction (all);
}

event cleanup (i = end, last) {
  free (evolving);
  free (layer);
  free (ul), ul = NULL;
  free (wl), wl = NULL;
}

#define radiation(ref) (sqrt (G*max(h[],0.)) - sqrt(G*max((ref) - zb[], 0.)))

#include "elevation.h"
#include "gauges.h"
