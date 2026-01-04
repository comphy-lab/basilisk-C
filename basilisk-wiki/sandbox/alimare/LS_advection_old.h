/**
# Advection of a level-set function 

We will use a RK3 time integration scheme for collocated variables.
With its own timestep (no restriction due to the presence of an embedded
boundary).
*/
#include "simple_discretization.h"
/**
For this function we need the following files.
*/
#include "LS_reinit.h"
#include "../ghigo/src/myquadratic.h" // best extrapolation so far

/**
If we solve the Navier-Stokes equations, we use ghigo's method to update the
variables in emerged cells.
*/
#if NS_emerged
#include "emerged_NS.h"
#endif


struct LSadv{
  scalar dist;
  scalar cs;
  face vector fs;
  scalar TS;
  scalar TL;
  vector vpcf;
  int itredist;
  double s_clean;
  double NB_width;
  scalar curve;
};

void advection_LS(struct LSadv p){

  scalar dist         = p.dist;
  scalar cs           = p.cs;
  face vector fs      = p.fs;
  scalar TS           = p.TS;
  scalar TL           = p.TL;
  vector vpcf         = p.vpcf;
  int    itredist     = p.itredist;
  double s_clean      = p.s_clean;
  double NB_width     = p.NB_width;
  scalar curve        = p.curve;

/**
Previous state of cs is saved into csm1
*/
  foreach(){
    csm1[]   = cs[];
  }

  boundary({csm1});
  restriction({csm1});
  face vector fsm1[];
  foreach_face()
    fsm1.x[] = fs.x[];
  boundary((scalar *){fsm1});
  restriction((scalar *){fsm1});

  RK3_WENO5(dist,vpcf,dt, NB_width);

  boundary ({dist});
  restriction({dist});

  LS_reinit(dist, it_max = itredist);

/**
We remove the overshoots that we might create with reinitialization.
*/
  foreach(){
    dist[] = clamp(dist[], -1.05*NB_width, 1.05*NB_width);
  }
  boundary({dist});
  restriction({dist});
/**
We update the volume and face fractions accordingly
*/
  LS2fractions(dist,cs,fs,s_clean);

#if CURVE_LS // mandatory in Gibbs Thomsmon cases
  curvature_LS(dist,curve);
#else
  curvature(cs,curve);
#endif
  boundary({curve});
  restriction({curve});

/**
Sometimes, a new cell becomes uncovered/covered because the interface has moved.
We must initialize the tracers field. We use Ghigo's methodology in
update_tracer() that can be found in this [page](../ghigo/myembed.h). The
basic idea is very similar to what can be found in the Dirichlet boundary
condition calculation in `embed.h`.
*/
#if 1
  foreach() {
    if (cs[] > 0. && csm1[] <= 0.) { // Emerged cells
      assert(cs[]!=1.); // cell shouldn't be full
      coord o = {0.,0.};
      TL[] = embed_extrapolate_ls (point, TL, cs, o, false);
    }
  }
  
  /**
  We update the boundary condition for *a*. */
  boundary({TL});
  restriction({TL});
  
#if NS_emerged
  init_emerged_NS(cs,csm1);
#endif
  invertcs(cs,fs);
  invertcs(csm1,fsm1);

  foreach() {
    if (cs[] > 0. && csm1[] <= 0.) { // Emerged cells
      assert(cs[]!=1.); // cell shouldn't be full
      coord o = {0.,0.};
      TS[] = embed_extrapolate_ls(point, TS,cs, o, false);
    }
  }
  boundary ({TS});
  invertcs(cs,fs);
  invertcs(csm1,fsm1);
#endif
}
