#include "fractions_ebm.h"

attribute {
  scalar angle;
  scalar cfg;
  scalar cet;
  vector mcl;
  vector hp;
}

event stability (i++) {
#if !EMBED
  if (CFL > 0.5)
    CFL = 0.5;
#else
  if (CFL > 0.5)
    CFL = 0.1;
#endif
}

static scalar * _interface = NULL;

// To avoid the advection of "f"
event vof (i++) {
  interfaces = _interface;
}

// need to disable the vof event
#include "vof.h"

// To avoid the advection of "f"
event vof (i++) {
  _interface = interfaces;
  interfaces = NULL;
}

event defaults (i = 0)
{
  for (scalar c in interfaces) {
    scalar cfg = c.cfg;
    if (!cfg.i)
      cfg = new scalar;
    c.cfg = cfg;

    scalar cet = c.cet;
    if (!cet.i)
      cet = new scalar;
    c.cet = cet;

    vector mcl = c.mcl;
    if (!mcl.x.i)
      mcl = new vector;
    c.mcl = mcl;

    vector hp = c.hp;
    if (!hp.x.i)
      hp = new vector;
    c.hp = hp;
  }
}

foreach_dimension()
static void sweep_ebm_x (scalar c, scalar cc, scalar * tcl)
{
  vector n[], ncs[];
  scalar alpha[], alphacs[], flux[];
  double cfl = 0.;

  scalar * tracers = c.tracers, * gfl = NULL, * tfluxl = NULL;
  if (tracers) {
    for (scalar t in tracers) {
      scalar gf = new scalar, flux = new scalar;
      gfl = list_append (gfl, gf);
      tfluxl = list_append (tfluxl, flux);
    }

    foreach() {
      scalar t, gf;
      for (t,gf in tracers,gfl)
	      gf[] = vof_concentration_gradient_x (point, c, t);
    }
  }

  reconstruction_cs (cs, fs, ncs, alphacs);
  reconstruction_emd (c, n, alpha, cs, ncs, alphacs);
  
  foreach_face(x, reduction (max:cfl)) {

    double un = uf.x[]*dt/(Delta*fm.x[] + SEPS), s = sign(un);
    int i = -(s + 1.)/2.;

    if (un*fm.x[]*s/(cm[] + SEPS) > cfl)
      cfl = un*fm.x[]*s/(cm[] + SEPS);

    double cf;
    if (cs[i] >= 1.)
      cf = (c[i] <= 0. || c[i] >= 1.) ? c[i] :
      rectangle_fraction ((coord){-s*n.x[i], n.y[i], n.z[i]}, alpha[i],
			  (coord){-0.5, -0.5, -0.5},
			  (coord){s*un - 0.5, 0.5, 0.5});
    else {
      if (c[i] <= 0.)
        cf = 0.;
      else if (c[i] >= cs[i])
        cf = 1.;
      else {
        assert (cs[i] > 0. && cs[i] < 1.);
        assert (c[i] > 0. && c[i] < cs[i]);
#if dimension == 2
        if (fs.x[] > 0. && fabs(un) > 0.) {
          coord m   = {1., 0., 0.};
          coord mcs = {-s*ncs.x[i], ncs.y[i], ncs.z[i]};
          double fv = uf.x[]*dt*s/Delta;
          assert (fv > 0. && fv < cs[i]);
          double unc = line_alpha_ebm (fv, cs[i], mcs, alphacs[i], m) + 0.5;
          assert (unc > fv - 1e-10);

          double vlq = rectangle_fraction_cs (cs[i], fs.x[], (coord){-s*ncs.x[i], ncs.y[i], ncs.z[i]}, alphacs[i],
                (coord){-s*n.x[i], n.y[i], n.z[i]}, alpha[i],
			          (coord){-0.5, -0.5, -0.5},
			          (coord){unc - 0.5, 0.5, 0.5});
          double vfluid = rectangle_fraction ((coord){-s*ncs.x[i], ncs.y[i], ncs.z[i]}, alphacs[i],
			          (coord){-0.5, -0.5, -0.5},
			          (coord){unc - 0.5, 0.5, 0.5});

          if (vfluid < vlq - 1e-10)
            fprintf (stdout, "Warning: t: %g vfluid - vlq = %g\n", t, (vfluid - vlq));
          cf = min (vlq/vfluid, 1.);
        }
        else
          cf = 0.;
#endif
      }
    }

    flux[] = cf*uf.x[];
    
    scalar t, gf, tflux;
    for (t,gf,tflux in tracers,gfl,tfluxl) {
      double cf1 = cf, ci = c[i];
      if (t.inverse)
	      cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
	      double ff = t[i]/ci + s*min(1., 1. - s*un)*gf[i]*Delta/2.;
	      tflux[] = ff*cf1*uf.x[];
      }
      else
	      tflux[] = 0.;
    }
  }
  delete (gfl); free (gfl);

  if (cfl > 0.1 + 1e-6)
    fprintf (stdout, 
	     "WARNING: CFL must be <= 0.1 for VOF (cfl - 0.1 = %g)\n", 
	     cfl - 0.1), fflush (stdout);

  foreach()
    if (cs[] > 0.) {
      c[] += dt*(flux[] - flux[1] + cc[]*(uf.x[1] - uf.x[]))/Delta;
      c[] = c[] < 1e-10 ? 0. : c[] > cs[] - 1e-10 ? cs[] : c[]; // filter
#if NO_1D_COMPRESSION
      for (t, tflux in tracers, tfluxl)
	t[] += dt*(tflux[] - tflux[1])/Delta;
#else // !NO_1D_COMPRESSION
      scalar t, tc, tflux;
      for (t, tc, tflux in tracers, tcl, tfluxl)
	t[] += dt*(tflux[] - tflux[1] + tc[]*(uf.x[1] - uf.x[]))/Delta;
#endif // !NO_1D_COMPRESSION
    }

  scalar cfg = c.cfg;
  scalar cfgn[];
  foreach()
    cfgn[] = cfg[];

  foreach()
    if (cfg[] == 1.) {
      if (is_three_phase (c[], cs[])) {
        coord m, mcs;
        foreach_dimension() {
          m.x = n.x[];
          mcs.x = ncs.x[];
        }
        double alphan = line_alpha_ebm (c[], cs[], mcs, alphacs[], m);
        coord clp = find_cl_pos (m, alphan, mcs, alphacs[]);
        if (fabs(clp.x) > 0.5 || fabs(clp.y) > 0.5)
          for (int k = 0; k < 2; k++) {
            int i = 2*k - 1;
            foreach_dimension()
              if (fs.x[k] > 0. && fs.x[k] < 1. && i*clp.x > 0.5 && is_three_phase (c[i], cs[i]))
                cfgn[] = (fabs(cfg[i]) == 2. ? - cfg[i] : nodata);
          }
      }
      else
        cfgn[] = (c[] >= cs[] ? 2. : -2.);
    }

  foreach()
    if (cs[] > 0. && cs[] < 1. && cfgn[] == nodata)
      for (int k = 0; k < 2; k++) {
        int i = 2*k - 1;
        foreach_dimension()
          if (fs.x[k] > 0. && fs.x[k] < 1. && fabs(cfg[i]) == 2.)
            cfgn[] = cfg[i];
      }

  foreach()
    if (cs[] > 0. && cs[] < 1. && fabs(cfgn[]) == 2.)
      for (int k = 0; k < 2; k++) {
        int i = 2*k - 1;
        foreach_dimension()
          if (fs.x[k] > 0. && fs.x[k] < 1. && (cfgn[] + cfgn[i]) == 0. && cfg[i] == 1.)
            cfgn[] = 1.;
      }

  foreach()
    cfg[] = cfgn[];

  delete (tfluxl); free (tfluxl);
}

/**
## Multi-dimensional advection

The multi-dimensional advection is performed by the event below. */

void vof_advection_ebm (scalar * interfaces, int i)
{
  for (scalar c in interfaces) {

    /**
    We first define the volume fraction field used to compute the
    divergent term in the one-dimensional advection equation above. We
    follow [Weymouth & Yue, 2010](/src/references.bib#weymouth2010) and use a
    step function which guarantees exact mass conservation for the
    multi-dimensional advection scheme (provided the advection velocity
    field is exactly non-divergent). */

    scalar cc[], * tcl = NULL, * tracers = c.tracers;    
    for (scalar t in tracers) {
#if !NO_1D_COMPRESSION
      scalar tc = new scalar;
      tcl = list_append (tcl, tc);
#endif // !NO_1D_COMPRESSION
#if TREE
      if (t.refine != vof_concentration_refine) {
	t.refine = t.prolongation = vof_concentration_refine;
	t.restriction = restriction_volume_average;
	t.dirty = true;
	t.c = c;
      }
#endif // TREE
    }
    foreach() {
#if !EMBED
      cc[] = (c[] > 0.5);
#else
      cc[] = (c[] > 0.5*cs[]);
#endif
#if !NO_1D_COMPRESSION
      scalar t, tc;
      for (t, tc in tracers, tcl) {
	if (t.inverse)
	  tc[] = c[] < 0.5 ? t[]/(1. - c[]) : 0.;
	else
	  tc[] = c[] > 0.5 ? t[]/c[] : 0.;
      }
#endif // !NO_1D_COMPRESSION
    }

    /**
    We then apply the one-dimensional advection scheme along each
    dimension. To try to minimise phase errors, we alternate dimensions
    according to the parity of the iteration index `i`. */

    void (* sweep[dimension]) (scalar, scalar, scalar *);
    int d = 0;
#if !EMBED
    foreach_dimension()
      sweep[d++] = sweep_x;
#else
    foreach_dimension()
      sweep[d++] = sweep_ebm_x;
#endif
    for (d = 0; d < dimension; d++) {
      sweep[(i + d) % dimension] (c, cc, tcl);
#if EMBED
      update_cfg (c, cs, fs, c.cfg);
      update_cet (c, cs, fs, c.cfg, c.cet, c.angle);
#endif

      assert (c.height.x.i);
#if !EMBED
      heights (c, c.height);
#else
      heights_ebm (c, cs, c.height, c.cet);
#endif
    }
    delete (tcl), free (tcl);
  }
}

event vof (i++)
  vof_advection_ebm (interfaces, i);
