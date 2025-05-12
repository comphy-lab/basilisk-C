/**
# Vertical remapping

This implements a simple vertical remapping to "$\sigma$-coordinates"
(equally-distributed by default).

We use the [PPR Library](https://github.com/dengwirda/PPR) of Engwirda
and Kelley to perform the remapping. The default settings are using
the Parabolic Piecewise Method without limiting. */

#include "ppr/ppr.h"

// int edge_meth = p1e_method, cell_meth = plm_method, cell_lim = null_limit;
int edge_meth = p3e_method, cell_meth = ppm_method, cell_lim = null_limit;
// int edge_meth = p5e_method, cell_meth = pqm_method, cell_lim = null_limit;

/**
The distribution of layers can be controlled using the *beta* array
which defines the ratio of the thickness of each layer to the total
depth $H$ (i.e. the relative thickness). By default all layers have
the same relative thickness. */

double * beta = NULL;

event defaults (i = 0)
{
  beta = malloc (nl*sizeof(double));
  for (int l = 0; l < nl; l++)
    beta[l] = 1./nl;
}

/**
The *vertical_remapping()* function takes a list of layer thicknesses
and the corresponding array of lists of tracers for each layer and
performs the remapping (defined by *beta*). */

trace
void vertical_remapping (scalar * hl, scalar ** tracers)
{
  int nvar = list_len(tracers[0]), ndof = 1, npos = nl + 1;
  foreach() {
    double H = 0.;
    for (scalar h in hl)
      H += h[];

    if (H > dry) {
      double zpos[npos], znew[npos];
      double fdat[nvar*nl], fnew[nvar*nl];
      zpos[0] = znew[0] = 0.;
      int l = 0;
      for (scalar h in hl) {
	zpos[l+1] = zpos[l] + h[];
	int i = nvar*l;
	for (scalar s in tracers[l])
	  fdat[i++] = s[];
	h[] = H*beta[l];
	znew[l+1] = znew[l] + h[];
	l++;
      }

      my_remap (&npos, &npos, &nvar, &ndof, zpos, znew, fdat, fnew,
		&edge_meth, &cell_meth, &cell_lim);

      l = 0;
      for (scalar h in hl) {
	int i = nvar*l;
	for (scalar s in tracers[l])
	  s[] = fnew[i++];
	l++;
      }
    }
  }
  
  scalar * list = list_copy (hl);
  for (int l = 0; l < nl; l++)
    for (scalar s in tracers[l])
      list = list_append (list, s);
  boundary (list);
  free (list);
}

/**
The remapping is applied at every timestep just before the vertical
viscosity i.e. just after horizontal advection in the [multilayer
solver](hydro.h). */

event viscous_term (i++) {
  if (nl > 1)
    vertical_remapping (hl, tracers);
}

/**
The *beta* array is freed at the end of the run. */

event cleanup (i = end)
{
  free (beta), beta = NULL;
}
