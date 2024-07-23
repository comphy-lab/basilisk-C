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
The default uniform layer distribution can be replaced with a
geometric progression for the layer thicknesses. This needs to be
called for example in the `init()` event. The `rmin` parameter
specifies the minimum layer thickness relative to the uniform layer
thickness (proportional to `1/nl`). If the `top` parameter is set to
`true` the minimum layer thickness is at the top (layer `nl - 1`),
otherwise it is at the bottom (layer 0). */

void geometric_beta (double rmin, bool top)
{
  if (rmin <= 0. || rmin >= 1. || nl < 2)
    return;
  double r = 1. + 2.*(1./rmin - 1.)/(nl - 1.);
  double hmin = (r - 1.)/(pow(r, nl) - 1.);
  for (int l = 0; l < nl; l++)
    beta[l] = hmin*pow(r, top ? nl - 1 - l : l);
}

/**
The *vertical_remapping()* function takes a (block) field of layer
thicknesses and the corresponding list of tracer fields and performs
the remapping (defined by *beta*). */

trace
void vertical_remapping (scalar h, scalar * tracers)
{
  int nvar = list_len(tracers), ndof = 1, npos = nl + 1;
  foreach() {
#if HALF
    double H0 = 0., H1 = 0., H;
    foreach_layer() {
      if (point.l < nl/2)
	H0 += h[];
      else
	H1 += h[];
    }
    H = H0 + H1;
#else
    double H = 0.;
    foreach_layer()
      H += h[];
#endif
    
    if (H > dry) {
      double zpos[npos], znew[npos];
      double fdat[nvar*nl], fnew[nvar*nl];
      zpos[0] = znew[0] = 0.;
      foreach_layer() {
	zpos[point.l+1] = zpos[point.l] + max(h[],dry);
	int i = nvar*point.l;
	for (scalar s in tracers) {
	  dimensional (fnew[i] = s[]);
	  fdat[i++] = s[];
	}
#if HALF
	if (point.l < nl/2)
	  h[] = 2.*H0*beta[point.l];
	else
	  h[] = 2.*H1*beta[point.l];	
#else
	h[] = H*beta[point.l];
#endif
	znew[point.l+1] = znew[point.l] + h[];
      }

      my_remap (&npos, &npos, &nvar, &ndof, zpos, znew, fdat, fnew,
		&edge_meth, &cell_meth, &cell_lim);

      foreach_layer() {
	int i = nvar*point.l;
	for (scalar s in tracers)
	  s[] = fnew[i++];
      }
    }
  }
}

/**
The remapping is applied at every timestep. */

event remap (i++) {
  if (nl > 1)
    vertical_remapping (h, tracers);
}

/**
The *beta* array is freed at the end of the run. */

event cleanup (i = end)
{
  free (beta), beta = NULL;
}
