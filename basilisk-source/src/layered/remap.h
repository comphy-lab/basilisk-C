/**
# Vertical remapping

This implements a simple vertical remapping to "$\sigma$-coordinates"
(equally-distributed by default).

We can optionally use the [PPR
Library](https://github.com/dengwirda/PPR) of Engwirda and Kelley to
perform the remapping. The default settings are using the Parabolic
Piecewise Method without limiting. 

By default we use the [C PPR implementation](remapc.h) without
limiting. */

#if USE_PPR
# include "ppr/ppr.h"
// int edge_meth = p1e_method, cell_meth = plm_method, cell_lim = null_limit;
int edge_meth = p3e_method, cell_meth = ppm_method;
// int edge_meth = p5e_method, cell_meth = pqm_method, cell_lim = null_limit;
#else // !USE_PPR
# include "remapc.h"
const bool null_limit = false, mono_limit = true;
#endif
int cell_lim = null_limit;

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
  const int npos = nl + 1;
#if USE_PPR
  int nvar = list_len(tracers), ndof = 1;
#endif
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
#if USE_PPR
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
#else // !USE_PPR
      zpos[0] = znew[0] = 0.;
      foreach_layer() {
	zpos[point.l+1] = zpos[point.l] + max(h[],dry);
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
      znew[npos -1] = zpos[npos - 1];
      double fdat[nl], fnew[nl];
      for (scalar s in tracers) {
        int j = 0;
        foreach_layer() {
          dimensional (fnew[j] = s[]);
          fdat[j++] = s[];
        }

        /**
        For the moment we simply use Neumann zero top and bottom
        boundary conditions for all fields. */

        remap_c (npos, npos, zpos, znew, fdat, fnew,
                 0, HUGE, 0,
                 0, HUGE, 0,
                 cell_lim);

        j = 0;
        foreach_layer()
	  s[] = fnew[j++];
      }
#endif // !USE_PPR
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
