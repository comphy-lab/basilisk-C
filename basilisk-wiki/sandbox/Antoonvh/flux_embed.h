/**
# A Cut-cell two-body heat-equation solver with net interfacial flux

Here we are inspired by the problem of an irradiated body that is in
thermal contact with its environment. The solver concerns the
continuous scalar field `s` (representing temperature),

$$\frac{\partial s}{\partial t} = \nabla\cdot \left( \alpha \nabla s\right).$$

whos derivatives are not continuous due to a jump in the (thermal)
diffusivity ($\alpha$) at the interface. Furthermore, we include a
finite (heat) flux that emerges at the embedded surface ($F$). Such
that on the interface we can write,

$$\kappa_1\frac{\partial s_1}{\partial n} - \kappa_2\frac{\partial s_2}{\partial n} = F,$$

where $\kappa_{1,2}$ is the heat conductivity of body 1 and 2,
respectively, $s_{1,2}$ is the scalar field at either side of the
interface and $n$ is the normal vector to the (embedded) interface
pointing towards the body associated with number `2`. The diffusivity
and conductivity are related by the volumetric heat capacity $\rho C_p$
via $\alpha = \frac{\kappa}{\rho C_p}$

## Implementation

The formulation involves an embedded boundary, the (implicit)
diffusion solver and the time loop.
 */
#include "embed.h"
#include "diffusion.h"
#include "run.h"
/**
   We choose to declare fields for the scalar field at either side of
   the interface and a scalar field which stores the local gradient at
   the interface for each field. Furthermore, default values are set
   for the thermal properties.
*/

scalar s1[], s2[], grad1[], grad2[];
double kappav1= 1, kappav2 = 1;
double rcp1 = 1., rcp2 = 1;

s1[embed] = neumann (grad1[]);
s2[embed] = neumann (grad2[]);
/**
On trees we need a modified version of the refinement and prolongation
methods for the scalar field that lives in the region of `cs[] =
0`. For clarity of presentation, they are implemented near the end of
this document. We also define a double-pointer type to workaround the
limitations of qcc, whilst maintaining the (highly desirable)
automatic boundary conditions.
 */
#if TREE
static inline void refine_embed_linear2 (Point point, scalar s);
static inline void restriction_volume_average2 (Point point, scalar s);
#endif

typedef double * doublep;

/**
A helper function is defined that interpolates a scalar value in the
vicinity of the boundary. This is adapted from [`embed.h`](/src/embed.h).
 */
foreach_dimension()
  double interpolate_nearby_x  (Point point, scalar s, scalar cs, face vector fs,
				coord n, coord p, double v[2], double d[2]) {
  foreach_dimension()
    n.x = - n.x;
  for (int i = 0; i < 2; i++) {
    d[i] = nodata;
    v[i] = nodata;
  }
  bool defined = true;
  foreach_dimension()
    if (defined && !fs.x[(n.x > 0.)])
      defined = false;
  if (defined)
    for (int l = 0; l <= 1; l++) {
      int i = (l + 1)*sign(n.x);
      d[l] = (i - p.x)/n.x;
      double y1 = p.y + d[l]*n.y;
      int j = y1 > 0.5 ? 1 : y1 < -0.5 ? -1 : 0;
      y1 -= j;
#if dimension == 2
      if (fs.x[i + (i < 0),j] && fs.y[i,j] && fs.y[i,j+1] &&
	  cs[i,j-1] && cs[i,j] && cs[i,j+1])
	v[l] = quadratic (y1, (s[i,j-1]), (s[i,j]), (s[i,j+1]));
#else // dimension == 3
      double z = p.z + d[l]*n.z;
      int k = z > 0.5 ? 1 : z < -0.5 ? -1 : 0;
      z -= k;
      bool defined = fs.x[i + (i < 0),j,k];
      for (int m = -1; m <= 1 && defined; m++)
	if (!fs.y[i,j,k+m] || !fs.y[i,j+1,k+m] ||
	    !fs.z[i,j+m,k] || !fs.z[i,j+m,k+1] ||
	    !cs[i,j+m,k-1] || !cs[i,j+m,k] || !cs[i,j+m,k+1])
	  defined = false;
      if (defined)
	// bi-quadratic interpolation
	v[l] =
	  quadratic (z,
		     quadratic (y1,
				(s[i,j-1,k-1]), (s[i,j,k-1]), (s[i,j+1,k-1])),
		     quadratic (y1,
				(s[i,j-1,k]),   (s[i,j,k]),   (s[i,j+1,k])),
		     quadratic (y1,
				(s[i,j-1,k+1]), (s[i,j,k+1]), (s[i,j+1,k+1])));
#endif // dimension == 3
      else
	break;
    }
  assert(v[0] != nodata);
  return v[0];
}
/**
The user needs to define the `total_embed_flux` function which returns
the `point`-local flux (F).
*/

double total_flux_embed (Point point, coord p, coord n);

/**
## Flux partitioning

We compute the gradients,

$$\kappa_1\frac{\partial s_1}{\partial n} - \kappa_2\frac{\partial s_2}{\partial n} = F,$$

Following the discretization,

$$-\kappa_1\frac{s_1(\mathbf{p} - \Delta d_0\mathbf{n}) -
s(\mathbf{p})}{\Delta d_0} - \kappa_2\frac{s_2(\mathbf{p}) + \Delta d_0
\mathbf{n} - s(\mathbf{p})}{\Delta d_0} = F,$$

where $\mathbf{p}$ is the center coordiante of the embedded boundary fragment,
$d_0$ is a distance parameter (order 1), and $s(\mathrm{p})$ is the
local surface value of $s$. We can compute the gradients after we
compute the surface value,

$$s(\mathrm{p}) = ...  .$$

The gradients are then stored in their respective fields.  
 */
void set_interface_val (scalar T1, scalar T2,
			scalar grad1, scalar grad2) {
  // The fractions for the `s2` scalar are precomputed
  scalar cs2[];
  face vector fs2[];
#if TREE
  cs2.prolongation = fraction_refine;
  foreach_dimension()
    fs2.x.prolongation = embed_face_fraction_refine_x;
#endif
  foreach() 
    cs2[] = 1. - cs[];
  foreach_face()
    fs2.x[] = 1. - fs.x[];
  foreach() {
    if (cs[] > 0. && cs[] < 1.) {
      doublep v1 = malloc (2*sizeof(double));
      doublep d1 = malloc (2*sizeof(double));
      doublep v2 = malloc (2*sizeof(double));
      doublep d2 = malloc (2*sizeof(double));
      coord p, n;
      embed_geometry (point, &p, &n);
      normalize (&n);
      double FLUXv = total_flux_embed (point, p, n);
#if dimension == 2
      foreach_dimension()
	if (fabs(n.x) >= fabs(n.y)) {
	  interpolate_nearby_x (point, s1, cs, fs, n, p, v1, d1);
	  foreach_dimension()
	    n.x = -n.x;
	  interpolate_nearby_x (point, s2, cs2, fs2, n, p, v2, d2);
	}
#else // dimension == 3
      if (fabs(n.x) >= fabs(n.y)) {
	if (fabs(n.x) >= fabs(n.z)) {
	  interpolate_nearby_x (point, s1, cs, fs, n, p, v1, d1);
	  foreach_dimension()
	    n.x = -n.x;
	  interpolate_nearby_x (point, s2, cs2, fs2, n, p, v2, d2);
	}
      }
      else if (fabs(n.y) >= fabs(n.z)) {
	interpolate_nearby_y (point, s1, cs, fs, n, p, v1, d1);
	foreach_dimension()
	  n.x = -n.x;
	interpolate_nearby_y (point, s2, cs2, fs2, n, p, v2, d2);
      } else {
	interpolate_nearby_z (point, s1, cs, fs, n, p, v1, d1);
	foreach_dimension()
	  n.x = -n.x;
	interpolate_nearby_z (point, s2, cs2, fs2, n, p, v2, d2);
      }
#endif // dimension == 3
      double Ts0 = (v1[0]*(kappav1)/(Delta*d1[0]) + v2[0]*(kappav2)/(Delta*d2[0]) + FLUXv)/
	((kappav1)/(Delta*d1[0]) + (kappav2)/(Delta*d2[0]));
      grad1[]  = -(v1[0] - Ts0)/(d1[0]*Delta);
      grad2[]  = -(v2[0] - Ts0)/(d2[0]*Delta);
      free (d1);free (v1);free (d2);free (v2);
    }
  }
}

event defaults (i++) {
  foreach()
   grad1[] = grad2[] =  nodata;
#if TREE
  s1.refine = s1.prolongation = refine_embed_linear;
  s2.refine = s2.prolongation = refine_embed_linear2;
  s1.restriction = s1.coarsen = restriction_volume_average;
  s2.restriction = s2.coarsen = restriction_volume_average2;
#endif
}

void invert_embed (scalar s, face vector f) {
  foreach()
    s[] = 1. - s[];
  foreach_face()
    f.x[] = 1. - f.x[];
}

/**
## Two-component diffusion solver

The diffusion-problem solver is called twice, one time for each
field. The embedded boundary fractions need to be inverted between the
calls.
 */

event tracer_diffusion (i++, last) {
  set_interface_val (s1, s2, grad1, grad2);
  scalar cp[];
  face vector kappa[];
  foreach()
    cp[] = rcp1*cm[];
  foreach_face()
    kappa.x[] = fs.x[] * kappav1;
  diffusion (s1, dt, kappa, theta = cp);
  invert_embed (cs, fs);
  foreach()
    cp[] = rcp2*cm[];
  foreach_face()
    kappa.x[] = fs.x[] * kappav2;
  diffusion (s2, dt, kappa, theta = cp);
  invert_embed (cs, fs);
}
/**
As promised, there are two attribites to be defined in order to do
consistent coarsening and refinement for the `s2` field. 
 */

static inline void refine_embed_linear2 (Point point, scalar s)
{
  foreach_child() {
    if (!(1. - cs[]))
      s[] = 0.;
    else {
      assert (1 - coarse(cs));
      int i = (child.x + 1)/2, j = (child.y + 1)/2;
#if dimension == 2
      if (1 - coarse(fs.x,i) && 1 - coarse(fs.y,0,j) &&
	  (1 - coarse(cs) == 1. || 1 - coarse(cs,child.x) == 1. ||
	   1 - coarse(cs,0,child.y) == 1. || 1 - coarse(cs,child.x,child.y) == 1.)) {
	assert (1 - coarse(cs,child.x) && 1 - coarse(cs,0,child.y));
	if (1 - coarse(fs.x,i,child.y) && 1 - coarse(fs.y,child.x,j)) {
	  // bilinear interpolation
	  assert (1 - coarse(cs,child.x,child.y));
	  s[] = (9.*coarse(s) + 
		     3.*(coarse(s,child.x) + coarse(s,0,child.y)) + 
		 coarse(s,child.x,child.y))/16.;
	}
	  else
	    // triangular interpolation	  
	    s[] = (2.*coarse(s) + coarse(s,child.x) + coarse(s,0,child.y))/4.;
	}
	else if (1 - coarse(cs,child.x,child.y) &&
	       ((1 - coarse(fs.x,i) && 1 - coarse(fs.y,child.x,j)) ||
		(1 - coarse(fs.y,0,j) && 1 - coarse(fs.x,i,child.y)))) {
	  // diagonal interpolation
	  s[] = (3.*coarse(s) + coarse(s,child.x,child.y))/4.;
      }
#else // dimension == 3
      int k = (child.z + 1)/2;
      if (1 - coarse(fs.x,i) > 0.25 && 1 - coarse(fs.y,0,j) > 0.25 &&
	  1 - coarse(fs.z,0,0,k) > 0.25 &&
	  (1 - coarse(cs) == 1. || 1 - coarse(cs,child.x) == 1. ||
	   1- coarse(cs,0,child.y) == 1. || 1-coarse(cs,child.x,child.y) == 1. ||
	   1-coarse(cs,0,0,child.z) == 1. || 1-coarse(cs,child.x,0,child.z) == 1. ||
	   1-coarse(cs,0,child.y,child.z) == 1. ||
	   1-coarse(cs,child.x,child.y,child.z) == 1.)) {
	assert (1-coarse(cs,child.x) && 1 - coarse(cs,0,child.y) &&
		1 - coarse(cs,0,0,child.z));
	if (1 - coarse(fs.x,i,child.y) && 1 - coarse(fs.y,child.x,j) &&
	    1 - coarse(fs.z,child.x,child.y,k) &&
	    1 - coarse(fs.z,child.x,0,k) && 1 - coarse(fs.z,0,child.y,k)) {
	  assert (1 - coarse(cs,child.x,child.y) && 1 - coarse(cs,child.x,0,child.z) &&
		  1 - coarse(cs,0,child.y,child.z) &&
		  1 - coarse(cs,child.x,child.y,child.z));
	  // bilinear interpolation
	  s[] = (27.*coarse(s) + 
		 9.*(coarse(s,child.x) + coarse(s,0,child.y) +
		     coarse(s,0,0,child.z)) + 
		 3.*(coarse(s,child.x,child.y) + coarse(s,child.x,0,child.z) +
		     coarse(s,0,child.y,child.z)) + 
		 coarse(s,child.x,child.y,child.z))/64.;
	}
	else
	  // tetrahedral interpolation
	  s[] = (coarse(s) + coarse(s,child.x) + coarse(s,0,child.y) +
		 coarse(s,0,0,child.z))/4.;
      }
      else if (1 - coarse(cs,child.x,child.y,child.z) &&
	       ((1 - coarse(fs.z,child.x,child.y,k) &&
		 ((1 - coarse(fs.x,i) && 1 - coarse(fs.y,child.x,j)) ||
		  (1 - coarse(fs.y,0,j) && 1 - coarse(fs.x,i,child.y))))
		||
		(1 - coarse(fs.z,0,0,k) &&
		 ((1 - coarse(fs.x,i,0,child.z) && 1- coarse(fs.y,child.x,j,child.z)) ||
		  (1 - coarse(fs.y,0,j,child.z) && 1 - coarse(fs.x,i,child.y,child.z))))
		||
		(1 - coarse(fs.z,child.x,0,k) &&
		 1 - coarse(fs.x,i) && 1- coarse(fs.y,child.x,j,child.z))
		||
		(1 -coarse(fs.z,0,child.y,k) &&
		 1 - coarse(fs.y,0,j) && 1 - coarse(fs.x,i,child.y,child.z))
		))
	// diagonal interpolation
	s[] = (3.*coarse(s) + coarse(s,child.x,child.y,child.z))/4.;
#endif // dimension == 3
      else {
	// Pathological cases, use 1D gradients.
	s[] = coarse(s);
	foreach_dimension() {
	  if (1 - coarse(fs.x,(child.x + 1)/2) && 1 - coarse(cs,child.x))
	    s[] += (coarse(s,child.x) - coarse(s))/4.;
	  else if (coarse(1 - fs.x,(- child.x + 1)/2) && 1 - coarse(cs,- child.x))
	    s[] -= (coarse(s,- child.x) - coarse(s))/4.;
	}
      }
    }
  }
}

static inline void restriction_volume_average2 (Point point, scalar s)
{
  double sum = 0.;
  foreach_child()
    sum += (1 - cm[])*s[];
  s[] = sum/(1 << dimension)/((1 - cm[]) + 1e-30);
}

/**
## Tests

* Diffusion of a Gaussian pulse for a one-solid problem
* *The* linear stationary solution
* Equilibration of a sphere with its environment
* Transient heating of a sphere

## Usage

* ...

## Todo

* Tests
* Usage
* Optimize by reducing the number of fields
* Think about the sense of the implicit/explicit mismatch
* Get a grip on the stability properties

*/
