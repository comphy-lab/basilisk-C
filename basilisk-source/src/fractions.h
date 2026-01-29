/**
# Volume fractions

These functions are used to maintain or define volume and surface
fractions either from an initial geometric definition or from an
existing volume fraction field. 

We will use basic geometric functions for square cut cells and the
"Mixed-Youngs-Centered" (MYC) normal approximation of Ruben
Scardovelli. */

#include "geometry.h"
#if dimension == 1
coord mycs (Point point, scalar c) {
  return (coord){sign(c[-1] - c[1])};
}
#elif dimension == 2
# include "myc2d.h"
#else // dimension == 3
# include "myc.h"
#endif

/**
By default the interface normal is computed using the MYC
approximation. This can be overloaded by redefining this macro. */

auto macro coord interface_normal (Point p, scalar c) {
  return mycs (p, c);
}

/**
## Coarsening and refinement of a volume fraction field 

On trees, we need to define how to coarsen (i.e. "restrict") or
refine (i.e. "prolongate") interface definitions (see [geometry.h]()
for a basic explanation of how interfaces are defined). */

#if TREE

void fraction_refine (Point point, scalar c)
{
  
  /**
  If the parent cell is empty or full, we just use the same value for
  the fine cell. */

  double cc = c[];
  if (cc <= 0. || cc >= 1.)
    foreach_child()
      c[] = cc;
  else {

    /**
    Otherwise, we reconstruct the interface in the parent cell. */

    coord n = mycs (point, c);
    double alpha = plane_alpha (cc, n);

    /**
    And compute the volume fraction in the quadrant of the coarse cell
    matching the fine cells. We use symmetries to simplify the
    combinations. */

    foreach_child() {
      static const coord a = {0.,0.,0.}, b = {.5,.5,.5};
      coord nc;
      foreach_dimension()
	nc.x = child.x*n.x;
      c[] = rectangle_fraction (nc, alpha, a, b);
    }
  }
}

/**
Finally, we also need to prolongate the reconstructed value of
$\alpha$. This is done with the simple formula below. We add an
attribute so that we can access the normal from the refinement
function. */

attribute {
  vector n;
}

static void alpha_refine (Point point, scalar alpha)
{
  vector n = alpha.n;
  double alphac = 2.*alpha[];
  coord m;
  foreach_dimension()
    m.x = n.x[];
  foreach_child() {
    alpha[] = alphac;
    foreach_dimension()
      alpha[] -= child.x*m.x/2.;
  }
}

#endif // TREE

/**
## Computing volume fractions from a "levelset" function 

Initialising a volume fraction field representing an interface is not
trivial since it involves the numerical evaluation of surface
integrals.

Here we define a function which allows the approximation of these
surface integrals in the case of an interface defined by a "levelset"
function $\Phi$ sampled on the *vertices* of the grid.

By convention the "inside" of the interface corresponds to $\Phi > 0$.

The function takes the vertex scalar field $\Phi$ as input and fills
`c` with the volume fraction and, optionally if it is given, `s`
with the surface fractions i.e. the fractions of the faces of the cell
which are inside the interface.

![Volume and surface fractions](/src/figures/fractions.svg) */

trace
void fractions (vertex scalar Phi, scalar c,
		face vector s = {0}, double val = 0.)
{
#if dimension > 1
  face vector as = automatic (s);
  
  /**
  We store the positions of the intersections of the surface with the
  edges of the cell in vector field `p`. In two dimensions, this field
  is just the transpose of the *line fractions* `s`, in 3D we need to
  allocate a new field. */
  
#if dimension == 3
  vector p[];
#else // dimension == 2
  vector p;
  p.x = as.y; p.y = as.x;
#endif
  
  /**
  ### Line fraction computation
  
  We start by computing the *line fractions* i.e. the (normalised)
  lengths of the edges of the cell within the surface. */

  foreach_edge() {

    /**
    If the values of $\Phi$ on the vertices of the edge have opposite
    signs, we know that the edge is cut by the interface. */

    if ((Phi[] - val)*(Phi[1] - val) < 0.) {

      /**
      In that case we can find an approximation of the interface position by
      simple linear interpolation. We also check the sign of one of the
      vertices to orient the interface properly. */

      p.x[] = (Phi[] - val)/(Phi[] - Phi[1]);
      if (Phi[] < val)
	p.x[] = 1. - p.x[];
    }

    /**
    If the values of $\Phi$ on the vertices of the edge have the same sign
    (or are zero), then the edge is either entirely outside or entirely
    inside the interface. We check the sign of both vertices to treat
    limit cases properly (when the interface intersects the edge exactly
    on one of the vertices). */

    else
      p.x[] = (Phi[] > val || Phi[1] > val);
  }

  /**
  ### Surface fraction computation 

  We can now compute the surface fractions. In 3D they will be
  computed for each face (in the z, x and y directions) and stored in
  the face field `s`. In 2D the surface fraction in the z-direction is
  the *volume fraction* `c`. */

#if dimension == 3

  /**
  In 3D we need to prevent boundary conditions, since this would
  impose vertex field BCs which are not (apparently) consistent for
  the edge intersection coordinates. This can probably be improved. */
  
  foreach_dimension()
    p.x.dirty = false;
  
  scalar s_x = as.x, s_y = as.y, s_z = as.z;
  foreach_face(z,x,y)
#else // dimension == 2
  scalar s_z = c;
  foreach()
#endif
  {

    /**
    We first compute the normal to the interface. This can be done easily
    using the line fractions. The idea is to compute the circulation of
    the normal along the boundary $\partial\Omega$ of the fraction of the
    cell $\Omega$ inside the interface. Since this is a closed curve, we
    have
    $$
    \oint_{\partial\Omega}\mathbf{n}\;dl = 0
    $$ 
    We can further decompose the integral into its parts along the edges
    of the square and the part along the interface. For the case pictured
    above, we get for one component (and similarly for the other)
    $$
    - s_x[] + \oint_{\Phi=0}n_x\;dl = 0
    $$
    If we now define the *average normal* to the interface as
    $$
    \overline{\mathbf{n}} = \oint_{\Phi=0}\mathbf{n}\;dl
    $$
    We have in the general case
    $$
    \overline{\mathbf{n}}_x = s_x[] - s_x[1,0]
    $$
    and
    $$
    |\overline{\mathbf{n}}| = \oint_{\Phi=0}\;dl
    $$ 
    Note also that this average normal is exact in the case of a linear
    interface. */

    coord n;
    double nn = 0.;
    foreach_dimension(2) {
      n.x = p.y[] - p.y[1];
      nn += fabs(n.x);
    }
    
    /**
    If the norm is zero, the cell is full or empty and the surface fraction
    is identical to one of the line fractions. */

    if (nn == 0.)
      s_z[] = p.x[];
    else {
    
      /**
      Otherwise we are in a cell containing the interface. We first
      normalise the normal. */

      foreach_dimension(2)
	n.x /= nn;

      /**
      To find the intercept $\alpha$, we look for edges which are cut by the
      interface, find the coordinate $a$ of the intersection and use it to
      derive $\alpha$. We take the average of $\alpha$ for all intersections. */
      
      double alpha = 0., ni = 0.;
      for (int i = 0; i <= 1; i++)
	foreach_dimension(2)
	  if (p.x[0,i] > 0. && p.x[0,i] < 1.) {
	    double a = sign(Phi[0,i] - val)*(p.x[0,i] - 0.5);
	    alpha += n.x*a + n.y*(i - 0.5);
	    ni++;
	  }

      /**
      Once we have $\mathbf{n}$ and $\alpha$, the (linear) interface
      is fully defined and we can compute the surface fraction using
      our pre-defined function. For marginal cases, the cell is full
      or empty (*ni == 0*) and we look at the line fractions to
      decide. */

      if (ni == 0)
	s_z[] = max (p.x[], p.y[]);
      else if (ni != 4)
	s_z[] = line_area (n.x, n.y, alpha/ni);
      else {
#if dimension == 3
	s_z[] = (p.x[] + p.x[0,1] + p.y[] + p.y[1] > 2.);
#else
	s_z[] = 0.;
#endif
      }
    }
  }
  
  /**
  ### Volume fraction computation

  To compute the volume fraction in 3D, we use the same approach. */
  
#if dimension == 3
  foreach() {

    /**
    Estimation of the average normal from the surface fractions. */
       
    coord n;
    double nn = 0.;
    foreach_dimension(3) {
      n.x = as.x[] - as.x[1];
      nn += fabs(n.x);
    }
    if (nn == 0.)
      c[] = as.x[];
    else {
      foreach_dimension(3)
	n.x /= nn;

      /**
      We compute the average value of *alpha* by looking at the
      intersections of the surface with the twelve edges of the
      cube. */
      
      double alpha = 0., ni = 0.;
      for (int i = 0; i <= 1; i++)
	for (int j = 0; j <= 1; j++)
	  foreach_dimension(3)
	    if (p.x[0,i,j] > 0. && p.x[0,i,j] < 1.) {
	      double a = sign(Phi[0,i,j] - val)*(p.x[0,i,j] - 0.5);
	      alpha += n.x*a + n.y*(i - 0.5) + n.z*(j - 0.5);
	      ni++;
	    }

      /**
      Finally we compute the volume fraction. */

      if (ni == 0)
	c[] = as.x[];
      else if (ni < 3 || ni > 6)
	c[] = 0.; // this is important for robustness of embedded boundaries
      else
	c[] = plane_volume (n, alpha/ni);
    }
  }
#endif // dimension == 3
#else  // dimension == 1
  if (s.x.i)
    foreach_face()
      s.x[] = Phi[] > 0.;
  foreach()
    if ((Phi[] - val)*(Phi[1] - val) < 0.) {
      c[] = (Phi[] - val)/(Phi[] - Phi[1]);
      if (Phi[] < val)
	c[] = 1. - c[];
    }
    else
      c[] = (Phi[] > val || Phi[1] > val);
#endif
}

/**
The convenience macros below can be used to define volume and surface
fraction fields directly from a function. */

macro fraction (scalar f, double func)
{
  {
    vertex scalar phi[];
    foreach_vertex()
      phi[] = func;
    fractions (phi, f);
  }
}

macro solid (scalar cs, face vector fs, double func)
{
  {
    vertex scalar phi[];
    foreach_vertex()
      phi[] = func;
    fractions (phi, cs, fs);
  }
}

/**
### Boolean operations

Implicit surface representations have the advantage of allowing simple
[constructive solid
geometry](https://en.wikipedia.org/wiki/Constructive_solid_geometry)
operations. */

#define intersection(a,b)   min(a,b)
#define union(a,b)          max(a,b)
#define difference(a,b)     min(a,-(b))

/**
## Interface reconstruction from volume fractions

The reconstruction of the interface geometry from the volume fraction
field requires computing an approximation to the interface normal.

### Youngs normal approximation 

This a simple, but relatively inaccurate way of approximating the
normal. It is simply a weighted average of centered volume fraction
gradients. We include it as an example but it is not used. */

coord youngs_normal (Point point, scalar c)
{
  coord n;
  double nn = 0.;
  assert (dimension == 2);
  foreach_dimension() {
    n.x = (c[-1,1] + 2.*c[-1,0] + c[-1,-1] -
	   c[+1,1] - 2.*c[+1,0] - c[+1,-1]);
    nn += fabs(n.x);
  }
  // normalize
  if (nn > 0.)
    foreach_dimension()
      n.x /= nn;
  else // this is a small fragment
    n.x = 1.;
  return n;
}

/**
### Normal approximation using MYC or face fractions
*/

coord facet_normal (Point point, scalar c, face vector s)
{
  if (s.x.i >= 0) { // compute normal from face fractions
    coord n;
    double nn = 0.;
    foreach_dimension() {
      n.x = s.x[] - s.x[1];
      nn += fabs(n.x);
    }
    if (nn > 0.)
      foreach_dimension()
	n.x /= nn;
    else
      foreach_dimension()
	n.x = 1./dimension;
    return n;
  }
  return interface_normal (point, c);
}

/**
### Interface reconstruction 

The reconstruction function takes a volume fraction field `c` and
returns the corresponding normal vector field `n` and intercept field
$\alpha$. */

trace
void reconstruction (const scalar c, vector n, scalar alpha)
{
  foreach() {

    /**
    If the cell is empty or full, we set $\mathbf{n}$ and $\alpha$ only to
    avoid using uninitialised values in `alpha_refine()`. */

    if (c[] <= 0. || c[] >= 1.) {
      alpha[] = 0.;
      foreach_dimension()
	n.x[] = 0.;
    }
    else {

      /**
      Otherwise, we compute the interface normal using the
      Mixed-Youngs-Centered scheme, copy the result into the normal field
      and compute the intercept $\alpha$ using our predefined function. */

      coord m = interface_normal (point, c);
      foreach_dimension()
	n.x[] = m.x;
      alpha[] = plane_alpha (c[], m);
    }
  }

#if TREE

  /**
  On a tree grid, for the normal to the interface, we don't use
  any interpolation from coarse to fine i.e. we use straight
  "injection". */

  foreach_dimension()
    n.x.refine = n.x.prolongation = refine_injection;

  /**
  We set our refinement function for *alpha*. */

  alpha.n = n;
  alpha.refine = alpha.prolongation = alpha_refine;
#endif
}

/**
## Interface output

This function "draws" interface facets in a file. The segment
endpoints are defined by pairs of coordinates. Each pair of endpoints
is separated from the next pair by a newline, so that the resulting
file is directly visualisable with gnuplot.

The input parameters are a volume fraction field `c`, an optional file
pointer `fp` (which defaults to stdout) and an optional face
vector field `s` containing the surface fractions.

If `s` is specified, the surface fractions are used to compute the
interface normals which leads to a continuous interface representation
in most cases. Otherwise the interface normals are approximated from
the volume fraction field, which results in a piecewise continuous
(i.e. geometric VOF) interface representation. */

trace
void output_facets (scalar c, FILE * fp = stdout, face vector s = {{-1}})
{
  foreach (serial)
    if (c[] > 1e-6 && c[] < 1. - 1e-6) {
      coord n = facet_normal (point, c, s);
      double alpha = plane_alpha (c[], n);
#if dimension == 1
      fprintf (fp, "%g\n", x + Delta*alpha/n.x);
#elif dimension == 2
      coord segment[2];
      if (facets (n, alpha, segment) == 2)
	fprintf (fp, "%g %g\n%g %g\n\n", 
		 x + segment[0].x*Delta, y + segment[0].y*Delta, 
		 x + segment[1].x*Delta, y + segment[1].y*Delta);
#else // dimension == 3
      coord v[12];
      int m = facets (n, alpha, v, 1.);
      for (int i = 0; i < m; i++)
	fprintf (fp, "%g %g %g\n",
		 x + v[i].x*Delta, y + v[i].y*Delta, z + v[i].z*Delta);
      if (m > 0)
	fputc ('\n', fp);
#endif
    }

  fflush (fp);
}

/**
## Interfacial area

This function returns the surface area of the interface as estimated
using its VOF reconstruction. */

trace
double interface_area (scalar c)
{
  double area = 0.;
  foreach (reduction(+:area))
    if (c[] > 1e-6 && c[] < 1. - 1e-6) {
      coord n = interface_normal (point, c), p;
      double alpha = plane_alpha (c[], n);
      area += pow(Delta, dimension - 1)*plane_area_center (n, alpha, &p);
    }
  return area;
}

/**
## Face fraction

We wish to calculate the fraction of face surface occupied by a phase defined
by a volume fraction field. This operation can be useful in different contexts,
for example the solution of the diffusion equation in a specific phase.

The surface fraction is computed as the intersection between the face of the
cell under investigation and the PLIC interface fragment. The problem boils down
to the calculation of the intersection between the interfacial plane and the
coordinate of the cell face ($x = x_o$). In 2D:

$$
  y = \frac{\alpha - m_x x_o}{m_y}
$$

after shifting the reference frame and handling degenerate cases ($m_y = 0$).
In 3D, this concept is extended by computing the area occupied by the phase.
The plane is considered either on the left or on the right side of the face, as
controlled using the boolean `right`. */

foreach_dimension()
static double interface_fraction_x (coord m, double alpha, bool right) {
#if dimension == 1
  return 1;
#elif dimension == 2
  if (fabs (m.y) < 1e-4) // degenerate case
    return right ? (m.x < 0. ? 1. : 0.)
                 : (m.x > 0. ? 1. : 0.);

  alpha += 0.5*(m.x + m.y);
  if (m.y < 0.) {
    alpha -= m.y;
    m.y *= -1;
  }
  double xo = right ? 1. : 0.;
  return clamp ((alpha - m.x*xo)/m.y, 0., 1.);
#elif dimension == 3
  if (fabs (m.y) < 1e-4 && fabs (m.z) < 1e-4) // degenerate case
    return right ? (m.x < 0. ? 1. : 0.)
                 : (m.x > 0. ? 1. : 0.);

  double n1 = m.y/(fabs (m.y) + fabs (m.z));
  double n2 = m.z/(fabs (m.y) + fabs (m.z));
  double j = right ? 0.5 : -0.5;
  alpha -= j*m.x;
  alpha /= (fabs(m.y) + fabs(m.z));
  return clamp (line_area (n1, n2, alpha), 0., 1.);
#endif
}

/**
This function calculates a single surface fraction value in a given face. The
PLIC interface representation cannot guarantee the continuity of the planar
segments across the faces. Therefore, we obtain a single value using a geometric
mean between the left and right sides. The tolerance defines the interfacial
cells and it can be modified.

Note that adaptation of the face fraction `s` is currently not performed through
a dedicated refinement procedure. Consequently, `s` must be recalculated after
each grid adaptation and update of the volume fraction. */

trace
void face_fraction (scalar f, face vector s, double tol = 1.e-6) {

  /**
  We calculate the interface normal field. */

  vector n[];
  foreach() {
    coord m = mycs (point, f);
    foreach_dimension()
      n.x[] = m.x;
  }

  foreach_face() {

    /**
    Case 1: the face is shared between a full and an empty cell. The surface
    fraction is null. */

    if (f[-1] < tol || f[] < tol)
      s.x[] = 0.;

    /**
    Case 2: if both cells are full, the surface fraction is unitary. */

    else if (f[-1] > 1. - tol && f[] > 1. - tol)
      s.x[] = 1.;

    /**
    Case 3: if both cells are interfacial, the face contains the interface, and
    we proceed with the calculation of the surface fraction. */

    else {
      double vleft = 1., vright = 1.;
      if (f[] < 1. - tol) {
        coord m = {0};
        m.x = n.x[];
#if dimension > 1
        m.y = n.y[];
#endif
#if dimension > 2
        m.z = n.z[];
#endif
        double alpha = plane_alpha (f[], m);
        vleft = interface_fraction_x (m, alpha, false);
      }
      if (f[-1] < 1. - tol) {
        coord m = {0};
        m.x = n.x[-1];
#if dimension > 1
        m.y = n.y[-1];
#endif
#if dimension > 2
        m.z = n.z[-1];
#endif
        double alpha = plane_alpha (f[-1], m);
        vright = interface_fraction_x (m, alpha, true);
      }
      s.x[] = sqrt (vleft*vright);
    }
  }
}

