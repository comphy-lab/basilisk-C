/**
This file makes it possible to impose a contact angle between 45 and 135 degrees for the case when the wall is on the left.
*/

#include "heights_angle.h"

float hoffman(float x){
        float y;
        y = acos(1 -(2*tanh(5.16*pow((x/(1 + 1.31*pow(x,0.99))),0.706))));
        return y;
}

float inverse_h(float x){
        float y0;
        float y1;
        float y2;
        y0 = (1-cos(x))/2;
        y1 = pow(atanh(y0),1.41643);
        y2 = pow(((-0.763359*y1) + 3.93893)/y1,1.0101);
        return y2;
        }
    
    
/**
# Volume fractions
x
These functions are used to maintain or define volume and surface
fractions either from an initial geometric definition or from an
existing volume fraction field. 

We will use basic geometric functions for square cut cells and the
"Mixed-Youngs-Centered" normal approximation of Ruben Scardovelli. */

#include "geometry.h"
#if dimension == 1
coord mycs (Point point, scalar c) {
  coord n = {1.};
  return n;
}
#elif dimension == 2
# include "myc2d.h"
#else // dimension == 3
# include "myc.h"
#endif

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

    coord a, b;
    foreach_dimension() {
      a.x = 0.; b.x = 0.5;
    }
    
    foreach_child() {
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

struct Fractions {
  vertex scalar Phi; // compulsory
  scalar c;          // compulsory
  face vector s;     // optional
};

trace
void fractions (struct Fractions a)
{
  vertex scalar Phi = a.Phi;
  scalar c = a.c;
  face vector s = automatic (a.s);

  /**
  We store the positions of the intersections of the surface with the
  edges of the cell in vector field `p`. In two dimensions, this field
  is just the transpose of the *line fractions* `s`, in 3D we need to
  allocate a new field. */
  
#if dimension == 3
  vector p[];
#else // dimension == 2
  vector p;
  p.x = s.y; p.y = s.x;
#endif
  
  /**
  ### Line fraction computation
  
  We start by computing the *line fractions* i.e. the (normalised)
  lengths of the edges of the cell within the surface. */

  foreach_edge() {

    /**
    If the values of $\Phi$ on the vertices of the edge have opposite
    signs, we know that the edge is cut by the interface. */

    if (Phi[]*Phi[1] < 0.) {

      /**
      In that case we can find an approximation of the interface position by
      simple linear interpolation. We also check the sign of one of the
      vertices to orient the interface properly. */

      p.x[] = Phi[]/(Phi[] - Phi[1]);
      if (Phi[] < 0.)
	p.x[] = 1. - p.x[];
    }

    /**
    If the values of $\Phi$ on the vertices of the edge have the same sign
    (or are zero), then the edge is either entirely outside or entirely
    inside the interface. We check the sign of both vertices to treat
    limit cases properly (when the interface intersects the edge exactly
    on one of the vertices). */

    else
      p.x[] = (Phi[] > 0. || Phi[1] > 0.);
  }

  /**
  ### Surface fraction computation 

  We can now compute the surface fractions. In 3D they will be
  computed for each face (in the z, x and y directions) and stored in
  the face field `s`. In 2D the surface fraction in the z-direction is
  the *volume fraction* `c`. The call to `boundary_flux()` defines
  consistent line fractions on trees. */

#if dimension == 3
  scalar s_x = s.x, s_y = s.y, s_z = s.z;
  foreach_face(z,x,y)
#else // dimension == 2
  boundary_flux ({s});
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
	    double a = sign(Phi[0,i])*(p.x[0,i] - 0.5);
	    alpha += n.x*a + n.y*(i - 0.5);
	    ni++;
	  }

      /**
      Once we have $\mathbf{n}$ and $\alpha$, the (linear) interface
      is fully defined and we can compute the surface fraction using
      our pre-defined function. For marginal cases, the cell is full
      or empty (*ni == 0*) and we look at the line fractions to
      decide. */

      s_z[] = ni ? line_area (n.x, n.y, alpha/ni) : max (p.x[], p.y[]);
    }
  }

  /**
  ### Volume fraction computation

  To compute the volume fraction in 3D, we use the same approach. */
  
#if dimension == 3
  boundary_flux ({s});
  foreach() {

    /**
    Estimation of the average normal from the surface fractions. */
       
    coord n;
    double nn = 0.;
    foreach_dimension(3) {
      n.x = s.x[] - s.x[1];
      nn += fabs(n.x);
    }
    if (nn == 0.)
      c[] = s.x[];
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
	      double a = sign(Phi[0,i,j])*(p.x[0,i,j] - 0.5);
	      alpha += n.x*a + n.y*(i - 0.5) + n.z*(j - 0.5);
	      ni++;
	    }

      /**
      Finally we compute the volume fraction. */
      
      c[] = ni ? plane_volume (n, alpha/ni) : s.x[];
    }
  }
#endif
  
  /**
  Finally we apply the boundary conditions. */

  boundary ({c});

int count = 0;
foreach_boundary(left){
  if (c[] < 1. && c[] > 0.){
  count += 1;
  }
  }  

  vector h[];
  //float theta = theta_0;


                             
  heights(c, h, theta);
   
  foreach_boundary(left)
	{
	
    	if (c[] < 1. && c[] > 0.)
      		{                             
                                
			if(theta <= pi/2. && c[1,-1] == 0)
                          {
				c[-1,-1] = 1.;
                        	double contactline;
				double xghost;
				double x1;
				x1 =((y + (Delta*height(h.y[])))/Delta) +  (((cos(theta)/sin(theta))/2.));
				contactline = x1 - floor(y/Delta);                   
				xghost = contactline + ((cos(theta)/sin(theta)));
				if(xghost < 1.)
				{
				c[-1, 1] = 0.; 
				c[-1, 0] = contactline + ((cos(theta)/sin(theta))/2.);


				}
				else if (xghost >= 1. && contactline >=1)
				{
				c[-1, 1] = ((sq(xghost-1))*tan(theta)/2.);
				c[-1, 0] = 1. - ((sq(1-contactline))*tan(theta)/2.);
				}
				else if (contactline >= 1.)
				{
				
				c[-1, 1] = ((sq(2. - contactline))*(tan(theta)/2));
				c[-1, 0] = 1.;
			        c[-1, 2] = (sq(xghost-2.))*tan(theta)/2.;
				c[0, 1]  = (sq(contactline - 1.))*tan(theta)/2.;
				

				} 
			}

			else if(theta > pi/2 && f[-1,-1] == 1)
			{
				c[-1, 1] = 0;
                        	double contactline;
				double xghost;
				double x1;
				double phi = pi - theta;
				double x0;
				x0 = ((y + (Delta*height(h.y[])))/Delta) - floor(y/Delta);
				x1 =((y + (Delta*height(h.y[])))/Delta) -  (((cos(phi)/sin(phi))/2.));
				contactline = x1 - floor(y/Delta);                   
				xghost = contactline - ((cos(phi)/sin(phi)));
				if (xghost> 0)
				{
				c[-1, -1] = 1;
				c[-1, 0] = 1 -(contactline + xghost)/2;
				}
				else if (xghost <= 0. && contactline>0.)
				{
				c[-1, 0] = ((sq(contactline)*tan(phi))/2.);
				c[-1, -1] = 1 - ((sq(contactline)*tan(phi))/2.);
				}
				else if (contactline <= 0.)
				{
				c[-1,0] = 0;
				c[-1, -1] = ((sq(contactline)*tan(phi))/2.);
				c[0, -1] = 1 - ((sq(contactline)*tan(phi))/2.);
				c[] = c[] - ((sq(x0)*tan(phi)) + (sq(x0)*tan(phi)));
				}
			}
      		 }
	}
       	


}


/**
The convenience macro below can be used to define a volume fraction
field directly from a function. */

#define fraction(f,func) do {			\
    vertex scalar phi[];			\
    foreach_vertex()				\
      phi[] = func;				\
    fractions (phi, f);				\
  } while(0)

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
### Interface reconstruction 

The reconstruction function takes a volume fraction field `c` and
returns the corresponding normal vector field `n` and intercept field
$\alpha$. */

trace
void reconstruction (scalar c, vector n, scalar alpha)
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

      coord m = mycs (point, c);
      // coord m = youngs_normal (point, c);
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

  /**
  Finally we apply the boundary conditions to define $\mathbf{n}$ and
  $\alpha$ everywhere (using the prolongation functions when necessary
  on tree grids). */

  boundary ({n, alpha});
	

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

struct OutputFacets {
  scalar c;
  FILE * fp;     // optional: default is stdout
  face vector s; // optional: default is none
};

trace
void output_facets (struct OutputFacets p)
{
  scalar c = p.c;
  face vector s = p.s;


  if (!p.fp) p.fp = stdout;
  foreach()
    if (c[] > 1e-6 && c[] < 1. - 1e-6) {
      coord n;
      if (!s.x.i){
      //compute normal from volume fraction
	n = mycs (point, c);
       }
      else {
	// compute normal from face fractions
	double nn = 0.;
	foreach_dimension() {
	  n.x = s.x[] - s.x[1];
	  nn += fabs(n.x);
	}
	assert (nn > 0.);
	foreach_dimension()
	  n.x /= nn;
      }

      double alpha = plane_alpha (c[], n);
#if dimension == 2      
      coord segment[2];
      if (facets (n, alpha, segment) == 2)
	fprintf (p.fp, "%g %g\n%g %g\n\n", 
		 x + segment[0].x*Delta, y + segment[0].y*Delta, 
		 x + segment[1].x*Delta, y + segment[1].y*Delta);
#else // dimension == 3
      coord v[12];
      int m = facets (n, alpha, v, 1.);
      for (int i = 0; i < m; i++)
	fprintf (p.fp, "%g %g %g\n",
		 x + v[i].x*Delta, y + v[i].y*Delta, z + v[i].z*Delta);
      if (m > 0)
	fputc ('\n', p.fp);
#endif
    }

  fflush (p.fp);
}

/**
## Interfacial area

This function returns the surface area of the interface as estimated
using its VOF reconstruction. */

trace
double interface_area (scalar c)
{
  double area = 0.;
  foreach()
    if (c[] > 1e-6 && c[] < 1. - 1e-6) {
      coord n = mycs (point, c), p;
      double alpha = plane_alpha (c[], n);
      area += pow(Delta, dimension - 1)*plane_area_center (n, alpha, &p);
    }
  return area;
}
