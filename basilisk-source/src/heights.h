/**
# Height-Functions

The "height-function" is a vector field which gives the distance,
along each coordinate axis, from the center of the cell to the closest
interface defined by a volume fraction field. This distance is
estimated using the "column integral" of the volume fraction in the
corresponding direction. This integral is not always defined (for
example because the interface is too far i.e. farther than 5.5 cells
in our implementation) in which case the value of the field is set to
*nodata*. See e.g. [Popinet, 2009](references.bib#popinet2009) for
more details on height functions.

We also store the "orientation" of the height function together with
its value by adding *HSHIFT* if the volume fraction is unity on the
"top" end. The function below applied to the value will return the
corresponding height and orientation.

The distance is normalised with the cell size so that the coordinates
of the interface are given by

~~~c
(x, y + Delta*height(h.y[])) or (x + Delta*height(h.x[]), y)
~~~
*/

#define HSHIFT 20.

static inline double height (double H) {
  return H > HSHIFT/2. ? H - HSHIFT : H < -HSHIFT/2. ? H + HSHIFT : H;
}

static inline int orientation (double H) {
  return fabs(H) > HSHIFT/2.;
}

/**
We make sure that two layers of ghost cells are defined on the
boundaries (the default is one layer). */

#define BGHOSTS 2

/**
## Half-column integration 

This helper function performs the integration on half a column, either
"downward" (*j = -1*) or "upward" (*j = 1*). */

static void half_column (Point point, scalar c, vector h, vector cs, int j)
{

  /**
  The 'state' of the height function can be: *complete* if both
  ends were found, zero or one if one end was found and between zero
  and one if only the interface was found. */

  const int complete = -1;

  foreach_dimension() {

    /**
     *S* is the state and *H* the (partial) value of the height
     function. If we are on the (first) downward integration (*j =
     -1*) we initialise *S* and *H* with the volume fraction in
     the current cell. */

    double S = c[], H = S, ci;

    /**
     On the upward integration (*j = 1*), we recover the state of the
     downward integration. Both the state and the (possibly partial)
     height value are encoded in a single number using a base 100
     shift for the state. */
    
    typedef struct { int s; double h; } HState;
    HState state = {0, 0};
    if (j == 1) {
      
      /**
      Check whether this is an inconsistent height. */

      if (h.x[] == 300.)
	state.s = complete, state.h = nodata;

      /**
      Otherwise, this is either a complete or a partial height. */
      
      else {
	int s = (h.x[] + HSHIFT/2.)/100.;
	state.h = h.x[] - 100.*s;
	state.s = s - 1;
      }

      /**
      If this is a complete height, we start a fresh upward
      integration. */
      
      if (state.s != complete)
	S = state.s, H = state.h;
    }
    
    /**
     We consider the four neighboring cells of the half column, the
     corresponding volume fraction *ci* is recovered either from the
     standard volume fraction field *c* (first two cells) or from the
     shifted field *cs* (last two cells). The construction of *cs* is
     explained in the next section. */
    
    for (int i = 1; i <= 4; i++) {
      ci = i <= 2 ? c[i*j] : cs.x[(i - 2)*j];
      H += ci;
      
      /**
       We then check whether the partial height is complete or not. */
      
      if (S > 0. && S < 1.) {
	S = ci;
	if (ci <= 0. || ci >= 1.) {
	  
	  /**
	   We just left an interfacial cell (*S*) and found a full or
	   empty cell (*ci*): this is a partial height and we can stop
	   the integration. If the cell is full (*ci = 1*) we shift
	   the origin of the height. */
	  
	  H -= i*ci;
	  break;
	}
      }
      
      /**
       If *S* is empty or full and *ci* is full or empty, we went
       right through he interface i.e. the height is complete and
       we can stop the integration. The origin is shifted
       appropriately and the orientation is encoded using the "HSHIFT
       trick". */
      
      else if (S >= 1. && ci <= 0.) {
	H = (H - 0.5)*j + (j == -1)*HSHIFT;
	S = complete;
	break;
      }
      else if (S <= 0. && ci >= 1.) {
	H = (i + 0.5 - H)*j + (j == 1)*HSHIFT;
	S = complete;
	break;
      }
      
      /**
       If *ci* is identical to *S* (which is empty or full), we
       check that *H* is an integer (i.e. its fractional value is
       zero), otherwise we are in the case where we found an
       interface but didn't go through it: this is an
       inconsistent height and we stop the integration. */
      
      else if (S == ci && trunc(H) != H)
	break;
    }

    /**
     We update the global state using the state *S* of the
     half-integration. */
    
    if (j == -1) {

      /**
       For the downward integration, we check that the partial heights
       (*S != complete*) are consistent: if the first cell is full
       or empty or if the last cell is interfacial, the partial
       height is marked as inconsistent. */
      
      if (S != complete && ((c[] <= 0. || c[] >= 1.) ||
			    (S > 0. && S < 1.)))
	h.x[] = 300.; // inconsistent
      else if (S == complete)
	h.x[] = H;
      else

	/**
	This is a partial height: we encode the state using a base 100
	shift. */
	
	h.x[] = H + 100.*(1. + (S >= 1.));
    }
    else { // j = 1

      /**
       For the upward integration, we update the current *state*
       using the state of the half-integration *S* only if the
       first downward integration returned a partial height, or if
       the upward integration returned a complete height with a
       smaller value than the (complete) height of the downward
       integration. */
	  
      if (state.s != complete ||
	  (S == complete && fabs(height(H)) < fabs(height(state.h))))
	state.s = S, state.h = H;
      
      /**
       Finally, we set the vector field *h* using the state and
       height. */
      
      if (state.s != complete)
	h.x[] = nodata;
      else
	h.x[] = (state.h > 1e10 ? nodata : state.h);
    }
  }
}

/**
## Column propagation

Once columns are computed on a local 9-cells-high stencil, we will
need to "propagate" these values upward or downward so that they are
accessible at distances of up to 5.5 cells from the interface. This is
important in 3D in particular where marginal (~45 degrees) cases may
require such high stencils to compute consistent HF curvatures. We do
this by selecting the smallest height in a 5-cells neighborhood along
each direction. */

static void column_propagation (vector h)
{
#if _OPENMP 
  foreach (serial) // not compatible with OpenMP
#else
  foreach()
#endif
    for (int i = -2; i <= 2; i++)
      foreach_dimension()
	if (fabs(height(h.x[i])) <= 3.5 &&
	    fabs(height(h.x[i]) + i) < fabs(height(h.x[])))
	  h.x[] = h.x[i] + i;
}

/**
## Multigrid implementation

The *heights()* function takes a volume fraction field *c* and returns
the height function vector field *h*. */

#if !TREE
trace
void heights (scalar c, vector h)
{

  /**
  We need a 9-points-high stencil (rather than the default
  5-points). To do this we store in *s* the volume fraction field *c*
  shifted by 2 grid points in the respective directions. We make sure
  that this field uses the same boundary conditions as *c*. */
  
  vector s[];
  foreach_dimension()
    for (int i = 0; i < nboundary; i++)
      s.x.boundary[i] = c.boundary[i];

  /**
  To compute the height function, we sum the volume fractions in a
  (half-)column starting at the current cell. We start by integrating
  downward (*j = -1*) and then integrate upward (*j = 1*). */
  
  for (int j = -1; j <= 1; j += 2) {

    /**
    We first build the shifted (by $\pm 2$) volume fraction field in each
    direction. */
    
    foreach()
      foreach_dimension()
        s.x[] = c[2*j];

    /**
    We sum the half-column, downward or upward. We (exceptionally)
    need to allow for stencil overflow. */

    foreach (overflow)
      half_column (point, c, h, s, j);
  }

  /**
  Finally we "propagate" values along columns. */
  
  column_propagation (h);
}

/**
## Tree implementation 

We first define the prolongation functions for heights. */

#else // TREE
foreach_dimension()
static void refine_h_x (Point point, scalar h)
{

  /**
  We try to prolongate columns from nearby non-prolongation cells. */

  bool complete = true;
  foreach_child() {
    for (int i = -2; i <= 2; i++)
      if (allocated(i) &&
	  !is_prolongation(neighbor(i)) && !is_boundary(neighbor(i)) &&
	  fabs(height(h[i])) <= 3.5 &&
	  fabs(height(h[i]) + i) < fabs(height(h[])))
	h[] = h[i] + i;
    if (h[] == nodata)
      complete = false;
  }
  if (complete)
    return;

  /**
  If some children have not been initialised, we first check that the
  (three in 2D, nine in 3D) coarse heights are defined and have
  compatible orientations. If not, the children heights are
  undefined. Otherwise, a (bi)quadratic fit of the coarse heights is
  used to compute the children heights. */

  int ori = orientation(h[]);
#if dimension == 2
  for (int i = -1; i <= 1; i++)
    if (h[0,i] == nodata || orientation(h[0,i]) != ori)
      return;

  double h0 = (30.*height(h[]) + height(h[0,1]) + height(h[0,-1]))/16.
    + HSHIFT*ori;
  double dh = (height(h[0,1]) - height(h[0,-1]))/4.;
  foreach_child()
    if (h[] == nodata)
      h[] = h0 + dh*child.y - child.x/2.;
#else // dimension == 3
  double H[3][3], H0 = height(h[]);
  for (int i = -1; i <= 1; i++)
    for (int j = -1; j <= 1; j++) {
      if (h[0,i,j] == nodata || orientation(h[0,i,j]) != ori)
	return;
      else
	H[i+1][j+1] = height(h[0,i,j]) - H0;
    }
      
  double h0 = 
    2.*H0 + (H[2][2] + H[2][0] + H[0][0] + H[0][2] +
	     30.*(H[2][1] + H[0][1] + H[1][0] + H[1][2]))/512.
    + HSHIFT*ori;
  double h1 = (H[2][2] + H[2][0] - H[0][0] - H[0][2] +
	       30.*(H[2][1] - H[0][1]))/128.;
  double h2 = (H[2][2] - H[2][0] - H[0][0] + H[0][2] +
	       30.*(H[1][2] - H[1][0]))/128.;
  double h3 = (H[0][0] + H[2][2] - H[0][2] - H[2][0])/32.;
  foreach_child()
    if (h[] == nodata)
      h[] = h0 + h1*child.y + h2*child.z + h3*child.y*child.z - child.x/2.;
#endif // dimension == 3
}

/**
The *heights()* function implementation is similar to the multigrid
case, but the construction of the shifted volume fraction field *cs*
is more complex. */

trace
void heights (scalar c, vector h)
{
  vector s[];
  foreach_dimension()
    for (int i = 0; i < nboundary; i++)
      s.x.boundary[i] = c.boundary[i];

  /**
  To compute the shifted field, we first need to *restrict* the volume
  fraction on all levels. */
  
  restriction ({c});
  for (int j = -1; j <= 1; j += 2) {

    /**
    We traverse the tree level by level, from coarse to fine.
    On the root cell the height function is undefined. */
  
    foreach_level(0)
      foreach_dimension()
        h.x[] = nodata;
  
    for (int l = 1; l <= depth(); l++) {

      /**
      We construct the ($\pm 2$) shifted field at this level. */
      
      foreach_level (l)
	foreach_dimension()
	  s.x[] = c[2*j];

      /**
      We then need to apply boundary conditions on the shifted
      field. This is more complex than for a constant resolution grid.

      We first construct the ($\pm 1$) shifted field for the
      immediately coarser level. This is done by copying the volume
      fraction field for pairs of adjacent cells. */
      
      foreach_level (l - 1)
	foreach_dimension() {
	  s.x[] = c[j];
	  s.x[j] = c[2*j];
        }

      /**
      We can now use this shifted coarse field (which matches the
      shifted fine field) to apply boundary conditions on coarse/fine
      prolongation halos. */
      
      foreach_halo (prolongation, l - 1)
	foreach_dimension()
	  c.prolongation (point, s.x);
      boundary_iterate (level, (scalar *){s}, l);

      /**
      We can now sum the half-column at this level, downward or upward
      according to *j*. */

      foreach_level (l)
        half_column (point, c, h, s, j);
    }
  }
    
  /**
  We fill the prolongation cells with "nodata". The restriction
  function does nothing as we have already defined *h* on all
  levels. */
  
  foreach_dimension() {
    h.x.prolongation = no_data;
    h.x.restriction = no_restriction;
    h.x.dirty = true;
  }

  /**
  Finally, we "propagate" values along columns. */
  
  column_propagation (h);
  
  /**
  Final prolongation cells will be filled with values obtained either
  from neighboring columns or by interpolation from coarser levels
  (see *refine_h_x()* above). */
  
  foreach_dimension()
    h.x.prolongation = refine_h_x;
}

#endif // TREE

/**
An attribute is added so that the height function field can be
associated to a (VOF) tracer. */

attribute {
  vector height;
}
