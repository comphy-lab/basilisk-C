/**
# Height-Functions

The contact angle boundary condition is imposed by calculating the heights function in the ghost cells
*/

#define HSHIFT 20.

static inline double height (double H) {
  return H > HSHIFT/2. ? H - HSHIFT : H < -HSHIFT/2. ? H + HSHIFT : H;
}

static inline int orientation (double H) {
  return fabs(H) > HSHIFT/2.;
}

#define BGHOSTS 2
static void half_column (Point point, scalar c, vector h, vector cs, int j)
{
  const int complete = -1;

  foreach_dimension() {
    double S = c[], H = S, ci, a;
    typedef struct { int s; double h; } HState;
    HState state = {0, 0};
    if (j == 1) {
      if (h.x[] == 300.)
	state.s = complete, state.h = nodata;      
      else {
	int s = (h.x[] + HSHIFT/2.)/100.;
	state.h = h.x[] - 100.*s;
	state.s = s - 1;
      }
      if (state.s != complete)
	S = state.s, H = state.h;
    }
    
    for (int i = 1; i <= 4; i++) {
      ci = i <= 2 ? c[i*j] : cs.x[(i - 2)*j];
      H += ci;
      
      
      if (S > 0. && S < 1.) {
	S = ci;
	if (ci <= 0. || ci >= 1.) {
	  
	  H -= i*ci;
	  break;
	}
      }
      
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
      
      else if (S == ci && modf(H, &a))
	break;
    }

    if (j == -1) {
      
      if (S != complete && ((c[] <= 0. || c[] >= 1.) ||
			    (S > 0. && S < 1.)))
	h.x[] = 300.; // inconsistent
      else if (S == complete)
	h.x[] = H;
      else

	
	h.x[] = H + 100.*(1. + (S >= 1.));
    }
    else { // j = 1
	  
      if (state.s != complete ||
	  (S == complete && fabs(height(H)) < fabs(height(state.h))))
	state.s = S, state.h = H;
      
      if (state.s != complete)
	h.x[] = nodata;
      else
	h.x[] = (state.h > 1e10 ? nodata : state.h);
    }
  }
}

static void column_propagation (vector h)
{
  foreach()
    for (int i = -2; i <= 2; i++)
      foreach_dimension()
	if (fabs(height(h.x[i])) <= 3.5 &&
	    fabs(height(h.x[i]) + i) < fabs(height(h.x[])))
	  h.x[] = h.x[i] + i;
  boundary ((scalar *){h});

  foreach_boundary(left)
	{
	h.y[-1,0] = height(h.y []) + (1/tan(theta_0));
	}
}


#if !TREE
trace
void heights (scalar c, vector h)
{
  
  vector cs[];
  foreach_dimension()
    for (int i = 0; i < nboundary; i++)
      cs.x.boundary[i] = c.boundary[i];

  for (int j = -1; j <= 1; j += 2) {
    
    foreach()
      foreach_dimension()
        cs.x[] = c[2*j];
    boundary ((scalar *){cs});

    
    foreach()
      half_column (point, c, h, cs, j);
  }
  boundary ((scalar *){h});
  
  column_propagation (h);
}

#else // TREE
foreach_dimension()
static void refine_h_x (Point point, scalar h)
{

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
    for (int j = -1; j <= 1; j++)
      if (h[0,i,j] == nodata || orientation(h[0,i,j]) != ori)
	return;
      else
	H[i+1][j+1] = height(h[0,i,j]) - H0;

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


trace
void heights (scalar c, vector h)
{
  vector cs[];
  foreach_dimension()
    for (int i = 0; i < nboundary; i++)
      cs.x.boundary[i] = c.boundary[i];
  
  restriction ({c});
  for (int j = -1; j <= 1; j += 2) {
    foreach_level(0)
      foreach_dimension()
        h.x[] = nodata;
  
    for (int l = 1; l <= depth(); l++) {
      
      foreach_level (l)
	foreach_dimension()
	  cs.x[] = c[2*j];
      
      foreach_level (l - 1)
	foreach_dimension() {
	  cs.x[] = c[j];
	  cs.x[j] = c[2*j];
        }
      
      foreach_halo (prolongation, l - 1)
	foreach_dimension()
	  c.prolongation (point, cs.x);
      boundary_iterate (level, (scalar *){cs}, l);

      foreach_level (l)
        half_column (point, c, h, cs, j);
    }
  }
    
  
  foreach_dimension() {
    h.x.prolongation = no_data;
    h.x.restriction = no_restriction;
  }
  boundary ((scalar *){h});

/**
## Contact angle boundary condition implementation

*/

  foreach_boundary(left)
	{
	h.y[-1,0] = height(h.y []) + (1/tan(theta_0));
	}

/**



*/

  foreach_dimension()
    h.x.prolongation = refine_h_x;
  column_propagation (h);
}

#endif // TREE