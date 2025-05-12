/**
# Discharge.h

We define several functions in order to impose a flux at
boundaries. The discharge function can work with adaptative
refinement. We first define the scalar river[] as a global scalar.*/

scalar river[];

/**
## Adaptation

When using an adaptive discretisation, we define the
refinement laws for the scalar river[]. The refine_injection function
give  to each child the parent value. The coarsen_max function fix the
value of the parent as the maximum value of all the children.*/

#if QUADTREE
static void coarsen_max (Point point, scalar s){
  s[] = max(max(fine(s,0,0),fine(s,1,0)),max(fine(s,0,1),fine(s,1,1)));
}
event defaults( i = 0 ){
 river.refine = river.prolongation = refine_injection;
 river.coarsen = coarsen_max;
}
#endif

/**
## Boundary structure "defb"

This structure stores all we need in order to define the boundary : a
double value which test the location of the scalar river[], two
doubles “xb" and "yb” to store the x and y coordinates of the
boundary, “dx” equal to the minimum spatial resolution, "zmin" which
store the minimum topography elevation on the boundary, "h0" the water
elevation on zmin, "lim" the limit of the boundary and an integer "kw"
which stores the relative location of the boundary (top, bottom, left
and right).*/

typedef struct {
  double value, xb, yb, dx, zmin, h0, lim;
  int kw;
} defb;

/**
## Structure definition

We allocate the right values to the "defb" structure*/

static defb defbf (int keyw, double value){
  defb b;
  b.dx = L0/N;
  b.value = value;
  b.kw = keyw;
  b.xb = X0 + 0.5 * b.dx;
  b.yb = Y0 + 0.5 * b.dx;
  b.zmin= 1e40;
  b.h0 = 0;
  if( keyw == right ) b.xb = X0 + L0 - 0.5 * b.dx;
  else if( keyw == top )  b.yb = X0 + L0 - 0.5 * b.dx;
  else if( keyw != left && keyw != bottom ){
    fprintf(stderr,"wrong keyword in discharge function\n");
    assert (false);
  }
  if (b.kw == bottom || b.kw == top)
    b.lim = X0 + L0;
  else
    b.lim = Y0 + L0;
  return b;
}

/**
## Altitude

The altitude function find the minimum value of the topography and the
water depth at this location.*/

static void altitude (defb * b) {
  double ztemp = b->zmin, htemp = 0, x = b->xb, y = b->yb;
  
# if _OPENMP
#   pragma omp parallel for reduction(min:ztemp)
# endif
  for (int j = 0; j < N-1; j++) { 
    if( b->kw == top || b->kw == bottom ) x = b->xb + j * b->dx;
    else y = b->yb + j * b->dx;
    Point point = locate (x,y);
    if (zb[] < ztemp && river[] == b->value) {
      ztemp = zb[];
      htemp = h[];
    }
  }
  b->zmin = ztemp;
  b->h0 = htemp;
}

/**
## Real flux

This function return the exact incoming flow on the boundary for an
imposed water depth.*/

static double realflux (defb b, double eta_imp) {
  double q = 0, deltax = b.dx , xp = b.xb, yp = b.yb, dtmax = 1;
  double bx;
  
  do {
    Point point = locate (xp, yp);
    if (river[] == b.value) {
#     if QUADTREE
      // Relocate at the center of the element only on the sliding coordinate
      if(b.kw == bottom || b.kw == top ) xp = x;
      else yp = y;
      // To not miss refined elements
      deltax = 0.75*Delta;
#     endif
	// Take care of the sign of velocity
	double ub;
	if (b.kw == bottom || b.kw == top)
	  ub = b.kw == bottom ? u.y[] : - u.y[];
	else
	  ub =  b.kw == left ? u.x[] : - u.x[];

	// computing the local numerical flux
	double hn = max (eta_imp - zb[], 0);
	if (h[] > dry || hn > dry) { 
	  double fh, fu;
	  kurganov (hn, max(h[],0), ub, ub, b.dx, &fh, &fu, &dtmax);
	  q += fh*b.dx;
	}
    }
    
    // Testing the boundary limit and advancing the local point
    if (b.kw == bottom || b.kw == top) {
      xp += deltax;
      bx = xp;
    }
    else {
      yp += deltax;
      bx = yp;
    }
  } while( bx < b.lim );
  
  return q;
}

/**
## Double false position method

Here, we find the water depth corresponding to the imposed inflow by
the [false position
method](http://en.wikipedia.org/wiki/False_position_method). This
method is a root-finding method of type "guess and check". In our
case, we want to solve the equation :

$$
flux(\eta) - q_{imp}=0
$$

From wikipedia : "The method proceeds by producing a sequence of
shrinking intervals [ak, bk] that all contain a root of f.(...) At
each step, the method compute ck which is the root of the secant line
through (ak, f(ak)) and (bk, f(bk)). If f(ak) and f(ck) have the same
sign, then we set ak+1 = ck and bk+1 = bk, otherwise we set ak+1 = ak
and bk+1 = ck. This process is repeated until the root is approximated
sufficiently well."*/

static double metfalsepos (defb bo, double prec, double q_imp) {
  double binf, bsup, qinf, newb, qsup;
  double eta0 = bo.zmin + bo.h0;
  /**
     First, we test the inflow at the precedent time step. */
  
  double newq = realflux (bo, eta0);
  if (fabs((newq - q_imp)/q_imp) <= prec)
    return eta0;
  
  /**
     Then, we find an inferior and a superior limit around the
     solution.*/
  
  // Find inferior limit
    if (newq > q_imp) {
    bsup = eta0;
    qsup = newq;
    binf = bo.zmin;
    qinf = 0;
  }
  else {
    // We must take care of the case h0=0
    if (bo.h0 <= dry) {
      bo.h0 = 1;
      bsup = bo.zmin;
    }
    else  bsup = eta0;
    do {
      // Find a "good" superior limit
      binf = bsup;
      bsup = bsup + bo.h0;
      qsup = realflux (bo, bsup);
    } while (qsup <= q_imp);    
    qinf = realflux (bo, binf);
  }
  // The false position method
  do {
    double alpha = (qsup - qinf)/(bsup - binf); // Define the slope
    newb = binf + (q_imp - qinf)/alpha; // new bound
    newq = realflux (bo, newb); // new total inflow
    if (newq > q_imp) { //test the new solution & replacing the good limit
      bsup = newb;
      qsup = newq;
    }
    else {
      binf = newb;
      qinf = newq;
    }
  } while (fabs((newq - q_imp)/q_imp) >= prec);

  return newb;
}

/**
## The discharge routine */

double discharge (double q_imp, int keyw, double value)
{
 
   // Inialisation of the Defb structure
  defb bo = defbf(keyw, value);
    
  // Fix zmin et h
  altitude (&bo);
  
  // If imposed inflow <= 0
  if (q_imp <= 0) return bo.zmin - 0.0001;
  
  // Return eta found by the false position method
  return metfalsepos (bo, 0.001, q_imp);
}

/**
# Hydrograph

This function return the inflow by linearly interpolating the data. */

double hydrograph (const char * name, double mult)
{
  FILE * fp;
  if ((fp = fopen (name , "r")) == NULL) {
    fprintf (stderr,"cannot open hydrograph data file.\n name = %s ",name);
    assert (false);
  }
  double time = 0, timea, q = 0, qa, alpha = 0;
  // We read the data at each call => Can definitively be optimised !
  do {
    qa = q;
    timea = time;
    if (fscanf (fp, "%lf \t %lf \n", &time, &q) == EOF) break;
    time *= mult;
    if (time - timea != 0) alpha = (q - qa)/(time - timea);
    else alpha = 0;
  } while (time < t);
  fclose (fp);
  /*
    If the solver time is sup to the larger time of the data, return the last value   */
  if (timea  >= t)  return q;
  // Else, return the interpolation
  else return alpha*(t - timea) + qa;
}
