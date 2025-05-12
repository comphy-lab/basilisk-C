/**
# Discharge.h

We define several functions in order to impose a flux at boundaries
*/

/**
## Boundary structure

This structure stores all what we need in order to define the boundary
: a scalar "river" which is equal to zero or 1, it defines the river
location, two doubles “xb" and "yb” to store the x and y coordinates
of the boundary, the double “dx” equal to the maximum spatial
resolution (this can be optimised) and an integer kw which stores the
relative location of the boundary (top, bottom, left and right).
*/

typedef struct {
  scalar river;
  double value;
  double xb;
  double yb;
  double dx;
  int kw;
} Defb;

/**
## Real flux

This function return the exact incoming flow on the boundary for an
imposed water height
 */

static double realflux (Defb b, double eta_imp)
{
  double q = 0, x = b.xb, y = b.yb, dtmax = 2000;
  scalar lim = b.river;
@if _QCCACC
  ACC(acc update device(N, X0, Y0, L0, G, CFL))
  ACC(acc kernels loop independent present(acc_dataarr) reduction(+:q))
  for (int j = 0; j < N-1; j++) {
    if (b.kw == bottom || b.kw == top) x = b.xb + j * b.dx;
    else y = b.yb + j * b.dx;
    Point point;
    locate2Dacc(x, y, &point);
    if (lim[] == b.value) {
      double ub;
      if (b.kw == bottom || b.kw == top)
	ub = b.kw == bottom ? u.y[] : - u.y[];
      else
	ub =  b.kw == left ? u.x[] : - u.x[];
      double hn = max (eta_imp - zb[], 0);
      if (h[] > dry || hn > dry) { 
	double fh, fu;
	kurganov (hn, h[], ub, ub, b.dx, &fh, &fu, &dtmax);
	q += fh*b.dx;
      }
    }         
  }
@else
  for (int j = 0; j < N-1; j++) {
    if (b.kw == bottom || b.kw == top) x = b.xb + j * b.dx;
    else y = b.yb + j * b.dx;
    Point point = locate (x, y);
    if (lim[] == b.value) {
      double ub;
      if (b.kw == bottom || b.kw == top)
	ub = b.kw == bottom ? u.y[] : - u.y[];
      else
	ub =  b.kw == left ? u.x[] : - u.x[];
      double hn = max (eta_imp - zb[], 0);
      if (h[] > dry || hn > dry) { 
	double fh, fu;
	kurganov (hn, h[], ub, ub, b.dx, &fh, &fu, &dtmax);
	q += fh*b.dx;
      }
    }         
  }
@endif
  return q;
}

/**
## Structure definition
*/

static void defbound (int keyw, scalar lim, double value, Defb * b)
{
  b->dx = L0/N;
  b->river = lim;
  b->value = value;
  b->kw = keyw;
  b->xb = X0 + 0.5 * b->dx;
  b->yb = Y0 + 0.5 * b->dx;
  if( keyw == right ) b->xb += L0 - b->dx;
  else if( keyw == top )  b->yb += L0 - b->dx;
  else if( keyw != left && keyw != bottom ) assert (false);
}

/**
## Altitude

The altitude function find the minimum value of the topography. It
stores this value and also the one of the water height. */

static void altitude (Defb b,  double * zmin, double * h0)
{
  double ztemp = *zmin, htemp = 0, x = b.xb, y = b.yb;
  scalar lim = b.river;

@if _QCCACC
  /* Note: this kernel will only run in scalar mode on the accelerator
     due to the simultaneous and complex reduction of two variables */ 
  ACC(acc update device(N, X0, Y0, L0, G, CFL))
  ACC(acc kernels loop present(acc_dataarr))
  for (int j = 0; j < N-1; j++) {
    if( b.kw == top || b.kw == bottom ) x = b.xb + j * b.dx;
    else y = b.yb + j * b.dx;
    Point point;
    locate2Dacc(x, y, &point);
    if (zb[] < ztemp && lim[] == b.value) {
      ztemp = zb[];
      htemp = h[];
    }
  }
@else
#if _OPENMP
#pragma omp parallel for reduction(min:ztemp)
#endif

  for (int j = 0; j < N-1; j++) {
    if( b.kw == top || b.kw == bottom ) x = b.xb + j * b.dx;
    else y = b.yb + j * b.dx;
    Point point = locate (x,y);
    if (zb[] < ztemp && lim[] == b.value) {
      ztemp = zb[];
      htemp = h[];
    }
  }
@endif
  *zmin = ztemp;
  *h0 = htemp;
}

/**
## False position method

Here, we find the water height corresponding to the imposed inflow by
the [false position
method](http://en.wikipedia.org/wiki/False_position_method). */

static double metfalsepos (Defb bo, double h0, double zmin,
			   double prec, double q_imp)
{
  double binf, bsup, qinf, newb, qsup;
  
  // inflow at n-1
  double newq = realflux (bo, h0 + zmin);
  if (fabs((newq - q_imp)/q_imp) <= prec)
    return h0 + zmin;
  
  // Find sup & inf boundaries
  if (newq > q_imp) {
    bsup = h0 + zmin;
    qsup = newq;
    binf = zmin;
    qinf = 0;
  }
  else {
    // We must take care of the case h0=0
    if (h0 <= dry) {
      h0 = 1;
      bsup = zmin;
    }
    else  bsup = h0 + zmin;
    do {
      // here we find a "good" superior limit
      binf = bsup;
      bsup = bsup + h0;
      qsup = realflux (bo, bsup);
    } while (qsup <= q_imp);    
    qinf = realflux (bo, binf);
  }

  // The false position method
  do {
    double alpha = (qsup - qinf)/(bsup - binf);
    newb = binf + (q_imp - qinf)/alpha;
    newq = realflux (bo, newb);
    if (newq > q_imp) {
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

double discharge (double q_imp, int keyw, scalar limit, double value)
{  
  
  // Inialisation of the Defb structure
  Defb bo;
  defbound (keyw, limit, value, &bo);
  
  // Finding zmin et h
  double h0 = 0, zmin = 200000;
  altitude (bo, &zmin, &h0);
  
  // If imposed inflow <= 0
  if (q_imp <= 0) return zmin - 0.1;

  // Return eta found by Newton met
  return metfalsepos (bo, h0, zmin, 0.001, q_imp);
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
