extern scalar f;

#if dimension == 2
foreach_dimension()
static double kappa_ebm_y (Point point, vector h)
{
  if (h.y[] == nodata)
    return nodata;

  int ori = orientation(h.y[]);
  double hl = nodata, hr = nodata, hc = h.y[];
  if (h.y[-1] != nodata && orientation(h.y[-1]) == ori)
    hl = h.y[-1];
  if (h.y[1] != nodata && orientation(h.y[1]) == ori)
    hr = h.y[1];
  if (hl == nodata && hr == nodata)
    return nodata;

  scalar cfg = f.cfg;
  if (cfg[] == 1.) {
    vector mcl = f.mcl;
    vector hp  = f.hp;
    double hx = nodata, hxx = nodata;
    double hxcl = - mcl.x[]/mcl.y[];
    double hcl = hc + hxcl*hp.x[];
    if (hl != nodata && cfg[-1] != 1.) {
      hx = (hcl - hl)/(1. + hp.x[]);
      hxx = 2.*(hxcl - hx)/(1. + hp.x[])/Delta;
    }
    else if (hr != nodata && cfg[1] != 1.) {
      hx = (hr - hcl)/(1. - hp.x[]);
      hxx = 2.*(hx - hxcl)/(1. - hp.x[])/Delta;
    }
    if (hx != nodata && hxx != nodata)
      return hxx/pow(1. + sq(hxcl), 3/2.);
  }

  if (hl == nodata || hr == nodata) {
    double hll = nodata, hrr = nodata;
    if (h.y[-2] != nodata && orientation(h.y[-2]) == ori)
      hll = h.y[-2];
    if (h.y[2] != nodata && orientation(h.y[2]) == ori)
      hrr = h.y[2];
    if (hl != nodata && hll != nodata) {
      hr = hc;
      hc = hl;
      hl = hll;
    }
    if (hr != nodata && hrr != nodata) {
      hl = hc;
      hc = hr;
      hr = hrr;
    }
  }

  if (hl == nodata || hr == nodata)
    return nodata;
  double hx = (hr - hl)/2.;
  double hxx = (hr + hl - 2.*hc)/Delta;
  return hxx/pow(1. + sq(hx), 3/2.);
}

foreach_dimension()
static coord normal_ebm_y (Point point, vector h)
{
  coord n = {nodata, nodata, nodata};

  scalar cfg = f.cfg;
  if (cfg[] == 1.) {
    vector mcl = f.mcl;
    n.x = mcl.x[]/(sq(mcl.x[]) + sq(mcl.y[]));
    n.y = mcl.y[]/(sq(mcl.x[]) + sq(mcl.y[]));
    return n;
  }

  if (h.y[] == nodata)
    return n;

  int ori = orientation(h.y[]);
  double hl = nodata, hr = nodata, hc = h.y[];
  if (h.y[-1] != nodata && orientation(h.y[-1]) == ori)
    hl = h.y[-1];
  if (h.y[1] != nodata && orientation(h.y[1]) == ori)
    hr = h.y[1];
  if (hl != nodata && hr != nodata)
    n.x = (hl - hr)/2.;
  else if (hl != nodata)
    n.x = hl - hc;
  else if (hr != nodata)
    n.x = hc - hr;
  else
    return n;

  double nn = (ori ? -1. : 1.)*sqrt(1. + sq(n.x));
  n.x /= nn;
  n.y = 1./nn;
  return n;
}
#else // dimension == 3
foreach_dimension()
static double kappa_ebm_z (Point point, vector h)
{
	return nodata;
}

foreach_dimension()
static coord normal2_ebm_z (Point point, vector h)
{
  return (coord){nodata, nodata, nodata};
}

foreach_dimension()
static coord normal_ebm_z (Point point, vector h) {
  return (coord){nodata, nodata, nodata};
}
#endif

static void half_column_ebm (Point point, scalar c, scalar cfg, vector h, vector cs, vector cfgs, int j)
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

    double S = c[], H = S, ci, cfgi, a;

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

    cfgi = cfg[];
    for (int i = 1; i <= 4; i++) {
      ci = i <= 2 ? c[i*j] : cs.x[(i - 2)*j];

      if (ci == nodata) {
        if (fabs(cfgi) == 2.)
          ci = (cfgi == 2. ? 1. : 0.);
        else
          break;
      }

      cfgi = i <= 2 ? cfg[i*j] : cfgs.x[(i - 2)*j];

      if (cfgi == 1.)
        break;

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
      
      else if (S == ci && modf(H, &a))
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

static void calculate_clh_ebm (Point point, scalar c, vector h)
{
  scalar cet = c.cet;
  vector mcl = c.mcl;

  coord ncl;
  foreach_dimension()
    ncl.x = mcl.x[];
  double alpha = plane_alpha (cet[], ncl);

  foreach_dimension()
    if (fabs(ncl.x) > 0. && fabs(alpha/ncl.x) <= 5.5)
      h.x[] = alpha/ncl.x + (ncl.x < 0.)*HSHIFT;
}

static void column_propagation_ebm (vector h)
{
  int change = 1;
  while (change) {
    change = 0;
    foreach (serial) // not compatible with OpenMP
      for (int i = -2; i <= 2; i++)
        foreach_dimension()
          if (fabs(height(h.x[i])) <= 3.5 &&
              fabs(height(h.x[i]) + i) < fabs(height(h.x[]))) {
              h.x[] = h.x[i] + i;
              change = 1;
          }
  }
}

#if !TREE
trace
void heights_ebm (scalar c, scalar cs, vector h, scalar cet)
{
  foreach()
    foreach_dimension()
      h.x[] = nodata;

  /**
  We need a 9-points-high stencil (rather than the default
  5-points). To do this we store in *s* the volume fraction field *c*
  shifted by 2 grid points in the respective directions. We make sure
  that this field uses the same boundary conditions as *c*. */

  scalar cfg = c.cfg;
  vector s[], sfg[];
  foreach_dimension()
    for (int i = 0; i < nboundary; i++) {
      s.x.boundary[i] = c.boundary[i];
      sfg.x.boundary[i] = cfg.boundary[i];
    }

  /**
  To compute the height function, we sum the volume fractions in a
  (half-)column starting at the current cell. We start by integrating
  downward (*j = -1*) and then integrate upward (*j = 1*). */
  
  for (int j = -1; j <= 1; j += 2) {

    /**
    We first build the shifted (by $\pm 2$) volume fraction field in each
    direction. */
    
    foreach()
      foreach_dimension() {
        s.x[] = cet[2*j];
        sfg.x[] = cfg[2*j];
      }

    /**
    We sum the half-column, downward or upward. We (exceptionally)
    need to allow for stencil overflow. */

    foreach (overflow)
      if (cfg[] != 1. && cet[] > 0. && cet[] < 1.)
        half_column_ebm (point, cet, cfg, h, s, sfg, j);
  }

  foreach()
    if (cfg[] == 1.)
      calculate_clh_ebm (point, c, h);

  /**
  Finally we "propagate" values along columns. */

  column_propagation_ebm (h);
}

#endif // TREE