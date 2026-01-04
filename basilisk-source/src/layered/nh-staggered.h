/**
# The vertically-staggered non-hydrostatic solver

See section 3.6.1 of [Popinet (2019)](/Bibliography#popinet2019). Note
that this uses the vertical momentum *qz* rather than vertical
velocity *w* and will not be compatible with [hydro.h](). */

#define NH 1

#include "poisson.h"

scalar * qzl = NULL, * phil = NULL;

mgstats mgp;
double * alpha = NULL, * gammap = NULL;

event defaults (i = 0)
{
  non_hydro = true;
  
  assert (qzl == NULL && phil == NULL &&
	  alpha == NULL && beta == NULL && gammap == NULL);
  assert (nl > 0);
  alpha = malloc (nl*sizeof(double));
  gammap = malloc (nl*sizeof(double));
  beta = malloc (nl*sizeof(double));  
  for (int l = 0; l < nl; l++) {
    scalar qz = new scalar;
    scalar phi = new scalar;
    qzl = list_append (qzl, qz);
    phil = list_append (phil, phi);
    beta[l] = 1./nl;
  }
  reset (qzl, 0.);
  reset (phil, 0.);

  if (!tracers)
    tracers = calloc (nl, sizeof(scalar *));

  int l = 0;
  for (scalar qz in qzl) {
    tracers[l] = list_append (tracers[l], qz);
    l++;
  }

  for (int l = 0; l < nl; l++)
    alpha[l] = gammap[l] = 1.;

#if 1
#if 0
  if (nl == 1) {
    alpha[0] = 3.;
    beta[0] = 1.;
  }
  else if (nl == 2) {
    // see section \pade2 of ale2.tm
    alpha[0] = 200./69., alpha[1] = 168070./109503.;
    beta[0] = 20./69., beta[1] = 49./69.;
  }
#if 0
  else if (nl == 3) {
    // Padé approximant: see dispersion3.tm
    alpha[0] = 147./50.;
    alpha[1] = 21884919./16356250.;
    alpha[2] = 2624395780083./1792301911300.;
    beta[0] = 7./50.;
    beta[1] = 19663./65425.;
    beta[2] = 14641./26170.;
  }
#else
  else if (nl == 3) {
    // Gnuplot fit: see dispersion3.tm
    alpha[0] = 3.483931265467430;
    alpha[1] = 1.916127288169944;
    alpha[2] = 1.710830494276016;
    beta[0] = 0.0569717048168878;
    beta[1] = 0.2367448039534642;
    beta[2] = 0.7062834912296478;
  }
#endif
  else if (nl == 4) {
    alpha[0] = 2.962285714285714;
    alpha[1] = 1.350414737345379;
    alpha[2] = 1.185189693944078;
    alpha[3] = 1.472922288188395;
    beta[0] = 0.0822857142857143;
    beta[1] = 0.1806217280518879;
    beta[2] = 0.2639407055948642;
    beta[3] = 0.4731518520675335;
  }
  else if (nl == 5) {
    alpha[0] = 2.974434611602751;
    alpha[1] = 1.364709853835999;
    alpha[2] = 1.176733678121684;
    alpha[3] = 1.143228438965934;
    alpha[4] = 1.490559988000380;
    beta[0] = 0.0540806293018682;
    beta[1] = 0.1209195157421197;
    beta[2] = 0.1769793771848123;
    beta[3] = 0.2310497798249295;
    beta[4] = 0.4169706979462702;
  }
  else {
    for (int l = 0; l < nl; l++)
      alpha[l] = 1., beta[l] = 1./nl;
    alpha[0] = 3.;
  }
#else
  if (nl == 1) {
    beta[0] = 1.;
    alpha[0] = 1.4;
  }
  else if (nl == 2) {
    // Gnuplot fit: see dispersion4.tm
#if 1
    beta[0] = 0.20052;
    alpha[0] = 1.57293;
    alpha[1] = 2.13295;
#else
    beta[0] = beta[1] = 0.5;
#if 0 // Keller matching
    alpha[0] = 3.569011576135351, alpha[1] = 17.93213572854291;
    gammap[0] = 1.314297124600639, gammap[1] = 8.11587147030185;
#else
    alpha[0] = 2.74711, alpha[1] = 12.6356;
    gammap[0] = 1., gammap[1] = 6.49546;
#endif
#endif
  }
  else if (nl == 3) {
#if 1 
    // Gnuplot fit: see dispersion4.tm
    beta[0] = 0.0640959;
    beta[1] = 0.231816;
    alpha[0] = 1.46978;
    alpha[1] = 2.04778;
    alpha[2] = 1.71597;
#else
    beta[0] = beta[1] = beta[2] = 1./3.;
#if 0 // Keller matching   
    alpha[0] = 618.4548736462094;
    alpha[1] = 4.337984496124031;
    alpha[2] = 9.782144219763449;
    gammap[0] = 33.44576523031203;
    gammap[1] = 4.741971207087486;
    gammap[2] = 6.894753476611884;
#else
    alpha[0] = 0.222769;
    alpha[1] = 3.27677;
    alpha[2] = 1.;
    gammap[0] = 0.0351922;
    gammap[1] = 0.355296;
    gammap[2] = 1.;
#endif
#endif
  }
  else if (nl == 4) {
    // Gnuplot fit: see dispersion4.tm
    beta[0]           = 0.0230779;
    beta[1]           = 0.0813287;
    beta[2]           = 0.252027;
    alpha[0]          = 1.4426;
    alpha[1]          = 2.03836;
    alpha[2]          = 1.69601;
    alpha[3]          = 1.6165;
  }
  else if (nl == 5) {
    // Gnuplot fit: see dispersion4.tm
    beta[0]           = 0.0170975;
    beta[1]           = 0.0468568;
    beta[2]           = 0.119315;
    beta[3]           = 0.260945;
    alpha[0]          = 1.34777;
    alpha[1]          = 1.76632;
    alpha[2]          = 1.56509;
    alpha[3]          = 1.4479;
    alpha[4]          = 1.52395;  
  }
#endif
  beta[nl-1] = 1.;
  for (int l = 0; l < nl - 1; l++)
    beta[nl-1] -= beta[l];
#endif
}

#define PHIT 6.

void correct_qz (double dt, (const) scalar phis)
{
  foreach() {
    scalar phi, qz, h;
    int l = 0;
    for (phi,qz,h in phil,qzl,hl) {
      double phit;
      if (l == nl - 1) {
	double phib, hb;
	if (l == 0) {
	  scalar s = hl[0];
	  phib = phi[], hb = s[];
	}
	else {
	  scalar s = phil[l-1];
	  phib = s[];
	  s = hl[l-1];
	  hb = s[];
	}
#if 0
	phit = ((2.*sq(h[])*phib - phi[]*sq(hb)
		 - 5.*h[]*phi[]*hb - 6.*sq(h[])*phi[])/
		(sq(hb) + 3.*h[]*hb + 2.*sq(h[])));
#else
	phit = (8.*phis[] + phib - PHIT*phi[])/3.; // top BC (Dirichlet)
#endif
      }
      else {
	scalar s = phil[l+1];
	phit = s[];
      }
      qz[] -= dt*alpha[l]*(phit - phi[]);
      l++;
    }
  }
}

#if 1
event viscous_term (i++)
{
  if (nu > 0.) {
#if 0
    correct_qz (dt);
#endif
    // fixme: ugly hack
    scalar lb = lambda_b, d = dut, u = u_b;
    lambda_b = dut = u_b = zeroc;
    foreach()
      // fixme: BCs should be different from those of horizontal velocity
      vertical_viscosity (point, hl, (scalar *) qzl, dt);
    lambda_b = lb; dut = d; u_b = u;
#if 0
    correct_qz (- dt);
#endif
  }
}
#endif

static void coeffs1 (Point point, scalar * hl,
		     double * a, double (* b)[2], double * c)
{
  int l = 0;
  double zr = zb[], zl = zb[-1];
  for (scalar h in hl) {
    b[l][0] = - (zr - zl - 2.*h[-1])/(2.*Delta);
    b[l][1] = - (zr - zl + 2.*h[])/(2.*Delta);
    if (l == nl - 1) {
      a[l] = (zr + h[] - zl - h[-1])/3./(2.*Delta);
      b[l][0] -= PHIT*a[l];
      b[l][1] -= PHIT*a[l];
      c[l] = 0.;
    }
    else {
      a[l] = 0.;
      c[l] = (zr + h[] - zl - h[-1])/(2.*Delta);
    }
    l++, zr += h[], zl += h[-1];
  }
}

// same as coeffs1 but without the [\phi\partial_x z] term
// (in blue in control.tm)
static void coeffs2 (Point point, scalar * hl,
		     double * a, double (* b)[2], double * c)
{
  int l = 0;
  for (scalar h in hl) {
    b[l][0] = h[-1]/Delta;
    b[l][1] = - h[]/Delta;
    if (l == nl - 1) {
      a[l] = 0.;
      b[l][0] -= PHIT*a[l];
      b[l][1] -= PHIT*a[l];
      c[l] = 0.;
    }
    else {
      a[l] = 0.;
      c[l] = 0.;
    }
    l++;
  }
}

static void matrix (Point point, scalar * phil, scalar * rhsl,
		    double * a, double * b, double * c, double * d)
{
  double al[nl], bl[nl][2], cl[nl];
  double ar[nl], br[nl][2], cr[nl];
#if 1
  coeffs1 (point, hl, al, bl, cl);
  coeffs1 (neighborp(1), hl, ar, br, cr);
#else
  coeffs2 (point, hl, al, bl, cl);
  coeffs2 (neighborp(1), hl, ar, br, cr);
#endif
  
  int l = 0;
  scalar phi, rhs, h;
  double zr = zb[1], zl = zb[-1];
  for (phi,rhs,h in phil,rhsl,hl) {
    d[l] = rhs[];
    
    c[l] = alpha[l];
    if (l == 0)
      a[l] = 0.; // bottom BC (Neumann)
    else {
      scalar hm = hl[l-1];
      a[l] = alpha[l-1]*h[]/hm[];
    }
    b[l] = - a[l] - c[l];
    if (l == nl - 1) {
      // top BC (Dirichlet)
      if (l == 0)
	b[l] -= (PHIT - 1.)*c[l]/3.;
      else {
	b[l] -= PHIT*c[l]/3.;
	a[l] += c[l]/3.;
      }
    }

#if 0 // no effect
    d[l] += 
      h[]*(h[-1]*phi[-1] - h[1]*phi[1])*(zr + h[1] - zl - h[-1])/(4.*sq(Delta));
    if (l > 0) {
      scalar hb = hl[l-1], phib = phil[l-1];
      d[l] -= h[]*(hb[-1]*phib[-1] - hb[1]*phib[1])*(zr - zl)/(4.*sq(Delta));
    }
#endif

    
    a[l] -= gammap[l]*h[]*(ar[l] - al[l])/Delta;
    b[l] -= gammap[l]*h[]*(br[l][0] - bl[l][1])/Delta;
    c[l] -= gammap[l]*h[]*(cr[l] - cl[l])/Delta;
    scalar phib = l > 0 ? phil[l-1] : phi;
    scalar phit = l < nl - 1 ? phil[l+1] : phi;
    d[l] += gammap[l]*h[]*
      ((ar[l]*phib[1] + br[l][1]*phi[1] + cr[l]*phit[1]) -
       (al[l]*phib[-1] + bl[l][0]*phi[-1] + cl[l]*phit[-1]))/Delta;
    l++, zr += h[1], zl += h[-1];
  }
}

static void matrix1 (Point point, scalar * phil, scalar * rhsl,
		     double * a, double * b, double * c, double * d)
{
  int l = 0;
  scalar phi, rhs, h;
  double zl = zb[-1], zc = zb[], zr = zb[1];
  for (phi,rhs,h in phil,rhsl,hl) {
    d[l] = rhs[];
#if 0 // h\partial_x(h\partial_x\phi), case (bbb)
    foreach_dimension()
      d[l] -= gammap[l]*h[]*((h[] + h[-1])*phi[-1] +
			     (h[] + h[1])*phi[1])/(2.*sq(Delta));
#else // h\partial_x\partial_x(h\phi), case (aaa)
    foreach_dimension()
      d[l] -= gammap[l]*h[]*(h[-1]*phi[-1] + h[1]*phi[1])/sq(Delta);
#endif
    
    c[l] = alpha[l];
    if (l == 0)
      a[l] = 0.; // bottom BC (Neumann)
    else
      a[l] = alpha[l-1];
#if 0 // \partial_x(h\partial_x\phi), case (bbb)
    // fixme: 1D only
    b[l] = - a[l] - c[l] - gammap[l]*h[]*(h[-1] + 2.*h[] + h[1])/(2.*sq(Delta));
#else // h\partial_x\partial_x(h\phi), case (aaa)
    b[l] = - a[l] - c[l] - gammap[l]*2.*dimension*sq(h[]/Delta);
#endif
    if (l == nl - 1) {
      // top BC (Dirichlet)
      if (l == 0)
	b[l] -= (PHIT - 1.)*c[l]/3.;
      else {
	b[l] -= PHIT*c[l]/3.;
	a[l] += c[l]/3.;
      }
    }
#if 0
    else {      
      /**
      gammap[l]*h[]*(((h[] + h[1])*(phi[1] - phi[]) -
		      (h[] + h[-1])*(phi[] - phi[-1]))/(2.*sq(Delta))
		     - ((zr + h[1] - zc - h[])*(s[1] + s[]) - 
                        (zr - zc)*(phi[] + phi[-1])
			-
			(zr + h[1] - zc - h[])*(s[1] + s[]) +
                        (zr - zc)*(phi[] + phi[-1]))/(2.*Delta));
      */
      b[l] += gammap[l]*h[]*(zr + h[1]/2. - 2.*(zc + h[]/2.) + zl + h[-1]/2.)/
	Delta;
      c[l] -= gammap[l]*h[]*(zr + h[1]/2. - 2.*(zc + h[]/2.) + zl + h[-1]/2.)/
	Delta;
      scalar s = phil[l+1];
      d[l] += gammap[l]*h[]*((zr + h[1]/2. - zc - h[]/2.)*(s[1] - phi[1]) -
			     (zc + h[]/2. - zl - h[-1]/2.)*(s[-1] - phi[-1]))/
	Delta;
    }
#else
#endif
    l++, zl += h[-1], zc += h[], zr += h[1];
  }
}

trace
static void relax_nh (scalar * phil, scalar * rhsl, int lev, void * data)
{
  foreach_level_or_leaf (lev) {
    double a[nl], b[nl], c[nl], d[nl]; // the tridiagonal system
    matrix (point, phil, rhsl, a, b, c, d);
    // Thomas algorithm
    for (int l = 1; l < nl; l++) {
      b[l] -= a[l]*c[l-1]/b[l-1];
      d[l] -= a[l]*d[l-1]/b[l-1];
    }
    scalar phi = phil[nl-1];
    phi[] = a[nl-1] = d[nl-1]/b[nl-1];
    for (int l = nl - 2; l >= 0; l--) {
      phi = phil[l];
      phi[] = a[l] = (d[l] - c[l]*a[l+1])/b[l];
    }
  }
}

static double residual_nh (scalar * phil, scalar * rhsl,
			   scalar * resl, void * data)
{
  (const) scalar phis = *((scalar *) data);
  double maxres = 0.;
  /* "naive" discretisation (only 1st order on trees) */
  foreach (reduction(max:maxres)) {
    double a[nl], b[nl], c[nl], d[nl]; // the tridiagonal system
    matrix (point, phil, rhsl, a, b, c, d);
    d[nl-1] -= 8.*alpha[nl-1]*phis[]/3.;
    int l = 0;
    scalar phi, res;
    for (phi,res in phil,resl) {
      res[] = d[l] - b[l]*phi[];
      if (l > 0) {
	scalar phim = phil[l-1];
	res[] -= a[l]*phim[];
      }
      if (l < nl - 1) {
	scalar phip = phil[l+1];
	res[] -= c[l]*phip[];
      }
      if (fabs (res[]) > maxres)
	maxres = fabs (res[]);
      l++;
    }
  }
  return maxres;
}

trace
static void relax_nh1 (scalar * phil, scalar * rhsl, int lev, void * data)
{
  foreach_level_or_leaf (lev) {
    double a[nl], b[nl], c[nl], d[nl]; // the tridiagonal system
    int l = 0;
    scalar phi, rhs, h;
    for (phi,rhs,h in phil,rhsl,hl) {
      d[l] = rhs[];
      foreach_dimension()
	d[l] -= gammap[l]*h[]*((h[] + h[-1])*phi[-1] +
			       (h[] + h[1])*phi[1])/(2.*sq(Delta));
      c[l] = alpha[l];
      if (l == 0)
	a[l] = 0.; // bottom BC (Neumann)
      else
	a[l] = alpha[l-1];
      // fixme: 1D only
      b[l] = - a[l] - c[l] - gammap[l]*h[]*(h[-1] + 2.*h[] + h[1])/sq(Delta);
      if (l == nl - 1) {
	// top BC (Dirichlet)
	if (l == 0)
	  b[l] -= (PHIT - 1.)*c[l]/3.;
	else {
	  b[l] -= PHIT*c[l]/3.;
	  a[l] += c[l]/3.;
	}	  
      }
      l++;
    }
    // Thomas algorithm
    for (int l = 1; l < nl; l++) {
      b[l] -= a[l]*c[l-1]/b[l-1];
      d[l] -= a[l]*d[l-1]/b[l-1];
    }
    phi = phil[nl-1];
    phi[] = a[nl-1] = d[nl-1]/b[nl-1];
    for (int l = nl - 2; l >= 0; l--) {
      phi = phil[l];
      phi[] = a[l] = (d[l] - c[l]*a[l+1])/b[l];
    }
  }
}

static double residual_nh2 (scalar * phil, scalar * rhsl,
			    scalar * resl, void * data)
{
  double maxres = 0.;
  foreach (reduction(max:maxres)) {
    scalar phi, rhs, res, h;
    int l = 0;
    double zl = zb[-1], z = zb[], zr = zb[1];
    for (phi,rhs,res,h in phil,rhsl,resl,hl) {
      res[] = rhs[];
      double phit, phib, hb, alphab;
      if (l == 0) {
#if 1
	phib = phi[]; // bottom BC (Neumann)
#else
	phib = phi[] - h[]*(zb[1] - zb[-1])*(phi[1] - phi[-1])/sq(2.*Delta);
#endif
	hb = h[];
	alphab = alpha[l];
      }
      else {
	scalar phi = phil[l-1], h = hl[l-1];
	phib = phi[];
	hb = h[];
	alphab = alpha[l-1];
      }
      if (l == nl - 1)
#if 0
	phit = ((2.*sq(h[])*phib - phi[]*sq(hb)
		 - 5.*h[]*phi[]*hb - 6.*sq(h[])*phi[])/
		(sq(hb) + 3.*h[]*hb + 2.*sq(h[])));
#else
	phit = (phib - PHIT*phi[])/3.; // top BC (Dirichlet)
#endif
      else {
	scalar s = phil[l+1];
	phit = s[];
#if 1 // not much effect
	res[] += h[]*((s[1] - s[-1])*(zr + h[1] - zl - h[-1]) -
		      (phi[1] - phi[-1])*(zr - zl))/sq(2.*Delta);
#endif
#if 0
	res[] += h[]*((zr + h[1]/2. - z - h[]/2.)*
		      (s[1] + s[] - phi[1] - phi[]) -
		      (z + h[]/2. - zl - h[-1]/2.)*
		      (s[] + s[-1] - phi[] - phi[-1]))/(2.*Delta);
#endif
      }
      res[] -= alpha[l]*(phit - phi[]) - alphab*h[]*(phi[] - phib)/hb;
#if 1 // almost no effect
      foreach_dimension()
	res[] += gammap[l]*h[]*((h[] + h[-1])*(phi[] - phi[-1]) -
				(h[] + h[1])*(phi[1] - phi[]))/(2.*sq(Delta));
#elif 0
      foreach_dimension()
	res[] += gammap[l]*sq(h[])*(face_gradient_x (phi, 0) -
				    face_gradient_x (phi, 1))/Delta;
#else // almost no effect
      foreach_dimension()
	res[] += gammap[l]*h[]*(2.*h[]*phi[] - h[-1]*phi[-1] - h[1]*phi[1])/
	sq(Delta);
#endif
      if (fabs (res[]) > maxres)
	maxres = fabs (res[]);
      l++, zl += h[-1], z += h[], zr += h[1];
    }
  }
  return maxres;
}

static double residual_nh3 (scalar * phil, scalar * rhsl,
			    scalar * resl, void * data)
{
  double maxres = 0.;
  foreach (reduction(max:maxres)) {
    double a[nl], b[nl], c[nl], d[nl]; // the tridiagonal system
    matrix (point, phil, rhsl, a, b, c, d);
    
    scalar phi, rhs, res, h;
    int l = 0;
    for (phi,rhs,res,h in phil,rhsl,resl,hl) {
      res[] = rhs[];
      double phit, phib, hb, alphab;
      if (l == 0) {
	phib = phi[]; // bottom BC (Neumann)
	hb = h[];
	alphab = alpha[l];
      }
      else {
	scalar phi = phil[l-1], h = hl[l-1];
	phib = phi[];
	hb = h[];
	alphab = alpha[l-1];
      }
      if (l == nl - 1)
	phit = (phib - PHIT*phi[])/3.; // top BC (Dirichlet)
      else {
	scalar s = phil[l+1];
	phit = s[];
      }
      res[] -= alpha[l]*(phit - phi[]) - alphab*h[]*(phi[] - phib)/hb;
      foreach_dimension()
	res[] += gammap[l]*h[]*((h[] + h[-1])*(phi[] - phi[-1]) -
				(h[] + h[1])*(phi[1] - phi[]))/(2.*sq(Delta));

#if 0
      double res2 = d[l] - b[l]*phi[];
      if (l > 0) {
	scalar phim = phil[l-1];
	res2 -= a[l]*phim[];
      }
      if (l < nl - 1) {
	scalar phip = phil[l+1];
	res2 -= c[l]*phip[];
      }

      fprintf (stderr, "res: %d %g %g\n", l, res[], res2);
#endif
      
      if (fabs (res[]) > maxres)
	maxres = fabs (res[]);
      l++;
    }
  }
  return maxres;
}

scalar phit = zeroc;

event pressure (i++)
{
  scalar * rhsl = list_clone (phil);
  foreach() {
    scalar rhs, h, qz;
    face vector uf;
    int l = 0;
    double zl = zb[-1], zr = zb[1];
    vector q;
    for (rhs,h,qz,uf,q in rhsl,hl,qzl,ufl,ql) {
      rhs[] = qz[] - h[]*(uf.x[] + uf.x[1])*(zr + h[1] - zl - h[-1])/(4.*Delta);
      if (l > 0) {
	scalar qzm = qzl[l-1], hm = hl[l-1];
	rhs[] -= h[]*qzm[]/hm[];
	vector ufm = ufl[l-1];
	rhs[] += h[]*(ufm.x[] + ufm.x[1])*(zr - zl)/(4.*Delta);
      }
      foreach_dimension()
	rhs[] += h[]*((h[] + h[1])*uf.x[1] - (h[] + h[-1])*uf.x[])/(2.*Delta);
      rhs[] /= dt;
      l++, zl += h[-1], zr += h[1];
    }
  }
  
  restriction (hl);
  mgp = mg_solve (phil, rhsl, residual_nh, relax_nh, &phit,
		  nrelax = 4, res = NULL, minlevel = 1, tolerance = TOLERANCE);
  delete (rhsl), free (rhsl);
  
  foreach_face() {
    double ac[nl], bc[nl][2], cc[nl];
    coeffs1 (point, hl, ac, bc, cc);

    scalar phi, h;
    vector uf, a;
    int l = 0;
    for (phi,uf,a,h in phil,ufl,al,hl) {
      scalar phit = l < nl - 1 ? phil[l+1] : phi;
      scalar phib = l > 0 ? phil[l-1] : phi;
      double ax = - fm.x[]*(ac[l]*phib[-1] + bc[l][0]*phi[-1] + cc[l]*phit[-1] +
			    ac[l]*phib[] + bc[l][1]*phi[] + cc[l]*phit[])/
	((h[] + h[-1])/2.);
      uf.x[] -= dt*ax;
      a.x[] -= ax;
      l++;
    }
  }

  correct_qz (dt, phit);
}

event cleanup (i = end, last) {
  free (alpha), alpha = NULL;
  free (gammap), gammap = NULL;
  delete (qzl), free (qzl), qzl = NULL;
  delete (phil), free (phil), phil = NULL;
}
