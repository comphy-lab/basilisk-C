/**
# Shallow-water flux computation at boundaries

In order to impose a given flow rate $Q_b$ through a boundary $b$, when
solving the [Saint-Venant equations](saint-venant.h), we need to
compute the water level $\eta_b$ corresponding to this flow rate. The
flow rate $Q$ is a complicated, non-linear function of the water
level, the bathymetry $z_b$ and the velocity normal to the boundary
$u_n$.

We first define a function which, given $\eta_b$ and the current
topography and velocity fields (defined by the [Saint-Venant
solver](saint-venant.h)), returns the flow rate through boundary *b*.

The extent of the boundary can also be limited to points for which
`limit[] = value`. */

struct Eta_b {
  // compulsory arguments
  double Q_b;
  bid b;
  // optional arguments
  scalar limit;
  double value;
  double prec; // precision (default 0.1%)
};

static double bflux (struct Eta_b p, double eta_b)
{
  double Q = 0.;
  scalar limit = p.limit;
  foreach_boundary (p.b, reduction(+:Q))
    if (limit.i < 0 || limit[] == p.value) {
      scalar u_n = ig ? u.x : u.y; // normal velocity component  
      double sign = - (ig + jg), ub = u_n[]*sign;
      double hn = max (eta_b - zb[], 0.), hp = max(h[],0.);
      if (hp > dry || hn > dry) {
	double fh, fu, dtmax;
	kurganov (hn, hp, ub, ub, 0., &fh, &fu, &dtmax);
	Q += fh*Delta; // fixme: metric
      }
    }
  return Q;
}

/**
To find the water level $\eta_b$ corresponding to the flow rate $Q_b$
we want to impose, we need to invert the function above i.e. find
$\eta_b$ such that
$$
Q[z_b,u_n](\eta_b) = Q_b
$$
We do this using the [false position
method](http://en.wikipedia.org/wiki/False_position_method). */

static double falsepos (struct Eta_b p,
			double binf, double qinf,
			double bsup, double qsup)
{
  int n = 0;
  double newb, newq;
  qinf -= p.Q_b;
  qsup -= p.Q_b;
  do {
    newb = (binf*qsup - bsup*qinf)/(qsup - qinf);
    newq = bflux (p, newb) - p.Q_b;
    if (newq > 0.)
      bsup = newb, qsup = newq;
    else
      binf = newb, qinf = newq;
    n++;
  } while (fabs(newq/p.Q_b) > p.prec && n < 100);

  if (n >= 100)
    fprintf (stderr, "src/discharge.h:%d: warning: eta_b(): convergence not reached\n", LINENO);
  
  return newb;
}

/**
## User interface

Given a target flux $Q_b$ and a boundary $b$ (optionally limited to
points for which `limit[] = value`), this function returns the
corresponding water level $\eta_b$. */

double eta_b (double Q_b, bid b,
	      scalar limit = {-1}, double value = 0, double prec = 0.001)
{
  double zmin = HUGE, etas = 0., hs = 0.;
  foreach_boundary (b, reduction(+:etas) reduction(+:hs) reduction(min:zmin))
    if (limit.i < 0 || limit[] == value) {
      if (zb[] < zmin)
	zmin = zb[];
      etas += Delta*h[]*eta[];
      hs += Delta*h[];
    }

  if (Q_b <= 0.)
    return zmin - 1.;

  /**
  We try to find good bounds on the solution. */
  
  double etasup = hs > 0. ? etas/hs : zmin;
  struct Eta_b p = { Q_b, b, limit, value, prec };
  double Qsup = bflux (p, etasup), etainf = zmin, Qinf = 0.;
  double h0 = etasup - zmin;
  if (h0 < dry)
    h0 = 1.;
  int n = 0;
  while (Qsup < p.Q_b && n++ < 100) {
    etainf = etasup, Qinf = Qsup;
    etasup += h0;
    Qsup = bflux (p, etasup);
  }
  if (n >= 100)
    fprintf (stderr, "src/discharge.h:%d: warning: eta_b() not converged\n", LINENO);
  return falsepos (p, etainf, Qinf, etasup, Qsup);
}
