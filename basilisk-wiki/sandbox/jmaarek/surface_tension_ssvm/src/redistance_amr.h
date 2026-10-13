/**
# Redistancing of a distance field — AMR-consistent variant

Drop-in replacement for src/redistance.h. Two changes:

1. Interfacial cells are detected from the *true leaf neighbours*
   (same level, finer children, or coarser parent-level leaf) instead of
   from `phi0[±1]`, which at a resolution boundary is a restricted
   average or a bilinear prolongation and can have the wrong sign.
   This makes the interfacial flag consistent on both sides of a
   coarse/fine face.

2. A sign-preservation limiter on every forward-Euler stage, so a cell
   can never be carried through zero in one stage. Because the RK3
   scheme is a convex combination of Euler stages (Shu & Osher), the
   full step then preserves sign(phi0) as well. (The order = 2 midpoint
   scheme is not SSP, so only order = 3 gets this guarantee.)
*/

static inline double minmod3 (double a, double b)
{
  if (a == b || a*b <= 0.)
    return 0.;
  return fabs(a) < fabs(b) ? a : b;
}

#define BGHOSTS 2

/**
Distance from the cell centre to the zero crossing of `phi0` in
direction `j` (±1), using the actual leaf neighbour across the face.
Returns -1 if there is no crossing. */

foreach_dimension()
static double crossing_x (Point point, scalar phi0, int j, double phixxmin)
{
#if TREE
  if (is_refined (neighbor(j))) {
    /* Finer leaves across the face: test the 2^(d-1) children that
       touch this face. Their centres are 3Δ/4 away along x. */
    double d = -1.;
    int ci = j > 0 ? 2 : -1;
#if dimension == 1
    double pc = fine(phi0, ci);
    if (phi0[]*pc < 0.)
      d = 0.75*Delta*phi0[]/(phi0[] - pc);
#elif dimension == 2
    for (int k = 0; k <= 1; k++) {
      double pc = fine(phi0, ci, k);
      if (phi0[]*pc < 0.) {
	double dk = 0.75*Delta*phi0[]/(phi0[] - pc);
	if (d < 0. || dk < d) d = dk;
      }
    }
#else
    for (int k = 0; k <= 1; k++)
      for (int l = 0; l <= 1; l++) {
	double pc = fine(phi0, ci, k, l);
	if (phi0[]*pc < 0.) {
	  double dk = 0.75*Delta*phi0[]/(phi0[] - pc);
	  if (d < 0. || dk < d) d = dk;
	}
      }
#endif
    return d;
  }
  if (is_prolongation (neighbor(j))) {
    /* Coarser leaf across the face (2:1 balance guarantees it is the
       parent-level neighbour). Its centre is 3Δ/2 away along x. */
    double pc = coarse(phi0, j);
    return phi0[]*pc < 0. ? 1.5*Delta*phi0[]/(phi0[] - pc) : -1.;
  }
#endif // TREE

  /* Same-level neighbour: original subcell estimate (Min, 2010). */
  if (phi0[]*phi0[j] >= 0.)
    return -1.;
  double dx = Delta;
  double phixx = minmod3 (phi0[2*j] + phi0[] - 2.*phi0[j],
			  phi0[1] + phi0[-1] - 2.*phi0[]);
  if (fabs(phixx) > phixxmin) {
    double D = sq(phixx/2. - phi0[] - phi0[j]) - 4.*phi0[]*phi0[j];
    dx *= 1/2. + (phi0[] - phi0[j] - sign2(phi0[] - phi0[j])*sqrt(D))/phixx;
  }
  else
    dx *= phi0[]/(phi0[] - phi0[j]);
  return dx;
}

static
void dphidt (scalar phi, scalar dphi, scalar phi0, const double cfl,
	     const double phixxmin)
{
  foreach() {
    double dt = cfl*Delta;

    /* ENO one-sided differences (unchanged). */
    coord gra[2];
    for (int i = 0, j = 1; i < 2; j = 1 - 2*++i)
      foreach_dimension() {
	double s1 = (phi[2*j] + phi[] - 2.*phi[j])/Delta;
	double s2 = (phi[1] + phi[-1] - 2.*phi[])/Delta;
	gra[i].x = j*((phi[j] - phi[])/Delta - minmod3(s1, s2)/2.);
      }

    /* Interfacial cells and subcell fix, from true leaf neighbours. */
    dphi[] = - sign2(phi0[]);
    bool interfacial = false;
    double size = HUGE;
    foreach_dimension()
      for (int i = 0, j = 1; i < 2; j = 1 - 2*++i) {
	double dx = crossing_x (point, phi0, j, phixxmin);
	if (dx >= 0.) {
	  interfacial = true;
	  if (dx != 0.) {
	    double sxx1 = phi[2 - 4*i] + phi[] - 2.*phi[1 - 2*i];
	    double sxx2 = phi[1] + phi[-1] - 2.*phi[];
	    gra[i].x = (2*i - 1)*(phi[]/dx +
				  dx*minmod3(sxx1, sxx2)/(2.*sq(Delta)));
	  }
	  else
	    gra[i].x = 0.;
	  size = min(size, dx);
	}
      }
    if (interfacial)
      dphi[] *= min(dt, fabs(size)/2.)/dt;

    /* Godunov Hamiltonian (unchanged). */
    double H_G = 0;
    if (phi0[] > 0)
      foreach_dimension() {
	double a = min(0., gra[0].x);
	double b = max(0., gra[1].x);
	H_G += max(sq(a), sq(b));
      }
    else
      foreach_dimension() {
	double a = max(0., gra[0].x);
	double b = min(0., gra[1].x);
	H_G += max(sq(a), sq(b));
      }
    dphi[] *= sqrt(H_G) - 1.;

    /* Sign-preservation limiter: one Euler stage of length dt may
       remove at most half of |phi| on the phi0 side. */
    if (dphi[]*phi0[] < 0.) {
      double lim = max(phi[]*sign2(phi0[]), 0.)/(2.*dt);
      if (fabs(dphi[]) > lim)
	dphi[] = sign2(dphi[])*lim;
    }
  }
}

trace
int redistance (scalar phi,
		int imax = 1,
		double cfl = 0.5,
		int order = 3,
		double eps = 1e-6,
		double band = HUGE,
		scalar resf = {-1},
		const double phixxmin = 1./HUGE)
{
  scalar phi0[];
  foreach()
    phi0[] = phi[];

  for (int i = 1; i <= imax; i++) {
    double maxres = 0.;
    if (order == 2) {
      scalar tmp1[];
      dphidt (phi, tmp1, phi0, cfl/2., phixxmin);
      foreach()
	tmp1[] = phi[] + Delta*cfl/2.*tmp1[];
      scalar tmp2[];
      dphidt (tmp1, tmp2, phi0, cfl, phixxmin);
      foreach (reduction(max:maxres)) {
	double res = tmp2[];
	if (resf.i >= 0) resf[] = res;
	if (fabs (res) > maxres && fabs(phi0[]) < band*Delta) maxres = fabs (res);
	phi[] += Delta*cfl*tmp2[];
      }
    }
    else {
      scalar tmp1[];
      dphidt (phi, tmp1, phi0, cfl, phixxmin);
      foreach()
	tmp1[] = phi[] + Delta*cfl*tmp1[];
      scalar tmp2[];
      dphidt (tmp1, tmp2, phi0, cfl, phixxmin);
      foreach()
	tmp1[] = (3.*phi[] + tmp1[] + Delta*cfl*tmp2[])/4.;
      dphidt (tmp1, tmp2, phi0, cfl, phixxmin);
      foreach (reduction(max:maxres)) {
	double res = 2./3.*((phi[] - tmp1[])/(Delta*cfl) - tmp2[]);
	if (resf.i >= 0) resf[] = res;
	if (fabs (res) > maxres && fabs(phi0[]) < band*Delta) maxres = fabs (res);
	phi[] = (phi[] + 2.*(tmp1[] + Delta*cfl*tmp2[]))/3.;
      }
    }
    if (maxres < eps)
      return i;
  }
  return imax;
}