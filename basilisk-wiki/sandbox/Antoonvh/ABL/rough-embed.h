/**
# Embedded law of the wall

The so-called "law of the wall" is a major experimental result on
wall-bounded turbulent flows. It may be exploited to include the main
effect of a rough surface in simulations without computing the
turbulent transfer directly. 

The variation of the "averaged" value of a velocity component $u_i$ at
a distance $d$ from a wall, inside a wall-bounded turbulent flow may
exibit the logaritmic law of a no-slip wall:

$$U(d) = \frac{{u}_{*}}{\kappa}
\mathrm{ln}\left(\frac{d}{z_0}\right),$$

with $\kappa$ the Von Karman constant and $z_0$ the roughness length
($d \gg z_0$). We can compute $u_*$ from,

$$u_* = \frac{\kappa}{\mathrm{ln}\left(\frac{d}{z_0}\right)}U(d).$$

This relates to the flux ($F_\phi$) of a dummy variable $\phi$:

$$F_\phi = \frac{u_* \kappa}{\mathrm{ln}(\frac{d}{z_{0,\phi}})}
\left( \phi (d) - \phi(0) \right),$$

where $z_{0,\phi}$ is an $\phi$-specific roughness length (by default,
it is equal to $z_0$. These scalars must be added to a list called
`rough_list`, their surface value and roughness length are defined as
scalar-field attributes.
**/
scalar * rough_list = NULL;

attribute {
  double rl;
  double s_val;
}

#ifndef Z_0 //z0 for momentum -> default for scalars 
#define Z_0 (L0/(10*(1 << grid->maxdepth)))
#endif

event defaults (i = 0) {
  foreach_dimension() {
    scalar s = u.x;
    rough_list = list_add (rough_list, u.x);
    s.rl = Z_0;
    s.s_val = 0; //no-slip;
  }
}

// A user interface for adding scalar fields
struct rough_s {
  scalar s;
  double s_val;
  double z0;
};

void add_rough_scalar (struct rough_s p) {
  scalar s = p.s;
  rough_list = list_add (rough_list, s);
  s.rl = p.z0 ? p.z0 : Z_0;
  s.s_val = p.s_val; 
}
/**
The code in this header file implements the closure, and applies the
resulting tendencies. It is supposed to be used with [embedded
boundaries](/src/embed.h) and the
[centered](/src/navier-stokes/centered.h) solver.

## Implementation 
 */
#if (dimension == 2)
#define VOL (sq(Delta))
#elif (dimension == 3)
#define VOL (cube(Delta))
#endif
double karman = 0.3679; //= exp(-1)
/**
Next, the "work horse": a function that computes $F_\phi$ for all
fields in `rough_list`. Much of the used method is not very different
from using `dirichlet()` conditions with [embed.h](/src/embed.h) and
we follow the wisdom presented there. Note that the interpolator mixes
the finite-volume representation of the scalar-field values with the
vaulues at cell-centered locations. The written implementation
considers a surface whose normal component is mostly in the x
direction. The function returns $u_*$.
 */
foreach_dimension() { 
  static inline double u_tau_x (Point point, scalar cs, coord n, coord p,
				scalar * list, double * flx) {// List and fluxes 
    foreach_dimension()
      n.x = - n.x;     // inward not outward
    double d[2] = {nodata, nodata}, v[2][list_len(list)], u_star = 0; //Distances and values;
    bool defined = true;
    foreach_dimension()
      if (defined && !fs.x[(n.x > 0.)])
	defined = false;
    if (defined) {
      for (int l = 0; l < 2; l++) {
	int i = (l + 1)*sign(n.x);
	d[l] = (i - p.x)/n.x;
	double y1 = p.y + d[l]*n.y;
	int j = y1 > 0.5 ? 1 : y1 < -0.5 ? -1 : 0;
	y1 -= j;
#if dimension == 2
	if (fs.x[i + (i < 0),j] && fs.y[i,j] && fs.y[i,j+1] &&
	    cs[i,j-1] && cs[i,j] && cs[i,j+1]) {
	  int c = 0;
	  for (scalar s in list) 
	    v[l][c++] = quadratic (y1, (s[i,j-1]), (s[i,j]), (s[i,j+1]));
	}
#else // dimension == 3
	double z = p.z + d[l]*n.z;
	int k = z > 0.5 ? 1 : z < -0.5 ? -1 : 0;
	z -= k;
	bool defined = fs.x[i + (i < 0),j,k];
	for (int m = -1; m <= 1 && defined; m++)
	  if (!fs.y[i,j,k+m] || !fs.y[i,j+1,k+m] ||
	      !fs.z[i,j+m,k] || !fs.z[i,j+m,k+1] ||
	      !cs[i,j+m,k-1] || !cs[i,j+m,k] || !cs[i,j+m,k+1])
	    defined = false;
	if (defined) {
	  // bi-quadratic interpolation
	  int c = 0;
	  for (scalar s in list) {
	    v[l][c++] =
	      quadratic (z,
			 quadratic (y1,
				    (s[i,j-1,k-1]), (s[i,j,k-1]), (s[i,j+1,k-1])),
			 quadratic (y1,
				    (s[i,j-1,k]),   (s[i,j,k]),   (s[i,j+1,k])),
			 quadratic (y1,
				    (s[i,j-1,k+1]), (s[i,j,k+1]), (s[i,j+1,k+1])));
	  }
	}
#endif // dimension == 3
	else {
	  d[l] = nodata; 
	  break;
	}
      }
      if (d[0] != nodata && d[1] != nodata) { // Two points
	double U[2] = {0., 0.};
	for (int l = 0; l < 2; l++) {
	  int dim = 0;
	  foreach_dimension() {
	    U[l] += sq(v[l][dim]);
	    dim++;
	  }
	}
	for (int l = 0; l < 2; l++) 
	  U[l] = U[l] > 0 ? sqrt(U[l]) : 0;
	u_star = karman*(U[0]/(log(Delta*d[0]/Z_0)) +
			 U[1]/(log(Delta*d[1]/Z_0)))/2.; //mean
	int c = 0;
	for (scalar s in list) {
	  flx[c] = -u_star*((v[0][c] - s.s_val)/(log(Delta*d[0]/s.rl)) +
			    (v[1][c] - s.s_val)/(log(Delta*d[1]/s.rl)))/2.; //mean
	  c++;
	}
      }
      else if (d[0] != nodata || d[1] != nodata) { //Single point interpolation
	int l;
	if (d[0] != nodata)
	  l = 0;
	else
	  l = 1;
	double U = 0.;
	int dim = 0;
	foreach_dimension() {
	  U += sq(v[l][dim]);
	  dim++;
	}
	U = U > 0 ? sqrt(U) : 0; 
	u_star = karman*U/(log(Delta*d[l]/Z_0));
	int c = 0;
	for (scalar s in list) {
	  flx[c] = -u_star*(v[l][c] - s.s_val)/(log(Delta*d[l]/s.rl));
	  c++;
	}
      }
    } else { //else Poor estimate based on the embedded cell itself :(
      d[0] = max (10*Z_0, fabs (p.x/n.x));
      double U = 0;
      foreach_dimension()
	U += sq(u.x[]);
      U = U > 0 ? sqrt(U) : 0.;
      u_star = karman*U/(log(Delta*d[0]/Z_0));
      int c = 0;
      for (scalar s in list) 
	flx[c++] = -u_star*(s[] - s.s_val)/(log(Delta*d[0]/s.rl));
    }
    return u_star;
  }
}
/**
   The function `rough_wall()` calles `u_tau_i()` depending on the dominant
   orientation of the boundary. The function returns $u_*$.
*/
static inline double u_tau (Point point, scalar cs, coord n, coord p,
			    scalar * list, double * flx) {
#if dimension == 2
  foreach_dimension()
    if (fabs(n.x) >= fabs(n.y))
      return u_tau_x (point, cs, n, p, list, flx);
#else // dimension == 3
  if (fabs(n.x) >= fabs(n.y)) {
    if (fabs(n.x) >= fabs(n.z))
      return u_tau_x (point, cs, n, p, list, flx);
  }
  else if (fabs(n.y) >= fabs(n.z))
    return u_tau_y (point, cs, n, p, list, flx);
  return u_tau_z (point, cs, n, p, list, flx);
#endif // dimension == 3
  return nodata;
}
/**
   The drag term is included just before the diffusive step. The total
   change ($m^{1+\mathrm{dimension}}s^{-1}= u_{*}$) for each boundary
   fragment is stored in the `flxdt` vector.
*/

event viscous_term (i++) {
  scalar * dsl = list_clone (rough_list); //change for each field
  foreach() { //Boundary cells.
    for (scalar ds in dsl)
      ds[] = 0.;
    double scs = 0;
    if (cs[] > 0. && cs[] < 1.) {
      foreach_neighbor(1)
	scs += cs[];
      double sf = 1/Delta; //Scale factor
      foreach_dimension()
	sf *= Delta;
      coord n = facet_normal (point, cs, fs), p;
      double alpha = plane_alpha (cs[], n);
      double area = sf*plane_area_center (n, alpha, &p);
      double flxs[list_len(rough_list)];
      u_tau (point, cs, n , p, rough_list, flxs);
      int c = 0;
      for (scalar ds in dsl)
	ds[] = dt*area*flxs[c++]/scs;
    }
  }
  boundary (dsl);
  /**
The drag is added to the cells *near* the boundary (smearing). For
additional stability and physical consistency, the new velocity must
be between the current value and zero...

Note that `VOL = dv()/cs[]`. 
  */
  foreach() {
    scalar s, ds;
    for (s, ds in rough_list, dsl) {
      double se = 0;
      foreach_neighbor(1)
	se += ds[];
      if (cs[] > 0) {
	double new = s[] > s.s_val ? max(s[] + se/VOL, s.s_val) : min (s[] + se/VOL , s.s_val);
	s[] = s[] > s.s_val ? min (new, s[]) : max (new, s[]);
      }
    }
  }
  boundary (rough_list);
}

/**
## References

Juan C. Bergmann, 1998, *A physical interpretation of von Karman’s
constant based on asymptotic considerations—A new value*. Journal of
the Atmospheric
Sciences. [doi](https://doi.org/10.1175/1520-0469(1998)055<3403:APIOVK>2.0.CO;2)
 */
