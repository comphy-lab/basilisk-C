/**
# A solver for the Green-Naghdi equations

<div class="message">
Note that the [multilayer solver](layered/hydro.h) provides the same
functionality and should be prefered for most applications.</div>

The Green-Naghdi equations (also known as the Serre or fully nonlinear
Boussinesq equations) can be seen as an extension of the [Saint-Venant
equations](saint-venant.h) to the next order $O(\mu^2)$ in term of the
*shallowness parameter*
$$
\mu = \frac{h_0^2}{L^2}
$$
with $h_0$ the typical depth and $L$ the typical horizontal scale. In
contrast to the Saint-Venant equations the Green-Naghdi equations have
*dispersive* wave solutions. 

A more detailed description of the context and numerical scheme
implemented here is given in [Popinet,
2015](/src/references.bib#popinet2015).

The solver is built by adding a source term to the momentum equation
of the [Saint-Venant solver](saint-venant.h). Following [Bonneton et
al, 2011](/src/references.bib#bonneton2011), this source term can be
written
$$
\partial_t \left( hu \right) + \ldots = h \left( \frac{g}{\alpha_d}
   \nabla \eta - D \right)
$$
where $D$ verifies
$$
\alpha_d h\mathcal{T} \left( D \right) + hD = b
$$
and
$$
b = \left[ \frac{g}{\alpha_d} \nabla \eta +\mathcal{Q}_1 \left( u \right)
\right]
$$
With $\mathcal{T}$ a linear operator to be defined below, as well as
$\mathcal{Q}_1 \left( u \right)$.

Before including the Saint-Venant solver, we need to overload the
default *update* function of the predictor-corrector scheme in order
to add our source term. */

#include "predictor-corrector.h"

static double update_green_naghdi (scalar * current, scalar * updates,
				   double dtmax);

event defaults (i = 0)
  update = update_green_naghdi;

#include "saint-venant.h"

/**
The linear system can be inverted with the multigrid Poisson
solver. We declare *D* as a global variable so that it can be re-used
as initial guess for the Poisson solution. The solver statistics will
be stored in *mgD*. The *breaking* parameter defines the slope above
which dispersive terms are turned off. The $\alpha_d$ parameter
controls the optimisation of the dispersion relation (see [Bonneton et
al, 2011](/src/references.bib#bonneton2011)). */

#include "poisson.h"

vector D[];
mgstats mgD;
double breaking = 1., alpha_d = 1.153;

/**
We first define some useful macros, following the notations in
[Bonneton et al, 2011](/src/references.bib#bonneton2011).

Simple centered-difference approximations of the first and second
derivatives of a field. */

#define dx(s)  ((s[1,0] - s[-1,0])/(2.*Delta))
#define dy(s)  ((s[0,1] - s[0,-1])/(2.*Delta))
#define d2x(s) ((s[1,0] + s[-1,0] - 2.*s[])/sq(Delta))
#define d2y(s) ((s[0,1] + s[0,-1] - 2.*s[])/sq(Delta))
#define d2xy(s) ((s[1,1] - s[1,-1] - s[-1,1] + s[-1,-1])/sq(2.*Delta))

/**
The definitions of the $\mathcal{R}_1$ and $\mathcal{R}_2$ operators
$$
\begin{array}{lll}
  \mathcal{R}_1 \left[ h, z_b \right] w & = & - \frac{1}{3 h} \nabla \left(
  h^3 w \right) - \frac{h}{2} w \nabla z_b\\
  & = & - h \left[ \frac{h^{}}{3} \nabla w + w \left( \nabla h + \frac{1}{2}
  \nabla z_b \right)\right]\\
  \mathcal{R}_2 \left[ h, z_b \right] w & = & \frac{1}{2 h} \nabla \left( h^2
  w \right) + w \nabla z_b\\
  & = & \frac{h}{2} \nabla w + w \nabla \left( z_b + h \right)
\end{array}
$$ */

#define R1(h,zb,w) (-h[]*(h[]/3.*dx(w) + w[]*(dx(h) + dx(zb)/2.)))
#define R2(h,zb,w) (h[]/2.*dx(w) + w[]*(dx(zb) + dx(h)))

/**
## Inversion of the linear system

To invert the linear system which defines $D$, we need to write
discretised versions of the residual and relaxation operators. The
functions take generic multilevel parameters and a user-defined
structure which contains solution-specific parameters, in our case a
list of the $h$, $zb$ and *wet* fields. */

static double residual_GN (scalar * a, scalar * r, scalar * resl, void * data)
{
  /**
  We first recover all the parameters from the generic pointers and
  rename them according to our notations. */
  
  scalar * list = (scalar *) data;
  scalar h = list[0], zb = list[1], wet = list[2];
  vector D = vector(a[0]), b = vector(r[0]), res = vector(resl[0]);

  /**
  The general form for $\mathcal{T}$ is
  $$
  \mathcal{T} \left[ h, z_b \right] w = \mathcal{R}_1 \left[ h, z_b
  \right]  \left( \nabla \cdot w \right) + \mathcal{R}_2 \left[ h, z_b
  \right]  \left( \nabla z_b \cdot w \right)
  $$
  which gives the linear problem
  $$
  \alpha_d h \mathcal{T} \left( D \right) + hD = b
  $$
  $$
  \alpha_d h \mathcal{R}_1 \left( \nabla \cdot D \right) + \alpha_d h
  \mathcal{R}_2 \left( \nabla z_b \cdot D \right) + hD = b
  $$
  $$
  - \alpha_d h^2  \left[ \frac{h}{3} \nabla \left( \nabla \cdot D \right)
  + \left( \nabla \cdot D \right)  \left( \nabla h + \frac{1}{2} \nabla z_b
  \right)\right] + \alpha_d h \left[ \frac{h}{2} \nabla \left( \nabla z_b \cdot D
  \right) + \left( \nabla z_b \cdot D \right) \nabla \eta \right] + hD = b
  $$
  Expanding the operators for the $x$-component gives
  $$
  \begin{aligned}
  - \frac{\alpha_d}{3} \partial_x \left( h^3 \partial_x D_x \right) + h \left[
  \alpha_d \left( \partial_x \eta \partial_x z_b + \frac{h}{2} \partial^2_x z_b
  \right) + 1 \right] D_x + & \\
  \alpha_d h \left[ \left( \frac{h}{2} \partial^2_{xy} z_b + \partial_x \eta
  \partial_y z_b \right) D_y + \frac{h}{2} \partial_y z_b \partial_x D_y -
  \frac{h^2}{3} \partial^2_{xy} D_y - h \partial_y D_y  \left( \partial_x h +
  \frac{1}{2} \partial_x z_b \right) \right] & = b_x
  \end{aligned}
  $$
  The $y$-component is obtained by rotation of the indices. */
 
  double maxres = 0.;
  foreach (reduction(max:maxres))
    foreach_dimension() {
      if (wet[-1] == 1 && wet[] == 1 && wet[1] == 1) {
	double hc = h[], dxh = dx(h), dxzb = dx(zb), dxeta = dxzb + dxh;
	double hl3 = (hc + h[-1])/2.; hl3 = cube(hl3);
	double hr3 = (hc + h[1])/2.;  hr3 = cube(hr3);
	
	/**
	Finally we translate the formula above to its discrete
	version. */
	
	res.x[] = b.x[] -
	  (-alpha_d/3.*(hr3*D.x[1] + hl3*D.x[-1] - 
			(hr3 + hl3)*D.x[])/sq(Delta) +
	   hc*(alpha_d*(dxeta*dxzb + hc/2.*d2x(zb)) + 1.)*D.x[] +
	   alpha_d*hc*((hc/2.*d2xy(zb) + dxeta*dy(zb))*D.y[] + 
		       hc/2.*dy(zb)*dx(D.y) - sq(hc)/3.*d2xy(D.y)
		       - hc*dy(D.y)*(dxh + dxzb/2.)));
      
	/**
	The function also need to return the maximum residual. */
	
	if (fabs (res.x[]) > maxres)
	  maxres = fabs (res.x[]);
      }
      else
	res.x[] = 0.;
    }
  
  /**
  The maximum residual is normalised by gravity i.e. the tolerance is
  the relative acceleration times the depth. */

  return maxres/G;
}

/**
The relaxation function is built by copying and pasting the
residual implementation above and inverting manually for $D_x$. */

static void relax_GN (scalar * a, scalar * r, int l, void * data)
{
  scalar * list = (scalar *) data;
  scalar h = list[0], zb = list[1], wet = list[2];
  vector D = vector(a[0]), b = vector(r[0]);
  
#if GAUSS_SEIDEL || _GPU
  for (int parity = 0; parity < 2; parity++)
    foreach_level_or_leaf (l)
      if (level == 0 || ((point.i + parity) % 2) != (point.j % 2))
#else
  foreach_level_or_leaf (l)
#endif
    foreach_dimension() {
      if (h[] > dry && wet[-1] == 1 && wet[] == 1 && wet[1] == 1) {
	double hc = h[], dxh = dx(h), dxzb = dx(zb), dxeta = dxzb + dxh;
	double hl3 = (hc + h[-1])/2.; hl3 = cube(hl3);
	double hr3 = (hc + h[1])/2.;  hr3 = cube(hr3);
	D.x[] = (b.x[] -
		 (-alpha_d/3.*(hr3*D.x[1] + hl3*D.x[-1])/sq(Delta) +
		  alpha_d*hc*((hc/2.*d2xy(zb) + dxeta*dy(zb))*D.y[] + 
			      hc/2.*dy(zb)*dx(D.y) - sq(hc)/3.*d2xy(D.y)
			      - hc*dy(D.y)*(dxh + dxzb/2.))))/
	  (alpha_d*(hr3 + hl3)/(3.*sq(Delta)) + 
	   hc*(alpha_d*(dxeta*dxzb + hc/2.*d2x(zb)) + 1.));
      }
      else
	D.x[] = 0.;
    }
}

/**
## Source term computation

To add the source term to the Saint-Venant system we overload the
default *update* function with this one. The function takes
a list of the current evolving scalar fields and fills the
corresponding *updates* with the source terms. */

static double update_green_naghdi (scalar * current, scalar * updates,
				   double dtmax)
{
  double dt = update_saint_venant (current, updates, dtmax);
  scalar h = current[0];
  vector u = vector(current[1]);

  /**
  We first compute the right-hand-side $b$. The general form for the
  $\mathcal{Q}_1$ operator is (eq. (9) of Bonneton et al, 2011).
  $$
  \mathcal{Q}_1 \left[ h, z_b \right] \left( u \right) = - 2 \mathcal{R}_1 
  \left( c \right) + \mathcal{R}_2  \left( d \right)
  $$
  with
  $$
  \begin{aligned}
  c & = \partial_1 u \cdot \partial_2 u^{\perp} + \left( \nabla \cdot u
  \right)^2\\
  & = - \partial_x u_x \partial_y u_y + \partial_x u_y \partial_y u_x +
  \left( \partial_x u_x + \partial_y u_y \right)^2\\
  d & = u \cdot \left( u \cdot \nabla \right) \nabla z_b\\
  & = u_x  \left( u_x \partial_x^2 z_b + u_y \partial_{xy}^2 z_b \right) +
  u_y  \left( u_y \partial_y^2 z_b + u_x \partial_{xy}^2 z_b \right)\\
  & = u^2_x \partial_x^2 z_b + u^2_y \partial_y^2 z_b + 2 u_x u_y
  \partial_{xy}^2 z_b
  \end{aligned}
  $$
  Note that we declare field *c* and *d* in a new scope, so that the
  corresponding memory is freed after we have computed $b$. */

  vector b[];
  {
    scalar c[], d[];
    foreach() {
      double dxux = dx(u.x), dyuy = dy(u.y);
      c[] = - dxux*dyuy + dx(u.y)*dy(u.x) + sq(dxux + dyuy);
      d[] = sq(u.x[])*d2x(zb) + sq(u.y[])*d2y(zb) + 2.*u.x[]*u.y[]*d2xy(zb);
    }

    /**
    We can now compute $b$
    $$
    b = \left[ \frac{g}{\alpha_d} \nabla \eta +\mathcal{Q}_1 \left( u \right)
    \right]
    $$ */

    foreach()
      foreach_dimension()
        b.x[] = h[]*(G/alpha_d*dx(eta) - 2.*R1(h,zb,c) + R2(h,zb,d));
  }

  /**
  We declare a new field which will track whether cells are completely
  wet. This is necessary to turn off dispersive terms close to the
  shore so that lake-at-rest balance is maintained. */

  scalar wet[];
  foreach()
    wet[] = h[] > dry;

  /**
  Finally we solve the linear problem with the multigrid solver using
  the relaxation and residual functions defined above. We need to
  restrict $h$ and $z_b$ as their values will be required for
  relaxation on coarse grids. */

  scalar * list = {h, zb, wet};
  restriction (list);
  mgD = mg_solve ((scalar *){D}, (scalar *){b},
		  residual_GN, relax_GN, list, mgD.nrelax);

  /**
  We can then compute the updates for $hu$. */

  vector dhu = vector(updates[1]);

  /**
  We only apply the Green-Naghdi source term when the slope of the
  free surface is smaller than *breaking*. The idea is to turn off
  dispersion in areas where the wave is "breaking" (i.e. in
  hydraulic jumps). We also turn off dispersive terms close to shore
  so that lake-at-rest balance is maintained. */

  foreach()
    if (fabs(dx(eta)) < breaking && fabs(dy(eta)) < breaking)
      foreach_dimension()
	if (wet[-1] == 1 && wet[] == 1 && wet[1] == 1)
	  dhu.x[] += h[]*(G/alpha_d*dx(eta) - D.x[]);

  return dt;
}
