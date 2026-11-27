/**
# Staggered non-hydrostatic extension

This adds the non-hydrostatic terms of the [vertically-Lagrangian
multilayer solver for free-surface flows](hydro.h) described in
[Popinet, 2020](/Bibliography#popinet2020). The corresponding system
of equations is
$$
\begin{aligned}
  \partial_t h_k + \mathbf{{\nabla}} \cdot \left( h \mathbf{u} \right)_k & =
  0,\\
  \partial_t \left( h \mathbf{u} \right)_k + \mathbf{{\nabla}} \cdot \left(
  h \mathbf{u}  \mathbf{u} \right)_k & = - gh_k  \mathbf{{\nabla}} (\eta)
  {\color{blue} - \mathbf{{\nabla}} (h \phi)_k + \left[ \phi
      \mathbf{{\nabla}} z \right]_k},\\
  {\color{blue} \partial_t (hw)_k + \mathbf{{\nabla}} \cdot
  \left( hw \mathbf{u} \right)_k} & {\color{blue} = - [\phi]_k,}\\
  {\color{blue} \mathbf{{\nabla}} \cdot \left( h \mathbf{u} \right)_k +
  \left[ w - \mathbf{u} \cdot \mathbf{{\nabla}} z \right]_k} &
  {\color{blue} = 0}.
\end{aligned}
$$

## Staggered fields

The vertical velocity $w_l$ is stored at the center of layer $l$.
The pressure $\phi_l$ is stored at its lower interface. Thus
$\phi_{l+1}$ is the pressure at the upper interface of the same
layer, and the atmospheric condition at the free surface is imposed by
the convention
$$
\phi_{nl} = 0.
$$

The vertical pressure acceleration and its dual centered difference are
$$
(G_v\phi)_l = \frac{\phi_{l+1}-\phi_l}{h_l},
\qquad
(D_v G_v\phi)_l = (G_v\phi)_l - (G_v\phi)_{l-1},
$$
with $(G_v\phi)_{-1}=0$. After the pressure solve, the vertical
velocity is corrected over the full timestep according to
$$
w_l^{n+1} = w_l^\star - \Delta t\,(G_v\phi)_l.
$$

## Mid-step continuity

Let $\theta=\theta_H$, and let stars denote the predictor values
before the non-hydrostatic pressure correction. On each vertical face,
define
$$
q_{\phi,l} = -pg_l,
\qquad
q_{\eta,l} = hf_l\,a_{baro}(\eta),
$$
where $pg_l$ is the pressure acceleration returned by `hpg`.
The midpoint velocities and fluxes have the pressure responses
$$
\begin{aligned}
w^{n+\theta} &= w^\star - \theta\Delta t\,G_v\phi,\\
Q^{n+\theta} &= Q^\star - \theta\Delta t\,q_\phi
                    + \theta^2\Delta t\,q_\eta .
\end{aligned}
$$
The last term contains $\theta^2$ because `implicit.h` stores the
$\theta$-weighted end-of-step flux in `hu`.

The discrete incompressibility condition is imposed at
$n+\theta$:
$$
D_v w^{n+\theta} + D_h Q^{n+\theta} = 0.
$$
Here
$$
(D_v w)_l = w_l - w_{l-1}, \qquad w_{-1}=0,
$$
and $D_h$ is the dual horizontal divergence documented below.
Substitution of the pressure responses gives
$$
S - \theta\Delta t\,D_vG_v\phi
  - \theta\Delta t\,D_hq_\phi
  + \theta^2\Delta t\,D_hq_\eta = 0,
\qquad
S = D_vw^\star + D_hQ^\star .
$$
After division by $\theta\Delta t$, the pressure right-hand side
and local operator are therefore
$$
rhs_\phi = -\frac{S}{\theta\Delta t},
\qquad
L_\phi = -D_vG_v\phi - D_hq_\phi + \theta D_hq_\eta .
$$
The full-step correction is still applied to `w`, while the horizontal
flux correction is $\theta\Delta t\,pg$. For
$\theta=1/2$, the linear constant-coefficient time integration
is centered; $\theta=1$ gives the first-order backward-Euler
form.
*/

#define NH 1
#include "layered/implicit.h"

scalar w, phi;
mgstats mgp;
double breaking = HUGE;

/**
## Setup
*/

event defaults (i = 0)
{
  hydrostatic = false;
  assert (nl > 0);
  w = new scalar[nl];
  phi = new scalar[nl];

  reset ({w, phi}, 0.);
  if (!linearised)
    tracers = list_append (tracers, w);
}

/**
## Viscous term

Vertical diffusion is added to the vertical component of velocity. The
free-surface condition is $\partial_z w=0$ by default.
*/

(const) scalar dwt[] = 0.;

event viscous_term (i++)
{
  if (nu > 0.)
    foreach()
      vertical_diffusion (point, h, w, dt, nu, dwt[], 0., 0.);
}

/**
## Geometric coefficients

For face offset `i = 0`, the face lies between cells `-1` and `0`;
for `i = 1`, it lies between cells `0` and `1`. Denote the two cells
by $L$ and $R$, and the lower and upper interface elevations
of layer $l$ by $z_b$ and $z_t$. Define
$$
s_b = z_{b,R}-z_{b,L}, \qquad s_t = z_{t,R}-z_{t,L}.
$$
In the formulas below, these differences are limited as in
`hpg`: each $s$ is replaced by $\Delta x\,slope\_limited(s/\Delta x)$.

With $\bar h=(h_L+h_R)/2$, the geometric coefficients are
$$
\begin{aligned}
C_b &= \frac{-h_L+s_b}{2\bar h\,\Delta x}, &
C_t &= \frac{-h_L-s_t}{2\bar h\,\Delta x},\\
D_b &= \frac{ h_R+s_b}{2\bar h\,\Delta x}, &
D_t &= \frac{ h_R-s_t}{2\bar h\,\Delta x}.
\end{aligned}
$$
These coefficients contain only interface geometry and use the
arithmetic normalization of `hpg`. The face thickness `hf` is applied
separately when forming a flux.
*/


/*
macro to work within a foreach_dimension().
*/

#define nh_face_coefficients(point, i, Cb, Ct, Db, Dt) do {           \
  double _nh_zL = zb[i-1];                                             \
  double _nh_zR = zb[i];                                               \
  foreach_layer() {                                                    \
    if (h[i-1] + h[i] > dry) {                                        \
      double _nh_sb =                                                 \
        Delta*slope_limited ((_nh_zR - _nh_zL)/Delta);                \
      double _nh_st =                                                 \
        Delta*slope_limited (((_nh_zR + h[i]) -                       \
                              (_nh_zL + h[i-1]))/Delta);              \
      double _nh_scale = 1./((h[i-1] + h[i])*Delta);                  \
                                                                      \
      Cb[point.l] = (- h[i-1] + _nh_sb)*_nh_scale;                    \
      Ct[point.l] = (- h[i-1] - _nh_st)*_nh_scale;                    \
      Db[point.l] = (  h[i]   + _nh_sb)*_nh_scale;                    \
      Dt[point.l] = (  h[i]   - _nh_st)*_nh_scale;                    \
      _nh_zL += h[i-1];                                               \
      _nh_zR += h[i];                                                 \
    }                                                                 \
    else {                                                            \
      /* The face is absent from the operator. */                \
      Cb[point.l] = 0.;                                               \
      Ct[point.l] = 0.;                                               \
      Db[point.l] = 0.;                                               \
      Dt[point.l] = 0.;                                               \
    }                                                                 \
  }                                                                   \
} while (0)


/**
On a one-dimensional Cartesian mesh, the face flux associated with the
pressure acceleration returned by `hpg` can be written
$$
q_{\phi,l} = hf_l\left(
 C_{b,l}\phi_{L,l} + D_{b,l}\phi_{R,l}
+C_{t,l}\phi_{L,l+1} + D_{t,l}\phi_{R,l+1}\right),
$$
with $\phi_{nl}=0$. The residual obtains the pressure flux by calling
`hpg`, so that its metric and wet-face treatment are identical to the
pressure correction. The relaxation uses this coefficient form to
assemble the current-column response without reevaluating `hpg`.


For left and right face flux arrays $q_L$ and $q_R$, the
horizontal contribution to layer-integrated incompressibility is
$$
\begin{aligned}
(D_hq)_l ={}& -C^R_{b,l}q_{R,l}-D^L_{b,l}q_{L,l}\\
&-C^R_{t,l-1}q_{R,l-1}-D^L_{t,l-1}q_{L,l-1},
\end{aligned}
$$
where the last two terms are omitted for $l=0$. The superscripts
indicate the face at which the coefficients are evaluated. This is the
dual of the horizontal pressure-gradient stencil: the same coefficients
are used in `hpg`, in this divergence and in the relaxation matrix.
*/


/**
## Residual

The combined horizontal flux is
$$
q = -q_\phi + \theta q_\eta ,
$$
so that the layer operator evaluated below is
$$
L_\phi = -D_vG_v\phi + D_hq,
\qquad res_\phi = rhs_\phi-L_\phi .
$$
With $g=\alpha_{eta}a_{baro}(\eta)$ and $\beta=1$ for a free
surface and $\beta=0$ for a rigid lid, the surface residual uses
$$
L_\eta = \beta\eta
 - \frac{g_R-g_L}{\Delta x\,cm}
 - \theta\Delta t^2\frac{\sum_l(q_{\phi,R,l}-q_{\phi,L,l})}
                         {\Delta x\,cm},
\qquad res_\eta = rhs_\eta-L_\eta .
$$
*/

static double residual_nh (scalar * pl, scalar * rhsl,
                    scalar * resl, void * data)
{
  scalar p = pl[0], eta = pl[1];
  scalar rhs = rhsl[0], surface_rhs = rhsl[1];
  scalar res = resl[0], surface_res = resl[1];

  face vector alpha = *((vector *)data);

  double maxres = 0.;

  foreach (reduction(max:maxres))  {

    foreach_layer() {
      double above = point.l < nl - 1 ? p[0,0,1] : 0.;
      double dv = (above - p[])/h[];

      if (point.l > 0)
        dv -= (p[] - p[0,0,-1])/h[0,0,-1];

      res[] = rhs[] + dv;
    }

    surface_res[] = surface_rhs[] - (rigid ? 0. : eta[]);

    foreach_dimension() {
      double CbL[nl], CtL[nl], DbL[nl], DtL[nl];
      double CbR[nl], CtR[nl], DbR[nl], DtR[nl];
      NOT_UNUSED(CbL); NOT_UNUSED(CtL); NOT_UNUSED(DbL); NOT_UNUSED(DtL);
      NOT_UNUSED(CbR); NOT_UNUSED(CtR); NOT_UNUSED(DbR); NOT_UNUSED(DtR);

      nh_face_coefficients (point, 0, CbL, CtL, DbL, DtL);
      nh_face_coefficients (point, 1, CbR, CtR, DbR, DtR);

      double axL = a_baro (eta, 0);
      double axR = a_baro (eta, 1);
      double divq = 0.;

      /*
      A face flux in layer l contributes to residual rows l and l + 1.
      Since res = rhs + dv - dh, these contributions have positive signs.
      */

      hpg (pg, p, 0, {
        double q = pg + theta_H*hf.x[]*axL;

        divq += pg;
        if (h[] > dry)
          res[] += DbL[point.l]*q;

        if (point.l + 1 < nl && h[0,0,1] > dry)
          (res[0,0,1] += DtL[point.l]*q);
      });

      hpg (pg, p, 1, {
        double q = pg + theta_H*hf.x[1]*axR;

        divq -= pg;
        if (h[] > dry)
          res[] += CbR[point.l]*q;

        if (point.l + 1 < nl && h[0,0,1] > dry)
          (res[0,0,1] += CtR[point.l]*q);
      });

      /*
      Accumulate the barotropic divergence and the pressure-flux
      contribution to the free-surface residual for this direction.
      */

      surface_res[] +=
        (alpha.x[1]*axR - alpha.x[]*axL)/(Delta*cm[])
        + theta_H*sq(dt)*divq/(Delta*cm[]);
    }

    for (int l = 0; l < nl; l++)
      maxres = max (maxres, fabs (res[0,0,l]));

    maxres = max (maxres, fabs (surface_res[]));
  }

  return maxres;
}


/**
## Tridiagonal coupled column relaxation

With neighboring columns fixed, the local coupled system has the form
$$
\begin{aligned}
b_\phi + Tp + c\eta &= rhs_\phi,\\
b_\eta + r^Tp + d\eta &= rhs_\eta ,
\end{aligned}
$$
where $T$ is the tridiagonal pressure-pressure matrix. One Thomas
factorization is used for the two right-hand sides
$$
Tp_0 = rhs_\phi-b_\phi, \qquad Tp_c=c.
$$
The local coupled solution is then
$$
p = p_0-\eta p_c,
\qquad
\eta = \frac{rhs_\eta-b_\eta-r^Tp_0}{d-r^Tp_c}.
$$
Only the current-column block is assembled in $T$. Contributions
from neighboring columns are contained in the baseline responses.
*/

static void nh_thomas_solve (const double * lower,
                    const double * diagonal,
                    const double * upper,
                    double * b, double * c)
{
  double d[nl];

  d[0] = diagonal[0];
  assert (fabs (d[0]) > 0.);

  for (int l = 1; l < nl; l++) {
    double multiplier = lower[l]/d[l-1];

    d[l] = diagonal[l] - multiplier*upper[l-1];
    assert (fabs (d[l]) > 0.);

    b[l] -= multiplier*b[l-1];
    c[l] -= multiplier*c[l-1];
  }

  b[nl-1] /= d[nl-1];
  c[nl-1] /= d[nl-1];

  for (int l = nl - 2; l >= 0; l--) {
    b[l] = (b[l] - upper[l]*b[l+1])/d[l];
    c[l] = (c[l] - upper[l]*c[l+1])/d[l];
  }
}

static void relax_nh (scalar * pl, scalar * rhsl, int lev,
                    void * data)
{
  scalar p = pl[0], eta = pl[1];
  scalar rhs = rhsl[0], surface_rhs = rhsl[1];

  face vector alpha = *((vector *)data);


foreach_level_or_leaf (lev) {
    double lower[nl], diagonal[nl], upper[nl];
    double base_phi[nl], p0[nl], pc[nl], rphi[nl];
    bool active[nl];

    for (int l = 0; l < nl; l++) {
      active[l] = h[0,0,l] > dry;
      p[0,0,l] = 0.;

      lower[l] = 0.;
      upper[l] = 0.;

      base_phi[l] = 0.;
      pc[l] = 0.;
      rphi[l] = 0.;

      if (active[l]) {
        diagonal[l] = 1./h[0,0,l];

        if (l + 1 < nl && h[0,0,l+1] > dry)
          upper[l] = - 1./h[0,0,l];

        if (l > 0 && h[0,0,l-1] > dry) {
          diagonal[l] += 1./h[0,0,l-1];
          lower[l] = - 1./h[0,0,l-1];
        }
      }
      else {
        /* Dry pressure unknown: p_l = 0. */
        diagonal[l] = 1.;
      }
    }

    eta[] = 0.;

    double base_eta = 0.;
    double d_eta = rigid ? 0. : 1.;

    foreach_dimension() {
      double CbL[nl], CtL[nl], DbL[nl], DtL[nl];
      double CbR[nl], CtR[nl], DbR[nl], DtR[nl];
      NOT_UNUSED(CbL); NOT_UNUSED(CtL); NOT_UNUSED(DbL); NOT_UNUSED(DtL);
      NOT_UNUSED(CbR); NOT_UNUSED(CtR); NOT_UNUSED(DbR); NOT_UNUSED(DtR);

      nh_face_coefficients (point, 0, CbL, CtL, DbL, DtL);
      nh_face_coefficients (point, 1, CbR, CtR, DbR, DtR);

      /*
      Evaluate the base and unit-eta barotropic responses. The current
      column pressure remains zero while neighboring columns are frozen.
      */

      double axL0 = a_baro (eta, 0);
      double axR0 = a_baro (eta, 1);

      eta[] = 1.;
      double axL1 = a_baro (eta, 0);
      double axR1 = a_baro (eta, 1);
      eta[] = 0.;

      double divq_base = 0.;

      for (int l = 0; l < nl; l++) {
        double aboveL = l + 1 < nl ? p[-1,0,l+1] : 0.;
        double aboveR = l + 1 < nl ? p[1,0,l+1]  : 0.;

        /* Pressure fluxes containing only neighboring columns. */

        double qpL = gmetric (0)*hf.x[0,0,l]*
          (CbL[l]*p[-1,0,l] + CtL[l]*aboveL);
        double qpR = gmetric (1)*hf.x[1,0,l]*
          (DbR[l]*p[1,0,l] + DtR[l]*aboveR);

        divq_base += qpR - qpL;

        double qL = theta_H*hf.x[0,0,l]*axL0 - qpL;
        double qR = theta_H*hf.x[1,0,l]*axR0 - qpR;

        double qeL = hf.x[0,0,l]*(axL1 - axL0);
        double qeR = hf.x[1,0,l]*(axR1 - axR0);

        /*
        A flux in layer l contributes to rows l and l + 1. This is
        the dual-divergence stencil without the *_below temporaries.
        */
        if (active[l]) {
          base_phi[l] -= CbR[l]*qR + DbL[l]*qL;
          pc[l] -= theta_H*(CbR[l]*qeR + DbL[l]*qeL);
        }

        if (l + 1 < nl && active[l+1]) {
          base_phi[l+1] -= CtR[l]*qR + DtL[l]*qL;
          pc[l+1] -= theta_H*(CtR[l]*qeR + DtL[l]*qeL);
        }
      }

      double divg0 = (alpha.x[1]*axR0 - alpha.x[]*axL0)
        /(Delta*cm[]);
      double divg1 = (alpha.x[1]*axR1 - alpha.x[]*axL1)
        /(Delta*cm[]);

      divq_base /= Delta*cm[];

      base_eta -= divg0 + theta_H*sq(dt)*divq_base;
      d_eta -= divg1 - divg0;

      double surface_factor = theta_H*sq(dt)/(Delta*cm[]);
      double gmL = gmetric (0), gmR = gmetric (1);

      for (int l = 0; l < nl; l++) {
        if (active[l])
          diagonal[l] += gmetric (1)*hf.x[1,0,l]*sq (CbR[l])
                      +  gmetric (0)*hf.x[0,0,l]*sq (DbL[l]);

        if (l + 1 < nl && active[l+1])
           diagonal[l+1] += gmetric (1)*hf.x[1,0,l]*sq (CtR[l])
                    + gmetric (0)*hf.x[0,0,l]*sq (DtL[l]);

        if (l + 1 < nl && active[l] && active[l+1]) {
           double coupling = gmetric (1)*hf.x[1,0,l]*CbR[l]*CtR[l]
                    + gmetric (0)*hf.x[0,0,l]*DbL[l]*DtL[l];

          upper[l] += coupling;
          lower[l+1] += coupling;
        }

        if (h[-1,0,l] + h[0,0,l] > dry) {
          double f = surface_factor*gmL*hf.x[0,0,l];

          if (active[l])
            rphi[l] += f*DbL[l];

          if (l + 1 < nl && active[l+1])
            rphi[l+1] += f*DtL[l];
        }

        if (h[0,0,l] + h[1,0,l] > dry) {
          double f = surface_factor*gmR*hf.x[1,0,l];
          if (active[l])
            rphi[l] -= f*CbR[l];

          if (l + 1 < nl && active[l+1])
            rphi[l+1] -= f*CtR[l];
        }
      }
    }

    for (int l = 0; l < nl; l++)
      p0[l] = active[l] ? rhs[0,0,l] - base_phi[l] : 0.;

    nh_thomas_solve (lower, diagonal, upper, p0, pc);

    double r_p0 = 0., r_pc = 0.;

    for (int l = 0; l < nl; l++) {
      r_p0 += rphi[l]*p0[l];
      r_pc += rphi[l]*pc[l];
    }

    double denominator = d_eta - r_pc;
    assert (fabs (denominator) > 0.);

    double new_eta =
      (surface_rhs[] - base_eta - r_p0)/denominator;

    for (int l = 0; l < nl; l++)
      p[0,0,l] = active[l] ? p0[l] - new_eta*pc[l] : 0.;
    eta[] = new_eta;
  }
}

/**
## Coupled pressure solution

The predictor defect is assembled using $D_v$ and $D_h$.
After the coupled solve, the non-hydrostatic pressure gradient is added
once to the face acceleration and to the midpoint flux. The summed
gradient `su` is retained for the second-pass free-surface correction
in `implicit.h`.
*/

event pressure (i++)
{

  scalar rhs = new scalar[nl];
  double h1 = 0., v1 = 0.;

  foreach (reduction(+:h1) reduction(+:v1)) {

    /*
    Accumulate the mid-step defect S = D_v w + D_h Q.
    Start with the vertical contribution.
    */

    foreach_layer() {
      if (h[] > dry) {
        double below =
          point.l > 0 && h[0,0,-1] > dry ? w[0,0,-1] : 0.;
        rhs[] = - (w[] - below);
      }
      else {
        rhs[] = 0.;
        w[] = 0.;
        phi[] = 0.;
      }
    }

    /* Add the horizontal divergence direction by direction. */

    foreach_dimension() {
      double CbL[nl], CtL[nl], DbL[nl], DtL[nl];
      double CbR[nl], CtR[nl], DbR[nl], DtR[nl];
      NOT_UNUSED(CbL); NOT_UNUSED(CtL); NOT_UNUSED(DbL); NOT_UNUSED(DtL);
      NOT_UNUSED(CbR); NOT_UNUSED(CtR); NOT_UNUSED(DbR); NOT_UNUSED(DtR);

      nh_face_coefficients (point, 0, CbL, CtL, DbL, DtL);
      nh_face_coefficients (point, 1, CbR, CtR, DbR, DtR);

      for (int l = 0; l < nl; l++) {
        if (h[0,0,l] > dry) {
          double dh = - CbR[l]*hu.x[1,0,l]
                      - DbL[l]*hu.x[0,0,l];

          if (l > 0)
            dh -= CtR[l-1]*hu.x[1,0,l-1]
                + DtL[l-1]*hu.x[0,0,l-1];

          rhs[0,0,l] -= dh;
        }
      }
    }

    /* rhs_phi = -S/(theta_H*dt). */

    foreach_layer() {
      if (h[] > dry)
        rhs[] /= theta_H*dt;
      else {
        rhs[] = 0.;
        w[] = 0.;
        phi[] = 0.;
      }
      h1 += dv()*h[];
      v1 += dv();
    }

  }

  /* All geometry used on coarse multigrid levels must be restricted. */

  restriction ({zb, h, hf, alpha_eta});

  double hbar = h1/v1;
  double scale = hbar/sq(dt);
  double tolerance = TOLERANCE*scale;

  scalar res;
  if (res_eta.i >= 0)
    res = new scalar[nl];

  mgp = mg_solve ({phi, eta_r}, {rhs, rhs_eta},
                  residual_nh, relax_nh, &alpha_eta,
                  res = res_eta.i >= 0 ? (scalar *){res, res_eta} : NULL,
                  nrelax = 4, minlevel = 1,
                  tolerance = tolerance);

  delete ({rhs});
  if (res_eta.i >= 0)
    delete ({res});

  face vector su[];
  foreach_face () {
    su.x[] = 0.;
    hpg (pg, phi, 0, {
      ha.x[] += pg;
      hu.x[] += theta_H*dt*pg;
      su.x[] -= pg;
    });
  }

  foreach() {
    double wmax = HUGE;
    if (breaking < HUGE) {
      wmax = 0.;
      foreach_layer()
	wmax += h[];
      wmax = wmax > 0. ? breaking*sqrt(G*wmax) : 0.;
    }
    foreach_layer() {
      if (h[] > dry) {
        bool above_active =
          point.l < nl - 1 && h[0,0,1] > dry;
        double above = above_active ? phi[0,0,1] : 0.;

        w[] -= dt*(above - phi[])/h[];
      }
      else {
        w[] = 0.;
        phi[] = 0.;
      }
      if (fabs(w[]) > wmax)
        w[] = (w[] > 0. ? 1. : -1.)*wmax;
    }
    foreach_dimension()
      rhs_eta[] += theta_H*sq(dt)*(su.x[1] - su.x[])/(Delta*cm[]);
  }
}

event cleanup (i = end, last) {
  delete ({w, phi});
}
