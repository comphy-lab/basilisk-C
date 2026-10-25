/**
# Anisotropic self-similar viscous operator for the 2D Cartesian Huppert heap.

<div class="message">
<div id="msg_logo"><img src="/img/warning.png"></div>
**Important note**

This file was generated in october 2026 by ChatGPT, but double-checked by the 
author of this Sandbox regarding the discretization of the self-similar viscous 
operator.
</div>

 This is based on Basilisk [`src/viscosity.h`](https://basilisk.fr/src/viscosity.h).  
 The multigrid interface and
 the relaxation/residual strategy are retained, but the stress coefficients
 are split because the transformed operator is:

$$ 
\boxed{ 
  \mathcal{D}_{\overline{u}_\xi} 
  =
  A_\xi \partial_\xi (2 \mu \partial_\xi \overline{u}_\xi) 
  + A_\eta \partial_\eta (\mu \partial_\eta \overline{u}_\xi) 
  + A_\eta \partial_\eta (\mu \partial_\xi \overline{u}_\eta) 
} 
$$

$$ 
\boxed{ 
  \mathcal{D}_{\overline{u}_\eta} 
  = 
  A_\eta \partial_\xi(\mu \partial_\xi \overline{u}_\eta) 
  + A_\eta \partial_\xi (\mu \partial_\eta \overline{u}_\xi) 
  + A_\eta \partial_\eta(2 \mu \partial_\eta \overline{u}_\eta) 
} 
$$

with:
$$ 
\boxed{ 
  \begin{aligned} 
    A_\xi(\tau)
    &= 
    \mathrm{e}^{3\tau/5},
    \\ 
    A_\eta(\tau)
    &=
    \mathrm{e}^{7\tau/5}. 
  \end{aligned} 
} 
$$

 Required coefficient fields (all face-centred):

  | On x-faces           | On y-faces           | 
  | ---------------------| -------------------- |
  | `mu.x`   = `Ax*mu.x` | `mu.y`   = `Ay*mu.y` |
  | `muAx.x` = `Ax*mu.x` | `muAx.y` = `Ax*mu.y` |
  | `muAy.x` = `Ay*mu.x` | `muAy.y` = `Ay*mu.y` |

 They are all defined in 
 [`selfsim_two-phase-generic_huppert.h`](https://basilisk.fr/sandbox/cailler/huppert_heap/selfsim_two-phase-generic_huppert.h).
 
 `mu` is the directional normal-stress coefficient used in the diagonal
 terms. `muAx` and `muAy` separately carry the $\xi-$ and $\eta-$derivative coefficients
 in the mixed shear stress. This header is deliberately restricted to 2D
 Cartesian simulations; it is not an AXI/SPHERISYM or 3D implementation.
*/

#include "poisson.h"

#if dimension != 2
#error "selfsim_viscosity_huppert.h is for 2D Cartesian simulations only"
#endif
#if AXI 
#error "selfsim_viscosity_huppert.h does not support AXI"
#endif

struct ViscosityHuppert {
  face vector mu;
  face vector muAx, muAy;
  scalar rho;
  double dt;
};

/* Relaxation for the coupled correction equations. */
static void relax_viscosity_huppert (scalar * a, scalar * b, int l,
                                     void * data)
{
  struct ViscosityHuppert * p = (struct ViscosityHuppert *) data;
  (const) face vector mu = p->mu;
  (const) face vector muAx = p->muAx, muAy = p->muAy;
  (const) scalar rho = p->rho;
  double dt = p->dt;
  vector u = vector(a[0]), r = vector(b[0]);

#if JACOBI
  vector w[];
#else
  vector w = u;
#endif

#if GAUSS_SEIDEL || _GPU
  vector ua[];
  foreach_level (l) {
    ua.x[] = u.x[];
    ua.y[] = u.y[];
  }
  boundary_level ((scalar *){ua}, l);
  for (int parity = 0; parity < 2; parity++)
    foreach_level_or_leaf (l) {
      if (level == 0 || ((point.i + parity) % 2) != (point.j % 2)) {
        double vu_top = (u.y[1,0] + ua.y[1,1])/4. -
                        (u.y[-1,0] + ua.y[-1,1])/4.;
        double vu_bottom = (ua.y[1,-1] + u.y[1,0])/4. -
                           (ua.y[-1,-1] + u.y[-1,0])/4.;
        double uu_right = (u.x[0,1] + ua.x[1,1])/4. -
                          (u.x[0,-1] + ua.x[1,-1])/4.;
        double uu_left = (ua.x[-1,1] + u.x[0,1])/4. -
                         (ua.x[-1,-1] + u.x[0,-1])/4.;
        w.x[] = (r.x[]*sq(Delta) + dt/rho[]*(
                    2.*mu.x[1]*u.x[1] + 2.*mu.x[]*u.x[-1]
                  + muAy.y[0,1]*u.x[0,1] + muAx.y[0,1]*vu_top
                  + muAy.y[]*u.x[0,-1] - muAx.y[]*vu_bottom)) /
          (sq(Delta) + dt/rho[]*(2.*mu.x[1] + 2.*mu.x[]
                                  + muAy.y[0,1] + muAy.y[]));
        w.y[] = (r.y[]*sq(Delta) + dt/rho[]*(
                    2.*mu.y[0,1]*u.y[0,1] + 2.*mu.y[]*u.y[0,-1]
                  + muAx.x[1]*u.y[1] + muAy.x[1]*uu_right
                  + muAx.x[]*u.y[-1] - muAy.x[]*uu_left)) /
          (sq(Delta) + dt/rho[]*(2.*mu.y[0,1] + 2.*mu.y[]
                                  + muAx.x[1] + muAx.x[]));
      }
    }
#else
  vector ua = u;
  foreach_level_or_leaf (l) {
    double vu_top = (u.y[1,0] + ua.y[1,1])/4. -
                    (u.y[-1,0] + ua.y[-1,1])/4.;
    double vu_bottom = (ua.y[1,-1] + u.y[1,0])/4. -
                       (ua.y[-1,-1] + u.y[-1,0])/4.;
    double uu_right = (u.x[0,1] + ua.x[1,1])/4. -
                      (u.x[0,-1] + ua.x[1,-1])/4.;
    double uu_left = (ua.x[-1,1] + u.x[0,1])/4. -
                     (ua.x[-1,-1] + u.x[0,-1])/4.;
    w.x[] = (r.x[]*sq(Delta) + dt/rho[]*(
                2.*mu.x[1]*u.x[1] + 2.*mu.x[]*u.x[-1]
              + muAy.y[0,1]*u.x[0,1] + muAx.y[0,1]*vu_top
              + muAy.y[]*u.x[0,-1] - muAx.y[]*vu_bottom)) /
      (sq(Delta) + dt/rho[]*(2.*mu.x[1] + 2.*mu.x[]
                              + muAy.y[0,1] + muAy.y[]));
    w.y[] = (r.y[]*sq(Delta) + dt/rho[]*(
                2.*mu.y[0,1]*u.y[0,1] + 2.*mu.y[]*u.y[0,-1]
              + muAx.x[1]*u.y[1] + muAy.x[1]*uu_right
              + muAx.x[]*u.y[-1] - muAy.x[]*uu_left)) /
      (sq(Delta) + dt/rho[]*(2.*mu.y[0,1] + 2.*mu.y[]
                              + muAx.x[1] + muAx.x[]));
  }
#endif

#if JACOBI
  foreach_level_or_leaf (l) {
    u.x[] = (u.x[] + 2.*w.x[])/3.;
    u.y[] = (u.y[] + 2.*w.y[])/3.;
  }
#endif

#if TRASH
  vector u1[];
  foreach_level_or_leaf (l) {
    u1.x[] = u.x[];
    u1.y[] = u.y[];
  }
  trash ({u});
  foreach_level_or_leaf (l) {
    u.x[] = u1.x[];
    u.y[] = u1.y[];
  }
#endif
}

/* Residual: res = r - u + dt/rho div(stress). */
static double residual_viscosity_huppert (scalar * a, scalar * b,
                                          scalar * resl, void * data)
{
  struct ViscosityHuppert * p = (struct ViscosityHuppert *) data;
  (const) face vector mu = p->mu;
  (const) face vector muAx = p->muAx, muAy = p->muAy;
  (const) scalar rho = p->rho;
  double dt = p->dt;
  vector u = vector(a[0]), r = vector(b[0]), res = vector(resl[0]);
  double maxres = 0.;

#if TREE
  /* Conservative flux-divergence form, matching the stock tree branch. */
  boundary ({u});
  face vector tauxx[], tauxy[], tauyx[], tauyy[];
  foreach_face (x) {
    tauxx.x[] = 2.*mu.x[]*(u.x[] - u.x[-1])/Delta;
    tauyx.x[] = (muAx.x[]*(u.y[] - u.y[-1]) +
                 muAy.x[]*((u.x[-1,1] + u.x[0,1])/4. -
                            (u.x[-1,-1] + u.x[0,-1])/4.))/Delta;
  }
  foreach_face (y) {
    tauxy.y[] = (muAy.y[]*(u.x[] - u.x[0,-1]) +
                 muAx.y[]*((u.y[1,-1] + u.y[1,0])/4. -
                            (u.y[-1,-1] + u.y[-1,0])/4.))/Delta;
    tauyy.y[] = 2.*mu.y[]*(u.y[] - u.y[0,-1])/Delta;
  }
  foreach (reduction(max:maxres)) {
    double dx = tauxx.x[1] - tauxx.x[] + tauxy.y[0,1] - tauxy.y[];
    double dy = tauyx.x[1] - tauyx.x[] + tauyy.y[0,1] - tauyy.y[];
    res.x[] = r.x[] - u.x[] + dt/rho[]*dx/Delta;
    res.y[] = r.y[] - u.y[] + dt/rho[]*dy/Delta;
    if (fabs(res.x[]) > maxres) maxres = fabs(res.x[]);
    if (fabs(res.y[]) > maxres) maxres = fabs(res.y[]);
  }
#else
  /* Naive uniform-grid stencil, matching the stock non-tree branch. */
  foreach (reduction(max:maxres)) {
    double dx = 2.*mu.x[1]*(u.x[1] - u.x[]) -
                2.*mu.x[]*(u.x[] - u.x[-1]) +
                muAy.y[0,1]*(u.x[0,1] - u.x[]) +
                muAx.y[0,1]*((u.y[1,0] + u.y[1,1])/4. -
                             (u.y[-1,0] + u.y[-1,1])/4.) -
                muAy.y[]*(u.x[] - u.x[0,-1]) -
                muAx.y[]*((u.y[1,-1] + u.y[1,0])/4. -
                           (u.y[-1,-1] + u.y[-1,0])/4.);
    double dy = 2.*mu.y[0,1]*(u.y[0,1] - u.y[]) -
                2.*mu.y[]*(u.y[] - u.y[0,-1]) +
                muAx.x[1]*(u.y[1] - u.y[]) +
                muAy.x[1]*((u.x[0,1] + u.x[1,1])/4. -
                           (u.x[0,-1] + u.x[1,-1])/4.) -
                muAx.x[]*(u.y[] - u.y[-1]) -
                muAy.x[]*((u.x[-1,1] + u.x[0,1])/4. -
                           (u.x[-1,-1] + u.x[0,-1])/4.);
    res.x[] = r.x[] - u.x[] + dt/rho[]*dx/sq(Delta);
    res.y[] = r.y[] - u.y[] + dt/rho[]*dy/sq(Delta);
    if (fabs(res.x[]) > maxres) maxres = fabs(res.x[]);
    if (fabs(res.y[]) > maxres) maxres = fabs(res.y[]);
  }
#endif
  return maxres;
}

/* Implicit treatment: same multigrid strategy as viscosity.h. */
mgstats viscosity_huppert (vector u, face vector mu,
                           face vector muAx, face vector muAy,
                           scalar rho, double dt,
                           int nrelax = 4, scalar * res = NULL)
{
  vector r[];
  foreach() {
    r.x[] = u.x[];
    r.y[] = u.y[];
  }
  restriction ({mu, muAx, muAy, rho});
  struct ViscosityHuppert p = {mu, muAx, muAy, rho, dt};
  return mg_solve ((scalar *){u}, (scalar *){r},
                   residual_viscosity_huppert,
                   relax_viscosity_huppert, &p, nrelax, res);
}

/* Explicit treatment: one residual evaluation, analogous to viscosity_explicit(). */
mgstats viscosity_huppert_explicit (vector u, face vector mu,
                                    face vector muAx, face vector muAy,
                                    scalar rho, double dt)
{
  vector r[];
  mgstats mg = {0};
  struct ViscosityHuppert p = {mu, muAx, muAy, rho, dt};
  mg.resb = residual_viscosity_huppert ((scalar *){u}, (scalar *){u},
                                        (scalar *){r}, &p);
  foreach() {
    u.x[] += r.x[];
    u.y[] += r.y[];
  }
  return mg;
}
