/**
# Soluble surfactant transport

This module solves the coupled transport of a soluble surfactant between
the bulk fluid and a moving fluid--fluid interface, in the hybrid
Volume-of-Fluid / Phase-Field framework of [Haouche *et al.*,
2026](#haouche2026) (extending the insoluble formulation of
[Farsoiya *et al.*, 2024](#farsoiya2024) to the soluble case).

Two scalar fields are evolved in time:

* `ci` : the *interfacial* surfactant concentration $f(\mathbf{x},t)$,
  a volumetric proxy for the surface concentration
  $\Gamma(\mathbf{x},t)$. It is localised in a thin layer of thickness
  $\epsilon$ around the interface through the phase-field
  $\phi(\mathbf{x},t)$.
* `cb` : the *bulk* surfactant concentration $F(\mathbf{x},t)$. It is
  transported in the exterior fluid ($\phi \to 0$) and exchanges mass
  with `ci` through the adsorption/desorption source term.

Both equations take the generalised elliptic form (see Appendix~B of
the paper)
$$
\nabla\!\cdot\!\bigl(\alpha\,\nabla c + \boldsymbol{\beta}\,c\bigr)
  + \lambda\,c = b,
$$
where $\boldsymbol{\beta}$ is an anti-diffusive flux aligned with the
interface normal (Jain's ACDI sharpening,
[Jain, 2022](#jain2022)). The system is solved with a custom multigrid
using the `h_relax` / `h_residual` pair defined below (the standard
Basilisk `diffusion.h` solver does not support the $\boldsymbol{\beta}$
term).

The phase-field $\phi$ (`pfield`) is reconstructed from the VoF
signed distance $\chi$ (`d2` in the code) via the tanh profile
$$
\phi = \tfrac{1}{2}\bigl(1 - \tanh(\chi / 2\epsilon)\bigr),
$$
and optionally advected/regularised by the ACDI equation at each time
step.

## Dependencies
*/

#include "redistance2.h"
#include "diffusion.h"
#include "tracer.h"

/**
## Interface thickness

The regularisation length $\epsilon$ is set to $0.75\,\Delta x$ at
the finest level, so that the tanh profile of $\phi$ resolves the
interface over $\sim 3$ cells.
*/

#define EPSILON ((L0/(1 << grid->maxdepth))*0.75)

/**
## Runtime switches

* `advect_diff_phase_field` --- when `1`, the ACDI
  advection--sharpening equation is solved for `pfield`. When `0`,
  `pfield` is only reconstructed from the signed distance `d2`.
  The flag is set automatically in [`properties2`](#properties2)
  depending on `reinit_skip_steps`.
* `reinit_skip_steps` --- number of time steps between two
  redistance operations on `d2`. Large values privilege ACDI
  transport; `reinit_skip_steps = 1` reduces to pure level-set
  reconstruction.
* `surfactant` --- master switch. When `0` the module is dormant
  (useful for clean-interface runs).
* `counter` --- internal time-step counter used by the redistance
  scheduler.
*/

bool advect_diff_phase_field = 0;
int reinit_skip_steps = 100;
bool surfactant = 1;
int counter = 0;

/**
## Fields

`ci` and `cb` are the interfacial and bulk concentrations. `d2` is the
signed distance function, `pfield` the phase-field, `gamma2` a
user-accessible surface density field. All of them are VoF tracers so
they are advected with the interface-capturing scheme in
[`tracer.h`](/src/tracer.h). */

scalar ci[], cb[], gamma2[], d2[], pfield[];
scalar * tracers = {ci, cb, d2, pfield};

/**
## Physical parameters

| Symbol       | Code name                  | Meaning                                                              |
|--------------|----------------------------|----------------------------------------------------------------------|
| $\zeta$      | `zeta`                     | ACDI sharpening velocity (set dynamically from $|\mathbf{u}|_{\max}$) |
| $D_f$        | `Di0`                      | Interfacial diffusion coefficient                                    |
| $D_F$        | `Db0`                      | Bulk diffusion coefficient                                           |
| $r_a$        | `ra`                       | Adsorption rate                                                      |
| $r_d$        | `rd`                       | Desorption rate                                                      |
| $f_\infty$   | `c_hat_inf`                | Saturation (maximum packing) concentration                           |
| $\varepsilon_{\text{reg}}$ | `varepsilon` | Numerical regularisation for $\log$ / division by zero             |
| $k_s$        | `sharpening_coefficient`   | ACDI interface-sharpening coefficient (default $2$)                  |
| -            | `inverse_pfield`           | Sign convention: $\phi$ or $1-\phi$ is the "inside" phase            |
*/

double zeta = 1. [*];
double Di0 = 0.01 [*], Db0 = 0.01 [*];
double ra = 1. [*], rd = 1. [*];
double c_hat_inf = 1. [*];
double varepsilon = 1.e-12;
double sharpening_coefficient = 2.;
double inverse_pfield = 1; //inverting the field according to problem

/**
## Generalised diffusion multigrid

We solve the linear system
$$
\nabla\!\cdot\!\bigl(D\,\nabla c + \boldsymbol{\beta}\,c\bigr)
  + \lambda\,c = b
$$
with a geometric multigrid. The coefficients are packed in
`HDiffusion`. The Jacobi relaxation `h_relax` and the residual
`h_residual` follow the second-order finite-volume discretisation
described in Section~3.2 and Appendix~A of the paper. */

struct HDiffusion {
    face vector D;
    face vector beta;
    scalar lambda;
};

static void h_relax(scalar *al, scalar *bl, int l, void *data) {
    scalar a = al[0], b = bl[0];
    struct HDiffusion *p = (struct HDiffusion *) data;
    face vector D = p->D, beta = p->beta;
    scalar lambda = p->lambda;

    scalar c = a;
    foreach_level_or_leaf(l) {
        double n = -sq(Delta)*b[], d = -lambda[]*sq(Delta);

        foreach_dimension() {
            n += D.x[1]*a[1] + D.x[]*a[-1] + 0.5*Delta*(beta.x[1]*a[1] - beta.x[]*a[-1]);
            d += D.x[1] + D.x[] - 0.5*Delta*(beta.x[1] - beta.x[]);
        }
        c[] = n / d;
    }
}

static double h_residual(scalar *al, scalar *bl, scalar *resl, void *data) {
    scalar a = al[0], b = bl[0], res = resl[0];
    struct HDiffusion *p = (struct HDiffusion *) data;
    face vector D = p->D, beta = p->beta;
    scalar lambda = p->lambda;
    double maxres = 0.;

#if TREE
    face vector g[];
    foreach_face() {
        g.x[] = D.x[]*face_gradient_x(a, 0) + beta.x[]*face_value(a, 0);
    }
    foreach(reduction(max:maxres)) {
        res[] = b[] - lambda[]*a[];
        // res[] = b[] + cm[]/dt*a[];

        foreach_dimension() {
            res[] -= (g.x[1] - g.x[]) / Delta;
        }
        if (fabs(res[]) > maxres) maxres = fabs(res[]);
    }
#else
    foreach(reduction(max:maxres)) {
        res[] = b[] - lambda[]*a[];
        foreach_dimension() {
            res[] -= (D.x[1]*face_gradient_x(a, 1) - D.x[0]*face_gradient_x(a, 0) + beta.x[1]*face_value(a, 1) - beta.x[0]*face_value(a, 0))/Delta;
        }
        if (fabs(res[] > maxres)) maxres = fabs(res[]);
    }
#endif

    return maxres;
}

/**
## Stability {#stability}

The time step is bounded by a CFL-like condition that aggregates the
ACDI diffusion, the anti-diffusive sharpening flux, and the physical
diffusivities $D_f$, $D_F$:
$$
\Delta t \leq 0.75
  \frac{\Delta x^2}{\zeta\,\epsilon},\,\,\,\Delta t \leq 0.75
  \frac{\Delta x^2}{2d\,\text{max}(D_f,\,D_F)},
$$
with $d$ the spatial dimension (factor $2d = 4$ in 2D/axi, $6$ in 3D).

$\zeta$ is set to $1.1\,|\mathbf{u}|_{\max}$ evaluated *at the
interface only* (weighted by $\delta_\phi = \phi(1-\phi)/\epsilon$),
which avoids stiffening $\Delta t$ with large bulk velocities far from
the interface. */

event stability(i++) {
    double maxvel = 1.e-6 [*];
    double Dmax = max(Di0, Db0);

    if (surfactant) {
        foreach_face(reduction(max:maxvel)) {
            if (u.x[] != 0.) {
                double deltas = (pfield[]*(1. - pfield[])) / EPSILON;
                if (deltas > 0.01) {
                    if (fabs(u.x[]) > maxvel) maxvel = fabs(u.x[]);
                }
            }
        }

        double deltaxmin = 1. [*];
        deltaxmin = L0 / (1 << grid->maxdepth);
        zeta = 1.1*maxvel;
        double dt = dtmax;

        if (advect_diff_phase_field) {
            dt = 0.75*sq(deltaxmin) / zeta / EPSILON;
            if (dt < dtmax) dtmax = dt;
        }

        dt = 0.75*sq(deltaxmin) / 4. / Dmax;
#if dimension == 3
        dt = 0.75*sq(deltaxmin) / 6. / Dmax;
#endif
    }
}

/**
## Safe clamp

Clips a value to $[0,1]$ with a small guard band to avoid
`log(0)` / division by zero in the ACDI formulation. */

double clamp2(double a) {
    if (a < 1.e-6) return 0;
    else if (a > 1. - 1.e-6) return 1.;
    else return a;
}

/**
## Phase-field reconstruction {#properties2}

Every `reinit_skip_steps` time steps, the signed distance `d2` is
recomputed from the volume fraction `f` by fast-marching
(`redistance`) and the phase-field is rebuilt from the tanh profile.

A sanity check on the amplitude of the signed distance rejects
pathological redistance results (e.g. early in the simulation when the
interface is not yet well-formed): in that case the redistance is
retried at the next time step.
*/

event properties2(i++) {
    if (reinit_skip_steps > 1) advect_diff_phase_field = 1;
    else advect_diff_phase_field = 0;

    if (counter >= 0 && counter%reinit_skip_steps == 0) {
        scalar d2temp[];
        double deltamin = L0 / (1 << grid->maxdepth);
        foreach() {
            d2temp[] = (2.*f[] - 1.)*deltamin*0.75;
        }
#if TREE
        restriction({d2temp});
#endif
        redistance(d2temp, imax = 0.5*(1 << grid->maxdepth));

        double d2max = statsf(d2temp).max;
        double d2min = statsf(d2temp).min;
        bool signed_distance_faulty = 0;

        if (d2max > 6. || d2min < -6.) {
            signed_distance_faulty = 1;
            counter = counter - 2;
        }

        if (!signed_distance_faulty) {
            foreach() {
                d2[] = inverse_pfield*d2temp[];
                pfield[] = 0.5*(1. - tanh((d2[] / 2. / EPSILON)));
                pfield[] = clamp2(pfield[]);
            }
            boundary({pfield});
        }
    }
    counter ++;
}

/**
## Coupled surfactant solve

This is the core event. Three linear systems are solved per time step:

1. **ACDI equation for the phase-field** (only when
   `advect_diff_phase_field = 1`):
   $$
   \frac{\partial\phi}{\partial t}
     = \nabla\!\cdot\!\bigl(\zeta\epsilon\,\nabla\phi
       - \tfrac{1}{4}\zeta\,(1-\tanh^2(\psi/2\epsilon))\,\mathbf{n}\bigr),
   $$
   with $\psi = \epsilon\log\!\bigl(\phi/(1-\phi)\bigr)$.

2. **Interfacial surfactant**
   (Eq.~(B.15) of the paper, implicit in time):
   $$
   \nabla\!\cdot\!\bigl(\alpha_f\nabla f^{n+1}
     + \boldsymbol{\beta}_f f^{n+1}\bigr) + \lambda_f f^{n+1} = b_f,
   $$
   $$
   \alpha_f = -D_f,\quad
   \boldsymbol{\beta}_f = -\mathbf{B}_f,\quad
   \lambda_f = \tfrac{1}{\Delta t} + r_d + r_a\,\tfrac{F^{n+1}}{\phi^{n+1}},\quad
   b_f = \tfrac{f^n}{\Delta t} + r_a\,\tfrac{F^{n+1}}{\phi^{n+1}} f_\infty.
   $$

3. **Bulk surfactant**
   (Eq.~(B.20) of the paper):
   $$
   \nabla\!\cdot\!\bigl(\alpha_F\nabla F^{n+1}
     + \boldsymbol{\beta}_F F^{n+1}\bigr) + \lambda_F F^{n+1} = b_F.
   $$

The anti-diffusive fluxes are aligned with the smoothed interface
normal $\mathbf{n} = \nabla\psi / |\nabla\psi|$ computed from the
ACDI potential $\psi$. The implicit source term in $\lambda_f$ ($\lambda_i$ in the code) /
$\lambda_F$ ($\lambda_b$ in the code) stabilises the stiff adsorption/desorption kinetics.
*/

event tracer_diffusion(i++) {
    if (surfactant) {

        scalar rhs[], psi[], r[], lambda[];
        face vector cflux[], geta[], beta[], D[];
        face vector betai[], betab[], Di[], Db[];
        scalar ri[], rb[], lambdai[], lambdab[];

        /**
        ### 1. ACDI phase-field sharpening
        */
        if (advect_diff_phase_field) {
            foreach() {
                pfield[] = clamp2(pfield[]);
                psi[] = EPSILON*log((pfield[] + varepsilon) / (1. - pfield[] + varepsilon));
                lambda[] = - cm[] / dt;
            }
            boundary({psi});

            foreach_face() {
                cflux.x[] = 0.;
                double psif = (psi[] + psi[-1])/2.;
                double gradpsi = (psi[] - psi[-1])/Delta;

                if (fabs(gradpsi) > varepsilon) {

#if dimension == 2
                    double psiup = 0.25*(psi[] + psi[-1] + psi[0,1] + psi[-1,1 ]);
                    double psidown = 0.25*(psi[] + psi[-1] + psi[0,-1] + psi[-1,- 1]);
                    double mag_grad_psi = sqrt(sq(psi[] - psi[-1]) + sq(psiup - psidown)) / Delta;
#endif

#if dimension == 3
                    double psiup = 0.25*(psi[0, 0, 0] + psi[-1, 0, 0] + psi[0, 1, 0] + psi[-1, 1, 0]);
                    double psidown = 0.25*(psi[0, 0, 0] + psi[-1, 0, 0] + psi[0, -1, 0] + psi[-1, -1, 0]);
                    double psifront = 0.25*(psi[0, 0, 0] + psi[-1, 0, 0] + psi[0, 0, 1] + psi[-1, 0, 1]);
                    double psiback = 0.25*(psi[0, 0, 0] + psi[-1, 0, 0] + psi[0, 0, -1] + psi[-1, 0, -1]);
                    double mag_grad_psi = sqrt(sq(psi[0, 0, 0] - psi[-1, 0, 0]) + sq(psiup - psidown) + sq(psifront - psiback)) / Delta;
#endif

                    cflux.x[] = fm.x[]*(1. - sq(tanh(psif/2./EPSILON)))*gradpsi/mag_grad_psi;
                }

                D.x[] = fm.x[]*zeta*EPSILON;
                beta.x[] = 0.;
            }

            foreach() {
#if dimension == 2
                rhs[] = - 0.25*zeta*(cflux.x[1] - cflux.x[] + (cflux.y[0, 1] - cflux.y[])/cm[]) / Delta*cm[];
                r[] = - cm[]*pfield[]/dt + 0.25*zeta*(cflux.x[1] - cflux.x[] + (cflux.y[0, 1] - cflux.y[])) / Delta;
#endif

#if dimension == 3
                rhs[] = - 0.25*zeta*(cflux.x[1] - cflux.x[] + cflux.y[0, 1] - cflux.y[] + cflux.z[0, 0, 1] - cflux.z[]) / Delta*cm[];
                r[] = - cm[]*pfield[]/dt + 0.25*zeta*(cflux.x[1, 0, 0] - cflux.x[0, 0, 0] + cflux.y[0, 1, 0] - cflux.y[0, 0, 0] + cflux.z[0, 0, 1] - cflux.z[0, 0, 0]) / Delta;
#endif
            }

            restriction({D, beta, cm, lambda});
            struct HDiffusion q1;
            q1.D = D;
            q1.beta = beta;
            q1.lambda = lambda;

            mg_solve({pfield}, {r}, h_residual, h_relax, &q1);
        }

        /**
        ### 2. Build cell-centred coefficients for `ci` and `cb`

        `deltas` is the regularised delta function
        $\delta_\phi = \phi(1-\phi)/\epsilon$ used to localise the
        adsorption term. Kinetics disabled by commenting out the
        full expression and using the *adsorption-only* variant
        (single uncommented line below each block).
        */
        foreach() {
            pfield[] = clamp2(pfield[]);
            double deltas = (pfield[]*(1. - pfield[]))/EPSILON;
            double phi = max(pfield[], 1.e-6);
            psi[] = EPSILON*log((clamp(pfield[],0.,1.) + varepsilon) / (1. - clamp(pfield[],0.,1.) + varepsilon));

#if dimension == 2
            double phiup   = 0.25*(pfield[]   + pfield[-1] + pfield[0,1] + pfield[-1,1]);
            double phidown = 0.25*(pfield[]   + pfield[-1] + pfield[0,-1]+ pfield[-1,-1]);
            double mag_grad_phi = sqrt(sq(pfield[] - pfield[-1]) + sq(phiup - phidown)) / Delta;
#endif

#if dimension == 3
            double phiup    = 0.25*(pfield[0,0,1] + pfield[-1,0,1] + pfield[0,1,1] + pfield[-1,1,1]);
            double phidown  = 0.25*(pfield[0,0,-1]+ pfield[-1,0,-1] + pfield[0,1,-1]+ pfield[-1,1,-1]);
            double phifront = 0.25*(pfield[0,0,0] + pfield[-1,0,0] + pfield[0,0,1] + pfield[-1,0,1]);
            double phiback  = 0.25*(pfield[0,0,0] + pfield[-1,0,0] + pfield[0,0,-1]+ pfield[-1,0,-1]);
            double mag_grad_phi = sqrt(sq(pfield[] - pfield[-1]) + sq(phiup - phidown) + sq(phifront - phiback)) / Delta;
#endif

            //lambdai[] = -cm[]*(1./dt); // to simulate only adsorption
            lambdai[] = -cm[]*(1./dt + ra*cb[]/phi + rd);

            //lambdab[] = -cm[]*(1./dt + ra*c_hat_inf*deltas/phi); // to simulate only adsorption
            lambdab[] = -cm[]*(1./dt + ra*(c_hat_inf*deltas - ci[])/phi);

            //ri[] = -cm[]*(ci[]/dt + ra*cb[]*c_hat_inf*deltas/phi); // to simulate only adsorption
            ri[] = -cm[]*(ci[]/dt + ra*cb[]*c_hat_inf*deltas/phi);

            rb[] = -cm[]*(cb[]/dt + rd*ci[]);
        }

        /**
        ### 3. Face-centred diffusion and anti-diffusive fluxes

        The anti-diffusive flux $\mathbf{B}$ is aligned with the
        interface normal computed from $\nabla\psi/|\nabla\psi|$:
        $$
        \mathbf{B}_f = -D_f\,k_s\,\frac{0.5 - \phi}{\epsilon}\,\mathbf{n},
        \qquad
        \mathbf{B}_F = -D_F\,\frac{1-\phi}{\epsilon}\,\mathbf{n}.
        $$
        */
        foreach_face() {
            double phif = (pfield[] + pfield[-1])/2.;
            Di.x[] = fm.x[]*Di0;
            Db.x[] = fm.x[]*Db0;
            double gradpsi = (psi[] - psi[-1])/Delta;
            betai.x[] = betab.x[] = 0.;

            if (fabs(gradpsi) > varepsilon) {
#if dimension == 2
                double psiup   = 0.5*(psi[0,1] + psi[-1,1]);
                double psidown = 0.5*(psi[0,-1] + psi[-1,-1]);
                double mag_grad_psi = sqrt(sq(psi[] - psi[-1]) + sq(psiup - psidown)/4.)/Delta;
#endif

#if dimension == 3
                double psiup    = 0.5*(psi[0,1,0] + psi[-1,1,0]);
                double psidown  = 0.5*(psi[0,-1,0] + psi[-1,-1,0]);
                double psifront = 0.5*(psi[0,0,1] + psi[-1,0,1]);
                double psiback  = 0.5*(psi[0,0,-1] + psi[-1,0,-1]);
                double mag_grad_psi = sqrt(sq(psi[] - psi[-1]) + sq(psiup - psidown)/4. + sq(psifront - psiback)/4.)/Delta;
#endif

                betai.x[] = -Di.x[]*sharpening_coefficient*(0.5 - phif)/EPSILON*gradpsi/mag_grad_psi;
                betab.x[] = -Db.x[]*(1. - phif)/EPSILON*gradpsi/mag_grad_psi;
            }
        }

        restriction({Di, Db, betai, betab, lambdai, lambdab, cm, ri, rb});

        /**
        ### 4. Solve the two decoupled linear systems

        The source terms `ri`, `rb` are treated as explicit (they use
        $c_b^n$ and $c_i^n$); the $\lambda_i$, $\lambda_b$ coefficients
        embed the implicit diagonal contributions from the kinetics.
        */
        struct HDiffusion qi;
        qi.D      = Di;
        qi.beta   = betai;
        qi.lambda = lambdai;

        mg_solve({ci}, {ri}, h_residual, h_relax, &qi);

        struct HDiffusion qb;
        qb.D      = Db;
        qb.beta   = betab;
        qb.lambda = lambdab;

        mg_solve({cb}, {rb}, h_residual, h_relax, &qb);
    }
}

/**
## Usage

## Test

* [Expanding circle](test_cases/expanding_circle.c)
* [Expanding sphere](test_cases/expanding_sphere.c)
* [Rotating circle](test_cases/rotating_circle.c)
* [Pure adsorption on a flat interface](test_cases/pure_adsorption_flat.c)
* [Pure adsorption on a circular interface](test_cases/pure_adsorption_circle.c)
* [Pure adsorption on a sphere interface](test_cases/pure_adsorption_sphere.c)
* [Pure desorption on a flat interface](test_cases/pure_desorption.c)
* [Soluble surfactants on a flat interface](test_cases/flat_interface.c)
* [Rising bubble in 2D-axi](test_cases/rising_bubble_axi.c)
* [Rising bubble in 3D](test_cases/rising_bubble_3D.c)

## References

~~~bib
@article{haouche2026,
  author  = {Haouche, I. and Reichert, B. and Baudoin, M. and Farsoiya, P. K.},
  title   = {A hybrid Volume of Fluid Phase-Field method for Direct
             Numerical Simulations of soluble surfactant-laden
             interfacial flows},
  journal = {Journal of Computational Physics},
  year    = {2026},
  note    = {Submitted}
}

@article{farsoiya2024,
  author  = {Farsoiya, P. K. and Popinet, S. and Stone, H. A. and
             Deike, L.},
  title   = {A coupled Volume-of-Fluid / Phase-Field method for
             insoluble surfactants},
  journal = {Journal of Fluid Mechanics},
  year    = {2024}
}

@article{jain2022,
  author  = {Jain, S. S.},
  title   = {Accurate conservative phase-field method for simulation
             of two-phase flows},
  journal = {Journal of Computational Physics},
  volume  = {469},
  pages   = {111529},
  year    = {2022}
}
~~~
*/