/**
# Vertical diffusion (delta form)

Variant of [diffusionH.h](diffusionH.h) that solves for the
*temperature change* instead of the absolute temperature.

The absolute form computes

  rhs_l = h_l s_l   (+ flux terms)

which for T = 20 °C and h = 20 m is a value around 400.
In single precision (float32, ~7 decimal digits) the rounding
error of each entry is of order 400 * 6e-8 ~ 2.4e-5, and these
rounding errors accumulate with a systematic sign over many
timesteps – this is the slow "cooling" drift seen on GPUs / in
single precision even though the scheme conserves energy
exactly in exact arithmetic.

The delta form computes the residual

  rhs_l = h_l*s_l - (A*s)_l + F_l

where A is the (exactly conservative) diffusion matrix. For a
steady state this residual is ~0, so all rounding errors are of
order ulp * |change| instead of ulp * |s|, and they accumulate
without a systematic sign. The matrix itself still conserves
energy exactly, so no drift accumulates.

Interface is identical to vertical_diffusion2().
*/
//@define double float

#define DIFFUSIONH_DELTA_VERSION "2026-09-22-a"
#if DIFF_DELTA
# pragma message ("diffusionH_delta.h version " DIFFUSIONH_DELTA_VERSION)
#endif

void vertical_diffusion_delta (Point point, scalar h, scalar s,
			       double dt, double D,
			       double dst, double dsb)
{
  double a[nl], b[nl], c[nl], rhs[nl], sold[nl];

  foreach_layer()
    sold[point.l] = s[];

  /**
  The lower, principal and upper diagonals $a$, $b$ and $c$
  (same as in diffusionH.h). */

  for (int l = 1; l < nl - 1; l++) {
    a[l] = - 2.*D*dt/(h[0,0,l-1] + h[0,0,l]);
    c[l] = - 2.*D*dt/(h[0,0,l] + h[0,0,l+1]);
    b[l] = h[0,0,l] - a[l] - c[l];
  }

  a[nl-1] = - 2.*D*dt/(h[0,0,nl-2] + h[0,0,nl-1]);
  b[nl-1] = h[0,0,nl-1] - a[nl-1];

  b[0] = h[] + 2.*dt*D*(1./(h[] + h[0,0,1]));
  c[0] = - 2.*dt*D*(1./(h[] + h[0,0,1]));

  /**
  Residual right-hand side: rhs = h·s − (A·s) + F, computed
  directly in delta form so that a steady state gives rhs ≈ 0
  and rounding errors stay of order ulp·|δs| instead of ulp·|s|. */

  if (nl > 1) {
    /* interior layers: A·s row = a·s[l-1] + b·s[l] + c·s[l+1],
       with b = h − a − c the residual becomes
       −a·(s[l−1]−s[l]) − c·(s[l+1]−s[l]) */
    for (int l = 1; l < nl - 1; l++)
      rhs[l] = - a[l]*(sold[l-1] - sold[l]) - c[l]*(sold[l+1] - sold[l]);

    /* top layer: row = a·s[nl-2] + b·s[nl-1], with b = h − a */
    rhs[nl-1] = - a[nl-1]*(sold[nl-2] - sold[nl-1]) + D*dt*dst;

    /* bottom layer: row = b·s[0] + c·s[1] */
    rhs[0] = - c[0]*(sold[1] - sold[0]) - D*dt*dsb;
  }
  else {
    /* nl == 1: flux-only right-hand side (see TODO below) */
    rhs[0] = (h[] - b[0] - c[0])*sold[0] + (- c[0]*h[] - D*dt)*dst;
    b[0] += c[0];
  }

  for (int l = 1; l < nl; l++) {
    b[l] -= a[l]*c[l-1]/b[l-1];
    rhs[l] -= a[l]*rhs[l-1]/b[l-1];
  }
  a[nl-1] = rhs[nl-1]/b[nl-1];
  s[0,0,nl-1] = sold[nl-1] + a[nl-1];
  for (int l = nl - 2; l >= 0; l--)
    s[0,0,l] = sold[l] + (a[l] = (rhs[l] - c[l]*a[l+1])/b[l]);
}

/* # TODO

- fix the case nl=1

**/
