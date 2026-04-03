/**
# Piecewise Polynomial Reconstruction

This is inspired by the [PPR Fortran90
library](https://github.com/dengwirda/PPR) of Darren Engwirda and
Maxwell Kelley, which is included in [/src/ppr]().

The initial version was written by Antoine Aubert.

It is used for [layer remapping](/src/layered/remap.h) within the
[multilayer framework](/src/layered/README).

The basic idea is to construct a second-order polynomial (a parabola)
fitting the integral over each layer and layer-interface values
obtained with fourth-order polynomial fits. Monotonic limiting can
also be applied.

A significant difference with the PPR library is that explicit boundary
conditions (Neumann or Robin/Navier) are imposed on the left and right
(or top and bottom) boundaries. Also, being written in C99, this code
runs transparently on GPUs.

## Neumann boundary conditions

We start with the simplest case i.e. the reconstruction of a parabola
with Neumann conditions on the left boundary and a Dirichlet condition
on the right boundary.

The polynomial is given by
$$
P(x) = \sum_i C_i x^i
$$
It must verify the following conditions
$$
\begin{aligned}
P'(0) &= df_b, \\
\int_0^1P(x)dx &= f \\
P(1) &= s_r
\end{aligned}
$$
which leads to the code below. */

static void remap_neumann (double s_r,
                           double f,
                           double df_b,
                           double C[3])
{
  C[0] = s_r - df_b - C[2];
  C[1] = df_b;
  C[2] = 3.*(s_r - df_b/2. - f)/2.;
}

/**
## Robin/Navier boundary conditions

This is a bit more complicated, with the following conditions
$$
\begin{aligned}
P(0) &= f_b + \lambda_b P'(0), \\
\int_0^1P(x)dx &= f \\
P(1) &= s_r
\end{aligned}
$$
which leads to the solution of the linear system coded below. */

static void remap_robin (double s_r,
                         double f,
                         double f_b, double lambda_b,
                         double C[3])
{
  double a = - lambda_b, s = - a + (a - 1.)/3. + 1/2. [0];
  double M[3][3] = {
    {    - a,    1./6.,     a/3.},
    {     1.,  - 2./3.,  - 1./3.},
    { a - 1.,    1./2,  1/2. - a}
  };
  double B[3] = { f, f_b, s_r };
  for (int i = 0; i < 3; i++) {
    C[i] = 0.;
    for (int j = 0; j < 3; j++)
      C[i] += M[i][j]*B[j];
    C[i] /= s;
  }
}

/**
## Fourth-order polynomial for interface values

The values between layers are approximated by fitting a fourth-order
polynomial to the two layers above and below the interface. This gives
the 4x4 linear system
$$
\int_{x_i}^{x_{i+1}} P(x') dx' = f_i
$$
with $k - 1 \leq i \leq k + 2$. The integration is normalized using
$$
x' = \frac{x}{x_{k+1} - x_k}
$$
*/

static double right_value (int k, int n,
                           const double x[n], const double f[n-1],
                           double f_b, double lambda_b, double df_b,
                           double f_t, double lambda_t, double df_t)
{
  double xk = x[k], dx = x[k+1] - xk;
  double M[4][4], B[4];
  for (int i = (k == 0); i < (k == n - 3 ? 3 : 4); i++) {
    double zeta1 = (x[i+k-1] - xk)/dx;
    double z1n = zeta1;
    double zeta0 = (x[i+k] - xk)/dx;
    double z0n = zeta0;
    for (int j = 0; j < 4; j++, z1n *= zeta1, z0n *= zeta0)
      M[i][j] = (z0n - z1n)/(j + 1);
    B[i] = (x[k+i] - x[k-1+i])/dx*f[k-1+i];
  }

  /**
  The left (or bottom) and right (or top) layers must use the left and
  right boundary conditions to close the system. 

  For the bottom layer we have: */
  
  if (k == 0) {

    /**
    If `df_b` is non-zero we apply a Neumann boundary condition, which
    changes one equation of the default linear system to */
    
    if (df_b) {
      M[0][0] = 0.;
      M[0][1] = 1.;
      M[0][2] = 0.;
      M[0][3] = 0.;
      B[0] = df_b*dx;
    }

    /**
    Otherwise we apply a Robin condition, as */
    
    else {
      M[0][0] = 1. [0];
      M[0][1] = - lambda_b/dx;
      M[0][2] = 0.;
      M[0][3] = 0.;
      B[0] = f_b;
    }
  }

  /**
  This is the equivalent for the top layer. The systems are different
  because at this location because $x'=1$ for the top layer and $x'=0$
  for the bottom layer. */
  
  if (k == n - 3) {
    double zetab = (x[k+2] - xk)/dx;
    if (df_t) {
      M[3][0] = 0.;
      M[3][1] = 1.;
      M[3][2] = 2.*zetab;
      M[3][3] = 3.*sq(zetab);
      B[3] = df_t*dx;
    }
    else {
      M[3][0] = 1.;
      M[3][1] = zetab - lambda_t/dx;
      M[3][2] = sq(zetab) - 2.*zetab*lambda_t/dx;
      M[3][3] = cube(zetab) - 3.*sq(zetab)*lambda_t/dx;
      B[3] = f_t;
    }
  }

  /**
  ### Solution of the linear system

  We invert the linear system to get the coefficients $C_i$ of the
  polynomial and return the value on the "right" side of the interval
  by computing $P(x'=1)$. */
  
  assert (smatrix_inverse (4, M, 1.e-30) > 1.e-20);

  double c[4] = {0.,0.,0.,0.};
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      c[i] += M[i][j]*B[j];
  double s_right = c[3];
  for (int i = 0; i < 3; i++)
    s_right += c[i];
  return s_right;
}

/**
## Polynomial reconstruction for "central" layers

We first define a simple monotonic limiter. */

static double minmodremap (double a, double b) {
  return a*b <= 0. ? 0. :
    fabs(a) < fabs(b) ? a : b;
}

/**
The reconstruction function take the left and right values, the layers
positions `x` and corresponding average values `f`, a limiting option
and returns the polynomial `C` on the interval $[x_k:x_{k+1}]$. */

static void remap_central (const double s_l, const double s_r,
                           const int n, const double x[n], const double f[n-1],
                           const int k,
                           double C[3],
                           const bool limiter)
{
  double xk = x[k], dx = x[k+1] - xk;
  double sl = s_l, sr = s_r;

  /**
  If limiting is used, we first limit the left and right values. */
  
  if (limiter) {
    if ((f[k+1] - f[k])*(f[k] - f[k-1]) < 0.) {
      sl = f[k];
      sr = f[k];
    }
    else {
      double sigma_left = 2.*(f[k] - f[k-1])/dx;
      double sigma_center = 2.*(f[k+1] - f[k-1])/(x[k+2] - x[k-1] + dx);
      double sigma_right = 2.*(f[k+1] - f[k])/dx;
      double sigma = minmodremap(sigma_center, minmodremap(sigma_left, sigma_right));
      if ((sl - f[k-1])*(f[k] - sl) < 0.)
        sl = f[k] - 1./2.*dx*sigma;
      if ((sr - f[k])*(f[k+1] - sr) < 0.)
        sr = f[k] + 1./2.*dx*sigma;
    }
  }

  /**
  The polynomial fit is given by the linear system
  $$
  \begin{aligned}
  P(0) &= s_l, \\
  \int_0^1P(x)dx &= f_k \\
  P(1) &= s_r
  \end{aligned}
  $$
  which has for solution */
  
  C[0] = sl;
  C[1] = 6.*f[k] - 2.*sr - 4.*sl;
  C[2] = 3.*(sr + sl - 2.*f[k]);

  /**
  If limiting is used, we check that the extrema of the polynomial do
  not over- or undershoot. For a second-order polynomial, the extrema (if
  $C_2\neq0$) is at
  $x=-\frac{C_1}{2C_2}$ which can be developed as
  $$
  x=-\frac{6f_k-2s_r-4s_l}{6(s_r+s_l-3f_k)}
  $$  
  If $x\in[0,0.5]$, we change $s_r$ such that
  $x=0$. The new value is then $s_r=3f_k-2s_l$.
  
  If $x\in[0.5,1]$, we change $s_l$ such that
  $x=1$. The new value is then $s_l=3f_k-2s_r$.
  */
  
  if (limiter && C[1]*C[2] < 0. && C[1]/C[2] > - 2.) {
    if (C[1]/C[2] > - 1.)
      sr = 3.*f[k] - 2.*sl;
    else
      sl = 3.*f[k] - 2.*sr;

    /**
    We update the polynomial coefficients using these new bounds. */
    
    C[0] = sl;
    C[1] = 6.*f[k] - 2.*sr - 4.*sl;
    C[2] = 3.*(sr + sl - 2.*f[k]);
  }
}

/**
## Remapping function

This is the remapping interface. It takes as arguments the number of
initial positions `npos`, the number of "remapped" positions `nnew`
and the corresponding arrays of initial and remapped positions `xpos`
and `xnew` respectively, as well as the initial averaged values on the
corresponding intervals `fdat`. The remapped averaged values will be
stored in `fnew`.

Note that the `xpos` and `xnew` must be stored in strict increasing
order i.e. negative interval values $x_{i+1} - x_i$ are not allowed.

The "bottom" and "top" boundary conditions are given by the
$(f_b,\lambda_b,df_b)$ parameters (resp. $(f_t,\lambda_t,df_t)$). If
$df_b$ (resp. $df_t$) is non-zero, then Neumann boundary conditions
are applied i.e.
$$
\partial_n P(0) = df_b
$$
where $n$ is the direction "normal" to the boundary i.e. $+x$ for the
bottom/left boundary and $-x$ for the top/right boundary.

If $df_b$ (resp. $df_t$) is zero, then Robin/Navier boundary conditions
are applied i.e.
$$
P(0) = f_b + \lambda_b \partial_n P(0)
$$

To impose a Neumann zero boundary condition one can set $\lambda_b$ to
`HUGE` and $f_b$ to zero.

If limiter is `true` then monotonic limiters are applied. In this case
the Neumann zero condition leads to a constant profile in the top or
bottom layer, which guarantees boundedness of the remapping. Note that
other boundary conditions will not guarantee boundedness. */

void remap_c (int npos, int nnew,
              const double xpos[npos], const double xnew[nnew],
              const double fdat[npos-1], double fnew[nnew-1],
              double f_b, double lambda_b, double df_b,
              double f_t, double lambda_t, double df_t,
              bool limiter)
{

  /**
  We can only remap on the same interval. It would be possible to
  remap on a smaller interval (a subset) but there is no need for now
  and this restriction makes things simpler. */
    
  assert (xnew[0] == xpos[0] && xnew[nnew-1] == xpos[npos - 1]);

  /**
  We go through each new interval. */
  
  double x = xnew[0], s_left = 0., C[3];
  fnew[0] = 0.;
  for (int inew = 0, k = -1; inew < nnew - 1 && x < xpos[npos - 1];) {

    /**
    ## Integration interval

    We first check if the starting point (x) of the integration
    interval is beyond the current interval [xpos[k]:xpos[k+1]]. */
    
    if (x >= xpos[k + 1]) {

      /**
      If this is the case, we need to consider the next interval and
      compute the corresponding polynomial coefficients. */
      
      k++;

      /**
      The top layer is a special case. */
      
      if (k == npos - 2) {

        /**
        If limiting is used together with Neumann zero fluxes, the
        value must be constant in the last layer. */
        
        if (limiter && lambda_t == HUGE) {
          C[0] = fdat[k];
          C[1] = 0.;
          C[2] = 0.;
        }

        /**
        Otherwise we use Neumann or Robin
        boundary conditions to compute the polynomial. */
        
        else {
          if (df_t)
            remap_neumann (s_left, fdat[k], - df_t*(xpos[k+1] - xpos[k]), C);
          else
            remap_robin (s_left, fdat[k], f_t, - lambda_t/(xpos[k+1] - xpos[k]), C);

          /**
          Since the functions above are written for the bottom layer, we
          need to apply symmetry conditions i.e. transform $x$ into $1 -
          x$. This gives the following polynomial coefficients. */
        
          C[0] += C[2] + C[1];
          C[1] = - C[1] - 2.*C[2];
        }
      }
      else {

        /**
        For all other layers, we need to compute the right value. If
        limiting is used together with Neumann zero fluxes we impose a
        constant profile in the bottom or top layer. */
        
        double s_right =
          limiter && lambda_b == HUGE && k == 0 ? fdat[0] :
          limiter && lambda_t == HUGE && k == npos - 3 ? fdat[npos - 2] :
          right_value (k, npos, xpos, fdat, f_b, lambda_b, df_b, f_t, lambda_t, df_t);

        /**
        We use boundary conditions for the bottom layer. */
        
        if (k == 0) {
          if (df_b)
            remap_neumann (s_right, fdat[0], df_b*(xpos[1] - xpos[0]), C);
          else
            remap_robin (s_right, fdat[0], f_b, lambda_b/(xpos[1] - xpos[0]), C);
        }

        /**
        Or central remapping, with optional limiting, using the left
        and right values for all the other layers. */
        
        else
          remap_central (s_left, s_right, npos, xpos, fdat, k, C, limiter);

        /**
        The new left value is just the old right value. */
        
        s_left = s_right;
      }
    }

    /**
    ## End of integration interval

    We determine the value of `xe`, the end of the integration
    interval, by comparing the end of the new interval, `xnew[inew+1]`
    with the end of the current interval `xpos[k+1]`. `jnew` records
    whether we need to change interval after the integration. */
    
    double xe;
    int jnew = inew;
    if (xnew[inew + 1] < xpos[k + 1]) {
      jnew = inew + 1;
      xe = xnew[jnew];
      fnew[jnew] = 0.;
    }
    else
      xe = xpos[k + 1];

    /**
    ## Integration

    We compute the integral of the current polynomial between `x` and
    `xe`. Since the polynomial is defined on the normalized interval
    [0:1], we first normalize the values of `x` and `xe`, in `a` and
    `b` respectively. We also need to de-normalize the integral with
    the ratio `dx1`. Note that it is particularly important to perform
    the integration in normalized space to minimize round-off errors
    which can cause instabilities when using single-precision on
    GPUs. */

    double xk = xpos[k], dx = xpos[k+1] - xk, dx1 = dx/(xnew[inew + 1] - xnew[inew]);
    double a = (x - xk)/dx, b = (xe - xk)/dx, an = a, bn = b;
    assert (a >= 0. && a <= b && b <= 1.); // xpos and xnew must be ordered properly
    for (int i = 0; i < 3; i++, an *= a, bn *= b)
      fnew[inew] += C[i]/(i + 1)*(bn - an)*dx1;

    /**
    ## Update of integration bounds

    The new integration interval starts at `xe` and concerns `jnew`. */

    x = xe;
    inew = jnew;
  }
}
