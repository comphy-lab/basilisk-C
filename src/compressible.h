/**
# Compressible gas dynamics

The Euler system of conservation laws for a compressible gas can be
written

$$
\partial_t\left(\begin{array}{c}
\rho \\
E    \\
w_x  \\
w_y  \\
\end{array}\right) 
+ 
\nabla_x \cdot\left(\begin{array}{c}
w_x \\
\frac{w_x}{\rho} ( E + p ) \\
\frac{w_x^2}{\rho} + p \\
\frac{w_y w_x}{\rho} \\
\end{array}\right)
+ 
\nabla_y \cdot\left(\begin{array}{c}
w_y \\
\frac{w_y}{\rho} ( E + p ) \\
\frac{w_y w_x}{\rho} \\
\frac{w_y^2}{\rho} + p \\
\end{array}\right)
= 0
$$

with $\rho$ the gas density, $E$ the total energy, $\mathbf{w}$ the
gas momentum and $p$ the pressure given by the equation of state

$$
p = (\gamma - 1)(E - \rho\mathbf{u}^2/2)
$$

with $\gamma$ the polytropic exponent. This system can be solved using
the generic solver for systems of conservation laws. */

#include "conservation.h"

/**
The conserved scalars are the gas density $\rho$ and the total energy
$E$. The only conserved vector is the momentum $\mathbf{w}$.  The constant
$\gamma$ is represented by *gammao* here, with a default value of 1.4. */

scalar rho[], E[];
vector w[];
scalar * scalars = {rho, E};
vector * vectors = {w};
double gammao = 1.4 ;

/**
The system is entirely defined by the `flux()` function called by the
generic solver for conservation laws. The parameter passed to the
function is the array `s` which contains the state variables for each
conserved field, in the order of their definition above (i.e. scalars
then vectors). */

void flux (const double * s, double * f, double e[2])
{

  /**
  We first recover each value ($\rho$, $E$, $w_x$ and $w_y$) and then
  compute the corresponding fluxes (`f[0]`, `f[1]`, `f[2]` and
  `f[3]`). */

  double rho = s[0], E = s[1], wn = s[2], w2 = 0.;
  for (int i = 2; i < 2 + dimension; i++)
    w2 += sq(s[i]);
  double un = wn/rho, p = (gammao - 1.)*(E - 0.5*w2/rho);

  f[0] = wn;
  f[1] = un*(E + p);
  f[2] = un*wn + p;
  for (int i = 3; i < 2 + dimension; i++)
    f[i] = un*s[i];

  /**
  The minimum and maximum eigenvalues for the Euler system are the
  characteristic speeds $u \pm \sqrt(\gamma p / \rho)$. */

  double c = sqrt(gammao*p/rho);
  e[0] = un - c; // min
  e[1] = un + c; // max
}
