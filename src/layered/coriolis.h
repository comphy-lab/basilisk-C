/**
# Coriolis/friction terms for the multilayer solver

This approximates
$$
\partial_t\mathbf{u} = \mathbf{B}\mathbf{u} + \mathbf{a}
$$
with
$$
  \mathbf{B} = \left( \begin{array}{cc}
    - K_0 & F_0\\
    - F_0 & - K_0
  \end{array} \right),
$$
and $K_0$ and $F_0$ the linear friction and Coriolis parameters
respectively.  The time-implicit discretisation of these terms can be
written
$$
  \frac{\mathbf{u}^{n + 1} -\mathbf{u}^n}{\Delta t} = \mathbf{B} [(1
  - \alpha_H) \mathbf{u}^n + \alpha_H \mathbf{u}^{n + 1}] + \mathbf{a}^n
$$
This then gives
$$
  \mathbf{u}^{n + 1}  (\mathbf{I}- \alpha_H \Delta t\mathbf{B}) =
  \mathbf{u}^n + (1 - \alpha_H) \Delta t\mathbf{B}\mathbf{u}^n + \Delta
  t\mathbf{a}^n
$$
The local $2 \times 2$ linear system is easily inverted analytically. The
final value is obtained by substracting the acceleration i.e.
$$
  \mathbf{u}^{\star} = \mathbf{u}^{n + 1} - \Delta t\mathbf{a}^n
$$

The `K0()` and/or `F0()` macros should be defined before including the
file. */

#ifndef K0
# define K0() 0.
#endif
#ifndef F0
# define F0() 0.
#endif
#ifndef alpha_H
# define alpha_H 0.5
#endif

event acceleration (i++)
{
  foreach()
    foreach_layer()
      if (h[] > dry) {
	coord b0 = { - K0(), - K0() }, b1 = { F0(), -F0() };
	coord m0 = { 1. - alpha_H*dt*b0.x, 1. - alpha_H*dt*b0.y };
	coord m1 = { - alpha_H*dt*b1.x, - alpha_H*dt*b1.y };
	double det = m0.x*m0.y - m1.x*m1.y;
        coord r, a;
	foreach_dimension() {
	  a.x = dt*(ha.x[] + ha.x[1])/(hf.x[] + hf.x[1] + dry);
	  r.x = u.x[] + (1. - alpha_H)*dt*(b0.x*u.x[] + b1.x*u.y[]) + a.x;
	}
	foreach_dimension()
	  u.x[] = (m0.y*r.x - m1.x*r.y)/det - a.x;
      }
}

#undef alpha_H
