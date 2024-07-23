/**
# Runge--Kutta time integrators
*/

static double update (scalar * ul, scalar * kl, double t, double dt,
		      void (* Lu) (scalar * ul, double t, scalar * kl),
		      scalar * dul, double w)
{
  scalar * u1l = list_clone (ul);
  foreach() {
    scalar u1, u, k;
    for (u1,u,k in u1l,ul,kl)
      u1[] = u[] + dt*k[];
  }
  Lu (u1l, t + dt, kl);
  foreach() {
    scalar du, k;
    for (du,k in dul,kl)
      du[] += w*k[];
  }
  delete (u1l), free (u1l);
  return w;
}

/**
The *runge_kutta()* function implements the classical first- (Euler),
second- and fourth-order Runge--Kutta time integrators for evolution
equations of the form
$$
\frac{\partial\mathbf{u}}{\partial t} = L(\mathbf{u}, t)
$$
with $\mathbf{u}$ a vector (i.e. list) of evolving fields and $L()$ a
generic, user-defined operator.

Given $\mathbf{u}$, the initial time *t*, a timestep *dt* and the
function $L()$ which should fill *kl* with the right-hand-side of the
evolution equation, the function below will return $\mathbf{u}$ at
time $t + dt$ using the Runge--Kutta scheme specified by *order*. */

void runge_kutta (scalar * ul, double t, double dt,
		  void (* Lu) (scalar * ul, double t, scalar * kl),
		  int order)
{
  scalar * dul = list_clone (ul);
  scalar * kl = list_clone (ul);
  Lu (ul, t, kl);
  foreach() {
    scalar du, k;
    for (du,k in dul,kl)
      du[] = k[];
  }

  double w = 1.;
  switch (order) {
  case 1: // Euler
    break;
  case 2:
    w += update (ul, kl, t, dt, Lu, dul, 1.);
    break;
  case 4:
    w += update (ul, kl, t, dt/2., Lu, dul, 2.);
    w += update (ul, kl, t, dt/2., Lu, dul, 2.);
    w += update (ul, kl, t, dt,    Lu, dul, 1.);
    break;
  default:
    assert (false); // not implemented
  }
  
  foreach() {
    scalar u, du;
    for (u,du in ul,dul)
      u[] += dt/w*du[];
  }

  delete (dul), free (dul);
  delete (kl), free (kl);
}
