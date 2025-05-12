/**
# Advection/diffusion of a soluble tracer

We consider the transport and diffusion of a tracer $c$ with different
solubilities in the two-phases described by
[two-phase.h](/src/two-phase.h).

The diffusion coefficients in the two phases are $D_1$ and $D_2$ and
the jump in concentration at the interface is given by
$$
c_1 = \alpha c_2
$$
The advection/diffusion equation for $c$ can then be written
$$
\partial_t c + \nabla\cdot(\mathbf{u} c) = 
   \nabla\cdot\left(D\nabla c - D \frac{c(\alpha - 1)}
                     {\alpha f + (1 - f)}\nabla f\right)
$$
with $f$ the volume fraction and
$$
D = \frac{D_1 D_2}{D_2 f + D_1 (1 - f)}
$$
the harmonic mean of the diffusion coefficients (see [Haroun et al.,
2010](#haroun2010)).

The diffusion coefficients and solubility are attributes of each
tracer.

The *stracers* list of soluble tracers must be defined by the calling
code. */

attribute {
  double D1, D2, alpha;  //D1 for f = 1, D2 for f = 0
  scalar phi1, phi2; // private
}

extern scalar * stracers;

/**
## Defaults

On trees we need to ensure conservation of the tracer when
refining/coarsening. */

#if TREE
event defaults (i = 0)
{
  for (scalar s in stracers) {
#if EMBED
    s.refine = s.prolongation = refine_embed_linear;
#else
    s.refine  = refine_linear;
#endif
    s.restriction = restriction_volume_average;
    s.dirty = true;
  }
}
#endif // TREE

/**
## Advection

To avoid numerical diffusion through the interface we use the [VOF
tracer transport scheme](/src/vof.h) for the temporary fields
$\phi_1$ and $\phi_2$, see section 3.2 of [Farsoiya et al.,
2021](#farsoiya2021). */

static scalar * phi_tracers = NULL;

event vof (i++)
{
  phi_tracers = f.tracers;
  for (scalar c in stracers) {
    scalar phi1 = new scalar, phi2 = new scalar;
    c.phi1 = phi1, c.phi2 = phi2;
    scalar_clone (phi1, c);
    scalar_clone (phi2, c);
    phi2.inverse = true;
    
    f.tracers = list_append (f.tracers, phi1);
    f.tracers = list_append (f.tracers, phi2);

    /**
    $\phi_1$ and $\phi_2$ are computed from $c$ as
    $$
    \phi_1 = c \frac{\alpha f}{\alpha f + (1 - f)}
    $$
    $$
    \phi_2 = c \frac{1 - f}{\alpha f + (1 - f)}
    $$
    */
		  
    foreach() {
      double a = c[]/(f[]*c.alpha + (1. - f[]));
      phi1[] = a*f[]*c.alpha;
      phi2[] = a*(1. - f[]);
    }
  }
}

/**
## Diffusion

We first define the relaxation and residual functions needed to solve
the implicit discrete system
$$
\frac{c^{n + 1} - c^n}{\Delta t} = 
\nabla\cdot (D \nabla c^{n + 1} + \beta c^{n + 1})
$$
see section 3.2 of [Farsoiya et al., 2021](#farsoiya2021).

Note that these functions are close to that in [poisson.h]() and
[diffusion.h]() but with the additional term $\nabla\cdot (\beta c^{n
+ 1})$. */

struct HDiffusion {
  face vector D;
  face vector beta;
};

static void h_relax (scalar * al, scalar * bl, int l, void * data)
{
  scalar a = al[0], b = bl[0];
  struct HDiffusion * p = (struct HDiffusion *) data;
  face vector D = p->D, beta = p->beta;

  scalar c = a;
  foreach_level_or_leaf (l) {
    double n = - sq(Delta)*b[], d = cm[]/dt*sq(Delta);
    foreach_dimension() {  
      n += D.x[1]*a[1] + D.x[]*a[-1] +
	Delta*(beta.x[1]*a[1] - beta.x[]*a[-1])/2.;
      d += D.x[1] + D.x[] - Delta*(beta.x[1] - beta.x[])/2.;
    }
    c[] = n/d;
  }
}

static double h_residual (scalar * al, scalar * bl, scalar * resl, void * data)
{
  scalar a = al[0], b = bl[0], res = resl[0];
  struct HDiffusion * p = (struct HDiffusion *) data;
  face vector D = p->D, beta = p->beta;
  double maxres = 0.;
#if TREE
  /* conservative coarse/fine discretisation (2nd order) */
  face vector g[];
  foreach_face()
    g.x[] = D.x[]*face_gradient_x (a, 0) + beta.x[]*face_value (a, 0);
  foreach (reduction(max:maxres)) {
    res[] = b[] + cm[]/dt*a[];
    foreach_dimension()
      res[] -= (g.x[1] - g.x[])/Delta;
    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);
  }
#else // !TREE
  /* "naive" discretisation (only 1st order on trees) */
  foreach (reduction(max:maxres)) {
    res[] = b[] + cm[]/dt*a[];
    foreach_dimension()
      res[] -= (D.x[1]*face_gradient_x (a, 1) -
		D.x[0]*face_gradient_x (a, 0) +
		beta.x[1]*face_value (a, 1) -
		beta.x[0]*face_value (a, 0))/Delta;  	  
    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);
  }
#endif // !TREE    
  return maxres;
}

event tracer_diffusion (i++)
{
  free (f.tracers);
  f.tracers = phi_tracers;
  for (scalar c in stracers) {

    /**
    The advected concentration is computed from $\phi_1$ and $\phi_2$ as
    $$
    c = \phi_1 + \phi_2
    $$
    and these fields are then discarded. */
    
    scalar phi1 = c.phi1, phi2 = c.phi2, r[];
    foreach() {
      c[] = phi1[] + phi2[];
      r[] = - cm[]*c[]/dt;
    }
    delete ({phi1, phi2});

    /**
    The diffusion equation for $c$ is then solved using the multigrid
    solver and the residual and relaxation functions defined above. */

    face vector D[], beta[];
    foreach_face() {
      double ff = (f[] + f[-1])/2.;
      D.x[] = fm.x[]*c.D1*c.D2/(c.D1*(1. - ff) + ff*c.D2);
      beta.x[] = - D.x[]*(c.alpha - 1.)/
	(ff*c.alpha + (1. - ff))*(f[] - f[-1])/Delta;
    }
  
    restriction ({D, beta, cm});
    struct HDiffusion q;
    q.D = D;
    q.beta = beta;
    mg_solve ({c}, {r}, h_residual, h_relax, &q);
  }
}

/**
## References

~~~bib
@article{haroun2010,
  title = {Volume of fluid method for interfacial reactive mass transfer: 
           application to stable liquid film},
  author = {Haroun, Y and Legendre, D and Raynal, L},
  journal = {Chemical Engineering Science},
  volume = {65},
  number = {10},
  pages = {2896--2909},
  year = {2010},
  doi = {10.1016/j.ces.2010.01.012}
}

@hal{farsoiya2021, hal-03227997}
~~~
*/
