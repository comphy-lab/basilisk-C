/**
# Double projection

This option for the [centered Navier--Stokes solver](centered.h) is
inspired by [Almgren et al., 2000](#almgren200). These authors first
recall that, while for the exact projection method, the definition of
the pressure is unambiguous, in the case of the *approximate*
projection method (used in [centered.h]()), several pressures can be
defined.

Exploiting the different properties of these different definitions may
be useful in some cases.

## Standard approximate projection

We first summarise the standard approximate projection scheme, as
implemented in [centered.h]().

1. Advection/viscosity
$$
 \frac{u^{\star} - u^n}{\Delta t} = (u \cdot \nabla u)^{n + 1 / 2} +
	\alpha \nabla \cdot (\mu \nabla D)^{\star}
$$
2. Acceleration
$$
  u^{\star}_f = \overline{u^{\star}} + \Delta ta_f
$$
3. Projection
$$
  u^{n + 1}_f = P (u^{\star}_f, p^{n + 1})
$$
4. Centered pressure gradient correction
$$
  g^{n + 1} = \overline{a_f - \alpha_f \nabla p^{n + 1}}
$$
$$
  u^{n + 1} = u^{\star} + \Delta tg^{n + 1}
$$

with $P(u,p)$ the projection operator and the overline designating 
either cell-center to cell-face averaging or reciprocally.

## Double approximate projection

As its name indicates this scheme adds an extra projection step and
defines two pressures: the standard one used to project the
face-centered velocity field $u_f$ (renamed $\mathbf{\delta p}$
below), and $p^{n+1}$ used to compute the cell-centered pressure
gradient $g$. This new pressure is obtained by projecting only the
*update* (i.e. its evolution in time) to the centered velocity
field. Note that in the case of an exact projection these two
projections are identical since the divergence of the velocity field
at the start of the timestep is zero.

The scheme can be summarised as:

1. Advection/viscosity
$$
 \frac{u^{\star} - u^n}{\Delta t} = (u \cdot \nabla u)^{n + 1 / 2} +
	\alpha \nabla \cdot (\mu \nabla D)^{\star} 
	\mathbf{+ g^n = A^{n+1/2} + g^n}
$$
$$
 \mathbf{A_f = \overline{A^{n+1/2}} + \Delta ta_f}
$$
2. Acceleration
$$
  u^{\star}_f = \overline{u^{\star}}
$$
3. Projection
$$
  u^{n + 1}_f = P (u^{\star}_f, \mathbf{\delta p})
$$
4. Approximate projection
$$
  u^{n + 1} = u^{\star} - \Delta t \overline{\alpha_f \nabla \delta p}
$$
5. Second projection
$$
\mathbf{P(A_f,p^n+1)}
$$
$$
  g^{n + 1} = \overline{a_f - \alpha_f \nabla p^{n + 1}}
$$

where the additions to the previous scheme are highlighted in bold.

Why is this useful? The new pressure does not feel the *history* of
divergence of the centered velocity field. This is useful in
particular when this history includes the noise induced by adaptive
mesh refinement.

The cost to pay is however significant since an extra (potentially
expensive) projection is required.

## Implementation

We need a (temporary) field to store the update $A_f$ and a
(permanent) field to store the projection pressure $\delta p$. */

face vector Af;
scalar dp[];

/**
We make heavy use of the [event inheritance mechanism](/Basilisk
C#event-inheritance). All the events below are first defined in
[centered.h]().

At the beginning of the timestep (i.e. before advection), we store the
(interpolated) value of the initial velocity field in $A_f$. */

event advection_term (i++)
{
  Af = new face vector;
  foreach_face()
    Af.x[] = - fm.x[]*face_value (u.x, 0);
}

/**
After the advection and diffusion terms have been added to $u$, we
recover the update by adding the new face-interpolated value of the
velocity field to the initial face velocity, and add the acceleration
i.e. we perform step 1 above:
$$
 A_f = \overline{A^{n+1/2}} + \Delta ta_f
$$
*/

face vector ab;

event acceleration (i++)
{
  foreach_face()
    Af.x[] += fm.x[]*(face_value (u.x, 0) + dt*a.x[]);

  /**
  We also add the centered gradient $\mathbf{g^n}$ to the centered
  velocity field. */
  
  correction (dt);

  /**
  Step 2 above (i.e. $u_f^{\star} = \overline{u^{\star}}$) is
  performed by the acceleration event of the [centered
  solver](centered.h#acceleration), but we need to reset the
  acceleration to zero. */

  ab = a;
  a = zerof;
}

/**
The projection step 3 is also performed by the [centered
solver](centered.h#projection). We want to store the resulting
pressure in $dp$ rather than $p$, so we swap the two fields, before
performing the projection. */

event projection (i++)
{
  scalar_clone (dp, p);
  swap (scalar, p, dp);
}

/**
Step 4 is performed by the [centered
solver](centered.h#projection). Step 5 is done at the end of the
timestep. We first restore the fields modified above, then perform the
second projection and compute the corresponding centered gradient
$g^{n+1}$. */

event end_timestep (i++)
{
  swap (scalar, p, dp);
  a = ab;
  // this could be optimised since we do not use Af
  mgp = project (Af, p, alpha, dt, mgp.nrelax);
  delete ((scalar *){Af});
  centered_gradient (p, g);
}

/**
## References

~~~bib
@article{almgren2000,
  title={Approximate projection methods: Part I. Inviscid analysis},
  author={Almgren, Ann S and Bell, John B and Crutchfield, William Y},
  journal={SIAM Journal on Scientific Computing},
  volume={22},
  number={4},
  pages={1139--1159},
  year={2000},
  publisher={SIAM},
  url={https://epubs.siam.org/doi/pdf/10.1137/S1064827599357024}
}
~~~
*/
