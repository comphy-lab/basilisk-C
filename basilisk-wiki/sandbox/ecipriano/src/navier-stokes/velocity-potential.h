/**
# Velocity Potential Model

We want to obtain a divergence-free velocity extrapolation from the
field velocity with phase change. A possible approach, proposed in
different works with small differences ([Scapin et al. 2020](#scapin2020volume),
[Malan et al. 2021](malan2021geometric), [Palmore et al. 2019](#palmore2019volume)),
consists in solving an additional Poisson equation that allows a velocity
potential $\phi$ to be computed:
$$
\nabla \cdot \left( \alpha \nabla \phi \right)
=
\dfrac{\dot{m}}{\Delta t} \left( \frac{1}{\rho_g} - \frac{1}{\rho_l} \right) \delta_\Gamma \\
$$
The stefan velocity can be obtained from the velocity potential as:
$$
\mathbf{u}^S = -\Delta t \alpha \nabla \phi
$$
We then calculate the extended velocity by subtracting the stefan velocity
from the field velocity. The resulting extended velocity field will be
divergence-free by construction:
$$
\mathbf{u}^E = \mathbf{u} - \mathbf{u}^S
$$
*/

/**
## Fields

We initialize the velocity potential field *ps*, the stefan velocity *ufs*,
and the divergence-free extended velocity *ufext*. */

scalar ps[];
face vector ufs[], ufext[];
mgstats mgpsf;

/**
The volume expansion term is declared in
[evaporation.h](/sandbox/ecipriano/src/evaporation.h). */

extern scalar stefanflow;

/**
## Helper functions

We define the function that performs the projection
of the stefan velocity onto a field with divergence
equal to the volume expansion term. */

trace
mgstats project_sv (face vector uf, scalar p,
     (const) face vector alpha = unityf,
     double dt = 1.,
     int nrelax = 4)
{
  scalar div[];
  foreach()
    div[] = stefanflow[]/dt;

  mgstats mgp = poisson (ps, div, alpha,
      tolerance = TOLERANCE/sq(dt), nrelax = nrelax);

  foreach_face()
    ufs.x[] = -dt*alpha.x[]*face_gradient_x (ps, 0);
  boundary((scalar*){ufs});

  return mgp;
}

/**
## Boundary conditions

It is important to impose for *ps* the same boundary conditions and the
same tolerance used in the Poisson equation for the pressure *p*. */

ps[right] = neumann (neumann_pressure(ghost));
ps[left]  = neumann (- neumann_pressure(0));

#if AXI
ufs.n[bottom] = 0.;
ufs.t[bottom] = dirichlet(0);
ps[top]    = neumann (neumann_pressure(ghost));
#else // !AXI
#  if dimension > 1
ps[top]    = neumann (neumann_pressure(ghost));
ps[bottom] = neumann (- neumann_pressure(0));
#  endif
#  if dimension > 2
ps[front]  = neumann (neumann_pressure(ghost));
ps[back]   = neumann (- neumann_pressure(0));
#  endif
#endif // !AXI

/**
## Extended velocity

We perform the projection of the stefan velocity
by overloading the event end_timestep, defined in
[navier-stokes/centered.h](/src/navier-stokes/centered.h). */

event end_timestep (i++, last)
{
  /**
  We set to zero the face-centered stefan velocity *ufs*
  and the velocity potential *ps*. */

  foreach_face()
    ufs.x[] = 0.;

  foreach()
    ps[] = 0.;

  /**
  We solve the Poisson equation using the multigrid solver. */

  mgpsf = project_sv (ufs, ps, alpha, dt, mgpsf.nrelax);

  /**
  We compute a divergence-free extended velocity by subtracting
  the stefan velocity from the field velocity. */

  foreach_face()
    ufext.x[] = uf.x[] - ufs.x[];
}

/**
## Notes
This approach works fine when the field velocity is larger than
the velocity due to the phase change. An example is a falling
droplet, or a droplet in forced convective conditions. Static
droplets, evaporating in reduced gravity conditions, can be
simulated using this method just for small vaporization rates
or for small density ratio values. If these conditions are not
met, the double pressure velocity coupling approach should be
preferred.

## References

~~~bib
@article{scapin2020volume,
  title={A volume-of-fluid method for interface-resolved simulations of phase-changing two-fluid flows},
  author={Scapin, Nicol{\`o} and Costa, Pedro and Brandt, Luca},
  journal={Journal of Computational Physics},
  volume={407},
  pages={109251},
  year={2020},
  publisher={Elsevier}
}

@article{malan2021geometric,
  title={A geometric VOF method for interface resolved phase change and conservative thermal energy advection},
  author={Malan, LC and Malan, Arnaud G and Zaleski, St{\'e}phane and Rousseau, PG},
  journal={Journal of Computational Physics},
  volume={426},
  pages={109920},
  year={2021},
  publisher={Elsevier}
}

@article{palmore2019volume,
  title={A volume of fluid framework for interface-resolved simulations of vaporizing liquid-gas flows},
  author={Palmore Jr, John and Desjardins, Olivier},
  journal={Journal of Computational Physics},
  volume={399},
  pages={108954},
  year={2019},
  publisher={Elsevier}
}
~~~

*/


