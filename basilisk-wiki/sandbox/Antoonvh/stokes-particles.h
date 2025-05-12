/**
# Stokes' drag law for spherical particles

The acceleration of a small particle $\mathbf{a}_p$ in a fluid medium
can be described by,

$$\mathbf{a}_p = \frac{\mathbf{F}_p}{m_p}$$

where $F_p$ are the forces acting on the particle and $m_p$ is the
mass of the particle. For spherical particles, we consider a drag
force due the the velocity difference of the particle and the fluid
medium.

$$\mathbf{F}_d = C_d r_p \mu\left(\mathbf{u}_f - \mathbf{v}_p\right),$$

where $C_d = 6\pi$, the Stokes drag coefficient, $r_p$ the particle
radius , $\mu$ is the fluid medium's dynamic viscosity and
$\left(\mathbf{u}_f - \mathbf{v}_p\right)$ is the volocity difference
between the undisturbed fluid flow and the particle. Furthermore, the
acceleration of gravity (`G`) is considered to act on the particle.

$$\mathbf{F}_p = \mathbf{F}_d + \mathbf{F}_G,$$

$$\mathbf{F}_G = G \times V_p \left(\rho_p - \rho_f\right),$$

with $V_p$ the volume of the particle. This gives a formulation for
one-way coupled partile-laden flow. It can be casted in an inertial
particle formalism.

## Implementation

It can be interesting to compute "`inertial-particles`" paths,
alongside [flow tracer-particles](tracer-particles.h). Therfore, it is
chosen to reuse the names of tracer-particle members for the mass
density (`u2.x`) and the radius (`u2.y`). Note if
[tracer-particles.h]() are also used, that headerfile should be
included first. Finally, if `u2.y` is not set `u2.z` can be used to
set the relaxation time scale $\tau_p$. By default it is computed from
using `u2.y` as $r_p$

$$\tau = \frac{\rho_p 4r_p^2}{18\mu}$$
*/
#ifndef ADD_PART_MEM
#define ADD_PART_MEM coord u; coord u2; long unsigned int tag;
#endif

#define PA_INP particle * pp; double dt;
#define INP    (PA_inp){p(), &p(), 2.*dt}

#include "inertial-particles.h"

				    /**
For two-way couling, the forces on the particles are stored in a
field.
				    */

extern vector u;        //Fluid medium flow
extern face vector mu;  //Fluid medium dynamic viscosity
extern scalar rho;      //Fluid medium density
coord G;                //Gravity acceleration vector (should it be external) ?
scalar * _automatics_ = NULL;

coord p_acc (PA_inp inp) {
  particle pa = inp.p;
  double xp = pa.x; double yp = pa.y; double zp = pa.z;  
  /**
### Implicity time discretization for $v_p$

Since the flow-adjustment timescale $\tau$ can be smaller than the
solver time step, it makes sense to use an implicit equation for
$v_p$,

$$\frac{v_p^{n+1}-v_p^n}{\mathtt{dt}} = \frac{1}{\tau}\left(u - v_p^{n+1}\right) + g',$$

which yields a weighted averaged between $v_p^n$ and
the ternimal velocity for $v_p^{n+1}$.

$$v_p^{n+1} =\frac{\mathtt{dt}^{-1}v_p^n  + \tau^{-1}u + g'}{\mathtt{dt}^{-1} + \tau^{-1}}.$$

Rather than computing the acceleration, the velocity is directly
updated in this function. 
   */
  double muc, tau, itau, rhof; //computed in `foreach_point`
  double rp = pa.u2.y;
  double rhop = pa.u2.x;
  double idt  = inp.dt > 0 ? 1./inp.dt: HUGE;
  foreach_point(xp, yp, zp, serial) {
    muc = is_constant (mu.x) ? constant(mu.x) : interpolate_linear (point, mu.x, xp, yp, zp);
    rhof = is_constant(rho) ? constant(rho) : interpolate_linear (point, rho, xp, yp, zp);
    tau  = rp ? (muc > 0 ? rhop * sq(2*rp)/(18*muc) : HUGE) : pa.u2.z ? pa.u2.z : HUGE;
    itau = tau > 0 ? 1/tau : HUGE;
    foreach_dimension()
      inp.pp->u.x = (idt*pa.u.x + itau*interpolate_linear(point, u.x, xp, yp ,zp))/(idt + itau);
  }
  // Gravity
  foreach_dimension() 
    inp.pp->u.x += rhop ? (G.x*(rhop - rhof)/rhop)/(idt + itau) : 0;
  // No additional acceleration
  return (coord){0., 0., 0.};
}

#if TWO_WAY
event defaults (i = 0) {  
  if (is_constant(a.x)) {
    a = new face vector;
    foreach_face()
      a.x[] = 0.;
    boundary ((scalar *){a});
  }
  foreach_dimension()
    _automatics_ = list_add (_automatics_ , u.x);
  if (!is_constant (rho))
    _automatics_ = list_add (_automatics_, rho);
  if (!is_constant (mu.x)) {
    foreach_dimension()
      _automatics_ = list_add (_automatics_, mu.x);
  }
}
// Fool(?) qcc, that we are not nesting foreach's...
void actually_add_particle (particle p, vector F) {
  double x = p.x, y = p.y, z = p.z;
  foreach_point(x,y,z, serial) {
    double muc = is_constant (mu.x) ? constant(mu.x) : interpolate_linear (point, mu.x, x, y, z);
    double rp = p.u2.y;
    double pref = 6*pi*muc*rp;
    double rhop = p.u2.x;
    foreach_dimension()
      F.x[] += pref*(p.u.x - interpolate_linear(point, u.x, x, y, z));
  }
}

void add_force (Particles P, vector F) {
  particle_boundary (P);
  double rhof = constant(rho);
  boundary ((scalar*){u});
  foreach_particle_in(P) {
    actually_add_particle (p(), F);
    //pref = 4./3.*pi*cube(rp)*(rhop-rhof);
  }
}

event acceleration (i++) {
  face vector av = a;
  vector Fp[];
  foreach() 
    foreach_dimension()
      Fp.x[] = 0;
  foreach_P_in_list (inertial_particles)
    add_force (P, Fp);  
  double rhof = constant(rho);
  boundary ((scalar*){Fp});
  foreach_face() 
    av.x[] += face_value(Fp.x, 0)/(dv()*rhof);
}
#endif

/**
## Tests

* [Sanity check for a sand granule in water](settling.c)
* [Particle in a Taylor-Green vortex](dagan_fig2a.c)

## Usage  

* [Particles in 2D turbulence](pc.c)
* [A particle-driven flow](ash.c)

## Todo

* ~~~Examplify with usages~~~
* ~~~Two-way coupling with the [centered](/src/navier-stokes/centered.h) solver.~~~
* ~~~Test if `inertial-particles` can indeed co-exist with `tracer-particles`~~~
* Critical tests  
*/
