/**
# Particle movements

This header file implements a generic particle-path solving
definition:

$$\frac{\mathrm{d}\mathbf{x}_p}{\mathrm{d}t} = \mathbf{v}_p.$$

Unlike [flow-tracer particles](tracer-particles.h), the intertial
particle velocity is governed by the definition of its acceleration
vector $\mathbf{a}_p$.

$$\frac{\mathrm{d}\mathbf{v}_p}{\mathrm{d}t} = \mathbf{a}_p.$$

For this purpose, a function can be specified following the prototype
for `p_acc()`. Furthermore, particles must atleast have some additonal
members to facilitate 2nd-order acurate advection *and* we aim to
achieve compatibility with [tracer-particles.h](), which should then
be included first.
*/
#include "run.h"
#ifndef ADD_PART_MEM 
#define ADD_PART_MEM coord u; long unsigned int tag; // u -> v_p
#endif

#include "particle.h"

#ifndef PA_INP
#define PA_INP //empty
#define INP (PA_inp){p()}
#endif

typedef struct PA_inp {
  particle p;
  PA_INP
} PA_inp;

coord p_acc (PA_inp inp);

Particles * inertial_particles = NULL;

extern scalar * _automatics_;
/**
## Time integration

A version of the second-order accurate [Verlet
method](https://en.wikipedia.org/wiki/Verlet_integration) is used.
(a.k.a. leaffrogging). Since `p_acc` may depend on `p.u`, we advance
the particles in two solver timesteps.

$$\mathbf{a}_{p,i} = \mathtt{p_{acc}}\left(\mathtt{p_i \ ADD\_ARGS}\right),$$

$$\mathbf{v}_{p, i + 1} = \mathbf{v}_{p, i - 1} +
\mathbf{a}_{p,i}\times \ 2 \mathtt{dt_i},$$

$$\mathbf{x}_{p, i + 2} = \mathbf{x}_{p, i} + \mathbf{v}_{p,i + 1}
\times \left( \mathtt{dt_i + dt_{i+1}}\right).$$

The following events take care of the time integration.
 */

double dtprev;

event defaults (t = 0);

event init (t = 0);

event set_dtmax (i++);

event velocity (i++);

void ip_step1 (Particles P, int i) {
  double Dt = i > 0 ? 2*dt : dt;
  boundary (_automatics_);
  foreach_particle_in(P) {
    coord a = p_acc (INP);
    foreach_dimension()
      p().u.x += a.x*Dt;
  }
}

void ip_step2 (Particles P) {
  foreach_particle_in(P) {
    foreach_dimension()
      p().x += p().u.x*(dt + dtprev);
  }
}

event intertial_particles_step1 (i += 2, last) {
  dtprev = dt;
  foreach_P_in_list(inertial_particles) {
    particle_boundary (P);
    ip_step1 (P, i);
  }
}

event intertial_particles_step2 (i = 1; i += 2) {
  foreach_P_in_list(inertial_particles) 
    ip_step2 (P);
}

event free_intertial_particles (t = end) {
  if (inertial_particles)
    free (inertial_particles);
  inertial_particles = NULL;
}

/**
## Utilities 

From [tracer-particles.h](), a bunch of functions are copied to
provide some user interface for particle initializaiton.
 */

Particles new_inertial_particles (long unsigned int n) {
  Particles p = new_particles (n);
  int l = 0, t = 0;
  if (inertial_particles != NULL) {
    while (pn[l] != terminate_int) {
      if (l == inertial_particles[t]) 
	t++;
      l++;
    }
  }
  inertial_particles = realloc (inertial_particles, (t + 1)*sizeof(Particles));
  inertial_particles[t] = p;
  return p;
}

void tag_ip_particles (Particles p) {
  long unsigned int offset = 0;
#if _MPI
  long unsigned int pntr[npe()], tag_start[npe()];
  MPI_Gather (&pn[p], 1, MPI_UNSIGNED_LONG,
	      pntr, 1, MPI_UNSIGNED_LONG,
	      0, MPI_COMM_WORLD);
  if (pid() == 0) {
    tag_start[0] = 0;
    for (int tr = 1; tr < npe(); tr++) 
      tag_start[tr] = tag_start[tr - 1] + pntr[tr - 1];
  }
  MPI_Scatter (tag_start, 1, MPI_UNSIGNED_LONG,
	       &offset, 1, MPI_UNSIGNED_LONG,
	       0, MPI_COMM_WORLD);
#endif
  foreach_particle_in(p)
    p().tag = _j_particle + offset;
}

void set_ip_particle_attributes (Particles p) {
  tag_ip_particles (p);
}

/**
### Inertial-particle (`ip`) initialization

Inertial-particles may be initialized by calling the following
functions.
 */

Particles init_ip_cells(void) {
  int np = 0;
  foreach(noauto)
    np++;
  Particles p = new_inertial_particles (np);
  place_in_cells(p);
  particle_boundary (p); //?
  set_ip_particle_attributes (p);
  return p;
}

Particles init_ip_square (int n = 10,
			  double xm = 0, double ym = 0,
			  double l = L0, double zp = 0)
{
  Particles p;
  if (pid() == 0) {
    p = new_inertial_particles (sq(n));
    place_in_square (p, n, xm, ym, l, zp);
  } else { 
    p = new_inertial_particles (0);
  }
  particle_boundary (p);
  set_ip_particle_attributes (p);
  return p;
}

Particles init_ip_circle (int n = 100,
			  double xm = 0, double ym = 0,
			  double l = L0/2., double zp = 0)
{
  Particles p;
  if (pid() == 0) {
    p = new_inertial_particles (n);
    place_in_circle (p, n, xm, ym, l, zp);
  }
  else
    p = new_inertial_particles (0);
  particle_boundary (p);
  set_ip_particle_attributes (p);
  return p;
}

/**
## Todo

* More Tests  
* ~~~Implement usages~~~

## Tests

* [Solving $\ddot{x} = -x$](tip.c)
* [Combine `tracer_particles` with `inertial_particles`](pc.c)

## Usage

* [Cosmological particles](cosmology.h)
* [Spherical particle-laden flow](stokes-particles.h)

*/
