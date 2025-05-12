/**
A copy of Antoon's [tracer-particle.h](TODO),
but adapted to work with Smoothed/Fluid Particle method.
I basically turn the particle advection terms.
*/
//#define Particle_Advection_ON 1

/**
# Langrangian particles

Tracer-particle coordinates $\mathbf{p_i}$ are advected by a vector field
$\mathbf{u}$: 

$$\frac{\mathrm{d}\mathbf{p_i}}{\mathrm{d}t} = \mathrm{u}|_{\mathbf{x} = \mathbf{p_i}}$$

## Definititions and types 

The generic `particle` class is used and tuned: For the
time-advancement scheme, additional storage for velocities and
locations is required. Furthermore, particles may be `tag`ged. The
tracer particles are stored in an array of groups of `particle`s. Note 
that `Particles` reference a list of `particle`s. 
 */
extern vector u;
extern scalar omega; // vorticity field

/*
# Hypothesis: Rigid Body Motion
Using variable viscosity, within _any_ point of the particle I've got some kind of RMB
Hence, I can interpolate the velocity at the particle's location,
making the rest of the treatment much more simple!
*/
#ifndef ADD_PART_MEM
//#define ADD_PART_MEM  coord u; coord u2; coord locn; long unsigned int tag; 
#define ADD_PART_MEM  coord u; coord u2; coord locn; long unsigned int tag; \
double theta; double thetan; double omega; double omega2; double r; \
double areaCenter; double areaTrans; double areaCis; 

#endif
#include "Antoonvh/particle.h"
Particles * tracer_particles = NULL;
bool P_RK2 = false;

/**
## Time integration  
 
The used particle-advection scheme is the RK-3 variant of Sanderse and
Veldman (2019) for $\alpha \neq \frac{2}{3} \pm \mathtt{tol}$:

$$ \begin{array}{c|ccc} 
0 & 0 & & \\ 
\alpha & \alpha & & \\ 
1 &1+\frac{1- \alpha}{\alpha (3\alpha -2)} & -\frac{1- \alpha}{\alpha
(3\alpha -2)} &\\ 
\hline
 & \frac{1}{2}-\frac{1}{6\alpha} &
\frac{1}{6\alpha(1-\alpha)} & \frac{2-3\alpha}{6(1-\alpha)} \\
\end{array} $$

see:  

B Sanderse and AEP Veldman, *Constraint-consistent Runge-Kutta
methods for one-dimensional incompressible multiphase flow*,
J. Comp. Phys. (2019) [link to
JCP.](https://www.sciencedirect.com/science/article/pii/S0021999119300683).

This scheme is used to advance the solution in two solver iterations
between $t_{n-2}$ and $t_n$. In case $\alpha \approx \frac{2}{3}$,
*or* `P_RK2 == true`, we skip the last stage and switch to a generic
second-order-accurate scheme.
 */
double tol = 3.e-2;
double dtf[2]; //timesteps
/**
The user is supposed to add particles in some relevant event, perhaps
called `init`. 
 */
event defaults (t = 0);

event init (t = 0);

event set_dtmax (i++);

event velocity (i++);
/**
By default, second-order-in-space "linear" interpolation is used to
approximate $\mathrm{u}|_{\mathbf{x} = \mathbf{p_i}}$. This could be
overloaded with [higher-order](higher-order.h) methods.

e.g for a 3rd order accurate method:

~~~literatec
...
#include "higher-order.h" //This is in /sandbox/Antoonvh/
#define interpolate_linear interpolate_quadratic
#include "tracer-particles.h" 
...
~~~

Velocities may found and stored in `p().u` or `p().u2`.
 */
void interpolate_vel (Particles p, int swtch) {
  if (swtch == 1) {
    foreach_particle_in(p) {
      if (locate (x, y, z).level < 0) {//freeze
	fprintf (stderr, "#pid: %d could not locate p.tag: %ld.\n",
		 pid(), p().tag); //This should not happen
	foreach_dimension()
	  p().u.x = 0;
	p().omega = 0;
      } else {
	foreach_dimension()
	  p().u.x =interpolate_linear (locate (x, y, z),
				       u.x, x, y, z);
	p().omega = OMEGA_EQUATION();
      }
    }
  } else {
    foreach_particle_in(p) {
      if (locate (x, y, z).level < 0) {//freeze
	fprintf (stderr, "#pid: %d could not locate p.tag: %ld.\n",
		 pid(), p().tag); //This should not happen
	foreach_dimension()
	  p().u2.x = 0;
	p().omega2 = 0;
      } else {
	foreach_dimension()
	  p().u2.x =interpolate_linear (locate (x, y, z),
					u.x, x, y, z);
	p().omega2 = OMEGA_EQUATION(); //interpolate_linear (locate (x, y, z),omega, x, y, z);
      }
    }
  }
}
/**
The three Runge-Kutta stages are defined below:
 */
void RK_step1 (Particles p, double dtf[2]) {
  interpolate_vel(p, 1);
  foreach_particle_in(p) { 
    foreach_dimension() {
      p().locn.x = p().x; //store the locations at t_n 
      p().x += dtf[0]*p().u.x;
    }
		p().thetan = p().theta;
		p().theta+=dtf[0]*p().omega; // simple Explicit Euler
  }
}

void RK_step2 (Particles p, double dtf[2]) {
  interpolate_vel(p, 2);
  double a1 = -1, a2 = 2., h = dtf[1] + dtf[0];
  if (dtf[1] != dtf[0] || P_RK2) {
    double c = dtf[0]/h;
    if (fabs (c - 2./3.) > tol && !P_RK2)
      a2 = (c - 1.)/(c*(3.*c - 2.));
    else // Raltson's 2nd order method
      a2 = 1./(2*c);
    a1 = 1 - a2;
  }
  foreach_particle_in(p) {
    foreach_dimension() 
      p().x = p().locn.x + h*(a1*p().u.x + a2*p().u2.x);
    p().theta = p().thetan + h*(a1*p().omega + a2*p().omega2);
  }
}

void RK_step3 (Particles p, double dtf[2]) {
  double h = dtf[1] + dtf[0];
  double c = dtf[0]/h;
  if (fabs(c - 2./3.) > tol) {  // RK-3
    double b1 = 0.5 - 1./(6.*c);
    double b2 = 1./(6.*c*(1. - c));
    double b3 = 1. - (b1 + b2);
    foreach_particle_in(p) {
      foreach_dimension() 
        p().x = p().locn.x + h*(b1*p().u.x + b2*p().u2.x +
				b3*interpolate_linear (locate (x, y, z),
				 u.x, x, y, z));
      p().theta = p().thetan + h*(b1*p().omega + b2*p().omega2 +
			b3*OMEGA_EQUATION());
			//b3*interpolate_linear (locate (x, y, z),
			//omega, x, y, z));
    }
  }
}
/**
The particles are advanced in time in various RK-stage events. Note
that the particle coordinates (`p().x`) only contain the proper
approximated location inbetween `step3` and `step1`. These locations
are otherwise stored in `p().locn`. 
*/
#ifdef Particle_Advection_ON // run this section...
event tracer_particles_step3 (i += 2, last);
event tracer_particles_step3 (i = 2; i += 2) {
  if (!P_RK2) {
    foreach_P_in_list (tracer_particles) {
      particle_boundary (P);
      RK_step3 (P, dtf);
    }
  }
}

event tracer_particles_step1 (i += 2, last) {
  dtf[0] = dt;
  foreach_P_in_list (tracer_particles) {
    particle_boundary (P);
    RK_step1 (P, dtf);
  }
}

event tracer_particles_step2 (i += 2, last);
event tracer_particles_step2 (i = 1; i += 2) {
  dtf[1] = dt;
  foreach_P_in_list (tracer_particles) {
    particle_boundary (P);
    RK_step2 (P, dtf);
    
  }
}
#else
event compute_velocity_without_advecting_location (i++) {
  foreach_P_in_list (tracer_particles) {
    particle_boundary (P);
  		interpolate_vel(P, 1);
  }
}
#endif

event free_tracers (t = end) {
  free (tracer_particles);
  tracer_particles = NULL;
}

/**
## Initializing tracer particles 

`n` New tracer particles can be declared by calling
`new_tracer_particles (n)`. This is a building-block of
tracer-particle initialization.
*/
Particles new_tracer_particles (long unsigned int n) {
  Particles p = new_particles (n);
  int l = 0, t = 0;
  if (tracer_particles != NULL) {
    while (pn[l] != terminate_int) {
      if (l == tracer_particles[t]) 
	t++;
      l++;
    }
  }
  tracer_particles = realloc (tracer_particles, (t + 1)*sizeof(Particles));
  tracer_particles[t] = p;
  return p;
}

void tag_particles (Particles p) {
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
/**
We need a ugly hack to access `ADD_PART_MEMBERS`, as they only become
available after the #define macros are processed.
 */
@define part p 
void set_a_particle_attributes (particle * p) {
  double x = p->x, y = p->y, z = p->z;
  foreach_dimension() {
    part->locn.x = p->x;
    part->u.x = interpolate (u.x, x, y, z);
    part->u2.x = 0.;
  }
}

void set_particle_attributes (Particles p) {
  interpolate_vel (p, 1);
  foreach_particle_in(p) {
    foreach_dimension() {
      p().locn.x = p().x;
      //p().u2.x = p().u.x; //Do or dont?
    }
  }
  tag_particles (p);
}

/**
### Tracer-particle (`tp`) initialization

Tracer-particles may be initialized by calling the following
functions.
 */

Particles init_tp_cells(void) {
  int np = 0;
  foreach(serial)
    np++;
  Particles p = new_tracer_particles (np);
  place_in_cells(p);
  particle_boundary (p); //?
  set_particle_attributes (p);
  return p;
}

Particles init_tp_square (int n = 10,
			  double xm = 0, double ym = 0,
			  double l = L0, double zp = 0)
{
  Particles p;
  if (pid() == 0) {
    p = new_tracer_particles (sq(n));
    place_in_square (p, n, xm, ym, l, zp);
  } else { 
    p = new_tracer_particles (0);
  }
  particle_boundary (p);
  set_particle_attributes (p);
  return p;
}

Particles init_tp_circle (int n = 100,
			  double xm = 0, double ym = 0,
			  double l = L0/2., double zp = 0)
{
  Particles p;
  if (pid() == 0) {
    p = new_tracer_particles (n);
    place_in_circle (p, n, xm, ym, l, zp);
  }
  else
    p = new_tracer_particles (0);
  particle_boundary (p);
  set_particle_attributes (p);
  return p;
}
