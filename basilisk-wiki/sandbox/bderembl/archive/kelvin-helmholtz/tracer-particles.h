/**
# Langrangian particles

Tracer particles can be advected by a vector field `u`. 

## how?  

The generic `particle` class is used. For the third-order accurate
time-advancement scheme, additional storage for velocities and
locations is required. Furthermore, particles may be `tag`ged.
 */
extern vector u;
#ifndef ADD_PART_MEM
#define ADD_PART_MEM  coord u; coord u2; coord locn; long unsigned int tag; 
#endif
#include "particle.h"
#include "run.h"
Particles * tracer_particles = NULL;
/**
`n` New tracer particles can be declared by calling
`new_tracer_particles (n)`.
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

/**
## Particle Advection  
 
The used advection scheme is the RK-3 variant of Sanderse and Veldman
(2019) for $\alpha \neq \frac{2}{3} \pm \mathtt{tol}$:

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

In case $\alpha \approx \frac{2}{3}$, we skip the last stage and
switch to a second-order-accurate scheme.
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
/**
Velocities may be interpolated and stored in `p.u` or `p.u2`.
 */
void interpolate_vel (Particles p, int swtch) {
  if (swtch == 1) {
    foreach_particle_in(p) {
      if (locate (x, y, z).level < 0) {//freeze
	fprintf (stderr, "#pid: %d could not locate p.tag: %ld.\n",
		 pid(), p().tag); //This should not happen
	foreach_dimension()
	  p().u.x = 0;
      } 
      foreach_dimension()
	p().u.x =interpolate_linear (locate (x, y, z),
				     (struct _interpolate)
				     {u.x, x, y, z});
    }
  } else {
    foreach_particle_in(p) {
      if (locate (x, y, z).level < 0) {//freeze
	fprintf (stderr, "#pid: %d could not locate p.tag: %ld.\n",
		 pid(), p().tag); //Tis should not happen
	foreach_dimension()
	  p().u2.x = 0;
      } 
      foreach_dimension()
	p().u2.x =interpolate_linear (locate (x, y, z),
				     (struct _interpolate)
				     {u.x, x, y, z});
    }
  }
}

void RK_step1 (Particles p, double dtf[2]) {
  interpolate_vel(p, 1);
  foreach_particle_in(p) { 
    foreach_dimension() {
      p().locn.x = p().x; //store the locations at t_n 
      p().x += dtf[0]*p().u.x;
    }
  }
}

void RK_step2 (Particles p, double dtf[2]) {
  interpolate_vel(p, 2);
  double a1 = -1, a2 = 2, h = dtf[1] + dtf[0];
  if (dtf[1] != dtf[0]) {
    double c = dtf[0]/h;
    if (fabs (c - 2./3.) > tol) 
      a2 = (c - 1.)/(c*(3.*c - 2.));
    else //Raltson's 2nd order method
      a2 = 1./(2*c);
    a1 = 1 - a2;
  }
  foreach_particle_in(p) {
    foreach_dimension() 
      p().x = p().locn.x + h*(a1*p().u.x + a2*p().u2.x);
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
				 (struct _interpolate){u.x, x, y, z}));
    }
  }
}

event advance_particles (i++, last) {
  foreach_P_in_list (tracer_particles) {
    particle_boundary (P);
    if (i%2 == 0) {
      if (i > 0) {
	RK_step3 (P, dtf);
	particle_boundary (P);
      }
      dtf[0] = dt;
      RK_step1 (P, dtf);
    } else { 
      dtf[1] = dt; 
      RK_step2 (P, dtf);
    }
  }
}

event free_tracers (t = end) {
  free (tracer_particles);
  tracer_particles = NULL;
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
    p().tag = j + offset;
}

/**
## Trace particle (`tp`) initialization
 */

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

Particles init_tp_cells(void) {
  int np = 0;
  foreach()
    np++;
  Particles p = new_tracer_particles (np);
  place_in_cells(p);
  particle_boundary (p); //?
  set_particle_attributes (p);
  return p;
}

Particles init_tp_square (struct Init_P inp) {
  if (!inp.n)
    inp.n = 10;
  if (!inp.l)
    inp.l = L0;
  Particles p;
  if (pid() == 0) {
    p = new_tracer_particles (sq(inp.n));
    place_in_square (p, inp);
  } else { 
    p = new_tracer_particles (0);
  }
  particle_boundary (p);
  set_particle_attributes (p);
  return p;
}

Particles init_tp_circle (struct Init_P inp) {
  Particles p;
   if (!inp.n)
     inp.n = 100;
  if (!inp.l)
    inp.l = L0/2.;
  if (pid() == 0) {
    p = new_tracer_particles (inp.n);
    place_in_circle (p, inp);
  } else { 
    p = new_tracer_particles (0);
  }
  particle_boundary (p);
  set_particle_attributes (p);
  return p;
}
