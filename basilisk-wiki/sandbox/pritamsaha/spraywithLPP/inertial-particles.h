/**
Taken from Antoon's sandbox
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
// Time integration
// A version of the second-order accurate Verlet method is used. (a.k.a. leaffrogging). Since p_acc may depend on p.u, we advance the particles in two solver timesteps.

// The following events take care of the time integration.

double dtprev;

event defaults (t = 0);

event init (t = 0);

event set_dtmax (i++);

event velocity (i++);

void ip_step1 (Particles P, int i) {
  double Dt = i > 0 ? 2*dt : dt;
  foreach_particle_in(P) {
    coord a = p_acc (INP);
    foreach_dimension()
      p().u.x += a.x*Dt;
  }
}

void ip_step2 (Particles P) {
  foreach_particle() {
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
// Utilities
// From tracer-particles.h, a bunch of functions are copied to provide some user interface for particle initializaiton.

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
/*
void extra_inertial_particles (inertial_particles IP, long unsigned int nextra) {
  //Particles p = new_particles (n);
  Particles P = inertial_particles[IP];
  change_plist_size(P, nextra); 
  int n0 = pn
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
*/

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
    p().tag = j + offset;
}

void set_ip_particle_attributes (Particles p) {
  tag_ip_particles (p);
}
// Inertial-particle (ip) initialization
// Inertial-particles may be initialized by calling the following functions.

Particles init_ip_cells(void) {
  int np = 0;
  foreach()
    np++;
  Particles p = new_inertial_particles (np);
  place_in_cells(p);
  particle_boundary (p); //?
  set_ip_particle_attributes (p);
  return p;
}

Particles init_ip_square (struct Init_P inp) {
  if (!inp.n)
    inp.n = 10;
  if (!inp.l)
    inp.l = L0;
  Particles p;
  if (pid() == 0) {
    p = new_inertial_particles (sq(inp.n));
    place_in_square (p, inp);
  } else { 
    p = new_inertial_particles (0);
  }
  particle_boundary (p);
  set_ip_particle_attributes (p);
  return p;
}

Particles init_ip_circle (struct Init_P inp) {
  Particles p;
   if (!inp.n)
     inp.n = 100;
  if (!inp.l)
    inp.l = L0/2.;
  if (pid() == 0) {
    p = new_inertial_particles (inp.n);
    place_in_circle (p, inp);
  } else { 
    p = new_inertial_particles (0);
  }
  particle_boundary (p);
  set_ip_particle_attributes (p);
  return p;
}
