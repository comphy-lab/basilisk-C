/**
# Particles and interfaces

Tracer-particle paths may cross a tracer interface [(represented by
volume fractions)](/src/vof.h), leading to inconsistencies. A *hack*
to circumvent this complicated problem is to *modify* estimated
particle positions based on the position of the reconstructed
interface. Such "corrections" could be small, and may have a huge
inpact on the resulting particle paths.

Using the code in this header file, the user can choose to restrict
particles to a specific side of the interface, to stay on the
interface or to be unbounded (regular tracer particles). This is
specified by an argument passed to `new_vof_tracer_particles()`. For
ten particles use:

~~~literatec
...
Particles Pin  = new_vof_tracer_particles (10, 1);  // f[] == 1 phase
Particles Pout = new_vof_tracer_particles (10, 0);  // f[] == 0 phase
Particles Pon  = new_vof_tracer_particles (10, 3);  // On the interface
Particles Pun  = new_vof_tracer_particles (10, -1); // Unbounded
...
~~~ 

These particles should then be positioned accordingly.

## Implementation

The `tracer-particles` code is extended to also include an array
(`Pphase`) that marks the phase of a particle.

The `reposition()` function takes the vof field and a list of
particles (referenced by `Particles P`) as the input. It returns the
number of particles that it repositioned. Note that only the particles
inside partial cells are repositioned.
*/

extern scalar f; // the vof field
//#define u uf     // Also used for VOF advection
#include "tracer-particles.h"
//#undef u
#include "fractions.h"

int * Pphase = NULL;

Particles new_vof_tracer_particles (long unsigned int n, int phase) {
  Particles P = new_tracer_particles (n);
  Pphase = realloc (Pphase, sizeof(int) * (P + 1));
  Pphase[P] = phase;
  return P;
}

event free_vof_tracers (t = end) {
  free (Pphase);
  Pphase = NULL;
}

int reposition (scalar f, Particles P) {
  int rep = 0;
  if  (Pphase[P] >= 0) { //Do something
    double mir = Pphase[P] < 2 ? 2. : 1.; //Mirror *or* on interface
    foreach_particle_in(P) {
      Point point = locate (x, y, z); 
      coord cc = {x,y,z}; //Cell centre
      if (f[] > 1e-6 && f[] < 1 - 1e-6) {
	coord n = interface_normal (point, f); //pointing outwards.
	double ff = Pphase[P] == 1 ? f[] : 1 - f[];
	if (Pphase[P] != 1) {
	  foreach_dimension()
	    n.x *= -1;
	}
	normalize (&n);
	double alpha = plane_alpha (ff, n);
	double ALP = 0;
	foreach_dimension()
	  ALP += n.x*(p().x - cc.x);
	ALP /= Delta;
	if (ALP > alpha || mir != 2) { //Reposition particle
	  foreach_dimension()
	    p().x -= mir*(ALP - alpha)*n.x*Delta; 
	  rep++;
	}
      }
    }
  }
  return rep;
}

event defaults (i = 0) 
  P_RK2 = true;        //Match the VOF advection order of accuracy

event tracer_particles_step1 (i += 2, last) {//After RK step 3 (i.e 2)
  foreach_P_in_list (tracer_particles) {
    
    particle_boundary (P);
    reposition (f, P);
  }
}
/**
## Todo
* ~~~Allow unbounded particles along side bounded counterparts.~~~
* ~~~Bind particles to the interface itself~~~

## Test

* [Test the `reposition()` function](trepos.c)
* [Reversed advection test cased](reversedp.c)

## Usage 

* [Droplet inpact on a pool](splash.c)
 */
