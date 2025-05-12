/**
# Reference particles with scalar-field data

It is possible to assign particles to cells, and find them.

The user interface is something like this:

~~~literatec
...
Particles parts = new_...
scalar p[];
assign_particles (parts, p);
foreach() {
 double ccx = x, ccy = y;
  foreach_particle_point(p) {
    printf ("%g %g %g %g", x, y, ccx, ccy);
  }
}
~~~

*Except*(!) qcc does not allow these nested loops and you must wrap it
in some function (see test). 
 */

attribute { 
  int plist; //offset by 1
}

#include "particle.h"

typedef union {
  int * p;
  double d;
} datic;

#define pointer_v(a) ((datic){.d = a}.p)
#define field_v(a) ((datic){.p = a}.d)

bool particle_in_cell (particle part, Point point) {
  coord cc = {x, y, z};
  foreach_dimension() {
    //  Particles can be *on* the rhs face 
    if (part.x >  cc.x + Delta/2. ||
	part.x <= cc.x - Delta/2.) 
      return false;
  }
  return true;
}

#if TREE
void p_refine (Point point, scalar s) {
  if (!s[]) { // No particles
    foreach_child()
      s[] = 0;
  } else { // distribute particles among children
    int * ind = pointer_v(s[]);
    particles parts = pl[s.plist - 1];
    foreach_child() {
      int n = 0, c = 0;
      int * indc = NULL;
      while (ind[n] >= 0) {
	if (particle_in_cell (parts[ind[n]], point)) {
	  indc = realloc (indc, (c + 2)*sizeof(int));
	  indc[c++] = ind[n];
	  indc[c] = -1;
	}
	n++;
      }
      if (indc == NULL)
	s[] = 0; 
      else
	s[] = field_v(indc);
    }
  }
}

void p_coarsen (Point point, scalar s) {
  s[] = 0;
  int nt = 0;
  foreach_child()
    if (s[]) {
      int * indc = pointer_v(s[]);
      int n = 0; while(indc[n++] >= 0);
      nt += n - 1;
    }
  if (nt) {
    int * indp = NULL;
    indp = malloc (sizeof(int)*nt + 1);
    nt = 0;
    foreach_child()
      if (s[]) {
	int * indc = pointer_v(s[]);
	int n = 0;
	while(indc[n] >= 0) 
	  indp[nt++] = indc[n++]; 
      }
    indp[nt] = -1;
    s[] = field_v(indp);
  }
}
#endif

void free_scalar_data (scalar s) {
  foreach_cell() {
    free(pointer_v(s[]));
  }
}

void assign_particles (Particles plist, scalar s) {
  if (s.plist > 0)
    free_scalar_data (s);
  s.plist = plist + 1;
  foreach()
    s[] = 0;
#if TREE
  s.refine = s.prolongation = p_refine;
  s.coarsen = s.restriction = p_coarsen;
#endif
  foreach_particle_in(plist) {
    Point point = locate (x, y, z);
    int * ind = NULL, n = 0;
    if (s[] == 0.) 
      n = 1;
    else {
      ind = pointer_v(s[]);
      while (ind[n++] >= 0);
    }
    ind = realloc (ind, (n + 1)*sizeof(int));
    ind[n - 1] = j;
    ind[n] = -1;
    s[] = field_v(ind);
  }
  boundary ({s}); // Find particles in ghosts and parent cells
}

// We need wrappers to help qcc
double value_p (scalar s, Point point) {
  return s[]; // Macro may not be expanded during `@def` preprocessing
}

int plist_s (scalar s) {
  return s.plist - 1; // Atribute may not exists at `@def` preprocessing
}

@def foreach_particle_point(s) {
  int l = plist_s(s);				
  if (value_p(s, point)) {					
    int * ind = pointer_v(value_p(s, point));			
    for (int n = 0; ind[n] >= 0; n++) {		
      int j = ind[n];
      PARTICLE_VARIABLES;			
      @						
	@def end_foreach_particle_point()
	}					
    }						
    }						
      @						
	
	/**
## Test

* [Find particles in a 3-point neighborhood](lp.c)

## Todo

* MPI
	 */
