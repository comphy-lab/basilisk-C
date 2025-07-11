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
//fix me: particles are in sibling ghost cells
void ref_prolongation (Point point, scalar s) {
  double c = s[];
  foreach_child() {
    if ((child.x + child.y) == 0)
      s[] = c;
    else
      s[] = 0;
  }
}

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
      if (c == 0)
	s[] = 0; 
      else
	s[] = field_v(indc);
    }
  }
}

void ref_restriction (Point point, scalar s) {
  int nt = 0;
  foreach_child()
    if (s[]) {
      int * indc = pointer_v(s[]);
      int n = 0;
      while(indc[n++] >= 0);
      nt += n - 1;
    }
  if (nt) {
    int * indp = NULL;
    if (s[]) 
      indp = pointer_v(s[]);
    indp = realloc (indp, sizeof(int)*(nt + 1));
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

void ref_coarsen (Point point, scalar s) {
  ref_restriction (point, s);
  foreach_child() {
    if (s[])
      free (pointer_v(s[]));
    s[] = 0;
    
  }
}

void free_scalar_data (scalar s) {
  foreach_cell_all() {
    fflush(stdout);
    if (!is_prolongation(cell) && s[] != 0) {
      fflush (stdout);
      free(pointer_v(s[]));
      pointer_v(s[]) = NULL;
    }
    s[] = 0;
    
  }
}

void assign_particles (Particles plist, scalar s) {
  if (s.plist > 0)
    free_scalar_data (s);
  s.plist = plist + 1;
  foreach()
    s[] = 0;
  // No box-boundary points
  foreach_dimension() {
    s[left] = 0;
    s[right] = 0;
  }
  s.prolongation = ref_prolongation;
  s.restriction = ref_restriction; 
#if TREE
  s.refine = p_refine;
  s.coarsen = ref_coarsen;
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
    ind[n - 1] = _j_particle;
    ind[n] = -1;
    s[] = field_v(ind);
  }
  multigrid_restriction ({s});
  boundary ({s}); // Find particles in ghosts and parent cells
}

// We need wrappers to help qcc
double value_p (scalar s, Point point) {
  return s[]; // Macro may not be expanded during `@def` preprocessing
}

int plist_s (scalar s) {
  return s.plist - 1; // Atribute may not exists at `@def` preprocessing
}

macro foreach_particle_point(scalar s, Point point) {
  int _l_particle = plist_s(s);
  NOT_UNUSED(_l_particle);
  if (value_p(s, point)) {
    int * ind = pointer_v(value_p(s, point));			
    for (int n = 0; ind[n] >= 0; n++) {		
      int _j_particle = ind[n];
      PARTICLE_VARIABLES;			
      {...}
    }					
  }
}
	
	/**
## Test

* [Find particles in a 3-point neighborhood](lp.c)

## Todo

* MPI
	 */
