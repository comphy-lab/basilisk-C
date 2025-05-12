/**
# A Particle class

This headerfile defines some useful definitions for generic particles
that may be tuned for a (single) particular purpose.

A particle is defined by a 3D position and possible
additional members that may be `#define`d via the hook (mass, size
etc.).
*/
typedef struct {
  double x;
  double y;
  double z;
#ifdef ADD_PART_MEM
  ADD_PART_MEM
#endif
} particle;
/**
In practice, more than one particle is stored in a list of
`particles`. It is defined as a list of such a type together with a
related `terminate_int`-terminated array that stores the numbers of
particles in the lists.
 */
typedef particle * particles; //omit double pointer types and mistakes
typedef int Particles;        //Maintaining an alias is hard
particles * pl = NULL;
long unsigned int * pn = NULL;
long unsigned int terminate_int = ULONG_MAX;
/**
`n` new particles can be declared using the `new_particles()`
function. All member-values are set to zero.
*/
#if _MPI //Trace allocated space
long unsigned int * pna;
void MPI_init_part (long unsigned int n) {
  int j = 1;
  if (pna != NULL) 
    while (pna[j++] != terminate_int);
  pna = realloc (pna, (j + 1)*sizeof(long unsigned int));
  pna[j - 1] = n + 1;
  pna[j] = terminate_int;
}
#endif

Particles new_particles (long unsigned int n) {
  assert (n >= 0 || n != terminate_int); //The cast already takes care
  int j = 1;
  if (pn != NULL) 
    while (pn[j++] != terminate_int);
  pn = realloc (pn, (j + 1)*sizeof(long unsigned int));
  pl = realloc (pl, j*sizeof(particles));
  pn[j - 1] = n;
  pn[j] = terminate_int;
  particle * pt = calloc (n + 1, sizeof(particle));
  pl[j - 1] = pt;
#if _MPI
  MPI_init_part (n); 
#endif
  return j - 1;
}
/**
the "undo-all" function can be used to prevent the all-important
memory leaks.
 */
void free_p (void) {
  int j = 0;
  while (pn[j++] != terminate_int) {
    free (pl[j - 1]);
    pl[j -1] = NULL;
  }
  free (pl);
  pl = NULL;
  free (pn);
  pn = NULL;
#if _MPI
  free (pna);
  pna = NULL;
#endif
}

event free_particles (t = end)
  free_p();
/**
## Iterators

A `foreach_particle()` iterator is defined that loops through `pn[i]`
particles in each `i`th entry of `pl`.
 */

@define p() pl[l][j] 

@def PARTICLE_VARIABLES
  double x = pl[l][j].x; NOT_UNUSED(x);
  double y = pl[l][j].y; NOT_UNUSED(y);
  double z = pl[l][j].z; NOT_UNUSED(z);
@

@def foreach_particle() {
  int l = 0;
  while (pn[l] != terminate_int) {
    for (int j = 0; j < pn[l]; j++) {
      PARTICLE_VARIABLES
@
@def end_foreach_particle()
    }
    l++;
  }
}
@
/**
Also, a `foreach_particle_in()` iterator is `@def`ined that loops over
the `particle`s in `PARTICLES`, which in turn is of the type
`Particles`...
 */

@def foreach_particle_in(PARTICLES) {
  int l = PARTICLES;
  for (int j = 0; j < pn[l]; j++) {
    PARTICLE_VARIABLES
@
@def end_foreach_particle_in()
  }
}
@
/**
Finally, an iterator exists that loops over a `Particles` in a
`Particles*`,  not the `partticle`s them selves.
 */
@def foreach_P_in_list(PARTICLES_LIST) {
  int lll = 0, ll = 0;
  while (pn[lll] != terminate_int) {
    if (lll == PARTICLES_LIST[ll]) {
	Particles P = lll;	
@
@def end_foreach_P_in_list()
        ll++;
    }
    lll++;
  }
}
@
/**
## MPI particle exchange 

With `_MPI`, particles may find themselves outside the realm of their
thread's domain. Here we implement a particle exchange function. If a
particle is not within any boundary, it is lost...
 */
#if _MPI
void update_mpi (Particles p) {
  int l = 0;
  while (pn[l] != terminate_int) {
    if (l == p) {
      int outt, in = 0, m = 0, out = 0;
      foreach_particle_in(p) 
	if (locate (p().x, p().y, p().z).level < 0) {
	  out++;
	}
      //get indices and outgoing data
      int ind[out];
      particle senddata[out];
      foreach_particle_in(p) { 
	if (locate (p().x, p().y, p().z).level < 0) {
	  ind[m] = j;
	  senddata[m++] = p();
	}
      }
      //remove the senddata from arrays (shrink)
      m = 0;
      int j = 0;
      while (j < pn[l] - out) {
	while (m < out ? j + m == ind[m] : 0)
	  m++;
	while (m < out ? j < pn[l] - out && j + m != ind[m] : j < pn[l] - out) {
	  pl[l][j]   = pl[l][j + m];
	  j++;
	}
      }
      // Gather lost particles among threads:
      // First, count all of them
      int outa[npe()], outat[npe()];
      outat[0] = 0;
      MPI_Allgather (&out, 1, MPI_INT, &outa[0], 1, MPI_INT, MPI_COMM_WORLD);
      //Compute displacements
      for (int j = 1; j < npe(); j++) 
	outat[j] = outa[j - 1] + outat[j - 1];
      outt = outat[npe() - 1] + outa[npe() - 1]; 
      // Allocate recieve buffer and gather
      particle recdata[outt];
    for (int j = 0; j < npe(); j++) {
      outat[j] *= sizeof(particle);
      outa[j]  *= sizeof(particle);
    }
    //send and recieve data
    MPI_Allgatherv (&senddata[0], outa[pid()], MPI_BYTE,
		    &recdata[0], outa, outat, MPI_BYTE,
		    MPI_COMM_WORLD); 
    //count new particles
    for (int j = 0; j < outt ; j++) 
      if (locate (recdata[j].x, recdata[j].y, recdata[j].z).level > 0)
	in++;
    int n_partn = pn[l] + in - out;
    //Manage the memory if required...
    if (n_partn > pna[l] || 2*(n_partn + 1) < pna[l]) {
      pna[l] = 2*(n_partn + 1);
      pl[l] = realloc (pl[l] , pna[l]*sizeof(particle));
    }
    //Collect new particles from `recdata`
    if (in > 0) {
      int indi[in];
      m = 0;
      for (int j = 0; j < outt; j++) 
	if (locate (recdata[j].x, recdata[j].y, recdata[j].z).level > 0)
	  indi[m++] = j;
      m = 0;
      for (j = pn[l] - out; j < n_partn; j++) {
	pl[l][j] = recdata[indi[m]];
	m++;
      }
    }
    //Update `c[l]p`
    pn[l] = n_partn;
    
    }
    l++;
  }
}
#endif

/**
## Utility functions

boundary conditions; peridoc box + MPI.
 */
void particle_boundary (Particles p) {
  coord mind = {X0, Y0, Z0}; 
  foreach_particle_in(p) { 
    foreach_dimension() {
      if (p().x < mind.x) 
	p().x += L0;
      else if (p().x > (mind.x + L0))
	p().x -= L0;
    }
  }
#if _MPI
  update_mpi (p);
#endif
}

/**
Place pre-allocated particles functions. 
 */
void place_in_cells (Particles p) {
  particles pp = pl[p];                 
  long unsigned int np = 0;
  foreach() {
    coord loc = {x, y, z};
    foreach_dimension()
      pp[np].x = loc.x;
    np++;
  }
}

struct Init_P {
  int n;
  double xm;
  double ym;
  double l;
  double zp;
};
				   
void place_in_square (Particles p, struct Init_P inp) {
  if (!inp.n)
    inp.n = 10;
  if (!inp.l)
    inp.l = L0;
  long unsigned int j = 0;
  particles pp = pl[p];
  double dx = inp.l/(double)inp.n;
  double x0 = inp.xm - inp.l/2. + dx/2;
  double y0 = inp.ym - inp.l/2. + dx/2;
  for (double xp = x0; xp < inp.xm + inp.l/2.; xp += dx) {
    for (double yp = y0; yp < inp.ym + inp.l/2.; yp += dx) {
      pp[j].x = xp;
      pp[j].y = yp;
      pp[j++].z = inp.zp;
    }
  }
}

void place_in_circle (Particles p, struct Init_P inp) {
  particles pp = pl[p];
  for (int j = 0; j < inp.n; j++) {
    double angle = noise()*pi;
    double R = sqrt(fabs(noise()))*inp.l;
    pp[j].x = R*sin(angle);
    pp[j].y = R*cos(angle);
    pp[j].z = inp.zp;
  }
}

/**
## Simple particles statistics

Average location, min, max and standard deviation. The correct
statistics are only available for `pid() == 0`).
 */

typedef struct {
  coord avg;
  coord min;
  coord max;
  coord stddev;
} pstats;

pstats statsp (Particles P) {
  coord avg, min, max, stddev;
  foreach_dimension(3) {
    avg.x = stddev.x = 0;
    min.x = HUGE;
    max.x = -HUGE;
  }
  long unsigned int np = pn[P];
#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, &np, 1, MPI_UNSIGNED_LONG,
		 MPI_SUM, MPI_COMM_WORLD);
#endif
  if (np) {
    foreach_dimension() { //reduction of coord members does not work
      double avgx = 0;
      double minx = HUGE;
      double maxx = -HUGE;
      foreach_particle_in(P reduction(+:avgx) reduction(max:maxx)
			  reduction (min:minx)) {
	avgx += p().x;
	if (p().x < minx)
	  minx = p().x;
	if (p().x > maxx)
	  maxx = p().x;
      }
      avg.x = avgx/(double)np;
      min.x = minx;
      max.x = maxx;
    }
    foreach_dimension() {
      double stddevx = 0;
      foreach_particle_in(P reduction(+:stddevx)) {
	stddevx += sq(avg.x - p().x);
      }
      stddev.x = sqrt(stddevx/(double)np);
    }
  }
  pstats s;
  s.max = max, s.min = min, s.avg = avg, s.stddev = stddev;
  return s;
}

  

