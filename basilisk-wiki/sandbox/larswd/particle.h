/**
# A Particle class

This headerfile defines some useful functions and types for generic
particles in terms of initialization and cleanup as well as iterators. 
It then includes the relevant sub-library based on solver type.

A particle is defined by a 3D position and possible additional members
that may be `#define`d via the hook (mass, size etc.).
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
`particle`s. The code maintains a list of all lists (`pl`) together
with a related `terminate_int`-terminated array (`pn`) that stores the
numbers of particles in each list and the (dynammically) allocated
space. Furthermore, a list can be referenced by an integer number ,
`Particles` (mind the case), indexing the how-manieth entry in `pn`
and `pl` a list of particles is.
 */

typedef particle * particles; //omit double pointer types and mistakes
typedef int Particles;        //An index reference
particles * pl = NULL;        
long unsigned int * pn = NULL, * pna =  NULL; // Particles and space
long unsigned int terminate_int = ULONG_MAX;

/**
`n` new particles can be declared using the `new_particles()`
function. All member-values are set to zero.
*/

Particles new_particles (long unsigned int n) {
  int j = 1;
  if (pn != NULL) 
    while (pn[j++] != terminate_int);
  pn = realloc (pn, (j + 1)*sizeof(long unsigned int));
  pna = realloc (pna, (j + 1)*sizeof(long unsigned int));
  pl = realloc (pl, j*sizeof(particle));
  pn[j] = terminate_int;
  pn[--j] = n;
  pna[j] = n;
  particle * pt = calloc (n, sizeof(particle));
  pl[j] = pt;
  return j;
}
/**
A list of particles can be added to a Particles list
*/

Particles * add_Particles_to_list (Particles p, Particles * plist) {
  int l = 0, t = 0;
  if (plist != NULL) {
    while (pn[l] != terminate_int) {
      if (l == plist[t]) 
	t++;
      l++;
    }
  }
  plist = realloc (plist, (t + 1)*sizeof(Particles));
  plist[t] = p;
  return plist;
}


/**
   A particle can be added to a list ...
*/
void change_plist_size (Particles Plist, int n) {
  long unsigned int n_new = n + pn[Plist];
  if (n_new > pna[Plist] || 2.5*(n_new + 1) < pna[Plist]) {
    pna[Plist] = 2*n_new;
    pl[Plist] = realloc (pl[Plist] , pna[Plist]*sizeof(particle));
  }
  pn[Plist] += n;
}

int add_particle (particle p, Particles Plist) {
#if _MPI
  if (!(locate (p.x, p.y, p.z).level > 0)) //Only place particle where it can be found
    return 0;
#endif
  
  change_plist_size (Plist, 1);
  pl[Plist][pn[Plist] - 1] = p; //append data
  return 1;
}
/**
... or removed (based on index or condition)
*/
void remove_particles_index (Particles Plist, int n, long unsigned int ind[n]) {
  int m = 0, j = 0;
  while (j < pn[Plist] - n) {
    while (m < n ? j + m == ind[m] : 0)
      m++;
    while (m < n ? j < pn[Plist] - n && j + m != ind[m] : j < pn[Plist] - n) {
      pl[Plist][j] = pl[Plist][j + m];
      j++;
    }
  }
  change_plist_size(Plist, -n);
}
/**
the "undo-everything" function can be used to prevent the
all-important memory leaks.
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
  free (pna);
  pna = NULL;
}

event free_particles (t = end)
  free_p();

/**
## Iterators

Two particle iterators are `@def`ined. Within these
`foreach_particle..` iterators, the coordinates `x, y, z` become
available and the `particle` data itself (`p()`).

* The `foreach_particle()` iterator loops over every particles.  

* The `foreach_particle_in(Particles P)` iterator loops over every
 particle in a single list referenced by `P`

Furthermore, the `foreach_P_in(Particles * Pl)` loops over lists of
`Particle`s. Inside the iterator, `Particles P` becomes available. 

For example, if you want to set the x-coordinate of all particles in
some list (say `tracer_particles`) to 1, you could do:

~~~literatec
...
Particles * tracer_particles;
...
foreach_P_in(tracer_particles) {
  foreach_particle_in(P) {
    p().x = 1.
  }
}
...
~~~

If all particles happen to be tracer_particles, this behaves identical to

~~~literatec
...
foreach_particle() {
  p().x = 1.
}
...
~~~

The Implementation makes use of `qcc`s excellent `foreach...`
functions. Meaning that these iterators can also do [simple
reductions](/Basilisk%20C#parallel-programming).
*/
#define p() pl[_l_particle][_j_particle] 

@def PARTICLE_VARIABLES
  double x = pl[_l_particle][_j_particle].x; NOT_UNUSED(x);
  double y = pl[_l_particle][_j_particle].y; NOT_UNUSED(y);
  double z = pl[_l_particle][_j_particle].z; NOT_UNUSED(z);
@

macro foreach_particle() {
  {
    int _l_particle = 0;
    while (pn[_l_particle] != terminate_int) {
      for (int _j_particle = 0; _j_particle < pn[_l_particle]; _j_particle++) {
	PARTICLE_VARIABLES;
	{...}
      }
      _l_particle++;
    }
  }
}

macro foreach_particle_in (int PARTICLES, Reduce reductions = None) {
  {
    int _l_particle = PARTICLES;
    for (int _j_particle = 0; _j_particle < pn[_l_particle]; _j_particle++) {
      PARTICLE_VARIABLES;
      {...}
    }
  }
}


macro foreach_P_in_list (int * PARTICLES_LIST) {
  {
    int _lll_particle = 0, _ll_particle = 0;
    while (pn[_lll_particle] != terminate_int) {
      if (_lll_particle == PARTICLES_LIST[_ll_particle]) {
	Particles P = _lll_particle; NOT_UNUSED(P);
	{...}
	_ll_particle++;
      }
      _lll_particle++;
    }
  }
}

#define remove_particles(p,CONDITION) do {		\
  int nr = 0;						\
  foreach_particle_in ((p)) {				\
    if ((CONDITION))					\
      nr++;						\
  }							\
  long unsigned int ind[nr];				\
  int i = 0;						\
  foreach_particle_in ((p)) {				\
    if ((CONDITION))					\
      ind[i++] = _j_particle;					\
  }							\
  remove_particles_index ((p), nr, ind);		\
  } while (0)


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
      int * ind = malloc (sizeof(int)*out);
      particle * senddata = malloc (sizeof(particle)*out);
      foreach_particle_in(p) { 
	if (locate (p().x, p().y, p().z).level < 0) {
	  ind[m] = _j_particle;
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
      particle * recdata = malloc (sizeof(particle)*outt);
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
    long unsigned int po = pn[l] - out;
    //Manage the memory if required...
    change_plist_size (l, in - out);
    //Collect new particles from `recdata`
    if (in > 0) {
      int indi[in];
      m = 0;
      for (int j = 0; j < outt; j++) 
	if (locate (recdata[j].x, recdata[j].y, recdata[j].z).level > 0)
	  indi[m++] = j;
      m = 0;
      for (j = po; j < pn[l]; j++) {
	pl[l][j] = recdata[indi[m]];
	m++;
      }
    }
    // clean the mess
    free (ind); free (senddata); free (recdata);
    }
    l++;
  }
  
}
#endif




// Some shared placement functions

void place_in_square (Particles p, int n = 10,
		      double xm = 0, double ym = 0,
		      double l = L0, double zp = 0)
{
  long unsigned int j = 0;
  particles pp = pl[p];
  double dx = l/(double)n;
  double x0 = xm - l/2. + dx/2;
  double y0 = ym - l/2. + dx/2;
  for (double xp = x0; xp < xm + l/2.; xp += dx) {
    for (double yp = y0; yp < ym + l/2.; yp += dx) {
      pp[j].x = xp;
      pp[j].y = yp;
      pp[j++].z = zp;
    }
  }
}

void place_in_circle (Particles p, int n = 10,
		      double xm = 0, double ym = 0,
		      double l = L0, double zp = 0) {
  particles pp = pl[p];
  for (int j = 0; j < n; j++) {
    double angle = noise()*pi;
    double R = sqrt(fabs(noise()))*l;
    pp[j].x = R*sin(angle) + xm;
    pp[j].y = R*cos(angle) + ym;
    pp[j].z = zp;
  }
}

// A shared particle statistics struct
typedef struct {
  coord avg;
  coord min;
  coord max;
  coord stddev;
} pstats;

/**
### Dump and restore particles
The `pdump()` and `prestore()` functions are implemented below.
 */

bool pdump (char * fname = "pdump",  //File name
	    Particles * td = NULL,   //Particle lists  Default: all
	    FILE * fp = NULL,        //File pointer    Default: Not used
	    bool Dmem = false,       //Member data     Default: false, only x,y,z
	    bool print = false) {    //Show dump data  Default: false
  // The file:
  char * file = fp ? NULL : fname;
  if (file) {
    if ((fp = fopen (file, "wb")) == NULL) {
      perror (file);
      exit (1);
    }
  }
  assert (fp);
  
  // Get particle lists data
  int j = 0, n = 0;        
  if (td) {
    foreach_P_in_list (td) 
      j++;
  } else  {// all
    while (pn[j++] != terminate_int);
    td = malloc (sizeof(int)*j);
    j = 0;
    while (pn[j] != terminate_int) { 
      td[j] = j;
      j++;
    }
  }

  // Nr of particles
  long unsigned int np[j];
  foreach_P_in_list (td) {
    np[n] = pn[P];
    n++;
  }
#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, np, j, MPI_UNSIGNED_LONG, MPI_SUM,
		 MPI_COMM_WORLD);
#endif

  // Print header
  if (pid() == 0) {
    fwrite (&j, sizeof(int), 1, fp);
    fwrite (np, sizeof(long unsigned int), j, fp);
    fwrite (&Dmem, sizeof(bool), 1, fp);
  }
  int Headersize = sizeof(int) + j*sizeof (long unsigned int) + sizeof(bool); 
  fseek (fp, Headersize, SEEK_SET); //offset

  if (print && pid() == 0) 
    fputs ("# P\tamount\n", stderr);

  // Print particle data
  for (int m = 0; m < j; m++) { //Nesting within foreach_P... did not work
    Particles P = td[m];
#if _MPI
    int size = Dmem ? sizeof(particle) : 3*sizeof(double);
    long unsigned int outa[npe()], outat[npe()  + 1];
    outat[0] = 0;
    MPI_Allgather (&pn[td[m]], 1, MPI_UNSIGNED_LONG,
		   &outa[0], 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    //Compute and set displacements
    for (int j = 1; j <= npe(); j++) 
      outat[j] = outa[j - 1] + outat[j - 1];
    fseek (fp, outat[pid()]*size, SEEK_CUR);
#endif
    foreach_particle_in (P) {
      if (Dmem) {
	fwrite (&p(), sizeof(particle), 1, fp);
      } else {
	double c[3] = {x, y, z};
	fwrite (c, sizeof(double), 3, fp);
      }
    }
    if (print && pid() == 0) 
      fprintf (stderr, "# %d\t%ld\n", m, np[m]);
#if _MPI
      fseek (fp, (np[m] - outat[pid() + 1])*size, SEEK_CUR);
#endif
  }
  fclose (fp);
  if (!td)
    free (td);
  return true; 
}

/**
These include statements separate out all functions which behave differently in 
the multilayer solver compared to the other solvers in basilisk C. The multilayer patch was developed by [larswd](basilisk.fr/sandbox/larswd/README). The [multilayer code](particle_multilayer.h) and the [code for the other solvers](particle_classic.h) are separated for ease of comparison.
*/

#if LAYERS
/**
Particles can be removed based on their position or another condition
 */
int remove_particle (particle p, Particles Plist) {
  long unsigned int * ind = malloc (sizeof(long int));
  int i = 0, a = 1; 
  foreach_particle_in(Plist) {
    foreach_dimension()
      if (p().x != p.x)
	      continue;
    if (p().z != p.z){
      continue;
    }
    if (i + 1 >= a) {
      ind = realloc (ind, a*2*sizeof(long int));
      a *= 2;
    }
    ind[i++] = _j_particle;
  }
  remove_particles_index (Plist, i, ind);
  free (ind); ind = NULL;
  return i;
}

/**
## Utility functions

### boundary conditions

peridoc box and `_MPI` boundaries.

Periodicity in vertical does not make sense
in the multilayer solver, and hence this function remains currently unchanged. 
 */
void particle_boundary (Particles P) {
  coord mind = {X0, Y0, 0}; 
  foreach_particle_in(P) {
    foreach_dimension() {
      if (p().x < mind.x) 
	      p().x += L0;
      else if (p().x > (mind.x + L0))
	      p().x -= L0;
    }
    Point point = locate(p().x,p().y,p().z);
    if (p().z > eta[]){
      p().z = eta[] - (p().z -eta[]);
    } else if (p().z < zb[]){
      p().z = zb[] - (p().z-zb[]);
    }
  }
#if _MPI
  update_mpi (P);
#endif
}


/**
### Place `particle`s

For already allocated `particle`s referenced by `Particles P`:

The foreach() and foreach_dimension()-iterators do not see z-axis in multilayer,
and hence the support for this axis must be added manually where relevant. 
*/
void place_in_cells (Particles P) {
  particles pp = pl[P];                 
  long unsigned int np = 0;
  foreach(serial) {
    double z = zb[];
    foreach_layer(){
      z = z + h[]/2;
      coord loc = {x, y,z};
      foreach_dimension()
        pp[np].x = loc.x;
    
      pp[np].z = loc.z;
      np++;
      z = z + h[]/2;
    }
  }
}

/**
### Simple particles statistics

The function computes, Average location, min, max and standard
deviation vectors. The correct statistics are only available for
`pid() == 0`.
 */

pstats statsp (Particles P) {
  coord avg = {0 ,0, 0}, min = {0,0,0},
    max = {0, 0, 0}, stddev = {0, 0, 0};
  foreach_dimension(3) {
    avg.x = stddev.x = 0;
    min.x = HUGE;
    max.x = -HUGE;
  }
  avg.z = stddev.z = 0;
  min.z = HUGE;
  max.z = -HUGE;
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
      foreach_particle_in(P, reduction(+:avgx) reduction(max:maxx) reduction (min:minx)) {
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
    double avgz =0; double minz = HUGE; double maxz = -HUGE;
    foreach_particle_in(P, reduction(+:avgz) reduction(max:maxz) reduction(min:minz)){
      avgz += p().z;
      if (p().z < minz){
        minz = p().z;
      }
      if (p().z > maxz){
        maxz = p().z; 
      }
    }
    avg.z = avgz;
    min.z = minz;
    max.z = maxz;
    
    foreach_dimension() {
      double stddevx = 0;
      foreach_particle_in(P, reduction(+:stddevx)) {
	      stddevx += sq(avg.x - p().x);
      }
      stddev.x = sqrt(stddevx/(double)np);
    }
    double stddevz = 0;
    foreach_particle_in(P, reduction(+:stddevz)){
      stddevz += sq(avg.z - p().z); 
    }
    stddev.z = sqrt(stddevz/(double)(np));
  }
  pstats s;
  s.max = max, s.min = min, s.avg = avg, s.stddev = stddev;
  return s;
}
/**
### Propability density function

Obtain a scalar PDF from particle locations 
 */
void particle_pdf (Particles P, scalar s) {
  
  foreach_layer(){
  foreach(){
    s[] = 0;
  }
  }
  particle_boundary (P);
  
  foreach_particle_in(P) {
    Point point = locate (x,y,z);
    double zc = zb[]; int k=0; int found = 0;
    while (!found){
      if (p().z <= zc + h[0,0,k]){
        found = 1;
      } else if (k == nl-1){
        found = 1;
      } else {
        k++;
        zc += h[0,0,k];
      }
    }
    if (point.level > 0)
      s[0,0,k]++;
  }

  foreach_layer(){
  foreach()
    s[] /= (pn[P]*dv());
  }
  boundary ({s});
}
/**
### Random step

Particles displace a certain distance (`step`) in a random direction
*/
void random_step (Particles P, double step) {
  foreach_particle_in(P) {
    double theta = noise()*pi;
    #if (dimension == 1)
      coord d = {sign(theta), 0, cos(theta)};
    #else
      double phi = acos(noise());
      coord d = {sin(phi)*cos(theta), sin(phi)*sin(theta), cos(phi)};
    #endif
    foreach_dimension()
      p().x += d.x*step;
    p().z += d.z*step;
  }
}
/**
### Length of particle loop

A function that computes the line length of particles places in a loop
(mind the ordering). Special care is taken to obtain 4th order
accuracy.
 */
double plength (Particles P) {
  double lr = 0;
  int np = pn[P];
  particles pli = pl[P];
  foreach_particle_in(P) {
    int il = _j_particle - 1 < 0 ? np - 1: _j_particle - 1;
    int ir = _j_particle + 1 >= np ? 0 : _j_particle + 1;
    coord pl = (coord){pli[il].x, pli[il].y, pli[il].z};
    coord pm = (coord){pli[_j_particle].x, pli[_j_particle].y, pli[_j_particle].z};
    coord pr = (coord){pli[ir].x, pli[ir].y, pli[ir].z};
    coord bv = {0, 0, 0}, a1 = {0,0,0};
    coord a2 = {0,0,0};
    double bl = 0, ka = 0;
    foreach_dimension() {
      a1.x = pm.x - pl.x;
      bv.x = pr.x - pl.x;
      bl += sq(bv.x);
    }
    a1.z = pm.z - pl.z;
    bv.z = pr.z - pl.z;
    bl  += sq(bv.z);
    bl = sqrt(bl);
    foreach_dimension()
      ka += a1.x*bv.x/bl;
    ka += a1.z*bv.z/bl;
    foreach_dimension()
      a2.x = a1.x - bv.x*ka/bl;
    a2.z = a1.z - bv.z*ka/bl;
    normalize (&a2);
    double al = 0, am = 0, ar = 0;
    foreach_dimension() {
      al -= a2.x*pl.x;
      am -= a2.x*pm.x;
      ar -= a2.x*pr.x;
    }
    al -= a2.z*pl.z;
    am -= a2.z*pm.z;
    ar -= a2.z*pr.z;

    double C = fabs(am - al);
    double xt = 0;
    double b = 0;
    foreach_dimension() {
      xt += sq(pl.x - pm.x);
      b += sq(pl.x - pr.x);
    }
    xt += sq(pl.z - pm.z);
    b  += sq(pl.z - pr.z);
    xt = sqrt (xt - sq(C));
    b = sqrt(b);
    double A =  C/((xt*(b - xt)));
    double xp1 = 0.5*b*(1 - 1/sqrt(3));
    double xp2 = 0.5*b*(1 + 1/sqrt(3));
    double l1 = b*sqrt(1. + sq(-2*A*xp1 + A*b));
    double l2 = b*sqrt(1. + sq(-2*A*xp2 + A*b));
    lr += (l1 + l2)/4;
  }
  return lr;
}
#else
/**
Particles can be removed based on their position or another condition
 */

 int remove_particle (particle p, Particles Plist) {
  long unsigned int * ind = malloc (sizeof(long int));
  int i = 0, a = 1; 
  foreach_particle_in(Plist) {
    foreach_dimension()
      if (p().x != p.x)
	      continue;
    if (i + 1 >= a) {
      ind = realloc (ind, a*2*sizeof(long int));
      a *= 2;
    }
    ind[i++] = _j_particle;
  }
  remove_particles_index (Plist, i, ind);
  free (ind); ind = NULL;
  return i;
}

/**
## Utility functions

### boundary conditions

peridoc box and `_MPI` boundaries.

Periodicity in vertical does not make sense
in the multilayer solver, and hence this function remains currently unchanged. 
 */

void particle_boundary (Particles P) {
  coord mind = {X0, Y0, Z0};
  foreach_particle_in(P) {
    foreach_dimension() {
      if (p().x < mind.x) 
	      p().x += L0;
      else if (p().x > (mind.x + L0))
	      p().x -= L0;
    }
  }
#if _MPI
  update_mpi (P);
#endif
}


/**
### Place `particle`s

For already allocated `particle`s referenced by `Particles P`:

The foreach() and foreach_dimension()-iterators do not see z-axis in multilayer,
and hence the support for this axis must be added manually where relevant. 
*/


void place_in_cells (Particles P) {
  particles pp = pl[P];                 
  int np = 0;
  foreach(noauto) {
    coord loc = {x, y,z};
    foreach_dimension()
      pp[np].x = loc.x;
    np++;
  }
}

/**
### Simple particles statistics

The function computes, Average location, min, max and standard
deviation vectors. The correct statistics are only available for
`pid() == 0`.
 */


pstats statsp (Particles P) {
  coord avg = {0 ,0, 0}, min = {0,0,0},
    max = {0, 0, 0}, stddev = {0, 0, 0};
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
      foreach_particle_in(P, reduction(+:avgx) reduction(max:maxx)
			  reduction (min:minx), serial) {
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
      foreach_particle_in(P, reduction(+:stddevx), serial) {
	stddevx += sq(avg.x - p().x);
      }
      stddev.x = sqrt(stddevx/(double)np);
    }
  }
  
  pstats s;
  s.max = max, s.min = min, s.avg = avg, s.stddev = stddev;
  return s;
}


/**
### Propability density function

Obtain a scalar PDF from particle locations 
 */

static inline void count_particle (double x, double y,
				   double z, scalar s) {
  foreach_point (x, y, z, serial)
    if (point.level > 0)
       s[]++;
}

void particle_pdf (Particles P, scalar s) {
  foreach()
    s[] = 0;
  particle_boundary (P);
  foreach_particle_in(P) 
    count_particle (x, y, z, s);
  foreach()
    s[] /= (pn[P]*dv());
  boundary ({s});
}


/**
### Random step

Particles displace a certain distance (`step`) in a random direction
*/

void random_step (Particles P, double step) {
  foreach_particle_in(P) {
    double theta = noise()*pi;
    #if (dimension == 1)
      coord d = {sign(theta)};
    #elif (dimension < 3)
      coord d = {sin(theta), cos(theta)};
    #else
      double phi = acos(noise());
      coord d = {sin(phi)*cos(theta), sin(phi)*sin(theta), cos(phi)};
    #endif
    foreach_dimension()
      p().x += d.x*step;
  }
}


/**
### Length of particle loop

A function that computes the line length of particles places in a loop
(mind the ordering). Special care is taken to obtain 4th order
accuracy.
 */
double plength (Particles P) {
  double lr = 0;
  int np = pn[P];
  particles pli = pl[P];
  foreach_particle_in(P) {
    int il = _j_particle - 1 < 0 ? np - 1: _j_particle - 1;
    int ir = _j_particle + 1 >= np ? 0 : _j_particle + 1;
    coord pl = (coord){pli[il].x, pli[il].y, pli[il].z};
    coord pm = (coord){pli[_j_particle].x, pli[_j_particle].y, pli[_j_particle].z};
    coord pr = (coord){pli[ir].x, pli[ir].y, pli[ir].z};
    coord bv = {0, 0, 0}, a1 = {0,0,0};
    coord a2 = {0,0,0};
    double bl = 0, ka = 0;
    foreach_dimension() {
      a1.x = pm.x - pl.x;
      bv.x = pr.x - pl.x;
      bl += sq(bv.x);
    }

    bl = sqrt(bl);
    foreach_dimension()
      ka += a1.x*bv.x/bl;

    foreach_dimension()
      a2.x = a1.x - bv.x*ka/bl;

    normalize (&a2);
    double al = 0, am = 0, ar = 0;
    foreach_dimension() {
      al -= a2.x*pl.x;
      am -= a2.x*pm.x;
      ar -= a2.x*pr.x;
    }
    double C = fabs(am - al);
    double xt = 0;
    double b = 0;
    foreach_dimension() {
      xt += sq(pl.x - pm.x);
      b += sq(pl.x - pr.x);
    }
    xt = sqrt (xt - sq(C));
    b = sqrt(b);
    double A =  C/((xt*(b - xt)));
    double xp1 = 0.5*b*(1 - 1/sqrt(3));
    double xp2 = 0.5*b*(1 + 1/sqrt(3));
    double l1 = b*sqrt(1. + sq(-2*A*xp1 + A*b));
    double l2 = b*sqrt(1. + sq(-2*A*xp2 + A*b));
    lr += (l1 + l2)/4;
  }
  return lr;
}

#endif
/**
#### Particle Restoration 

The restore function forgets all existing particles and does not
assign particles to any list, nor does it know about the names of the
particle lists references. See also [this test](test_prestore.c). */
int prestore (char * fname = "pdump",  //File name
	      FILE * fp = NULL,        //File pointer    Default: Not used
	      bool print = false)      //Show dump data  Default: false
{
  // Open file
  if (fp) fname = NULL;
  else if (fname && (fp = fopen (fname, "r")) == NULL)
    return 0;
  assert (fp);

  //read nr. of lists and nr. of particles
  int j;
  bool Dmem;
  fread (&j, sizeof(int), 1, fp);
  long unsigned int np[j];
  fread (np, sizeof(long unsigned int), j, fp);
  fread (&Dmem, sizeof(bool), 1, fp);
  
  // Print some reference data
  if (print) {
    if (pid() == 0) {
      if (Dmem)
	fputs ("# Restoring members...\n", stderr);
      fputs ("# P\tamount\n", stderr);
      for (int m = 0; m < j; m++)
	fprintf (stderr, "# %d\t%ld\n", m, np[m]);
    }
  }
  
  // Remove existing particles
  if (pl != NULL)
    free_p();

  // Allocate and read
  Particles P;
  for (int m = 0; m < j; m++) {
    if (pid() == 0) { // This could be more balanced
      P = new_particles (np[m]);
      if (Dmem) 
	fread (pl[m], sizeof(particle), np[m], fp);
      else {
	foreach_particle_in (P) {
	  double c[3];
	  fread (&c, sizeof(double), 3, fp);
	  p().x = c[0]; p().y = c[1]; p().z = c[2];
	}
      }
    } else // slaves do not read data
      P = new_particles (0);
     // Apply boundary conditions.
    particle_boundary (P);
  }
  if (fname)
    fclose (fp);
  return j;
}


/**
## todo

* ~~~Dynamically add and delete single `particle`s and `Particles` lists~~~
* Sort particles along grid-iterator curve ([radix.h]())
* More flexible attributes/members for `particle`s. 
* ~~~[Tie particles to Basilisk fields](tie_coord_to_cell.c)~~~
* ~~~dump and restore particles~~~


## Tests

* [Dump and restore particles with MPI](test_prestore.c)
* [Length of a loop](tlengt.c)
* [Adding and removing particles](inject_particles.c)

## Usage

* [flow tracer particles](tracer-particles.h)
* [Inertial particles](inertial-particles.h)
* [Brownian motion](brownian.c)

## See also

* [Bview draw function for particles](scatter2.h)
* [Reference particles using scalar data](particle_reference.h)
 */

