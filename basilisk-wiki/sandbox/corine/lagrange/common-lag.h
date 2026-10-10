/**

# Utilities for Lagrange solver coupled with VOF method

## Precompiler identifier, variables and constants
*/

#if dimension == 1
#define normcoord(v) fabs(v.x)
#define normsq(v) sq(v.x)
#define prodscal(u,v) u.x*v.x
#elif dimension == 2
#define normcoord(v) (sqrt(sq(v.x) + sq(v.y)))
#define normsq(v) (sq(v.x) + sq(v.y)
#define prodscal(u,v) u.x*v.x + u.y*v.y
#else // dimension == 3
#define normcoord(v) (sqrt(sq(v.x) + sq(v.y) + sq(v.z)))
#define normsq(v) (sq(v.x) + sq(v.y) + sq(v.z))
#define prodscal(u,v) u.x*v.x + u.y*v.y + u.z*v.z 
#endif

// Constants
const double Cd = 0.45;
const double Cdpara = 0.5;
const double Cdperp = 0.5;
const double distgasvel = 1.5;

/**
time discretisation.\
With EXAKT_KONSTRELAXTIME a constant relaxation time is assumed,
so that the differential equation for the particles (position and velocity) can be solved analytically.
To compare with, an first order EXPLICIT sheme is available.*/

enum {EXAKT_KONSTRELAXTIME, EXPLICIT, RUNGE_KUTTA} timediscr;

/**
treatment of the point-droplets at bottom wall.*/

enum {ELASTICDEFLECTION} dropatwall;

/**
Models for the drag force for spherical droplets.*/

enum {SPHERICAL_DROP, SOLID_PART} dragforce;


/**
Two way coupling: source term from particles in Lagrange framework for momentum equation. 
There is the possibility to neglect the source with NOCOUPLING (one way coupling).
With GAUSSDISTRIBUTION we adopt the procedure of [Tomar et al.](https://doi.org/10.1016/j.compfluid.2010.06.018) 
with a gaussian distribution of particle momentum change during the time step around particle position.
As particle trajectory covers several cells, and extend over several droplet diameter, due to the Lagrange time step, 
the source term is spread over the whole trajectory.
*/

enum {NOCOUPLING, GAUSSDISTRIBUTION} twowaycoupling;


/**
Methods to reconstruct particles coming back to VOF:\

- REFINEFRACTION: refine region using high resolution and initialise droplets with fraction function. Accuracy of mass conservation controlled by number of cells per droplet`s diameter.

- FRACTIONMASSCONS: refine region using low resolution and correct the volume fractions obtained to conserve droplet`s mass.
*/

enum {REFINEFRACTION, FRACTIONMASSCONS} partreconstruction;


/** 
## derived datatype and further variables      
*/

/** derived datatype to store the neigbors of a cell where the source term needs to be evaluated.
This neighbors can be adressed with the indices _l, _m, _n with values {-1,0,1}.
Default values are: min_x = -1, max_x = 1, min_y = -1, max_y = 1, min_z = -1, max_z = 1 to adress all 27 neighbors.
*/

typedef struct {
  int min_x,min_y,min_z,max_x,max_y,max_z;
} ind_neighbors;

/** Derived datatype to store the relative indices of the neighbors of *point* containing the values to interpolate with.
With point.i, point.j, point.k the indices at point.level refinement,
there are 8 sets of indices for the trilinear interpolation in interpolate_vel, with i-shift.m_x or i+shift.p_x in x-direction,
j-shift.m_y or j+shift.p_y in y-direction, k-shift.m_z or k+shift.p_z in z-direction.

note: we do not combine the datatypes indices and ind_neighbors,
since they have a different meaning and their attributes have different ranges of values.
*/

typedef struct {
  int m_x, m_y, m_z, p_x, p_y, p_z;
} indices;

/**
The primary variables are a list of Nd particles */
#ifndef ADD_PART_MEM 
#define ADD_PART_MEM double dp,vol; coord oldpos; coord u; coord uf; coord cf; int Ninterval; coord spot[51]; ind_neighbors indngb[51]; \
int ibig[51]; int jbig[51]; int kbig[51]; int lbig[51]; bool backtoVOF; long unsigned int tag; Point point;
#endif

/**
*double x,y,z*;             particle position\
*double dp,vol*;            Particle diameter and volume\
*coord oldpos*;             old particle position\
*coord u*;                  particle velocity\
*coord uf*;                 velocity of the undisturbed surrounding gasphase\
*coord cf*;                 acceleration due to drag\
*int Ninterval*;            number of intervals between old and new particle position (momentum source term)\
*ind_neighbors indngb[11]*  relative indices of big neighbors for Ninterval+1 (maximum 11) positions between old and new particle position (momentum source term)\
*int ibig[11]*;             indices of big cell for Ninterval+1 (maximum 11) positions between old and new particle position (momentum source term)\
*int jbig[11]*; \
*int kbig[11]*; \
*int lbig[11]*; \
*bool backtoVOF*;           Initialised to false, true when (We>12 not yet implemented) or particle is near to VOF field.\
*long unsigned int tag;*    not used so far\
*Point point*;              cell including particle position
*/

/** Antoonvh sandbox [particle.h](https://basilisk.fr/sandbox/Antoonvh/particle.h) */
#include "lagrange/particle.h"
/** Antoonvh sandbox [adapt_wavelet2.h](https://basilisk.fr/sandbox/Antoonvh/adapt2.h), struct Adapt2 is removed */
#include "lagrange/adapt_wavelet2.h"

/** */

double dtlag = 1.;    // time step for Lagrange solver
bool lagstep = false;

int partexist = 0;
int partexistold = 0;

double fcut = 1e-3;   // minimum VOF value taken into account to identify cohesive droplets

int Nd = 0;           // Number of particles: set/change this number exclusively over functions for particle generation in common-lag.h
                      // to guarantee that Nd fits to the number of particles stored in particles.
long unsigned int Ndold = 0;    // Number of particles at previous time step.         

particle * myparticles = NULL;  

coord gravity = {0., 0., 0.};   // constant gravity over the computation domain, 0 per default.

vector acoupling[];        // momentum source term summed over the list of particles.
                           // added to gravity in event acceleration to obtain the acceleration *a*.

/**
## functions
 */

Particles * inertial_particles = NULL;

/** new_inertial_particles can be found in Antoon van Hooft sandbox  ([inertial-particles.h](https://basilisk.fr/sandbox/Antoonvh/inertial-particles.h)).
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

/** Several functions inspired by update_mpi in ([particle.h](https://basilisk.fr/sandbox/Antoonvh/particle.h)).
*/
#if _MPI
// checks if the big region around a given position located in the cell *Point* has local leaf cells in the processor domain.
bool check_bigcells(Point point, ind_neighbors indngb)
{
  bool reset = false;
  foreach_neighbor(1) {
    if (point.level <= depth() && _l>=indngb.min_x && _l<=indngb.max_x && _m>=indngb.min_y && _m<=indngb.max_y && _n>=indngb.min_z && _n<=indngb.max_z &&
        point.i>=2 && point.i<=pow(2,point.level)+1 && point.j>=2 && point.j<=pow(2,point.level)+1 && point.k>=2 && point.k<=pow(2,point.level)+1)
      if (allocated(0) && is_local(cell)) {
        reset = true;
      }
  }  
  return reset;
}

/** Collect all particles, and assign the particles to the list of particles for a given processor,
if one of the spots on particle trajectory is located in the processor domain 
and it was not in the list before (locate provides a negative level for particle position).
Used to generate particles for traverse_coupling_neighbors.
*/

trace
void update_mpi_locatespot_all (Particles p) {
  bool spotlocated;
  int l = 0;
  while (pn[l] != terminate_int) {
    if (l == p) {
      int newt, in = 0, m = 0, new = 0;
      foreach_particle_in(p) {
        new++;
      }
      // get new particle data
      particle * senddata = malloc (sizeof(particle)*new);
      foreach_particle_in(p) { 
        //ind[m] = _j_particle; particles are not removed, no need to store the index
        senddata[m++] = p();
      }

      // Gather new particles among threads: npe() is the number of threads. First, count all of them
      int newa[npe()], newat[npe()];
      newat[0] = 0;
      
      // new is the number of new particles per thread
      // newa stores these numbers in an array, and distributes the array to new threads
     
      MPI_Allgather (&new, 1, MPI_INT, &newa[0], 1, MPI_INT, MPI_COMM_WORLD);
      // Compute displacements
     
      for (int j = 1; j < npe(); j++) 
        newat[j] = newa[j - 1] + newat[j - 1];
      newt = newat[npe() - 1] + newa[npe() - 1]; 
      // Allocate receive buffer and gather
     
      particle * recdata = malloc (sizeof(particle)*newt);
      for (int j = 0; j < npe(); j++) {
        newat[j] *= sizeof(particle);
        newa[j]  *= sizeof(particle);
      }
      // send and receive data: each thread sends his lost particles
      // newa[pid()]: number of elements in send buffer (this is *new* also)
      // newa: integer array containing the number of elements to be received from each process
      // newat: integer array, entry i specifies the displacement at which to place the incoming data from process i
     
      MPI_Allgatherv (&senddata[0], newa[pid()], MPI_BYTE,
		      &recdata[0], newa, newat, MPI_BYTE,
		      MPI_COMM_WORLD);

      // count new particles
      for (int j = 0; j < newt ; j++) {
        spotlocated = false;
        for (int Ntraj = 0; Ntraj <= recdata[j].Ninterval; Ntraj++) {
          if (locate (recdata[j].spot[Ntraj].x, recdata[j].spot[Ntraj].y, recdata[j].spot[Ntraj].z).level >= 0)
            spotlocated = true;
        }
        // particles are in (new for the thread), if one of the big cell of the 3diam region is local, but the cell including particle position is not local.
        
        if (spotlocated && locate (recdata[j].x, recdata[j].y, recdata[j].z).level < 0) {
          in++;
        }
      }
      long unsigned int po = pn[l];   
      // Manage the memory if required...
      change_plist_size (l, in);
      // Collect new particles from `recdata`*/
      if (in > 0) {
        int indi[in];
        m = 0;
        for (int j = 0; j < newt; j++) {
          spotlocated = false;
          for (int Ntraj = 0; Ntraj <= recdata[j].Ninterval; Ntraj++) {
            if (locate (recdata[j].spot[Ntraj].x, recdata[j].spot[Ntraj].y, recdata[j].spot[Ntraj].z).level >= 0)
              spotlocated = true;
          }
          if (spotlocated && locate (recdata[j].x, recdata[j].y, recdata[j].z).level < 0) {
            indi[m++] = j;
          }
        }      
        m = 0;
        for (int j = po; j < pn[l]; j++) {
          pl[l][j] = recdata[indi[m]];
          m++;
        }
      }
      // clean the mess
      free (senddata); free (recdata);
    }
    l++;
  }
}

/** Collect every particles, and assign the particles to the list of particles for a given processor,
if the big region around at least one spot contains local leaf cells for that processor (check_bigcells is true for that spot),
but the spot is not *located* in the processor (locate for the position of the spot provides a level<0).
Used to generate particles for traverse_coupling_neighbors.
*/

trace
void update_mpi_spots3diamregion_all (Particles p) {
  bool bigcell_local, spotlocated;
  int l = 0;
  while (pn[l] != terminate_int) {
    if (l == p) {
      int newt, in = 0, m = 0, new = 0;
      foreach_particle_in(p) {
        new++;
      }
      particle * senddata = malloc (sizeof(particle)*new);
      foreach_particle_in(p) { 
        senddata[m++] = p();
      }
      int newa[npe()], newat[npe()];
      newat[0] = 0;
      MPI_Allgather (&new, 1, MPI_INT, &newa[0], 1, MPI_INT, MPI_COMM_WORLD);
      for (int j = 1; j < npe(); j++) 
        newat[j] = newa[j - 1] + newat[j - 1];
      newt = newat[npe() - 1] + newa[npe() - 1]; 
      particle * recdata = malloc (sizeof(particle)*newt);
      for (int j = 0; j < npe(); j++) {
        newat[j] *= sizeof(particle);
        newa[j]  *= sizeof(particle);
      }

      MPI_Allgatherv (&senddata[0], newa[pid()], MPI_BYTE,
		      &recdata[0], newa, newat, MPI_BYTE,
		      MPI_COMM_WORLD);

      for (int j = 0; j < newt ; j++) {
       Point ptparent;
        int newbigcells = 0;
        for (int Ntraj = 0; Ntraj <= recdata[j].Ninterval; Ntraj++) {
          if (recdata[j].lbig[Ntraj] >= 0) {
            ptparent.i = recdata[j].ibig[Ntraj];
            ptparent.j = recdata[j].jbig[Ntraj];
            ptparent.k = recdata[j].kbig[Ntraj];
            ptparent.level = recdata[j].lbig[Ntraj];
            bigcell_local = false;
            if (check_bigcells(ptparent,recdata[j].indngb[Ntraj]))
              bigcell_local = true;
            spotlocated = true;
            if (locate (recdata[j].spot[Ntraj].x, recdata[j].spot[Ntraj].y, recdata[j].spot[Ntraj].z).level < 0)
              spotlocated = false;
            if (bigcell_local && !spotlocated)
              newbigcells += 1;
          }
        }
        if (newbigcells > 0) {
          in++;
        }
      }
      long unsigned int po = pn[l];
      change_plist_size (l, in);

      if (in > 0) {
        int indi[in];
        m = 0;
        for (int j = 0; j < newt; j++) {
          Point ptparent;
          int newbigcells = 0;
          for (int Ntraj = 0; Ntraj <= recdata[j].Ninterval; Ntraj++) {
            if (recdata[j].lbig[Ntraj] >= 0) {
              ptparent.i = recdata[j].ibig[Ntraj];
              ptparent.j = recdata[j].jbig[Ntraj];
              ptparent.k = recdata[j].kbig[Ntraj];
              ptparent.level = recdata[j].lbig[Ntraj];
              bigcell_local = false;
              if (check_bigcells(ptparent,recdata[j].indngb[Ntraj]))
                bigcell_local = true;
              spotlocated = true;
              if (locate (recdata[j].spot[Ntraj].x, recdata[j].spot[Ntraj].y, recdata[j].spot[Ntraj].z).level < 0)
                spotlocated = false;
              if (bigcell_local && !spotlocated)
                newbigcells += 1;
            }
          }
          if (newbigcells > 0) {
            indi[m++] = j;
          }
        }      
        m = 0;
        for (int j = po; j < pn[l]; j++) {
          pl[l][j] = recdata[indi[m]];
          m++;
        }
      }
      free (senddata); free (recdata);
    }
    l++;
  }
}

/** Collect all the new particles transfered from VOF to the lagrange solver.
Add a particle to the list for a given processor if the big region around particle's position is local for the processor
(i.e. if it has local leafs), but particle position is not inside the processor domain (tested with locate).
In the event droplet_to_particle, this function is needed to reset the gas velocities to the undisturbed gas velocity after removal of the droplet,
for traverse_resetgasvel_neighbors.
*/

trace
int update_mpi_3diamregion (Particles p, long unsigned int Ndold) {
  int l = 0;
  int newin = 0;
  while (pn[l] != terminate_int) {
    if (l == p) {
      int newt, in = 0, m = 0, new = 0;
      foreach_particle_in(p) {
        if (_j_particle > Ndold-1) {
          new++;
        }
      }
      particle * senddata = malloc (sizeof(particle)*new);
      foreach_particle_in(p) { 
        if (_j_particle > Ndold-1) {
          senddata[m++] = p();
        }
      }
      int newa[npe()], newat[npe()];
      newat[0] = 0;
      MPI_Allgather (&new, 1, MPI_INT, &newa[0], 1, MPI_INT, MPI_COMM_WORLD);
      for (int j = 1; j < npe(); j++) 
        newat[j] = newa[j - 1] + newat[j - 1];
      newt = newat[npe() - 1] + newa[npe() - 1]; 
      particle * recdata = malloc (sizeof(particle)*newt);
      for (int j = 0; j < npe(); j++) {
        newat[j] *= sizeof(particle);
        newa[j]  *= sizeof(particle);
      }

      MPI_Allgatherv (&senddata[0], newa[pid()], MPI_BYTE,
		      &recdata[0], newa, newat, MPI_BYTE,
		      MPI_COMM_WORLD); 

      for (int j = 0; j < newt ; j++) {
        Point ptparent;
        ptparent.i = recdata[j].ibig[0];
        ptparent.j = recdata[j].jbig[0];
        ptparent.k = recdata[j].kbig[0];
        ptparent.level = recdata[j].lbig[0];
        if (check_bigcells(ptparent,recdata[j].indngb[0]) && locate (recdata[j].x, recdata[j].y, recdata[j].z).level < 0) {
          in++;
        }
      }
      long unsigned int po = pn[l];
      newin = in;
      change_plist_size (l, in);

      if (in > 0) {
        int indi[in];
        m = 0;
        for (int j = 0; j < newt; j++) {
          Point ptparent;
          ptparent.i = recdata[j].ibig[0];
          ptparent.j = recdata[j].jbig[0];
          ptparent.k = recdata[j].kbig[0];
          ptparent.level = recdata[j].lbig[0];
          if (check_bigcells(ptparent,recdata[j].indngb[0]) && locate (recdata[j].x, recdata[j].y, recdata[j].z).level < 0) {
            indi[m++] = j;
         }
        }      
        m = 0;
        for (int j = po; j < pn[l]; j++) {
          pl[l][j] = recdata[indi[m]];
          m++;
        }
      }
      free (senddata); free (recdata);
    }
    l++;
  }
  return newin;
}

/** For event particles_out in the inputfile:
processor 0 writes out all the particles (to the out file), so it has to collect them.
If each processor writes out its big number of particles to out, they could get mixed up.
*/

trace
int update_mpi_all (Particles p) {
  int l = 0;
  int newin = 0;
  while (pn[l] != terminate_int) {
    if (l == p) {
      int newt, in = 0, m = 0, new = 0;
      foreach_particle_in(p) {
        new++;
      }
      particle * senddata = malloc (sizeof(particle)*new);
      foreach_particle_in(p) { 
        senddata[m++] = p();
      }
      int newa[npe()], newat[npe()];
      newat[0] = 0;
      MPI_Allgather (&new, 1, MPI_INT, &newa[0], 1, MPI_INT, MPI_COMM_WORLD);
      for (int j = 1; j < npe(); j++) 
        newat[j] = newa[j - 1] + newat[j - 1];
      newt = newat[npe() - 1] + newa[npe() - 1]; 
      particle * recdata = malloc (sizeof(particle)*newt);
      for (int j = 0; j < npe(); j++) {
        newat[j] *= sizeof(particle);
        newa[j]  *= sizeof(particle);
      }

      MPI_Allgatherv (&senddata[0], newa[pid()], MPI_BYTE,
		      &recdata[0], newa, newat, MPI_BYTE,
		      MPI_COMM_WORLD); 

      for (int j = 0; j < npe(); j++) {
        newat[j] /= sizeof(particle);
        newa[j]  /= sizeof(particle);
      }
      if (pid() > 0) {
        for (int j = 0; j < newat[pid()] ; j++) {
          in++;
        }
      }
      if (pid() < npe() - 1) {
        for (int j = newat[pid()]+newa[pid()]; j < newt ; j++) {
          in++;
        }
      }
      long unsigned int po = pn[l];
      newin = in;
      change_plist_size (l, in);

      if (in > 0) {
        int indi[in];
        m = 0;
        if (pid() > 0) {
          for (int j = 0; j < newat[pid()] ; j++) {
            indi[m++] = j;
          }
        }
        if (pid() < npe() - 1) {
          for (int j = newat[pid()]+newa[pid()]; j < newt ; j++) {
            indi[m++] = j;
          }
        }  
        m = 0;
        for (int j = po; j < pn[l]; j++) {
          pl[l][j] = recdata[indi[m]];
          m++;
        }
      }
      free (senddata); free (recdata);
    }
    l++;
  }
  return newin;
}

#endif

/** Collect every particles that return from Lagrange to DNS, so that each processor knows all of them.
This is because each processor has to reach the functions refine and fraction.
*/

trace
long unsigned int collect_backtoVOFparticles (Particles p0, Particles p1, long unsigned int Ndold) {
  long unsigned int newin = 0;
#if _MPI
  int l = 0;
  int newt = 0;
  particle * recdata = malloc (0.);
  while (pn[l] != terminate_int) {
    if (l == p0) {
      int m = 0, new = 0;
      foreach_particle_in(p0) {
        if (p().backtoVOF && _j_particle < Ndold) {
          new++;
        }
      }
      newin = new;
      particle * senddata = malloc (sizeof(particle)*new);
      foreach_particle_in(p0) { 
        if (p().backtoVOF && _j_particle < Ndold) {
          senddata[m++] = p();
        }
      }
      int newa[npe()], newat[npe()];
      newat[0] = 0;
      MPI_Allgather (&new, 1, MPI_INT, &newa[0], 1, MPI_INT, MPI_COMM_WORLD);
      for (int j = 1; j < npe(); j++) 
        newat[j] = newa[j - 1] + newat[j - 1];
      newt = newat[npe() - 1] + newa[npe() - 1]; 
      if (newt > 0) {
        recdata = realloc (recdata , newt*sizeof(particle));
        for (int j = 0; j < npe(); j++) {
          newat[j] *= sizeof(particle);
          newa[j]  *= sizeof(particle);
        }

        MPI_Allgatherv (&senddata[0], newa[pid()], MPI_BYTE,
            &recdata[0], newa, newat, MPI_BYTE,
            MPI_COMM_WORLD); 
        free (senddata);
      }
    }
    if (l == p1) {
      int n_partn = newt; 
      if (newt > 0) {
        if (n_partn > pna[l] || 2*(n_partn + 1) < pna[l]) {
          pna[l] = 2*(n_partn + 1);
          pl[l] = realloc (pl[l] , pna[l]*sizeof(particle));
        }
        for (int j = 0; j < n_partn; j++) {
          pl[l][j] = recdata[j];
        }
      }
      pn[l] = n_partn;
      free (recdata);
    }
    l++;
  }
#else
    int m = 0;
    if (Ndold > 0) {
      for (int j=0; j < Ndold; j++) {
        if (pl[0][j].backtoVOF) {
          m++;
        }
      }
    }
    pl[1] = realloc (pl[1] , m*sizeof(particle));
    if (m > 0) {
      m = 0;
      for (int j=0; j < Ndold; j++) {
        if (pl[0][j].backtoVOF) {
          pl[1][m++] = pl[0][j];
        }
      }
    }
    pn[1] = m;
    newin = m;
#endif  
  return newin;
}

/**

## Own functions for lagrange.h

drag coefficient for spherical droplets (p. 112 in the [book](https://doi.org/10.1017/S0022112079221290)).
The drag coefficient $c_d$ is given as function of the Reynolds number $Re_p=\dfrac{d_p\rho_g\left|\mathbf{u}_{rel}\right|}{\mu_g}$ 
for spherical droplets in different intervals:
$$	c_d=\frac{24}{Re_p}+4.5 \text{ for } Re_p < 0.006 $$
$$	c_d=\frac{24}{Re_p}\left(1+0.1315\cdot Re_p^{0.82-0.05w}\right) \text{ for } 0.006 < Re_p < 26 $$
$$	c_d=\frac{24}{Re_p}\left(1+0.1935\cdot Re_p^{0.6305}\right) \text{ for } 26 < Re_p < 259.7889 $$
and further:
$$	c_d=\text{exp}\left( 3.78430-2.58857w+0.35874w^2\right) \text{ for } 259.7889 < Re_p < 1526.087 $$
$$	c_d=\text{exp}\left( -5.65768+5.88495w-2.14025w^2+0.24154w^3\right) \text{ for } 1526.087 < Re_p < 1.200229\cdot 10^4 $$
$$	c_d=\text{exp}\left( -4.41659+1.46675w-1.4644w^2\right) \text{ for } 1.200229\cdot 10^4 < Re_p < 4.410367\cdot 10^4 $$
$$	c_d=\text{exp}\left(-9.99092-3.64016w-0.35598w^2\right) \text{ for } 4.410367\cdot 10^4 < Re_p < 3.384249\cdot 10^5 $$
$$	c_d= 29.78-5.3w \text{ for } 3.384249\cdot 10^5 < Re_p < 4.032325\cdot 10^5 $$
$$	c_d= -0.49+0.1w \text{ for } 4.032325\cdot 10^5 < Re_p < 10^6 $$
$$	c_d= 0.19-8\cdot 10^4/Re_p \text{ for } Re_p > 10^6 $$
with $w=log_{10}(Re_p)$.
*/

double cd_sphere(double Re)
{
    double w = 0.;
    if (Re >= 259.7889 && Re <= 1.e+06) 
      w = log10(Re);
    
    if (Re <= 259.7889){
        if (Re > 26.)
            return (24./Re)*(1. + 0.1935*pow(Re,0.6305));
        else if (Re > 0.006)
            return (24./Re)*(1. + 0.1315*pow(Re,0.82-0.05*w));
        else if (Re > 0.)
            return 24./Re + 4.5;
        else 
            return 0.;
    }
    else if (Re <= 1526.087)
        return exp( 3.78430 - 2.58857*w + 0.35874*w*w); 
    else if (Re < 1.200229e+04)
        return exp(-5.65768 + 5.88495*w - 2.14025*w*w + 0.24154*w*w*w);
    else if (Re <= 4.410367e+04)
        return exp(-4.41659 + 1.46675*w - 0.14644*w*w);
    else if (Re < 3.384249e+05)
        return exp(-9.99092 + 3.64016*w - 0.35598*w*w);
    else if (Re <= 4.032325e+05)
        return 29.78 - 5.3*w;                        
    else if (Re < 1.e+06)
        return -0.49 + 0.1*w;                         
    else
        return  0.19 - 8.e+04/Re;                            
}

bool refinecond (int cellsperdiam, double x, double y, double z, double Delta, int level, Particles p)
{
  bool cond = false;
  double invlog2 = 1./log(2.);
  for (int j=0; j < pn[p]; j++) {  
    int levelref = ceil (invlog2*log(sqrt(3.)*cellsperdiam*L0/pl[p][j].dp));
    cond = cond || ( sqrt(sq(x-pl[p][j].x) + sq(y-pl[p][j].y) + sq(z-pl[p][j].z)) - 0.5*pl[p][j].dp < 0.5*Delta && level < levelref );
  }
  return cond;
}

void add_droplets (int n, double v[n], coord droploc[n], coord dropvel[n], Particles Plist) {
  // pid=0 collects the droplets
  if (pid() == 0) {
    change_plist_size (Plist, n);
    for (int j = 0; j < n; j++) {
      pl[Plist][j].vol = v[j];
      pl[Plist][j].dp = pow(6.*v[j]/pi,1./3.);
      foreach_dimension(){
        pl[Plist][j].x = droploc[j].x;
        pl[Plist][j].u.x = dropvel[j].x;  
      }
    }
  }
  return;
}

/** Functions with Point as an arguments: not sure if it is the elegant way!*/

Point return_parent_ifpointexists(Point point)
{
  return parent;
}

Point return_parent(Point point)
{
  Point ptparent;
  ptparent.level = -1;

  if (point.i >= 0 && point.i < (1 << point.level) + 2*GHOSTS &&
      point.j >= 0 && point.j < (1 << point.level) + 2*GHOSTS &&
      point.k >= 0 && point.k < (1 << point.level) + 2*GHOSTS)
    if (allocated(0))
      ptparent = parent;

  return ptparent;
}

double return_value(Point point, scalar s)
{
  return s[];
}

void set_value(Point point, scalar s, double value)
{
  s[] = value;
  return;
}

bool test_leaf(Point point)
{
  if (allocated(0)) {
    if (is_leaf(cell))
      return true;
    else
      return false;
  }
  else
    return false;  
}

double interpolate_vof (scalar v, double xp = 0., double yp = 0., double zp = 0.)
{
  double val = nodata;
  foreach_point (xp, yp, zp, reduction (min:val))
    val = ( point.level >= 0 ? 
            (interfacial (point, v) ? interpolate_linear (point, v, xp, yp, zp) : return_value(point, f)) 
            : interpolate_linear (point, v, xp, yp, zp));
  return val;
}

/** 

### For event droplet_to_particle  */

/** Inspired from interpolate_linear in cartesian-common.h, this function interpolates the scalar field p.v at position p.x,p.y,p.z,
but not from the values at nearest cell centers like in interpolate_linear.
Instead, the indices of the cells with the values used to interpolate are *shifted* from the indices i,j,k of *point* (cell containing the position p.x,p.y,p.z)
by integer values stored in the derived type indices *shift*.
*/

double interpolate_vel (Point point, scalar v, double xp = 0., double yp = 0., double zp = 0., indices shift)
{
#if dimension == 1
  int im = shift.m_x;
  int ip = shift.p_x;
  x = (im ==0 && ip == 0) ? 0. : (xp - x)/(((double)(ip+im))*Delta) + im/((double)(ip+im));
  /* linear interpolation */
  return v[-im]*(1. - x) + v[ip]*x;
#elif dimension == 2
  int im = shift.m_x, jm = shift.m_y;
  int ip = shift.p_x, jp = shift.p_y;
  x = (im ==0 && ip == 0) ? 0. : (xp - x)/(((double)(ip+im))*Delta) + im/((double)(ip+im));
  y = (jm ==0 && jp == 0) ? 0. : (yp - y)/(((double)(jp+jm))*Delta) + jm/((double)(jp+jm));
  /* bilinear interpolation */
  return ((v[-im,-jm]*(1. - x) + v[ip,-jm]*x)*(1. - y) + 
	  (v[-im,jp]*(1. - x) + v[ip,jp]*x)*y);
#else
  // With the definition of the shift indices in set_undistflowvelMPI, xm=xp etc
  // Nevertheless we keep the possibility of setting different indices,
  // to handle correctly non-spherical droplets in the future.
  int im = shift.m_x, jm = shift.m_y, km = shift.m_z;
  int ip = shift.p_x, jp = shift.p_y, kp = shift.p_z;
  x = (im ==0 && ip == 0) ? 0. : (xp - x)/(((double)(ip+im))*Delta) + im/((double)(ip+im));
  y = (jm ==0 && jp == 0) ? 0. : (yp - y)/(((double)(jp+jm))*Delta) + jm/((double)(jp+jm));
  z = (km ==0 && kp == 0) ? 0. : (zp - z)/(((double)(kp+km))*Delta) + km/((double)(kp+km));

  return (((v[-im,-jm,-km]*(1. - x) + 
        v[ip,-jm,-km]*x)*(1. - y) + 
        (v[-im,jp,-km]*(1. - x) + 
        v[ip,jp,-km]*x)*y)*(1. - z) +
        ((v[-im,-jm,kp]*(1. - x) +
         v[ip,-jm,kp]*x)*(1. - y) + 
          (v[-im,jp,kp]*(1. - x) + 
          v[ip,jp,kp]*x)*y)*z);
  
#endif  
}

/** Determines the velocity of the undisturbed (by the droplet) flow velocity around the droplet transformed to a Lagrange particle.
First, determine iteratively the first parent cell (i,j,k,level) of the leaf cell containing droplets center-of-mass, 
for which 8 neighbor cells defined by 6 *shifts* respective to i,j,k (3 space directions, in positive and negative direction),
inlude the 8 intersections of the surface of a sphere of radius distgasvel*diam
with 4 lines through droplets center of mass, and supported by the vectors (1,1,1), (1,1,-1), (1,-1,1), (1,-1,-1),
It is required that the shifts do not exceed the value of 2.
Take the velocity values in the 8 cells to perform a trilinear interpolation at droplet position.
The neighbor cells are not further than the 2nd neighbors, so that no MPI exchange is necessary, since the 2nd neighbors of each cell 
containing local leafs is known for a given processor.
*/

coord set_undistflowvelMPI(Point point, vector u, coord droploc, double diam, coord dropvel){

  coord uflow;
  indices shift;
  double distance = distgasvel*diam/sqrt(3.);
  int levelinterpol = point.level;
    
  printf ("pid:%d, set_undistflowvelMPI, vel before: %g %g %g\n",pid(),u.x[],u.y[],u.z[]);

  int maxshift = 0;
  coord p = {x,y,z};
  foreach_dimension(){
    shift.p_x = round((distance-(p.x-droploc.x))*pow(2,point.level)/L0);
    if (shift.p_x > maxshift) maxshift = shift.p_x;
    shift.m_x = round((distance+(p.x-droploc.x))*pow(2,point.level)/L0);
    if (shift.m_x > maxshift) maxshift = shift.m_x;
  }

  while (maxshift > 2) {
    levelinterpol--;
    Point parentpoint = return_parent_ifpointexists(point);
    point = parentpoint;
    maxshift = 0;
    coord p = {x,y,z};
    foreach_dimension(){
      shift.p_x = round((distance-(p.x-droploc.x))*pow(2,point.level)/L0);
      if (shift.p_x > maxshift) maxshift = shift.p_x;
      shift.m_x = round((distance+(p.x-droploc.x))*pow(2,point.level)/L0);
      if (shift.m_x > maxshift) maxshift = shift.m_x;
    }
  }

  double p1 = droploc.x;
  double p2 = droploc.y;
  double p3 = droploc.z;   
  foreach_dimension()
    uflow.x = interpolate_vel (point, u.x, p1, p2, p3, shift);

  printf ("pid:%d, set_undistflowvelMPI, interpolated vel: %g %g %g\n",pid(),uflow.x,uflow.y,uflow.z);

  return uflow;
}

/** This function checks if the range of the gaussian distribution for which the acceleration *acoupling* is evaluated,
-this is a spherical region with radius *radregion* around position partloc- is inside a given cell defined by *point* *from one side*.
More precisely and for example, if partloc.x > x(cell center), we check if partloc.x - radregion > x - 0.5*Delta,
the left boundary (partloc.x - radregion) of the spherical range is inside the cell (x - 0.5*Delta is the position of the *left* face).
If it is not the case for this cell *point*, the parent cell is returned to be checked.
Position partloc is always inside the cell, since the first cell in this recursive search is the leaf cell containing partloc.
*/

Point parent_to_try(Point point, double radregion, coord partloc)
{
  if ((partloc.x > x ? partloc.x - radregion > x - 0.5*Delta : partloc.x + radregion < x + 0.5*Delta) &&
      (partloc.y > y ? partloc.y - radregion > y - 0.5*Delta : partloc.y + radregion < y + 0.5*Delta) &&
      (partloc.z > z ? partloc.z - radregion > z - 0.5*Delta : partloc.z + radregion < z + 0.5*Delta)){
    point.level = -1; // here the level of the input parameter point is changed, and is returned! But hopefully not changed outside the function!
    return point;
  }
  else
    return parent;
}

/** Determines the neighbors of *point* including the desired range radregion=3*sigma around partloc where we want to evaluate the acceleration,
reset the gas velocity after removal of the droplet, or determine the volume fraction in this big region to evaluate the proximity to the fluid phase.
*/

ind_neighbors big_neighbors(Point point, double radregion, coord partloc)
{
  ind_neighbors indngb = {-1,-1,-1,1,1,1};
  coord pos = {x,y,z};
  foreach_dimension(){
    if (partloc.x + radregion < pos.x + 0.5*Delta) indngb.max_x = 0;
    if (partloc.x - radregion > pos.x - 0.5*Delta) indngb.min_x = 0;
  } 
  return indngb;
}

/** volume fraction of fluid in the big region defined by the big cell *point* and its neighbors *indngb* around a droplet.
All these cells define the smallest box including the sphere centered at droplets/particles center of mass,
with a radius radregion.
*/

double VOFintheBox(Point point, ind_neighbors indngb, double dropvol)
{
  int counter = 0;
  double fdensity = 0.;
  foreach_neighbor(1){
    if (_l>=indngb.min_x && _l<=indngb.max_x && _m>=indngb.min_y && _m<=indngb.max_y && _n>=indngb.min_z && _n<=indngb.max_z 
          && point.i>=2 && point.i<=pow(2,point.level)+1 && point.j>=2 && point.j<=pow(2,point.level)+1 && point.k>=2 && point.k<=pow(2,point.level)+1)
      {
        counter++;      
        fdensity += return_value(point,f);
      }
  }
  return ((counter == 0) ? 0. : (fdensity*dv()-dropvol)/(dv()*counter));
}

/** Find recursively the smallest parent cell of point containing the sphere 
(centered at particle position *partloc* with radius *radregion*) *from one side* in each space direction.
*/

Point det_bigcell (Point point, coord partloc, double radregion)
{
  bool parentfound = false;
  Point nextparent;

  while (!parentfound && point.level > -1) {
    nextparent = parent_to_try(point,radregion,partloc);
    if (nextparent.level == -1){
      parentfound = true;
    }
    else{
      point = nextparent;
    } 
  }

  return point;
}

/** Set volume fraction 0 < f < fcut to zero,
and add removed value to the volume of the corresponding particle with index j in neighbor cells,
to ensure mass conservation.
*/

double setValZeroIfSmall(Point point, int j){
  double smallvol = 0.;
  if (f[] > 0 && f[] <= fcut){
    smallvol = f[]*dv();
    f[] = 0.;
  }
  return smallvol;
}

/** Check neighbor leaf cells of a droplet cell for volume fraction values 0 < f < fcut and remove them. */

double setSmallNgbrZero(Point point, int j){
  
  double smallvol = 0.;
  if (is_coarse())
    foreach_child()
      foreach_neighbor(2)
        if (allocated(0))
          if(is_leaf(cell))
            smallvol += setValZeroIfSmall(point,j);
   
  foreach_neighbor(1){
    if (is_leaf(cell)){
      smallvol += setValZeroIfSmall(point,j);
    }
    else{
      if (!is_refined(cell)){
        Point parentpoint=return_parent_ifpointexists(point);
        if (test_leaf(parentpoint))   
            smallvol += setValZeroIfSmall(point,j);
      }
    }
  }
  return smallvol;
}

/** See traverse_resetgasvel_neighbors: here we reset the velocities in the leaf cells found in traverse_resetgasvel.
This values at cell centers are interpolated from cells at a distance distgasvel*diam from posrel.
Similar to set_undistflowvelMPI, we don't want to go further than the second neighbors of a parent cell to avoid MPI exchange.
*/

coord reset_gasvelMPI (Point point, vector ua, coord posrel, double radregion)
{
  coord uflow, p = {x,y,z};
  int levelinterpol = point.level;
  indices shift;

  if (normcoord(posrel) == 0.) {
    // This happens when the droplet has values for the volume fraction in one single cell.
    shift.p_x = shift.m_x = shift.p_y = shift.m_y = shift.p_z = shift.m_z = 2;
  }
  else{
    int maxshift = 0;
    foreach_dimension(){
      shift.p_x = round((radregion*(fabs(posrel.x)/normcoord(posrel))-posrel.x)*pow(2,point.level)/L0);
      if (shift.p_x > maxshift) maxshift = shift.p_x;
      shift.m_x = round((radregion*(fabs(posrel.x)/normcoord(posrel))+posrel.x)*pow(2,point.level)/L0);
      if (shift.m_x > maxshift) maxshift = shift.m_x;
    }

    while (maxshift > 2) {
      levelinterpol--;
      Point parentpoint = return_parent_ifpointexists(point);
      point = parentpoint;
      maxshift = 0;
      foreach_dimension(){
        // TODO: limit xm,xp to the boundary of computation domain
        shift.p_x = round((radregion*(fabs(posrel.x)/normcoord(posrel))-posrel.x)*pow(2,point.level)/L0);
        if (shift.p_x > maxshift) maxshift = shift.p_x;
        shift.m_x = round((radregion*(fabs(posrel.x)/normcoord(posrel))+posrel.x)*pow(2,point.level)/L0);
        if (shift.m_x > maxshift) maxshift = shift.m_x;
      }
    }
  }  

  // Do not wonder why the interpolated velocity is the same than the original when shift.p_x=shift.p_y=shift.p_z=0
  // or shift.m_x=shift.m_y=shift.m_z=0 that's normal!
  double p1 = p.x;
  double p2 = p.y;
  double p3 = p.z;
  foreach_dimension()
    uflow.x = interpolate_vel (point, u.x, p1, p2, p3, shift);

  return uflow;
}

/** See traverse_resetgasvel_neighbors: here we identify the leaf cells recursively.
*/

void traverse_resetgasvel(Point point, double radregion, coord partloc, vector ua)
{
  if (is_leaf(cell))
  {
    if (is_local(cell)) {
      // We interpolate the velocity at cell center from 8 values located at a radregion distgasvel*(droplets diameter) from droplet's center-of-mass.
      coord p = {x,y,z};
      coord posrel;
      foreach_dimension()
        posrel.x = p.x - partloc.x;
      if (normcoord(posrel) < radregion) {
        coord u_new = reset_gasvelMPI(point,ua,posrel,radregion);
        foreach_dimension()
          u.x[] = u_new.x;
      }
    }
  }
  else{
    if (is_coarse())
      if (is_local(cell)) {
        foreach_child(){
          traverse_resetgasvel(point,radregion,partloc,ua);
      }
    }
  }
}

/** Reset velocity of surrounding gas where the droplet was, and where the disturbed region around the old droplet is.
This region is the same sphere centered at particle center of mass *partloc* and radius *radregion* used to calculate the undisturbed flow velocity *uf*,
stored in the particle datatype. First we traverse all the neighbors of point.
*/

trace
void traverse_resetgasvel_neighbors(Point point, vector ua, coord partloc, double radregion, ind_neighbors indngb)
{     
  foreach_neighbor(1){
    if (point.level <= depth() && _l>=indngb.min_x && _l<=indngb.max_x && _m>=indngb.min_y && _m<=indngb.max_y && _n>=indngb.min_z && _n<=indngb.max_z &&
        point.i>=2 && point.i<=pow(2,point.level)+1 && point.j>=2 && point.j<=pow(2,point.level)+1 && point.k>=2 && point.k<=pow(2,point.level)+1){    
      if (allocated(0))
        if (is_local(cell)) {
          traverse_resetgasvel(point,radregion,partloc,ua);
        }
    }
  }
}

/**
### For event update_flowvel_coupling_force:
 */

void set_coupling(Point point, double sigma, coord partloc, coord partforce, coord p)
{
  double dist,twosigmasq;
  coord distvec;
  foreach_dimension()
    distvec.x = partloc.x - p.x;
  dist = normcoord(distvec);
  if (dist < 3.*sigma ){
    twosigmasq = 2.*sq(sigma);
    foreach_dimension()
      acoupling.x[] += (partforce.x/rho2) * exp(- sq(dist)/twosigmasq) / pow(pi*twosigmasq,1.5);
// Comment: leaf cells inside bigregions are set by only one processor, there is no risk to set them twice.
// Moreover, since regions from different spots, or different particles maybe, can ovelapp, different contributions to acoupling are summed up. 
  }
}

/** Evaluate the acceleration at *point's* cell faces
with a gaussian distribution around position *partloc*, the width is *sigma*.
*/

void traverse_coupling(Point point, double sigma, coord partloc, coord partforce)
{
  if (is_leaf(cell))
  {
    if (is_local(cell)) {
      coord pos = {x,y,z};
// We use sigma=dp.
      set_coupling(point,sigma,partloc,partforce,pos); 
    }
  }
  else{
    if (is_coarse())
      if (is_local(cell)) {
        foreach_child(){
          traverse_coupling(point,sigma,partloc,partforce);
      }
    }
  }
}

/** Evaluate the acceleration for every neighbors of *point* including the desired spherical range 3*sigma around partloc.
*/

trace
void traverse_coupling_neighbors(Point point, double sigma, coord partloc, coord partforce, ind_neighbors indngb)
{
  foreach_neighbor(1){
    if (point.level <= depth() && _l>=indngb.min_x && _l<=indngb.max_x && _m>=indngb.min_y && _m<=indngb.max_y && _n>=indngb.min_z && _n<=indngb.max_z &&
        point.i>=2 && point.i<=pow(2,point.level)+1 && point.j>=2 && point.j<=pow(2,point.level)+1 && point.k>=2 && point.k<=pow(2,point.level)+1){    
      if (allocated(0))
        if (is_local(cell)) {   
          traverse_coupling(point,sigma,partloc,partforce);
        }
    }
  }
}
