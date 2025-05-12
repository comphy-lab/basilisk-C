/**
# Lagrangian particle advection

`n_part` Partciles can be seeded in an array of `coord`s called
`loc`. The formulations on this page aim to advect it in a time loop
by an external vector field `u`. There is an option for a
third-order-accurate time-marhcing scheme, 2nd order is the default.
*/
#include "run.h"
long unsigned int n_part;  //Number of particles (foreach MPI thread) 
coord * loc;               //coordinates of particles
bool P_RK3 = false;        //Switch for RK-3 scheme
/**
## Particle methods

First, some macros are defined to create a `foreach_particle()`
iterator. It will give access to the positions (`x,y,z`), their
writeable counterparts (`p_x()`, etc), velocity components (`k1_x()`,
etc) and other stuff. It can even do simple reductions (well done
`qcc`!).
*/

  @define ArrIP_x (j*dimension + cind.x)
  @define p_x()  loc[j].x
  @define pn_x() loc_n[j].x
  @define k1_x() k1[j*dimension + cind.x]
  @define k2_x() k2[j*dimension + cind.x]
  @define ArrIP_y (j*dimension + cind.y)
  @define p_y()  loc[j].y
  @define pn_y() loc_n[j].y
  @define k1_y() k1[j*dimension + cind.y]
  @define k2_y() k2[j*dimension + cind.y]
  @define ArrIP_z (j*dimension + cind.z)
  @define p_z()  loc[j].z
  @define pn_z() loc_n[j].z
  @define k1_z() k1[j*dimension + cind.z]
  @define k2_z() k2[j*dimension + cind.z]


@def PARTICLE_VARIABLES  
  bool LN = false;
  if (loc_n != NULL)
    LN = true;
  double x = loc[j].x; NOT_UNUSED(x);
  double xn = 0;
  if (LN)
    xn = loc_n[j].x;
  NOT_UNUSED(xn);
  double y = loc[j].y; NOT_UNUSED(y);
  double yn = 0;
  if (LN)
    yn = loc_n[j].y;
  NOT_UNUSED(yn);
  double z = loc[j].z; NOT_UNUSED(z);
  double zn = 0;
  if (LN)
    zn = loc_n[j].z;
  NOT_UNUSED(zn);
@
/**
In order to facilitate the higher-order advection, some storage arrays
are declared. 
 */
long unsigned int n_part_a;                //allocated size of arrays
coord * loc_n;                             //Location storage
double *k1, *k2, dtf[2];                   //Velocities and timesteps
bool part_linear = true, start = true;     //Linear interpolation and RK3 first iteration 
double tol = 1.e-2;                        //RK3: Tolerance on alpha != 2/3  

macro foreach_particle() {
  for (int j = 0; j < n_part; j++) {
    PARTICLE_VARIABLES;
    {...}
  }
}

typedef struct {
  int x, y, z;
} coordi;
coordi cind = {0,1,2};
/**
The above data arrays should be little concern to the user as the code
tries to manage these scratch arrays in the background.
 */
event init (t = 0) {
  if (loc_n == NULL) {//We have not already been here
    n_part_a = n_part + 1;
#if _MPI
    n_part_a *= 2;
    loc = realloc (loc, n_part_a*sizeof(coord));
#endif
    loc_n = malloc (n_part_a*sizeof(coord));
    k1 = malloc (n_part_a*dimension*sizeof(double));
    k2 = malloc (n_part_a*dimension*sizeof(double));
  }
}

event set_dtmax (i++);

void part_boundaries () {
  coord mind = {X0, Y0, Z0}; 
  foreach_particle() { 
    foreach_dimension() {
      if (p_x() < mind.x) 
	p_x() += L0;
      else if (p_x() > (mind.x + L0))
	p_x() -= L0;
    }
  }
}

#if _MPI
/**
   `Interpolate_array` is convenient. We re-implement it, assuming all
   local particles.
*/
void update_mpi (int step); //particle exchange function prototype
#endif
trace
void interpolate_array_local (scalar * list, coord * a, int n, double * v, bool linear) {
  int j = 0;
  for (int i = 0; i < n; i++)
    foreach_point (a[i].x, a[i].y, a[i].z,serial)
      for (scalar s in list)
	v[j++] =interpolate_linear (point,s, a[i].x, a[i].y, a[i].z);
}
#define interpolate_array interpolate_array_local

/**
## Advection scheme

To facilitate the general case of time-dependend flow, particles are
advected in two iterations between $t_{n-2}$ and $t_n$. The
formulation follows a variable two-stage Runge-Kutta scheme:

$$ \begin{array}{c|cc}
0 & 0 & \\ 
\alpha & \alpha &  \\ 
\hline
  & 1 - \frac{1}{2\alpha} & \frac{1}{2\alpha}
\end{array}
$$

Which can be extended to the RK-3 scheme of Sanderse and Veldman
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

B. Sanderse and AEP Veldman, *Constraint-consistent Runge-Kutta
methods for one-dimensional incompressible multiphase flow*,
J. Comp. Phys. (2019) [link to
JCP.](https://www.sciencedirect.com/science/article/pii/S0021999119300683).
*/

void RK_step1 (coord * loc, coord * loc_n, double * k1, double dtf[2]) {
  interpolate_array ((scalar*){u}, loc, n_part, k1, part_linear);
  foreach_particle() {
    foreach_dimension() {
      pn_x() = p_x(); //store the locations at t_n 
      p_x() += dtf[0]*k1_x();
    }
  }
}

void RK_step2 (coord * loc, coord * loc_n, double * k1, double * k2, double dtf[2]) {
   interpolate_array ((scalar*){u}, loc, n_part, k2, part_linear);
  double a1 = -1, a2 = 2, h = dtf[1] + dtf[0];
  if (dtf[1] != dtf[0] || !P_RK3) {
    double c = dtf[0]/h;
    if (fabs (c - 2./3.) > tol && P_RK3) 
      a2 = (c - 1.)/(c*(3.*c - 2.));
    else //Raltson's 2nd order method
      a2 = 1./(2*c);
    a1 = 1 - a2;
  }
  foreach_particle() 
    foreach_dimension() 
      p_x() = pn_x() + h*(a1*k1_x() + a2*k2_x());
}

void RK_step3 (coord * loc, coord * loc_n, double * k1, double * k2, double dtf[2]) {
  double h = dtf[1] + dtf[0];
  double c = dtf[0]/h;
  if (fabs(c - 2./3.) > tol) {// RK-3
    double b1 = 0.5 - 1./(6.*c);
    double b2 = 1./(6.*c*(1. - c));
    double b3 = 1. - (b1 + b2);
    double V[dimension*n_part];
    interpolate_array ((scalar*){u}, loc, n_part, V, true);
    foreach_particle() {
      foreach_dimension() { 
        p_x() = pn_x() + h*(b1*k1_x() + b2*k2_x() + b3*V[ArrIP_x]);
      }
    }
  }
}
/**
Particle advection is performed in the `advance_particles`
event. Varying between even and uneven iterations.
 */
event advance_particles (i++, last) {
  part_boundaries();
  if (i%2 == 0) {
    if (i > 0 && P_RK3) {
#if _MPI
      update_mpi(3);
#endif
      RK_step3 (loc, loc_n, k1, k2, dtf);
      part_boundaries();
    }
#if _MPI
    update_mpi(1);
#endif
    dtf[0] = dt;
    RK_step1 (loc, loc_n, k1, dtf);
  } else { 
#if _MPI
    update_mpi(2);
#endif
    dtf[1] = dt; 
    RK_step2 (loc, loc_n, k1, k2, dtf);
  }
}

event free_particles (t = end, last) {
  free (loc);
  free (loc_n);
  free (k1);
  free (k2);
  loc_n = NULL;
}

/**
## User functions for particle seeding
   
   The following function initializes a particle at the centre of each
   grid cell.
*/
void init_particles_in_cells(){
  n_part = 0;
  foreach(serial)
    n_part++;
  loc = malloc (n_part*sizeof(coord));
  int n = 0;
  foreach(serial) {
    coord cc = {x, y, z};
    foreach_dimension()
      loc[n].x = cc.x;
    n++;
  }
}
/**
   Initialize particles in a 2D `l`$\times$`l` grid centered at {`xo,
   yo`} with `n`$\times$`n` particles:
*/
void init_particles_2D_square_grid (int n, double xo, double yo, double l){
  n_part = 0;
  if (pid() == 0) {
    n_part = sq(n);
    loc = malloc (n_part*sizeof(coord));
    int i = 0;
    for (int j = 0; j < n; j++) {
      for (int k = 0; k < n; k++) {
	loc[i].x = xo - l/2. + (double)j*(l/((double)n - 1.));
	loc[i].y = yo - l/2. + (double)k*(l/((double)n - 1.));
	i++;
      }
    }
  } else
    loc = malloc (sizeof(coord));
}
/**
   The following function places `n` particles randomly in a circle with 
   radius `R` at location {`xo, yo`}. 
*/
void init_particles_random_circle (int n, double xo, double yo, double R) {
  n_part = 0;
  if (pid() == 0) {
    n_part = n;
    loc = malloc (n_part*sizeof(coord));
    int j = 0;
    while (j < n_part) {
      double a = noise();
      double b = noise();
      if (sq(a) + sq(b) < R) {
	loc[j].x = a + xo;
	loc[j].y = b + yo;
	j++;
      }
    }
  } else
    loc = malloc (sizeof(coord));
}


/**
## Visualization 

   For visualization with bview, the `scatter()` function can be used
   when `BVIEW` $\neq 0$.

use something like:

~~~literatec
...
#include "view.h"
#define BVIEW 1
#include "particles.h"
...
~~~
*/ 
#if BVIEW
#include "scatter.h"
#endif
/**
## Implementation of the MPI particle exchange:  

The ugliest bit is saved for last. Local partciles outside the MPI
domain are comunicated globally. In turn, each thread selects those
inside their domain. This facilitates advection at large CFL numbers,
grid adaptation and rebalancing.
 */

#if _MPI
void update_mpi (int step) {
  int outt, out = 0, in = 0, m = 0;
  int psize = step < 2 ? 1 : step < 3 ? 3 : 4;
  psize *= dimension;
  //Count the number of outgoing particles per thread
  foreach_particle() 
    if (locate (x, y, z).level < 0)
      out++;
  //get indices and outgoing data
  int ind[out];
  double senddata[out*psize];
  foreach_particle() { 
    if (locate (x, y, z).level < 0) {
      ind[m] = j;
      int c = 0;
      foreach_dimension()
	senddata [m*psize + c++] = loc[j].x;
      if (step > 1) {
	foreach_dimension() {
	  senddata [m*psize + c++] = loc_n[j].x;
	  senddata [m*psize + c++] = k1[j*dimension + cind.x];
	}
      }
      if (step > 2) {
	foreach_dimension() {
	  senddata [m*psize + c++] = k2[j*dimension + cind.x];
	}
      }
      m++;
    }
  }
  //remove the senddata from arrays (shrink)
  m = 0;
  int j = 0;
  fflush (stdout);
  while (j < n_part - out) {
    while (m < out ? j + m == ind[m] : 0)
      m++;
    while (m < out ? j < n_part - out && j + m != ind[m] : j < n_part - out) {
      loc[j]   = loc[j + m];
      if (step > 1) {
	loc_n[j]   = loc_n[j + m];
	foreach_dimension()
	  k1[j*dimension + cind.x]   =  k1[(j + m)*dimension + cind.x];
      }
      if (step > 2) {
	foreach_dimension()
	  k2[j*dimension + cind.x]   =  k2[(j + m)*dimension + cind.x];
      }
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
  //compute total
  outt = outat[npe() - 1] + outa[npe() - 1]; 
  // Allocate recieve buffer and gather
  double recdata[outt*psize];
  for (int j = 0; j < npe(); j++) {
    outat[j] *= psize;
    outa[j]  *= psize;
  }
  //send and recieve data
  MPI_Allgatherv (senddata, outa[pid()], MPI_DOUBLE,
		  recdata, outa, outat, MPI_DOUBLE,
		  MPI_COMM_WORLD); 
  //count new particles
  for (int j = 0; j < outt ; j++) { 
    coord a;
    foreach_dimension()
      a.x = recdata[j*psize + cind.x];
    if (locate (a.x, a.y, a.z).level > 0)
      in++;
  }
  int n_partn = n_part + in - out;
  //Manage the memory if required...
  if (n_partn > n_part_a || 2*(n_partn + 1) < n_part_a ) {
    n_part_a = 2*(n_partn + 1);
    loc = realloc   (loc  , n_part_a*sizeof(coord));
    loc_n = realloc (loc_n, n_part_a*sizeof(coord));
    k1 = realloc    (k1, n_part_a*sizeof(double)*dimension);
    k2 = realloc    (k2, n_part_a*sizeof(double)*dimension);
  }
  //Collect new particles from `recdata`
  if (in > 0) {
    int indi[in];
    m = 0;
    for (int j = 0; j < outt; j++) {
      coord a;
      foreach_dimension()
	a.x = recdata[j*psize + cind.x];
      if (locate (a.x, a.y, a.z).level > 0) {
	indi[m++] = j;
      }
    }
    m = 0;
    for (j = n_part - out; j < n_partn; j++) {
      int c = 0;
      foreach_dimension()
	loc[j].x = recdata[indi[m]*psize + c++] ;
      if (step > 1) {
	foreach_dimension() {
	  loc_n[j].x = recdata [indi[m]*psize + c++];
	  k1[j*dimension + cind.x] = recdata [indi[m]*psize + c++]; 
	}
      }
      if (step > 2) {
	foreach_dimension() {
	  k2[j*dimension + cind.x] = recdata [indi[m]*psize + c++];
	}
      }
      m++;
    }
  }
  //Update `n_part`
  n_part = n_partn;
}
#endif
/**
## Tests 
* [A quality test for the advection schemes](parttest.c)
* [Quantative tests for the advection schemes](tp.c)
* [Test for the MPI particle exchange](test_manyparticles.c)
* [Particles, vof and tracer field in a vortex comparison](reversed.c)

## Usage   
* [All pages using `particles.h`](http://basilisk.fr/_search?patterns=particles.h)
* [Tag a portion of a fluid](splash.c)
* [Settling of volcanic ash](ash.c)
* [Flow in a thermosyphon](tube.c)
* [Planetary core convection](core.c)
* [LES of isotropic turbulence](isotropicLES.c)
* [Laminar mixing in 2D](laminarmixing.c)
* [Axisymmetric mixing](coffee.c)

## Todo 
  
* Tag and trace particles with a constant unique number  
* Sort particles along grid iterator curve  
* Implement a proper Basilisk `particles` coordinates field  
 */
