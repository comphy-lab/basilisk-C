/**
tHree
*/

#include "navier-stokes/centered.h"
#include "view.h"
/**
I Will define here if I want a fancier/more complex way of computing angular velocity. */

double NP = 2.; // number of agents I want to work with...
double HEIGHT =128.0;
double THETA = M_PI/2.*1;
double RADIUS = 1.0;
double B = 10.0; // reorientation owing to gravity
double STERIC_GAP = 5.0; // interaction GAP, before an action is taken,
												 // to avoid steric interaction.
												 //
												 //
double BINARY_DISTANCE = 10.;
double BINARY_ANGLE    = 0.;
double BUOYANCY = 0.;

#include "fpicella/src/periodic-shift-treatment.h"

#define ADD_PART_MEM  coord u; coord u2; coord locn; long unsigned int tag; \
double theta; double thetan; double omega; double omega2; double r; \
coord B; double Thrust;

/** A pusher? A puller?
		draw your microswimmer! .*/
#define alpha +2.0
#define beta  +2.0

/**
I Will define here if I want a fancier/more complex way of computing angular velocity. */
//#define OMEGA_EQUATION() cos(p().theta)/(2.*B) + interpolate_linear(locate (x, y, z), omega, x, y, z)/2.  // FP 20250308 
#define OMEGA_EQUATION() interpolate_linear(locate (x, y, z), omega, x, y, z)/2.  // FP 20250308 
//#define OMEGA_EQUATION() 0.
#include "fpicella/src/tracer-particles-FP.h"
#include "fpicella/src/driver-tHree-smooth.h"

#include "Antoonvh/scatter2.h"

/**
	Vorticity field. */
scalar omega[];

/**
	Variable associated to particles that are modelling microswimmers */
Particles microswimmers;

FILE * output_file; // global output file pointer

scalar un[]; // reference velocity field, for stopping the simulation.

int main()
{
	
  /**
  Space and time are dimensionless. This is necessary to be able to
  use the 'mu = fm' trick. */

  size (64. [0]);
  //DT = HUGE [0];
  DT = 1.0 [0];
  
  origin (-L0/2., -L0/2.);

  periodic (right);
  //periodic (top);
  
  stokes = true;
  //TOLERANCE = 1e-3;

  output_file = fopen ("Output.dat", "w");

  
//  for (N = 32; N <= 256; N *= 2)
//		for(RADIUS = 0.05; RADIUS <= 0.45; RADIUS += 0.05)
	N=64;
	BINARY_DISTANCE = L0/2.;
		for(BINARY_DISTANCE = L0/2.*0.99; BINARY_DISTANCE>=RADIUS; BINARY_DISTANCE*=0.9){
			for(BINARY_ANGLE = -M_PI/2.; BINARY_ANGLE<=+M_PI/2.; BINARY_ANGLE+=M_PI*0.1){
				for(BUOYANCY = -1.0; BUOYANCY <= -1.0 ; BUOYANCY+=1.0){
    			run();
				}
			}
		}
	fclose(output_file);
}

#define WIDTH 0.5
#define EPS 1e-14

event init (t = 0) {
/**
  Boundary condition on top-bottom boundaries.*/
	u.n[top] = dirichlet(0.);
	u.t[top] = dirichlet(0.);
	u.n[bottom] = dirichlet(0.);
	u.t[bottom] = dirichlet(0.);

/**
	Initialize microswimmer particles. 
	When using init_tp_square, the total number
	of particles is NP*NP.*/
	//microswimmers = init_tp_square(NP);
	microswimmers = init_tp_circle(NP);
	foreach_particle_in(microswimmers){
/** This is where I set a non-zero buoyancy force.*/
		p().B.y=BUOYANCY;
		p().r = RADIUS;
		p().Thrust = 4.;
		p().theta = THETA+noise()*M_PI*0.; 
		/**
		First particle is always centered wrt the domain.*/
		if(_j_particle==0){
			p().x = 0.;
			p().y = 0.;
		}
		if(_j_particle==1){
			p().x = BINARY_DISTANCE*cos(BINARY_ANGLE);
			p().y = BINARY_DISTANCE*sin(BINARY_ANGLE);
		}
			
	}

  /**
  Viscosity is  unity. */
  
  mu = fm;

  ///**
  //The channel geometry is defined using Constructive Solid Geometry. */  
  mask (y > +HEIGHT/2. ? top : none);
  mask (y < -HEIGHT/2. ? bottom : none);

  //mask (x > +HEIGHT/2. ? right : none);
  //mask (x < -HEIGHT/2. ? left : none);
	
	u.n[top]    = dirichlet(0.);
	u.t[top]    = dirichlet(0.);
	u.n[bottom] = dirichlet(0.);
	u.t[bottom] = dirichlet(0.);
	u.t[right]  = dirichlet(0.); 
	u.n[right]  = dirichlet(0.);
	u.n[left]   = dirichlet(0.);
	u.t[left]   = dirichlet(0.);
	
	//* Initialize reference velocity field.*/
	  foreach()
    un[] = u.y[];
}


event acceleration (i++){
	compute_microswimmer_forcing_smooth(microswimmers);
	a = av; // assign acceleration term. That's capital!
}

event properties(i++){
	vorticity(u,omega);
	compute_bodyPlot(microswimmers);
}

/** Microswimmer's kinematics, is turned off.*/

/** When computation has conveged, write the output file.*/
event logfile ( i += 1; i <= 1000) {
  double du = change (u.y, un);
  if (i > 5 && du < 1e-3){ // since I'm looking for a steady solution, this must be quite small...
		foreach_particle_in(microswimmers){
  		fprintf(output_file,"%+3.2e %+04i %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+04i\n",
                L0, N, RADIUS, alpha, beta, p().Thrust, BUOYANCY, B, p().x,p().y, p().u.x,p().u.y,OMEGA_EQUATION(),_j_particle);
  		fflush(output_file);
		}
    return 1; /* stop */
  }
}

/**
~~~pythonplot streamwise velocities as a function of position. Black line is the isovalue for zero swimming velocity.
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import scipy.interpolate as interpolate

from matplotlib.colors import TwoSlopeNorm

cmap = plt.get_cmap('viridis')

data = np.loadtxt("Output.dat")

fig, ax = plt.subplots(1,3,figsize=(8,3))

satellite = data[(data[:,13]== 1) & (data[:,6]==-1),:]; # extract only data of the second (satellite) particle


# Extract x, y, z for interpolation
x = satellite[:, 8]
y = satellite[:, 9]

### ZOOM FACTOR
Z = 1.0;

# Create grid
xi = np.linspace(x.min()*Z, x.max()*Z, 50)
yi = np.linspace(y.min()*Z, y.max()*Z, 50)
xi, yi = np.meshgrid(xi, yi)

### PLOT FIRST ONE, on the left
z = satellite[:, 11]  # Streamwise velocity or scalar field

## Normalize z by value at maximum x
#max_x_index = x.argmax()
#z_max_x = z[max_x_index]
#if z_max_x != 0:
#    z = z / z_max_x
#else:
#    print("Warning: Zero velocity at max x.")

# Interpolate z-values onto the grid
zi = interpolate.griddata((x, y), z, (xi, yi), method='linear')

# Plot filled contours
contourf = ax[0].contourf(xi, yi, zi, levels=20, cmap=cmap)
cbar = plt.colorbar(contourf, ax=ax[0])

# Add isoline for zero
contours = ax[0].contour(xi, yi, zi, levels=[0.0], colors='white', linewidths=1)
ax[0].clabel(contours, inline=True, fontsize=8)

# Optionally overlay original points
#ax.scatter(x, y, color='k', s=5, alpha=0.1)

ax[0].set_xlabel(r'$x$')
ax[0].set_ylabel(r'$y$')
ax[0].set_title(r'$u_{streamwise}$')

#ax.legend()
ax[0].set_aspect('equal')  # Make axes equal


### PLOT SECOND ONE, in the center
z = satellite[:, 10]  # Spanwise velocity or scalar field

# Interpolate z-values onto the grid
zi = interpolate.griddata((x, y), z, (xi, yi), method='cubic')

# Plot filled contours
contourf = ax[1].contourf(xi, yi, zi, levels=20, cmap=cmap)
cbar = plt.colorbar(contourf, ax=ax[1])

# Add isoline for zero
contours = ax[1].contour(xi, yi, zi, levels=[-0.025,0.025], colors='white', linewidths=1)
ax[1].clabel(contours, inline=True, fontsize=8)

# Optionally overlay original points
#ax.scatter(x, y, color='k', s=5, alpha=0.1)

ax[1].set_xlabel(r'$x$')
ax[1].set_ylabel(r'$y$')
ax[1].set_title(r'$u_{spanwise}$')
ax[1].set_aspect('equal')  # Make axes equal

### PLOT SECOND ONE, in the center
z = satellite[:, 12]  # Spanwise velocity or scalar field

# Interpolate z-values onto the grid
zi = interpolate.griddata((x, y), z, (xi, yi), method='cubic')

# Plot filled contours
contourf = ax[2].contourf(xi, yi, zi, levels=20, cmap=cmap)
cbar = plt.colorbar(contourf, ax=ax[2])

# Add isoline for zero
contours = ax[2].contour(xi, yi, zi, levels=[-0.025,0.025], colors='white', linewidths=1)
ax[2].clabel(contours, inline=True, fontsize=8)

# Optionally overlay original points
#ax.scatter(x, y, color='k', s=5, alpha=0.1)

ax[2].set_xlabel(r'$x$')
ax[2].set_ylabel(r'$y$')
ax[2].set_title(r'$\omega$')
ax[2].set_aspect('equal')  # Make axes equal



plt.tight_layout();
plt.savefig('matplotlib_output.svg')
~~~
*/
