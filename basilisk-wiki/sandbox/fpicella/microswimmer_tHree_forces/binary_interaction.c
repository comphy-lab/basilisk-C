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
		for(BINARY_DISTANCE = L0/2.; BINARY_DISTANCE>=RADIUS; BINARY_DISTANCE*=0.9){
			for(BUOYANCY = -1; BUOYANCY <= +1. ; BUOYANCY+=1.0){
    		run();
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
/** Random initialization .*/
		p().x = -BINARY_DISTANCE/2. +_j_particle*p().r*BINARY_DISTANCE;
		p().y = 0.;
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
  if (i > 5 && du < 1e-5){ // since I'm looking for a steady solution, this must be quite small...
		foreach_particle_in(microswimmers){
  		fprintf(output_file,"%+3.2e %+04i %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+04i\n",
                L0, N, RADIUS, alpha, beta, p().Thrust, BUOYANCY, B, BINARY_DISTANCE, p().u.x,p().u.y,OMEGA_EQUATION(),_j_particle);
  		fflush(output_file);
		}
    return 1; /* stop */
  }
}

/**
# Binary interactions between tHree microswimmers
~~~gnuplot Binary interactions between two tHree microswimmers
set grid
set xlabel "d/R"
set title "binary interactions, spanwise aligned, theta = pi/2. (poiting y direction)"

# Plot data from file
plot 'Output.dat' using ( $13==0 && $7==0 ? ($9/$3) : 1/0 ):( $13==0  && $7==0 ? +$11 : 1/0 ) \
     with points pt 7 ps 2 lc rgb "green" title "u_y",\
     'Output.dat' using ( $13==0 && $7==0 ? ($9/$3) : 1/0 ):( $13==0  && $7==0 ? +$10 : 1/0 ) \
     with points pt 7 ps 2 lc rgb "blue" title "u_x ",\
     'Output.dat' using ( $13==0 && $7==0 ? ($9/$3) : 1/0 ):( $13==0  && $7==0 ? +$12 : 1/0 ) \
     with points pt 7 ps 2 lc rgb "red" title "omega"
~~~
*/


/**
~~~pythonplot velocities as a function of binary interaction, no buoyancy
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import scipy.interpolate as interpolate

cmap = plt.get_cmap('plasma')

data = np.loadtxt("Output.dat")

fig, ax = plt.subplots(figsize=(4,3))

ax.plot(data[(data[:, 12] == 0) & (data[:, 6] == 0), 8] / data[(data[:, 12] == 0) & (data[:, 6] == 0), 2], data[(data[:, 12] == 0) & (data[:, 6] == 0), 9], 'o', color=cmap(0.3), markersize=8, label="$u_x$")
ax.plot(data[(data[:, 12] == 0) & (data[:, 6] == 0), 8] / data[(data[:, 12] == 0) & (data[:, 6] == 0), 2], data[(data[:, 12] == 0) & (data[:, 6] == 0), 10], 'o', color=cmap(0.6), markersize=8, label="$u_y$")
ax.plot(data[(data[:, 12] == 0) & (data[:, 6] == 0), 8] / data[(data[:, 12] == 0) & (data[:, 6] == 0), 2], data[(data[:, 12] == 0) & (data[:, 6] == 0), 11], 'o', color=cmap(0.9), markersize=8, label="$\omega$")

plt.xlabel('$d/R$');
ax.legend()
plt.tight_layout();
plt.savefig('matplotlib_output.svg')
~~~
*/

/**
~~~pythonplot Swimming velocity, influenced by body distance and buoyancy.
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import scipy.interpolate as interpolate

cmap = plt.get_cmap('plasma')

data = np.loadtxt("Output.dat")

fig, ax = plt.subplots(figsize=(4,3))

ax.plot(data[(data[:, 12] == 0) & (data[:, 6] == -1), 8] / data[(data[:, 12] == 0) & (data[:, 6] == -1), 2], data[(data[:, 12] == 0) & (data[:, 6] == -1), 10], 'o', color=cmap(0.0), markersize=8, label="$-0.25$")
ax.plot(data[(data[:, 12] == 0) & (data[:, 6] == +0), 8] / data[(data[:, 12] == 0) & (data[:, 6] == +0), 2], data[(data[:, 12] == 0) & (data[:, 6] == +0), 10], 'o', color=cmap(0.3), markersize=8, label="$+0.0$")
ax.plot(data[(data[:, 12] == 0) & (data[:, 6] == +1), 8] / data[(data[:, 12] == 0) & (data[:, 6] == +1), 2], data[(data[:, 12] == 0) & (data[:, 6] == +1), 10], 'o', color=cmap(0.6), markersize=8, label="$+0.25$")

plt.xlabel('$d/R$');
plt.ylabel('$U_{swim}$');

# Add a dashed line at y=0
plt.axhline(y=0, color='black', linestyle=':', linewidth=1)

ax.legend(title = "buoyancy/thrust")
plt.tight_layout();
plt.savefig('velocity_distance_buoyancy.svg')
~~~
*/

/**
- Positive values, the agent is swimming towards the top (same as direction)
- Negative values, swimmer is _sinking_ (velocity sign is opposite than direction)


### So what?
- When particle come close enough, they swimming velocity is reduced.
- If negative buoyancy is added, they could start sinking when close enough!
*/

/**
~~~pythonplot Normalized swimming velocity, influenced by body distance and buoyancy.
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import scipy.interpolate as interpolate

cmap = plt.get_cmap('plasma')

data = np.loadtxt("Output.dat")

fig, ax = plt.subplots(figsize=(4,3))

# A fancy function to make it simpler later...

def get_normalized_y_at_max_x(data, mask, x_col, y_col, norm_col):
    """
    Returns x values and y values normalized so that y at max x becomes 1.

    Parameters:
        data     : numpy array with the data.
        mask     : boolean array to filter rows.
        x_col    : column index for the numerator of x-values.
        y_col    : column index for y-values.
        norm_col : column index for the denominator of x-values.

    Returns:
        x_vals                : array of x values (not normalized)
        y_vals_normalized     : y values scaled so that y at max x = 1
    """
    x_vals = data[mask, x_col] / data[mask, norm_col]
    y_vals = data[mask, y_col]

    if len(x_vals) == 0:
        return x_vals, y_vals  # empty

    max_x_index = x_vals.argmax()
    y_at_max_x = y_vals[max_x_index]

    if y_at_max_x == 0:
        y_vals_normalized = y_vals  # avoid division by zero
    else:
        y_vals_normalized = y_vals / y_at_max_x

    return x_vals, y_vals_normalized

""" negative buoyancy """
mask = (data[:, 12] == 0) & (data[:, 6] == -1)
x_vals, y_vals_norm = get_normalized_y_at_max_x(data, mask, x_col=8, y_col=10, norm_col=2)
ax.plot(x_vals, y_vals_norm, 'o', color=cmap(0.0), label="$-0.25$")
""" neutral buoyancy """
mask = (data[:, 12] == 0) & (data[:, 6] == +0)
x_vals, y_vals_norm = get_normalized_y_at_max_x(data, mask, x_col=8, y_col=10, norm_col=2)
ax.plot(x_vals, y_vals_norm, 'o', color=cmap(0.3), label="$+0.00$")
""" positive buoyancy """
mask = (data[:, 12] == 0) & (data[:, 6] == +1)
x_vals, y_vals_norm = get_normalized_y_at_max_x(data, mask, x_col=8, y_col=10, norm_col=2)
ax.plot(x_vals, y_vals_norm, 'o', color=cmap(0.6), label="$+0.25$")


plt.xlabel('$d/R$');
plt.ylabel('$U_{swim}/U_{swim}(d/R=\infty)$');

plt.axhline(y=0, color='black', linestyle=':', linewidth=1)

ax.legend(title = "buoyancy/thrust", framealpha =1 )
plt.tight_layout();
plt.savefig('normalized_velocity_distance_buoyancy.svg')
~~~
*/
/**
### Reudction in swimming velocity  _adding buoyancy to the tHree model_:
 - I normalize the swimming velocity wrt the velocity that would be observed in an unbounded configuration ($U_{swim}(d/R = \infty)$)
 - I reduction (up to a change of sign!!!) of swimming velocity is enhanced when using negative buoyancy.
 - This could explain the _Rayleigh-Taylor_-like instability when *local* concentration is high enough
 - whilst the swimming velocity of an isolated tHree microswimmer is positive, even with negative buoyancy.
*/
