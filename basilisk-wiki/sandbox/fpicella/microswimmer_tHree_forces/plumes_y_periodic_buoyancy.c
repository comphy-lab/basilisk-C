/**
# Plumes without a surface accumulation.
Swimming velocity of a multiple, non force-free (buoyant) particle.
Placed between two plates.
(like an infinite vertical pipe...)
![32 heart-shaped microswimmers in a y periodic domain + BUOYANCY (non-force-free!)](plumes_y_periodic_buoyancy/plumes_y_periodic_buoyancy.mp4)(loop)
*/

#include "navier-stokes/centered.h"
#include "view.h"
/**
I Will define here if I want a fancier/more complex way of computing angular velocity. */

double NP = 32.; // number of agents I want to work with...
double HEIGHT =64.0;
double THETA = M_PI/2.*1;
double RADIUS = 1.0;
double THRUST = 4.;
double BUOYANCY = -4*0.025*1.;
double B = 1.0; // reorientation owing to gravity
double STERIC_GAP = 4.0; // interaction GAP, before an action is taken,
												 // to avoid steric interaction.

#include "fpicella/src/periodic-shift-treatment.h"

#define ADD_PART_MEM  coord u; coord u2; coord locn; long unsigned int tag; \
double theta; double thetan; double omega; double omega2; double r; \
coord B; double Thrust;

/** A pusher? A puller?
		draw your microswimmer! .*/
#define alpha +2
#define beta  +2

//#define Particle_Advection_ON 1
/** In this implementation, I will use a custom (and simpler...)
		timestepper for advecting the particle's location.
		In particular, simple Explicit Euler is used.
		It is way simpler than Antoon's Runge Kutta, but it does the job!*/
/** Most important, it couples well with the anti
		penetration alghoritm that I've employed.
		Here I just work on particle's location, no on velocity. Why?
		My suggestion (it's a model!) is that we're in Stokes regime.
		You can not _bounce_ when quasi-steady ;) */

/**
I Will define here if I want a fancier/more complex way of computing angular velocity. */
#define OMEGA_EQUATION() cos(p().theta)/(2.*B) + interpolate_linear(locate (x, y, z), omega, x, y, z)/2.  // FP 20250308 
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

  //periodic (right);
  periodic (top);
  
  stokes = true;
  //TOLERANCE = 1e-3;

  //output_file = fopen ("Output.dat", "w");

  
//  for (N = 32; N <= 256; N *= 2)
//		for(RADIUS = 0.05; RADIUS <= 0.45; RADIUS += 0.05)
	N=64;
    	run();

	//fclose(output_file);
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
		p().Thrust = THRUST;
		p().theta = THETA;//+noise()*M_PI*2.; 
/** Random initialization .*/
		p().x =     L0/2.*noise();
		p().y = HEIGHT/2.*noise();
		if(_j_particle==0){
			p().x = 0.;
			p().y = 0.;
		}
	}

  /**
  Viscosity is  unity. */
  
  mu = fm;

  ///**
  //The channel geometry is defined using Constructive Solid Geometry. */  
  //mask (y > +HEIGHT/2. ? top : none);
  //mask (y < -HEIGHT/2. ? bottom : none);

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


event log_make_video (i += 1) {
  squares ("omega", linear = true);
  scatter (microswimmers, pc = {256, 256, 256});
  isoline ("bodyPlot", 0.25, lc = {256,256,256}, lw = 2);
  vectors (u = "u", scale = 0.25, lc = {256,256,256}, level = 8);
  box();
  save ("plumes_y_periodic_buoyancy.mp4");
}

/** End event, to properly unallocate variables.*/
event end (t = 10000);

event acceleration (i++){
	compute_microswimmer_forcing_smooth(microswimmers);
	//compute_repulsion(microswimmers,4.0,false,false,false,false,6.); // look at the driver-tHree file
																																		// to know what argument is meant
																																		// to do ;)
	a = av; // assign acceleration term. That's capital!
}

event properties(i++){
	vorticity(u,omega);
	compute_bodyPlot(microswimmers);
}

/** Microswimmer's kinematics, using simple Forward Euler.*/
event end_timestep(i++){
	foreach_particle(){
	/** Standard behaviour.*/
		foreach_dimension()
			p().x += p().u.x*dt; // The most simple timestepper possible:
													 // Explicit Euler.
	/** Angular kinematics.*/
		p().omega = OMEGA_EQUATION();
		p().theta += p().omega*dt;

//		/** Detect if particle is in contact with top or bottom boundary.*/
//		/** Hence, that's simply a matter of y-aligned dynamics.*/
//		if(p().y+p().r*STERIC_GAP*2.>HEIGHT/2.) // impacting with top boundary
//			p().y = HEIGHT/2.-p().r*STERIC_GAP*2.;
//		if(p().y-p().r*STERIC_GAP<-HEIGHT/2.) // impacting with bottom boundary
//			p().y = -HEIGHT/2.+p().r*STERIC_GAP;
//
		/** Detect if particle is in contact with top or bottom boundary.*/
		/** Hence, that's simply a matter of y-aligned dynamics.*/
		/** IN THE X DIRECTION NOW! */
		if(p().x+p().r*STERIC_GAP*2.>HEIGHT/2.) // impacting with top boundary
			p().x = HEIGHT/2.-p().r*STERIC_GAP*2.;
		if(p().x-p().r*STERIC_GAP<-HEIGHT/2.) // impacting with bottom boundary
			p().x = -HEIGHT/2.+p().r*STERIC_GAP;



  	/** Detect particle-particle contact (GAP<0) and correct kinematics..*/
  	for (int _k_particle = 0; _k_particle < pn[_l_particle]; _k_particle++) {
  	  if(_k_particle != _j_particle){
  	    coord DISTANCE; // relative distance
  	    foreach_dimension()
  	    DISTANCE.x = pl[_l_particle][_k_particle].x - pl[_l_particle][_j_particle].x;
  	    double distance = sqrt(sq(DISTANCE.x)+sq(DISTANCE.y));
        double angle = atan2(DISTANCE.y/distance,DISTANCE.x/distance);
				double GAP = distance - STERIC_GAP*p().r;
				/** What if the GAP < 0?
				p() are getting too close, action required!*/
			  if(GAP<0.){
					/** Correct position to avoid penetration.*/
					/** Correction is simply proportional to the amount
							of penetration (GAP<0) that is observed.*/
					p().x += GAP*cos(angle);
					p().y += GAP*sin(angle);
					/** Since penetration is symmetric between the two
							particles, the new outcome GAP will be symmetrical
							distributed to the two particles as well.*/
			  }
  	  }
  	}
	}
}


/**
~~~pythonplot 

"""
Analyze a sequence of particle_*.dat files:
 - Compute ⟨vy⟩ (average vertical velocity) at each timestep
 - Compute ⟨ρ⟩ (average local KDE-based density) at each timestep
 - Plot and save both as a function of time
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import glob
from tqdm import tqdm

# ────────────── Settings ──────────────
#folder        = "plumes_y_periodic_buoyancy"  # Update if needed
#file_pattern  = os.path.join(folder, "particle_*.dat")
#output_plot   = os.path.join(folder, "summary_vy_density.svg")
file_pattern = "particle_*.dat"
output_plot  = "summary_vy_density.svg"

# Domain extents
x_min, x_max = -64, 64
y_min, y_max = -64, 64

# Grid resolution
x_bins = 256
y_bins = 256
dx = (x_max - x_min) / x_bins
dy = (y_max - y_min) / y_bins

# KDE Gaussian kernel
kernel_radius_physical = 1.0
sigma_physical = kernel_radius_physical / 0.3
### IMPORTANT PARAMETER TO WORK OUT THE KERNEL RADIUS!!!!
ker_rad_x = 10*int(np.ceil(kernel_radius_physical / dx))
ker_rad_y = 10*int(np.ceil(kernel_radius_physical / dy))
ker_size_x = 2 * ker_rad_x + 1
ker_size_y = 2 * ker_rad_y + 1
xk = np.linspace(-ker_rad_x * dx, ker_rad_x * dx, ker_size_x)
yk = np.linspace(-ker_rad_y * dy, ker_rad_y * dy, ker_size_y)
Xk, Yk = np.meshgrid(xk, yk)
gaussian_kernel = np.exp(-(Xk**2 + Yk**2) / (2 * sigma_physical**2)).astype(np.float32)

# ────────────── Load files ──────────────
particle_files = sorted(glob.glob(file_pattern))
if not particle_files:
    raise FileNotFoundError(f"No files matching {file_pattern}")

# Load all data
all_particle_data = [np.loadtxt(f) for f in particle_files]
num_timesteps = all_particle_data[0].shape[0]

# ────────────── Time series containers ──────────────
avg_vy_per_timestep = []
avg_kde_per_timestep = []
timesteps = []

# ────────────── Main loop ──────────────
print("Computing time series of ⟨vy⟩ and ⟨ρ⟩ (KDE)...")
for t_idx in tqdm(range(num_timesteps)):
    presence_map = np.zeros((y_bins, x_bins), dtype=np.float32)
    x_list, y_list, vy_list = [], [], []

    for pdata in all_particle_data:
        x, y  = pdata[t_idx, 1], pdata[t_idx, 2]
        vy    = pdata[t_idx, 5]

        if not (x_min <= x <= x_max and y_min <= y <= y_max):
            continue

        x_list.append(x)
        y_list.append(y)
        vy_list.append(vy)

        xi = int((x - x_min) / dx)
        yi = int((y - y_min) / dy)

        x0, y0 = xi - ker_rad_x, yi - ker_rad_y
        x1, y1 = x0 + ker_size_x, y0 + ker_size_y

        x_start = max(x0, 0)
        y_start = max(y0, 0)
        x_end   = min(x1, x_bins)
        y_end   = min(y1, y_bins)

        kx_start = x_start - x0
        ky_start = y_start - y0
        kx_end   = kx_start + (x_end - x_start)
        ky_end   = ky_start + (y_end - y_start)

        presence_map[y_start:y_end, x_start:x_end] += \
            gaussian_kernel[ky_start:ky_end, kx_start:kx_end]

    kde_values = []
    for x, y in zip(x_list, y_list):
        xi = int((x - x_min) / dx)
        yi = int((y - y_min) / dy)
        if 0 <= xi < x_bins and 0 <= yi < y_bins:
            kde_values.append(presence_map[yi, xi])

    kde_values = np.array(kde_values)
    vy_list = np.array(vy_list)

    avg_vy = np.mean(vy_list) if vy_list.size > 0 else 0.0
    avg_kde = np.mean(kde_values) if kde_values.size > 0 else 0.0

    avg_vy_per_timestep.append(avg_vy)
    avg_kde_per_timestep.append(avg_kde)
    timesteps.append(t_idx)

# ────────────── Parameters ──────────────
sliding_window = 500  # <-- Adjust this to change the smoothing window

# ────────────── Sliding average function ──────────────
def sliding_average(data, window):
    if len(data) < window:
        return np.array(data)
    return np.convolve(data, np.ones(window)/window, mode='valid')

# Compute sliding averages
avg_vy_smoothed  = sliding_average(avg_vy_per_timestep, sliding_window)
avg_kde_smoothed = sliding_average(avg_kde_per_timestep, sliding_window)
smoothed_timesteps = timesteps[sliding_window - 1:]  # Align with smoothed data

# ────────────── Plot time series + smoothed ──────────────
plt.rc('text',usetex=True);
plt.rc('font',family='serif');

fig, ax = plt.subplots(2, 1, figsize=(5, 4), sharex=True)

# Plot ⟨vy⟩
#ax[0].plot(timesteps, avg_vy_per_timestep, label="⟨vy⟩ raw", color='crimson', alpha=0.4)
#ax[0].plot(smoothed_timesteps, avg_vy_smoothed, label=f"⟨vy⟩ (avg {sliding_window})", color='crimson', linewidth=2)
ax[0].plot(timesteps, avg_vy_per_timestep, color='crimson', alpha=0.25)
ax[0].plot(smoothed_timesteps, avg_vy_smoothed, color='crimson', linewidth=2)
ax[0].set_ylabel(r"$\left<U_y\right>, \overline{\left< U_y \right>}$")
#ax[0].legend()
ax[0].grid(True)

# Plot ⟨ρ⟩
#ax[1].plot(timesteps, avg_kde_per_timestep, label="⟨ρ⟩ raw", color='teal', alpha=0.4)
#ax[1].plot(smoothed_timesteps, avg_kde_smoothed, label=f"⟨ρ⟩ (avg {sliding_window})", color='teal', linewidth=2)
ax[1].plot(timesteps, avg_kde_per_timestep, color='teal', alpha=0.25)
ax[1].plot(smoothed_timesteps, avg_kde_smoothed, color='teal', linewidth=2)
ax[1].set_ylabel(r"$\left<\rho\right>, \overline{\left<\rho\right>} $")
ax[1].set_xlabel("Time")
#ax[1].legend()
ax[1].grid(True)

plt.tight_layout()
plt.savefig(output_plot)
print(f"Saved smoothed time series plot to: {output_plot}")
~~~
*/
