/**
Single air bubble, rising in a quiescent fluid, inmpacting a free-interface.

# Can a bubble change local particulate density in a solution of passive particles ?

Some technical points
 - two-vof method, to avoid coalescence
 - non inertial (black) and inertial (red) particles. 

Here I've got a bubble of diameter $1.0mm$, rising in a pool with a free surface $50 mm$ above.

As a rule of thumb, characteristic time here is $T \approx \sqrt{\dfrac{dBubble}{g}}$.

Total simulation time represents about $0.5s$

![An air bubble (blue interface) rising in a water pool containing tracer particles. Particle concentration around the bubble evolves in time, especially when the bubble reaches a free-interface (green). Video in the bubble's moving reference frame.](axi_particles_no_coalescence_quantitative/movie_follow-dBubble_+1.000e+00+5.000e+01.mp4)(width="75%")

# Temporal evolution of tracer particles near a rising bubble.
~~~gnuplot Bubble position and near-bubble particle concentration in time
set grid
set xlabel "t"
set key top left


file = 'particles_near_cg-dBubble_+1.000e+00-interface_location_+5.000e+01-radius_+5.680e-01.dat'

plot file u 1:2 w l title 'yBubble', \
     file u 1:5 w l title '# particles'

~~~

### What I see:

 - After a transient, the bubble rises at almost constant speed, up until reaching the free surface ($y = 50, t\approx 30$)
 - The number of particles _close_ to the interface rapidly raises in a transient, reaching a plateau when constant bubble speed is attained.
 - A steady particle concentration pattern is obtained around the bubble.
 - When the bubble gets closer up to an alt next to the free-surface, flow does not stop immediately, because we're already in a mildly inertial regime. 
 - Passive particles then continues to flow, generating still another concentration parttern around the bubble.
 - This provides also a new transient in particle concentration around the bubble. 
 
 
### The code
*/



#include "grid/quadtree.h"
#include "axi.h"
#include "navier-stokes/centered.h"
//// This block acts like a _two-phase.h_, but without advecting f field.
#include "vof.h"
scalar f[], f1[], f2[], * interfaces = {f1,f2}; // f1 will be for the bubble, f2 for the outer gas
// please note, I have deliberately avoided having f as a _interface_ field
#include "two-phase-generic.h"
//#include "navier-stokes/conserving.h" // it is not compatible with double vof.
//// end of block acting as _two-phase.h_
#include "tension.h"
#include "view.h"

#include "Antoonvh/tracer-particles.h"
#include "Antoonvh/stokes-particles.h"
#include "Antoonvh/scatter2.h"
Particles particles00, particles01, particles02;

//////////////////////////////////////////////////////////////////
/**
Domain size, final simulation time, timestepping for post-processing purposes.*/
// BUBBLE AVERAGE DIAMETER IS SET ALWAYS TO 1!
double box_size = 64.;
double endTime = 60.;
double saveStepTime= 0.25;
double interface_location; // is defined within the main()
double domain_shift = -1.;
// Bubble is seeded at location 0.
double phi_threshold = - 0.025; // to initialize seeded particle location far enough from interface.
double particle_density = 16.; // how many particles per unit line.
// For post-processing, how close from bubble center of gravity to count for particles
double distance_from_cg = 0.5*(1.136); // threshold radius, dimensionless length



/**
The only physical parameter I can change in my problem. Average bubble size (i.e. bubble volume...)
injected within the system.*/
double dBubble; // global variable, what is the present bubble diameter for which I'm doing my test?
/**
Phisical parameters of the problem. Air and water.*/
double rhoOut = 1000.0;
double nuOut  = 1.0e-6;
double sigma  = 0.072;
double gCons  = 9.81;

int MAXLEVEL = 9;
int MINLEVEL = 4;

long unsigned int Nparticles;

// // // // // // PSEUDO-LUBRICATION // // //
#include "heights.h"
double hcut  = 1.0;    // force cutoff  // range over which repulsion acts
double href  = 0.1;   // reference gap // raises the whole curve
double Klub  = 1.0;    // strength      // global multiplier
int    plub  = 8;      // exponent      // makes the force rise faster as h gets small

static inline bool vof_interfacial_cell (double c)
{
  return (c > 1e-6 && c < 1. - 1e-6);
}

static inline double Pi_lub (double h)
{
  double hmin = 0.5*L0/(1 << MAXLEVEL); // regularisation ~ half finest cell
  if (h <= 0. || h >= hcut)
    return 0.;

  h = max(h, hmin);

  // smooth cutoff to avoid a force discontinuity at h = hcut
  double w = 1. - h/hcut;
  return Klub*pow(href/h, plub)*sq(w);
}
// // // // // // // // // // // // // // //
// // // INITIAL BUBBLE SHAPE // // // // //
/**
# First-order analytical bubble shape (Bond number controlled)

We approximate the static bubble shape as a first-order perturbation of a
sphere in the Bond number

$$
Bo = \frac{\Delta \rho g R^2}{\gamma}.
$$

- `Bo = 0` → exact sphere  
- valid for moderate Bond numbers  
- does not capture necking or detachment
*/

static inline double P2 (double mu)
{
  return 0.5*(3.*sq(mu) - 1.);
}

static inline double bubble_phi_axi_bary_unit (double x, double y,
                                               double xBary,
                                               double Bo)
{
  const double eps = Bo/6.;
  const double vf = 1. + (3./5.)*sq(eps) + (2./35.)*eps*sq(eps);

  // radius corrected so that the bubble volume equals that of a sphere of diameter 1
  const double Rs = 0.5*pow(vf, -1./3.);

  const double X = x - xBary;
  const double R = sqrt(sq(X) + sq(y));

  if (R < 1e-14)
    return (Rs*(1. + eps));

  return (Rs*(1. + eps*P2(X/R)) - R);
}

/*
Function to track bubble location in time.*/

double interface_x_center (scalar f)
{
  double sx = 0., sl = 0.;

  foreach (reduction(+:sx) reduction(+:sl)) {
    if (f[] > 0. && f[] < 1.) {
      coord n = interface_normal (point, f);
      double alpha = plane_alpha (f[], n);

      coord p[2];
      int m = facets (n, alpha, p);

      if (m == 2) {
        foreach_dimension()
          p[0].x = x + p[0].x*Delta,
          p[1].x = x + p[1].x*Delta;

        double l = sqrt(sq(p[1].x - p[0].x) +
                        sq(p[1].y - p[0].y));

        double xm = 0.5*(p[0].x + p[1].x);

        sx += xm*l;
        sl += l;
      }
    }
  }

  return sl > 0. ? sx/sl : nodata;
}

//////////////////////////////////////////////////////////////////
int main()
{
  size (box_size);

  origin(domain_shift,0.); // domain shifted so not to interfere with bubble.

//#if !TREE
  N = 1 << (MAXLEVEL-2);
//#endif

/** The case is made dimensionless using the same approach as the one described in the [bubble example](https://basilisk.fr/src/examples/bubble.c). */

#define RHOR 1000.
#define MUR 100.

/**
For my convenience, I define Ga and Bo as macros... */
#define Ga (sqrt(gCons*(dBubble)*(dBubble)*(dBubble)*1e-9)/(nuOut))
#define Bo (rhoOut*(gCons)*(dBubble)*(dBubble)*1e-6/(sigma))

/**
I the only variable parameter here is the bubble diameter dBubb.
I provide it naively within a list .
As a convention, here dBubble is provided in millimiter.
This is the only _dimensional_ parameter, on which I will work out all the other
dimensionless parameters (Bond and Galillei) that drives the problem.
*/
  //double dBubble_list[] = {0.5,1.0,1.5,2.0,2.5};
  double dBubble_list[] = {1.0};
  int n = sizeof(dBubble_list)/sizeof(dBubble_list[0]);

  for (int k = 0; k < n; k++) { // loop on all bubble sizes
    dBubble = dBubble_list[k];
    fprintf(stderr,"### ### ### \n");
    fprintf(stderr, "iteration %04d/%04d, dBubble[mm] %+06.5e \n",k,n,dBubble);
    fprintf(stderr,"### ### ### \n");
/* Determine Bond and Galillei number for the simulation.*/

    //double Ga = sqrt(gCons*dBubble*dBubble*dBubble*1e-9)/nuOut;
    //double Bo = rhoOut*gCons*dBubble*dBubble*1e-6/sigma;

    fprintf(stderr, "dBubble = %+06.3e mm\n", dBubble);
    fprintf(stderr, "Ga      = %+06.3e\n", Ga);
    fprintf(stderr, "Bo      = %+06.3e\n", Bo);

    rho2 = 1. [0];    // liquid
    rho1 = 1./RHOR; // gas
    mu2 = 1./Ga;      // liquid
    mu1 = 1./(MUR*Ga);// gas
											//
    f1.sigma = 1./Bo;
    f2.sigma = 1./Bo;

// Variable interface location
  	double interface_location_list[] = {50.0};
  	int nn = sizeof(interface_location_list)/sizeof(interface_location_list[0]);
  	for (int kk = 0; kk < nn; kk++) { // loop on all bubble sizes
    	interface_location = interface_location_list[kk];
    	// How many particles to track?
    	// I try to have the same _density_ (here 50) regardless on the 
    	// domain size (i.e. the part of the fluid that is below the free-interface).
//    	Nparticles = (long unsigned int) llround(particle_density*(interface_location-domain_shift-3.14159/8));
//    fprintf(stderr,"I'm tracking %06lu particles \n",Nparticles);
    	run();
		}
  }
}

event init (t = 0)
{
// Remove the part of the domain that is not useful for the computation...
//  mask (y > L0/2 ? top : none); // mask does not work (apparently...) in MPI.
#if TREE
    // refine bubble location
  refine (bubble_phi_axi_bary_unit (x, y, 0., Bo) < 0. && level < MAXLEVEL);
    // refine interface location
  refine (-x+interface_location && level < MAXLEVEL);
#endif

  //fraction (f1, - (sq(x) + sq(y) - sq(0.5))); // spherical bubble
  fraction (f1, bubble_phi_axi_bary_unit(x,y,0.,Bo)); // non-spherical...
  fraction (f2, x-interface_location); // free-interface
	particles00 = new_tracer_particles(0);
// particles distributed homogeneously within the fluid domain
// First, I fill the whole domain (and even further)
// then I cut-out particles that sits in areas that are not of my interest.
  Nparticles = (long unsigned int) llround(particle_density*box_size);
  // NB, within init_ip square, this is a _linear_ density
  fprintf(stderr,"I'm tracking %06lu particles \n",Nparticles);
  particles01 = init_ip_square (Nparticles,
                                xm=box_size/2.+domain_shift,
                                ym=0.,
                                l=box_size);
  particles02 = init_ip_square (Nparticles,
                                xm=box_size/2.+domain_shift,
                                ym=0.,
                                l=box_size);
  // remove particles that are beyond the interface location
  remove_particles(particles01, x>interface_location);
  remove_particles(particles02, x>interface_location);

  // remove particles that are far enough from the axis
	remove_particles(particles01, y>1.0);
	remove_particles(particles02, y>1.0);
	// remove particles that are within the bubble, with a small threshold
	remove_particles(particles01, bubble_phi_axi_bary_unit(x,y,0.,Bo)>-0.01);
	remove_particles(particles02, bubble_phi_axi_bary_unit(x,y,0.,Bo)>-0.01);

	// Add particle tag
  set_particle_attributes(particles01);
  set_particle_attributes(particles02);
  // Consider variable inertia
  foreach_particle_in(particles01)
    p().u2.z = 1e-5;
  foreach_particle_in(particles02)
    p().u2.z = 1e-2;

// How many particles I'm tracking?

	int n01 = pn[particles01];
	int n02 = pn[particles02];

#if _MPI
	MPI_Allreduce(MPI_IN_PLACE, &n01, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &n02, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif

	fprintf(stderr,"I'm tracking n01 %05d and n02 %05d particles \n",n01,n02);

}
/*
This is not the most elegant way to do this,
but I will still employ a _common_ f[] field to account for
f1[] and f2[], together.
I need this so to easily assign variable viscosity to the field.
*/
event properties(i++){
	foreach()
  	f[] = max(f1[], f2[]);
	boundary ({f});
}

/*
Employing naive two-tracer VOF to avoid coalescence, 
reduced.h approach is not valid anymore.
For this reason, I try to work out a _custom_ one, working on fields f1 and f2.*/
coord G = {-1.,0.,0.}, Z = {0.,0.,0.};

event acceleration (i++)
{
  scalar phi1 = f1.phi;
  scalar phi2 = f2.phi;
  coord G1;
  foreach_dimension()
    G1.x = (rho2 - rho1)*G.x;

  // --- gravity contribution, same logic as your current code ---
  if (phi1.i)
    position (f1, phi1, G1, Z, add = true);
  else {
    phi1 = new scalar;
    position (f1, phi1, G1, Z, add = false);
    f1.phi = phi1;
  }

  if (phi2.i)
    position (f2, phi2, G1, Z, add = true);
  else {
    phi2 = new scalar;
    position (f2, phi2, G1, Z, add = false);
    f2.phi = phi2;
  }

  // --- lubrication-like repulsion contribution ---
	/**
	The local gap $h$ between interfaces is estimated along $x$ using height
	functions:
	$$x_b = x + \Delta\,h_1,\quad x_s = x + \Delta\,h_2,\quad h = x_s - x_b.$$
	Here $h$ has the same units as $x$, $L_0$ and $\Delta$ (dimensionless, scaled
	by the bubble diameter), so $h \sim \Delta$ corresponds to a one-cell gap.
	We model a short-range lubrication-like repulsion as an interfacial (surface)
	force by adding a potential $\phi$ to each tracer, following Basilisk’s
	CSF framework:
	$$\mathbf{a}_{\text{rep}} = \alpha\,\phi\,\nabla f.$$
	The potential is defined as
	$$\phi = \Pi(h), \quad
	\Pi(h) = K_{\text{lub}}\left(\frac{h_{\text{ref}}}{h}\right)^p
	\left(1 - \frac{h}{h_{\text{cut}}}\right)^2 \ \text{for } h < h_{\text{cut}},$$
	and $\Pi(h)=0$ otherwise, with $h_{\text{ref}}, h_{\text{cut}} = O(\Delta)$.
	*/
  vector h1[], h2[];
  heights (f1, h1);
  heights (f2, h2);

  foreach() {
    double Pi = 0.;

    // This assumes the thin film is mainly aligned with x,
    // which matches your axisymmetric bubble/free-surface setup.
    if (h1.x[] != nodata && h2.x[] != nodata) {
      double xb = x + Delta*height(h1.x[]);
      double xs = x + Delta*height(h2.x[]);
      double gap = xs - xb;   // free surface should be to the right of the bubble
      Pi = Pi_lub (gap);
    }

    if (Pi > 0.) {
      if (vof_interfacial_cell(f1[]))
        phi1[] += Pi;
      if (vof_interfacial_cell(f2[]))
        phi2[] += Pi;
    }
  }

  boundary ({phi1, phi2});
}

event logfile (i++){
  double xb = interface_x_center(f1);
	fprintf(stderr,"%05d %+06.5e %+06.5e %+06.5e %+06.5e \n", i, t,dBubble,xb,interface_location);
}

#if TREE
event adapt (i++) {
  adapt_wavelet ({f1,f2, u},(double[]) {0.001,0.001,1e-2,1e-2,1e-2},
     maxlevel = MAXLEVEL, minlevel = MINLEVEL);
}
#endif

//event movie (t += saveStepTime; t <= endTime)
//{
//  clear();
//	view(fov=20, tx = 0., ty = -0.4, psi=-pi/2., width = 1200, height = 1200);
//	box();
//  //squares ("u.x", spread = -1, linear = true);
//	cells(lc = {0.75, 0.75, 0.75});
//  draw_vof("f1", lc = {0, 0, 1}, lw = 2);
//  draw_vof("f2", lc = {0, 1, 0}, lw = 2);
//  scatter(particles01, s = 1);
//  // mirror part of the plot, for ease of view...
//  mirror(n={0,1,0}){
//    draw_vof("f1", lc = {0, 0, 1}, lw = 2);
//    draw_vof("f2", lc = {0, 1, 0}, lw = 2);
//    //scatter(particles02, s = 3, pc = {1, 0, 0});
//	  vectors ("u", scale = 0.03, lc = {0,0,0}, lw = 1);
//  };
//	char moviename[256];
//	snprintf(moviename, sizeof(moviename),
//	         "movie-dBubble_%+06.3e%+06.3e.mp4",
//	         dBubble, interface_location);
//	save (moviename);
//}

event movie_follow_bubble (t += saveStepTime; t <= endTime) {
  double xb = interface_x_center(f1);

  view(fov = 2.5, tx = 0., ty = -0.5, psi = -pi/2.,
       width = 1200, height = 1200);
  clear();
  translate (x = -xb+L0/2., y = 0.) {
		box();
		cells(lc = {0.75, 0.75, 0.75});
  	draw_vof("f1", lc = {0, 0, 1}, lw = 2);
  	draw_vof("f2", lc = {0, 1, 0}, lw = 2);
  	scatter(particles01, s = 1.5);
  	// mirror part of the plot, for ease of view...
  	mirror(n={0,1,0}){
  	  draw_vof("f1", lc = {0, 0, 1}, lw = 2);
  	  draw_vof("f2", lc = {0, 1, 0}, lw = 2);
  	  //scatter(particles02, s = 3, pc = {1, 0, 0});
		  vectors ("u", scale = 0.03, lc = {0,0,0}, lw = 1);
  	};
  }

	char moviename[256];
	snprintf(moviename, sizeof(moviename),
	         "movie_follow-dBubble_%+06.3e%+06.3e.mp4",
	         dBubble, interface_location);
	save (moviename);
}


event interface_tracking (t += saveStepTime; t <= endTime) {
  if(pid() == 0){
    char name[256];

    snprintf(name, sizeof(name),
             "interface_f1-dBubble_%+06.3e%+06.3e-t_%08.3f.dat",
             dBubble, interface_location, t);
    FILE * fp1 = fopen(name, "w");
    output_facets(f1, fp1);
    fclose(fp1);

    snprintf(name, sizeof(name),
             "interface_f2-dBubble_%+06.3e%+06.3e-t_%08.3f.dat",
             dBubble, interface_location, t);
    FILE * fp2 = fopen(name, "w");
    output_facets(f2, fp2);
    fclose(fp2);
  }
}

event particle_tracking (t += saveStepTime; t <= endTime) {
  if(pid() == 0){
    char name[256];

    snprintf(name, sizeof(name),
             "particles01-dBubble_%+06.3e%+06.3e-t_%08.3f.dat",
             dBubble, interface_location, t);
    FILE * fp1 = fopen(name, "w");

    fprintf(fp1,
            "# %14s %14s %14s %14s %14s %14s %14s %10s\n",
            "t", "x", "y", "ux", "uy", "u2x", "u2y", "tag");

    foreach_particle_in (particles01)
      fprintf(fp1,
              "%+06.5e %+06.5e %+06.5e %+06.5e %+06.5e %+06.5e %+06.5e %010lu\n",
              t,
              p().locn.x, p().locn.y,
              p().u.x, p().u.y,
              p().u2.x, p().u2.y,
              p().tag);

    fclose(fp1);

    snprintf(name, sizeof(name),
             "particles02-dBubble_%+06.3e%+06.3e-t_%08.3f.dat",
             dBubble, interface_location, t);
    FILE * fp2 = fopen(name, "w");

    fprintf(fp2,
            "# %14s %14s %14s %14s %14s %14s %14s %10s\n",
            "t", "x", "y", "ux", "uy", "u2x", "u2y", "tag");

    foreach_particle_in (particles02)
      fprintf(fp2,
              "%+06.5e %+06.5e %+06.5e %+06.5e %+06.5e %+06.5e %+06.5e %010lu\n",
              t,
              p().locn.x, p().locn.y,
              p().u.x, p().u.y,
              p().u2.x, p().u2.y,
              p().tag);

    fclose(fp2);
  }
}
/*
Simple runtime post-processing.
Count the number of particles that are closer than a given distance from the 
bubble center (xcg,0) in time.*/
event count_particles_near_interface (i++; t <= endTime)
{
  double xcg = interface_x_center (f1);
  double ycg = 0.;

  int n01 = 0, n02 = 0;

  foreach_particle_in (particles01) {
    double r = sqrt(sq(p().x - xcg) + sq(p().y - ycg));
    if (r <= distance_from_cg)
      n01++;
  }

  foreach_particle_in (particles02) {
    double r = sqrt(sq(p().x - xcg) + sq(p().y - ycg));
    if (r <= distance_from_cg)
      n02++;
  }

  if (pid() == 0) {
    static FILE * fp = NULL;

    if (fp == NULL) {
      char name[256];
      snprintf(name, sizeof(name),
               "particles_near_cg-dBubble_%+06.3e-interface_location_%+06.3e-radius_%+06.3e.dat",
               dBubble, interface_location, distance_from_cg);

      fp = fopen(name, "w");
      fprintf(fp, "# %14s %14s %14s %14s %14s %14s\n",
              "t", "xcg", "ycg", "radius", "n_particles01", "n_particles02");
    }

    fprintf(fp, "%+06.5e %+06.5e %+06.5e %+06.5e %05d %05d\n",
            t, xcg, ycg, distance_from_cg, n01, n02);

    fflush(fp);
  }
}