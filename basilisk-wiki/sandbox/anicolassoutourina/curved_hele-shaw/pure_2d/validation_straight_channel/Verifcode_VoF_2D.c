#include "grid/quadtree.h" // 2D adaptable grid
#include "navier-stokes/centered.h" // Navier-Stokes solver
#include "navier-stokes/perfs.h" // Output of performance statistics on the various steps of solvers in a file named perfs, allows to see where code goes slow
#include "two-phase.h" // VoF solver, coupled with above Navier-Stokes solver
#include "navier-stokes/conserving.h" // Momentum-conserving advection of velocity for the VoF solver
#include "tension.h" // Surface tension, switch to integral.h for better Laplace pressure conservation, but only works in 2D and in Levelset, no VoF or coupled VoF-Levelset (cf Jun)
#include "maxruntime.h" // To have the simulation end cleanly 5 mins before the max allocated cluster time, instead of brutal cutoff (not sure if it works with srun ?)
#include "view.h" // Visualisation tools
#include "curvature.h" // necessary to get the bool interfacial function

/* ---------------------------------------------------------------------------- */
/* ----------------- Defining parameters of the simulation -------------------- */
/* ---------------------------------------------------------------------------- */

const bool Restore = false ; // Whether or not to restore from a "dump" file, that must be located in the working directory of the executable.

/* Geometric parameters of the channel (cf Gauri)) */
const double R_c = 0.75 ; // Radius of curvature (in mm here). If not specified in simu results, assume 0.75mm
const double W_channel = 1. ; // Width of the channel (in mm here). If not specified in simu results, assume 1mm
const double h_channel = 0.05 ; // Thickness of the microfluidic channel (in mm). Only used in this 2D simulation to have unified Re definition

const double BOXSIZE = (2*R_c + W_channel) ; // Smallest fitting of the curved channel in a 2D square box
const double BUBBLERADIUS = 0.1 ; // Radius of the bubble, in mm as well. If not specified in simu results, assume 0.1mm (100µm)

const double bubble_init_x = 0., bubble_init_y = -0.5*BOXSIZE + 2*BUBBLERADIUS; // (x,y) position of the bubble center at t = 0.

const double U_inlet = 0.1 ; // speed at the inlet, (in mm/s). FIXME : understand why simulations run faster for u_in = 0.1 compared to u_in = 1. , since all results are the same (Re equivalence works)

/* Defining parameters */
const double Re = 0.1 ; // Reynolds number, cf calculations on drafts
const double Ca = pow(0.894*(0.05/(2*0.1))*2., 3) ; // Capillary number, defined from Katie's Bretherton, pow is bitchy so have to initialize by hand

const double Mu_ratio = 1000. ; // Ratio mu_liq / mu_air
const double Rho_ratio = 1000. ; // Ratio rho_liq / rho_air

const double t_max = 15. ; // max time the simulation will run for

const int MAXLEVEL = 10 ; // Max grid refinement at 2^MAXLEVEL for timestep adaptation
const int INITLEVEL = 9 ; // Initial grid refinement

const double n_frames = 5. ; // nb-1 of image frames we want to output (10 makes a nice-ish movie, but probs unnecessary, if bubble deformation is not checked)
const double n_points = 250. ; // nb-1 of data points that will be logged for the plots
const double n_dumps = 4. ; // nb-1 of dumps for recovery

/* ---------------------------------------------------------------------------- */
/* -------------------------- Boundary conditions ----------------------------- */
/* ---------------------------------------------------------------------------- */

/* General boundary conditions */

// For validation, symmetry on all walls

/* Inlet and outlet boundary conditions */

u.n[bottom] = dirichlet(U_inlet) ; // for validation
u.n[top] = neumann(0.);
p[top] = dirichlet(0.) ; // for validation

/* ---------------------------------------------------------------------------- */
/* ---------------------- Initialisation and main loop ------------------------ */
/* ---------------------------------------------------------------------------- */


/* Main loop */
int main(){
  size(BOXSIZE); // Domain size
  origin(-BOXSIZE/2,-BOXSIZE/2); // Recenter domain
  init_grid(1 << INITLEVEL);
  
  /* Dimensionless definition of physical parameters of the fluids */
  rho1 = 1., mu1 = (rho1*U_inlet*pow(h_channel,2))/(2*BUBBLERADIUS*Re); // Liquid params, viscosity defined with Re as defined in Gauri's report
  rho2 = rho1/Rho_ratio , mu2 = mu1/Mu_ratio; // Air params
  f.sigma = (mu1*U_inlet)/Ca; // Surface tension, defined with Ca
  
  fprintf(stderr, "Re = %f, Ca = %f, mu1 = %f, sigma = %f",Re,Ca,mu1,f.sigma);
  
  TOLERANCE = 1e-4 [*]; // Makes solver more accurate, but more expensive : from https://basilisk.fr/src/examples/bubble.c, "important to minimise mass conservation errors for very long simulations"
  
  run();
}

/* Initialisation */
event init (t = 0){
  if (Restore){
    restore();
  }
  else{
    fraction (f, pow(x-bubble_init_x,2) + pow(y-bubble_init_y,2) - pow(BUBBLERADIUS,2));// VoF initialisation of the bubble
  }
}

/* ---------------------------------------------------------------------------- */
/* ---------- Bubble center tracking w/ interface reconstruction -------------- */
/* ---------------------------------------------------------------------------- */

/* Tracking of bubble location in time, from https://basilisk.fr/sandbox/fpicella/rising_bubble_to_droplet/axi_particles_no_coalescence_quantitative.c#interface_x_center, with interface reconstruction that only works in 2D (FIXME : actually could work in 3D, but would need to force functions used here to behave as 2D by rewriting them in a custom header, and use antoonvh's slice functions or equivalents to feed 2D slices of the 3D simulation to these 2D-forced functions, TODO) */

/* Helper function to extract x-position of bubble from the VoF color function w/ interface reconstruction inside cells */
/* Getting x-position of bubble from the VoF color function */
double interface_x_center (scalar f){
  double sx = 0., sl = 0.;

  foreach (reduction(+:sx) reduction(+:sl)) {
    if (interfacial(point, f)){
      coord n = interface_normal (point, f);
      double alpha = plane_alpha (f[], n);

      coord p[2];
      int m = facets (n, alpha, p);

      if (m == 2) {
        foreach_dimension(){
          p[0].x = x + p[0].x*Delta,
          p[1].x = x + p[1].x*Delta;
	}
	
        double l = sqrt(sq(p[1].x - p[0].x) + sq(p[1].y - p[0].y));

        double xm = 0.5*(p[0].x + p[1].x);

        sx += xm*l;
        sl += l;
      }
    }
  }
  return sl > 0. ? sx/sl : nodata;
}

/* ANS : same helper function but for y-position */
double interface_y_center (scalar f){
  double sy = 0., sl = 0.;

  foreach (reduction(+:sy) reduction(+:sl)) {
    if (interfacial(point, f)){
      coord n = interface_normal (point, f);
      double alpha = plane_alpha (f[], n);

      coord p[2];
      int m = facets (n, alpha, p);

      if (m == 2) {
        foreach_dimension(){
          p[0].y = y + p[0].y*Delta,
          p[1].y = y + p[1].y*Delta;
	}
	
        double l = sqrt(sq(p[1].x - p[0].x) + sq(p[1].y - p[0].y));

        double ym = 0.5*(p[0].y + p[1].y);

        sy += ym*l;
        sl += l;
      }
    }
  }
  return sl > 0. ? sy/sl : nodata;
}

// For the Bretherton forces, interface_x_center needs to be called once for the whole of each iteration, otherwise each interface cell calls it, O(n*n^2) bs

// Tracking of bubble center w/ interface reconstruction, from fpicella's https://basilisk.fr/sandbox/fpicella/rising_bubble_to_droplet/axi_particles_no_coalescence_quantitative.c#interface_x_center
double x_bubble_slice = bubble_init_x, y_bubble_slice = bubble_init_y ;

event bubble_center(i++){
  x_bubble_slice = interface_x_center(f);
  y_bubble_slice = interface_y_center(f);
  //fprintf (stderr, "Bubble center at x = %f , y = %f \n", x_bubble_slice, y_bubble_slice);
}

/* ---------------------------------------------------------------------------- */
/* ------------- Logging of all stability events' resulting dts --------------- */
/* ---------------------------------------------------------------------------- */

event log_stability(i++){

  // "tension.h" stability event, copied to log its dtmax
  double tension_dtmax = 0. ;
  double amin = HUGE, amax = -HUGE, dmin = HUGE;
  foreach_face (reduction(min:amin) reduction(max:amax) reduction(min:dmin))
    if (fm.x[] > 0.) {
      if (alpha.x[]/fm.x[] > amax) amax = alpha.x[]/fm.x[];
      if (alpha.x[]/fm.x[] < amin) amin = alpha.x[]/fm.x[];
      if (Delta < dmin) dmin = Delta;
    }
  double rhom = (1./amin + 1./amax)/2.;

  double tension_dt = sqrt (rhom*cube(dmin)/(pi*f.sigma));
  if (tension_dt < dtmax){
    tension_dtmax = tension_dt;
  }
  
  // "vof.h" 's stability event is just CFL < 0.5 ? CFL = 0.5 : nothing, so not too sure how to extract a dtmax from that ?
  
  // "navier-stokes/conserving.h" stability event, copied to log its dtmax
  double conserving_dtmax = timestep (uf, dtmax);
  
  // "navier-stokes/centered.h" stability event, copied to log its dtmax
  double centered_dtmax = dtnext (stokes ? dtmax : timestep (uf, dtmax));
  
  fprintf (stderr, "Log of the dtmax : tension.h dtmax = %f , conserving.h dtmax = %f , centered.h dtmax = %f ,\n", tension_dtmax, conserving_dtmax, centered_dtmax);
}

/* ---------------------------------------------------------------------------- */
/* ------------ Logging of bubble position, and image production -------------- */
/* ---------------------------------------------------------------------------- */

/* Writing of the position data into a centralized logfile */
char logfilename[256] = "bubble_positions.dat"; // Logfile initialization

//event logposition (i++){ // When testing on local machine
event logposition (t += (t_max/n_points)){

  // FIXME : each parallel process calls the .dat file, so #Nprocs repeats for the whole file. Annoying but not code-breaking, not sure if fixable from Basilisk ? https://stackoverflow.com/questions/23586682/how-to-use-printf-in-multiple-threads feels lower-level than what is done in Basilisk, should still be tried. Otherwise running "awk '!seen[$0]++' bubble_positions.dat" removes all duplicate lines in the logfile
  FILE * logfile = fopen(logfilename, "a+"); // Open logfile to write at the end of it
  
  fprintf(logfile, "%d %f %f %f\n", i, t, x_bubble_slice, y_bubble_slice) ; // Written into the logfile : iteration n°, time, x and y of bubble center with 2D slice and interface reconstruction
  fclose(logfile);
}

/* Image frames to view full movement of the bubble, and use Python's OpenCV to measure bubble deformation */
char bubbleview[80] ; // Declaration of file name for image frames
int view_nb = 0 ; // Frame counter for frame names

//event XY_bubbleview (i++){
event XY_bubbleview (t += (t_max/n_points)){
  view(quat={0.,0.,0.,1.}, bg = {1.,1.,1.}); // Shows XY plane
  squares ("f", min = 0., max = 1., linear = true);
  if (view_nb < 10.){ // Creates name of output file while making sure they are properly ordered in a terminal
    sprintf (bubbleview, "volfrac_XY-00%d.ppm", view_nb);
  }
  else if (view_nb < 100.){
    sprintf (bubbleview, "volfrac_XY-0%d.ppm", view_nb);
  }
  else {
    sprintf (bubbleview, "volfrac_XY-%d.ppm", view_nb);
  }
  save (bubbleview);
  clear();
  
  view_nb ++ ;
}

/* Generation of .ppm frames of the y velocity component, FIXME : rewrite BC for Neumann outlet to add velocity streamlines back */

char framename[80] ; // Init of file name for image frames
int frame_nb = 0 ; // Frame counter for frame names

//event movie (i++){ // For testing only on 5 iterations, runs locally on 4 cores just fine
event movie (t += (t_max/n_frames)) {

  // Velocity in y direction
  view(quat={0.,0.,0.,1.}, bg = {1.,1.,1.}); // Quaternion is Theta, e_x,e_y,e_z ? Here shows XY plane
  //draw_vof ("cs", filled = -1, fc = {1,1,1}); // Only draws following profiles for the non-solid cells
  squares ("u.y", min = -1.5*U_inlet, max = 1.5*U_inlet, linear = true); // Draws y velocity component inside squares
  colorbar(pos = {0.8,-0.25}, label = "u_y");
  box(); // Adds coordinates of box
  sprintf (framename, "y_velocities-%d.ppm", frame_nb); // Creates name of output file
  save (framename);
  clear();
  
  // Density images
  static FILE * densityfile = fopen ("density.ppm", "w"); // Creates name of output file
  output_ppm (rho, densityfile, min = rho2, max = rho1, linear=true); // FIXME : Maybe cleaner way to do this w/ bview, but squares and co. functions don't know what "rho" is. Advantage here is that all frames are saved into one .ppm function, and Imagemagick's "animate" can then show its evolution
  
  // Pressure images
  static FILE * pressurefile = fopen ("pressure.ppm","w");
  output_ppm(p, pressurefile, linear = true);
    
  // y-Acceleration (forces) images
  static FILE * accfile = fopen ("y-acceleration.ppm","w");
  output_ppm(a.y, accfile, linear = true);
  
  // Frame counter
  frame_nb ++ ;
}

// Dumping of the whole simulation fields in case of crash, to be able to restart somewhat close to the end
//event recovery_dump (i++){
event recovery_dump (t += (t_max/n_dumps)){
  dump(); // By default, dumps all fields in a file named "dump"
}

// Force end of simulation
//event end (i = 2){
event end (t=t_max){
  return 1 ;
}