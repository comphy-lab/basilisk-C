#include "grid/multigrid3D.h" // 3D non-cubic grid
#include "embed.h" // Allows definition of solid cells and associated boundary conditions
#include "navier-stokes/centered.h" // Navier-Stokes solver
#include "navier-stokes/perfs.h" // Output of performance statistics on the various steps of solvers in a file named perfs, allows to see where code goes slow
#include "two-phase.h" // VoF solver, coupled with above Navier-Stokes solver
#include "navier-stokes/conserving.h" // Momentum-conserving advection of velocity for the VoF solver
#include "tension.h" // Surface tension, switch to integral.h for better Laplace pressure conservation, but only works in 2D and in Levelset, no VoF or coupled VoF-Levelset (cf Jun)
#include "maxruntime.h" // To have the simulation end cleanly 5 mins before the max allocated cluster time, instead of brutal cutoff (not sure if it works with srun ?)
#include "view.h" // Visualisation tools
#include "ANS_2D_Slice_bubble_tracking.h" // Tools for bubble position tracking via 2D interface reconstruction, cf fpicilla's algorithm : https://basilisk.fr/sandbox/fpicella/rising_bubble_to_droplet/axi_particles_no_coalescence_quantitative.c#interface_x_center

/* ---------------------------------------------------------------------------- */
/* ----------------- Defining parameters of the simulation -------------------- */
/* ---------------------------------------------------------------------------- */

const bool Restore = false ; // Whether or not to restore from a "dump" file, that must be located in the working directory of the executable.

/* Geometric parameters of the channel (cf Gauri)) */
const double R_c = 0.75 ; // Radius of curvature (in mm here)
const double W_channel = 1. ; // Width of the channel (in mm here)
const double h_channel = 0.05 ; // Thickness of the microfluidic channel (in mm)

const double BOXSIZE = (2*R_c + W_channel) ; // Smallest fitting of the curved channel in a cubic box
const double BOXSIZE_real = 1. ; // With weirdness due to "dimensions(nx,ny,nz)"
const double BUBBLERADIUS = 0.1 ; // Radius of the bubble, in mm as well. If not specified in simu results, assume 0.1mm (100µm)
const double BUBBLEHEIGHT = 0.5*h_channel ; // Total height of the bubble initially
const double z_size = 0.64 ;

const double x_offset = 0.5, y_offset = 0.5, z_offset = 0.39*z_size ; // FIXME : dimensions work strangely. All X and Y lenths have to be divided by 2.5, BOXSIZE
const double bubble_init_x = (-R_c/2.5) + x_offset ; // WARNING : cannot divide by any variable while initializing constants
const double bubble_init_y = ((-0.5*BOXSIZE + 2*BUBBLERADIUS)/2.5) + y_offset ; // (x,y) position of the bubble center at t = 0.

const double U_inlet = 0.04 ; // speed at the inlet, (in mm/s) Have to divide by 2.5 from 0.1 value since it is a value in the XY plane for a non cubic domain. FIXME : understand why simulations run faster for u_in = 0.1 compared to u_in = 1. , since all results are the same (Re equivalence works)

/* Defining parameters */
const double Re = 0.1 ; // Reynolds number, cf calculations on drafts
const double Ca = 0.1 ; // Capillary number, cf calculations on drafts

const double Mu_ratio = 1000. ; // Ratio mu_liq / mu_air
const double Rho_ratio = 1000. ; // Ratio rho_liq / rho_air

const double t_max = 9. ; // max time the simulation will run for (in s)

const int MAXLEVEL = 9 ; // Max grid refinement at 2^MAXLEVEL for timestep adaptation
const int INITLEVEL = 9 ; // Initial grid refinement

const double n_frames = 5. ; // nb-1 of image frames we want to output (10 makes a nice-ish movie, but probs unnecessary, if bubble deformation is not checked)
const double n_points = 250. ; // nb-1 of points for log plots and deformation images
const double n_dumps = 4. ; // nb-1 of dumps for recovery

/* ---------------------------------------------------------------------------- */
/* -------------------------- Boundary conditions ----------------------------- */
/* ---------------------------------------------------------------------------- */

/* General boundary conditions */
u.n[embed] = dirichlet(0.); // 0-velocity boundary conditions for the solid embedded walls
u.t[embed] = dirichlet(0.);
u.r[embed] = dirichlet(0.); // Lol, be careful not to forget this one for 3D simulations, everything aborts while running otherwise.

u.t[left] = dirichlet(0.); // No-slip condition for the simulation box's walls
u.t[right] = dirichlet(0.);
u.t[front] = dirichlet(0.);
u.t[back] = dirichlet(0.);
u.t[top] = dirichlet(0.);

/* Inlet and outlet boundary conditions */

u.n[bottom] = (z-z_offset) > 0.5*h_channel ? dirichlet(0.) : ( (z-z_offset) < -0.5*h_channel ? dirichlet(0.) : ((x-x_offset) > (R_c - 0.5*W_channel)/2.5 ? dirichlet(-U_inlet) : ( (x-x_offset) < (-R_c + 0.5*W_channel)/2.5 ? neumann(0.) : dirichlet(0.) )) );

p[bottom] = (z-z_offset) > 0.5*h_channel ? dirichlet(0.) : ( (z-z_offset) < -0.5*h_channel ? dirichlet(0.) : ( (x-x_offset) > (R_c - 0.5*W_channel)/2.5 ? neumann(0.) : ( (x-x_offset) < (-R_c + 0.5*W_channel)/2.5 ? dirichlet(0.) : neumann(0.) ) ) );

/* ---------------------------------------------------------------------------- */
/* ---------------------- Initialisation and main loop ------------------------ */
/* ---------------------------------------------------------------------------- */

double bubble_init(double x, double y, double z){ // Helper function to initialize the shape of the bubble : here, cylinder of radius BUBBLERADIUS and height 50% of channel height w/ rounded sides
  if (fabs(z-z_offset) > 0.5*BUBBLEHEIGHT){
    return (fabs(z-z_offset) - 0.5*BUBBLEHEIGHT);
  }
  else{
    double effective_radius = BUBBLERADIUS/2.5 - 0.5*BUBBLEHEIGHT + sqrt(pow(0.5*BUBBLEHEIGHT,2) - pow(z-z_offset,2)); // Effective radius, to have a "cylindrical bubble" with semicircles on its sides, closed torus
    return (pow(x-bubble_init_x,2) + pow(y-bubble_init_y,2) - pow(effective_radius,2)) ;
  }
}

double curve(double x, double y, double z){ // Helper function to define the channel shape via solid embedded cells
  if (fabs(z-z_offset) > 0.5*h_channel){
    return (0.5*h_channel - fabs(z-z_offset)) ;
  }
  else if ((y-y_offset) < 0.){
    return (fabs(x-x_offset) - (R_c - 0.5*W_channel)/2.5) ; // entry-exit channels
  }
  else {
    if (pow(x-x_offset,2) + pow(y-y_offset,2) > pow(R_c/2.5,2)){
      return (pow((R_c + 0.5*W_channel)/2.5,2) - (pow(x-x_offset,2)+pow(y-y_offset,2))) ; // extrados
    }
    else {
      return ((pow(x-x_offset,2)+pow(y-y_offset,2)) - pow((R_c - 0.5*W_channel)/2.5,2)) ; // intrados
    }
  }
}

/* Main loop */
int main(){
  dimensions(nx = BOXSIZE, ny = BOXSIZE, nz = z_size);
  init_grid(1 << INITLEVEL);
  
  /* Dimensionless definition of physical parameters of the fluids */
  rho1 = 1., mu1 = (rho1*U_inlet*pow(h_channel,2))/(2*BUBBLERADIUS*Re); // Liquid params, viscosity defined with Re as defined in Gauri's report
  rho2 = rho1/Rho_ratio , mu2 = mu1/Mu_ratio; // Air params
  f.sigma = (mu1*U_inlet)/Ca; // Surface tension, defined with Ca
  
  fprintf(stderr, "Re = %f, Ca = %f, mu1 = %f, sigma = %f",Re,Ca,mu1,f.sigma);
  
  TOLERANCE = 1e-4 [*]; // Makes solver more accurate, but more expensive : from https://basilisk.fr/src/examples/bubble.c, "important to minimise mass conservation errors for very long simulations" TODO : Jun suggested forcefully decreasing TOL to a lower value if convergence problems with poisson.h
  
  run();
}

/* Initialisation */
event init (t = 0){
  if (Restore){
    restore();
  }
  else{
    fraction(f, bubble_init(x,y,z));// VoF initialisation of the bubble
    //fractions_cleanup (cs, fs); // https://groups.google.com/g/basilisk-fr/c/SiBCu4txrmA, not sure if useful
    solid(cs, fs, curve(x,y,z)); // Solid initialisation of the curved Hele-Shaw channel
  }
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

/* ANS : Helper function to define the position of the centerline of the channel, to have a measure of the difference between the bubble's center position and the geometric centerline of the channel : Radial distance w/ origin on the centerline : centerline_distance > 0. means bubble center on the outside of the centerline */
/* Take x,y position (real ones, with correct dimensions and offsets) and return the distance to the closest centerline of the channel */
double centerline_distance (double x_in, double y_in){
  if (y_in < 0.){ // If bubble is in the straight entry/exit channels, distance to centerline is straightforward-ish : x - (+- R_c)
    return x_in < 0. ? -R_c -x_in : x_in - R_c ; // Weird calculation, but it's to always have "centerline_distance > 0. means bubble center on the outside of the centerline"
  }
  else { // If bubble is in the curved channel, since perfect half-disc, distance to centerline is : sqrt(x^2 + y^2) - R_c
    return sqrt(pow(x_in,2) + pow(y_in,2)) - R_c ;
  }
}

/* Writing of the position data into a centralized logfile (TODO : probably other metrics to process, get them in there as well : e.g. difference in velocities of the bubble compared to the fluid) */
char logfilename[256] = "bubble_positions.dat"; // Logfile initialization

//event logposition (i++){ // When testing on local machine for a few iterations
event logposition (t += (t_max/n_points)){
  
  // Tracking of bubble center via slice in 2D, with interface reconstruction, from ANS_2D_Slice_bubble_tracking.h
  double x_bubble_slice = 2.5*(ANS_interface_x_center(f,MAXLEVEL, BOXSIZE_real, x_offset, y_offset, z_offset) - x_offset) ;
  double y_bubble_slice = 2.5*(ANS_interface_y_center(f,MAXLEVEL, BOXSIZE_real, x_offset, y_offset, z_offset) - y_offset) ; // 1. instead of BOXSIZE for domain size since coords are weird ?

  // FIXME : each parallel process calls the .dat file, so #Nprocs repeats for the whole file. Annoying but not code-breaking, not sure if fixable from Basilisk ? https://stackoverflow.com/questions/23586682/how-to-use-printf-in-multiple-threads feels lower-level than what is done in Basilisk, should still be tried. Otherwise running "awk '!seen[$0]++' bubble_positions.dat" removes all duplicate lines in the logfile
  FILE * logfile = fopen(logfilename, "a+"); // Open logfile to write at the end of it
  
  fprintf(logfile, "%d %f %f %f %f\n", i, t, x_bubble_slice, y_bubble_slice, centerline_distance(x_bubble_slice,y_bubble_slice)) ; // Written into the logfile : iteration n°, time, x and y of bubble center with 2D slice and interface reconstruction, distance to centerline with this x and y
  fclose(logfile);
}

/* Image frames to view full movement of the bubble, and use Python's OpenCV to measure bubble deformation */
char bubbleview[80] ; // Declaration of file name for image frames
int view_nb = 0 ; // Frame counter for frame names

//event XY_bubbleview (i++){
event XY_bubbleview (t += (t_max/n_points)){
  view(quat={0.,0.,0.,1.}, tx = -BOXSIZE_real/2, ty = -BOXSIZE_real/2, bg = {1.,1.,1.});
  squares ("f", alpha= z_offset, min = 0., max = 1., linear = false);
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

/* Generation of .ppm frames of the y velocity component */
char framename[80] ; // Declaration of file name for image frames
int frame_nb = 0 ; // Frame counter for frame names

//event movie (i++){ // For testing only on a few iterations
event movie(t+=(t_max/n_frames)){

  // Velocity in y direction
  view(quat={0.,0.,0.,1.}, tx = -BOXSIZE_real/2, ty = -BOXSIZE_real/2, bg = {1.,1.,1.}); // Quaternion is Theta, e_x,e_y,e_z ? Here shows XY plane
  squares ("u.y", alpha = z_offset, min = -1.5*U_inlet, max = 1.5*U_inlet, linear = false); // Draws y velocity component inside squares
  colorbar(pos = {0.8,-0.25}, label = "u_y"); // TODO : Add max/min value labels ?
  box(); // Adds coordinates of box
  sprintf (framename, "uy_XY-%d.ppm", frame_nb); // Creates name of output file
  save (framename);
  clear();
  
  view(theta=-1.5708,tx = -BOXSIZE_real/2, ty = -BOXSIZE_real/2, bg = {1.,1.,1.}); // Allows visualization on the side (i.e. in YZ plane)
  squares ("u.y", alpha = bubble_init_x, min = -1.5*U_inlet, max = 1.5*U_inlet, n = {1.,0.,0.}, linear = false); // Draws y velocity component inside squares in the middle of the entry channel, for a YZ plane
  colorbar(pos = {0.8,-0.25}, label = "u_y"); // TODO : Add max/min value labels ?
  box(); // Adds coordinates of box
  sprintf (framename, "uy_YZ-%d.ppm", frame_nb); // Creates name of output file
  save (framename);
  clear();
  
  // Color function i.e. air-liquid interface visualization : XY unnecessary with XY_bubbleview event
  view(theta=-1.5708,tx = -BOXSIZE_real/2, ty = -BOXSIZE_real/2,  bg = {1.,1.,1.}); // Allows visualization on the side (i.e. in YZ plane)
  squares ("f", alpha = bubble_init_x, min = 0., max = 1., n = {1.,0.,0.}, linear = false); // Draws VoF color function in the middle of the entry channel, for a YZ plane
  box(); // Adds coordinates of box
  sprintf (framename, "volfrac_YZ_entry-%d.ppm", frame_nb); // Creates name of output file
  save (framename);
  clear();
  
  // Pressure
  view(quat={0.,0.,0.,1.}, tx = -BOXSIZE_real/2, ty = -BOXSIZE_real/2, bg = {1.,1.,1.}); // Quaternion is Theta, e_x,e_y,e_z ? Here shows XY plane
  squares ("p", linear = false);
  colorbar(pos = {0.8,-0.25}, label = "P"); // TODO : Add max/min value labels ?
  box(); // Adds coordinates of box
  sprintf (framename, "P_XY-%d.ppm", frame_nb); // Creates name of output file
  save (framename);
  clear();

  // Frame counter
  frame_nb ++ ;
}

// Dumping of the whole simulation fields in case of crash, to be able to restart somewhat close to the end
//event recovery_dump (i++){
event recovery_dump (t += (t_max/n_dumps)){
  dump(); // By default, dumps all fields in a file named "dump"
}

// Force end of simulation
//event end (i = 5){
event end (t=t_max){
  return 1 ;
}