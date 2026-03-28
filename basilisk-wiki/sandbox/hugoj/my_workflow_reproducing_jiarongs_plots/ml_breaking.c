/**
  Field scale wave breaking (multilayer solver)

No stratification


USAGE

  * Compilation and run (with mpirun)
  make ml_breaking.tst

  * HPC source file generation
  make _ml_breaking.c

  * Running on HPC

  -> first compile on hpc with
  mpicc -std=c99 -O2 _*.c

  -> then run (in a slurm file)
  srun ./ml_breaking
TODO: more details here

HOW TO CREATE A RESTART
  ncks -d time,N out.nc restart.nc

  with N = index of timestep in the file "out.nc"
*/


#include "grid/multigrid.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
#include "bderembl/libs/extra.h"      // parameters from namlist
#include "bderembl/libs/netcdf_bas.h"


#define g_ 9.81
#include "spectrum.h" // Initial conditions generation


/*

TO DO : 
- add dimensions so that they can be saved in the netcdf
- add ability to restart from netcdf file

*/

char my_name[40] = "ml_breaking";

/*
DEFAULT PARAMETERS

Dimensions : [Length, Time]
*/
char namlist[80] = "namelist.toml";    // file name of namlist
char file_out[20] = "out.nc";         // file name of output
char file_restart[20] = "restart.nc"; // file name of restart
// -> Initial conditions
double P = 0.2 [1, -1];     // energy level (estimated so that kpHs is reasonable)
int coeff_kpL0 = 10 [];     // kpL0 = coeff_kpL0 * pi
int N_mode = 32 [];         // Number of modes in wavenumber space
int N_power = 5 [];         // directional spreading coeff
int F_shape = 0 [];         // shape of the initial spectrum
// -> Domain definition
int N_grid = 5 [];       // 2^N_grid : number of x and y gridpoints
double L = 200.0 [1];       // domain size
int N_layer = 2 [];         // number of layers
double kp = PI*10/200.0 [-1];// peak wave number
double h0 = 1.0 [1];        // depth of water
// -> Runtime parameters
int restart = 0;            // 1: restart, 0: no restart
double tend = 2.0;          // end time of simulation
// -> saving outputs
double dtout = 2.0;        // dt for output in netcdf
int pad = 4;                // number of 0-padding for ouput files
int nout = 1;               // number of the outfile
char fileout[100];          // name of outfile
// -> physical properties
double nu0 = 0.;            // Viscosity for vertical diffusion
int RANDOM = 2;             // For random number generator
double thetaH = 0.5;        // theta_h for dumping fast barotropic modes


// diag
//double *dudz, *eps, *u_profile;
double *u_profile;
double dt_mean = 0.4;


static FILE * fp;


int main(int argc, char *argv[])  
{
   
  
  // #if _MPI
  //   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  //   MPI_Comm_size(MPI_COMM_WORLD, &msize);
  // #endif
  // Building a 'params' array with all parameters from the namlist
  params = array_new();
  add_param("P", &P, "double");
  add_param("coeff_kpL0", &coeff_kpL0, "int");
  add_param("N_mode", &N_mode, "int");
  add_param("N_power", &N_power, "int");
  add_param("F_shape", &F_shape, "int");
  add_param("N_grid", &N_grid, "int");
  add_param("L", &L, "double");
  add_param("N_layer", &N_layer, "int");
  add_param("h0", &h0, "double");
  add_param("tend", &tend, "double");
  add_param("nu0", &nu0, "double");
  add_param("RANDOM", &RANDOM, "int");
  add_param("thetaH", &thetaH, "double");
  add_param("restart", &restart, "int");
  add_param("dtout", &dtout, "double");
  kp = PI * coeff_kpL0 / L; // kpL=coeff x pi peak wavelength
  
  add_param("dt_mean", &dt_mean, "double");

  // Search for the configuration file with a given path or read params.in
  if (argc == 2)
    strcpy(file_param, argv[1]);
  else
    strcpy(file_param, namlist);
  read_params(file_param);
  
  
  // Settings solver values from namlist values
  L0 = L;
  nu = nu0;
  N = 1 << N_grid; // 1*2^N_grid
  nl = N_layer;
  G = g_;
  theta_H = thetaH;
  CFL_H = 1; // Smaller time step

  // Boundary condition
  origin (-L0/2., -L0/2.);
  periodic (right);
  periodic (top);
  
  // allocate diag
  // dudz = (double*)malloc(nl * sizeof(double));
  // eps = (double*)malloc(nl * sizeof(double));
  u_profile = (double*)calloc(nl, sizeof(double));
  fp  = fopen("u_profile.dat","w"); // reset file
  fclose(fp);
    
  fprintf (stderr, "Read in parameters!\n");
  run();
}

event init(i =  0) {
  geometric_beta (1./3., true); // Varying layer thickness
  if (restart!=1) {
    // Generate linealy spaced kx, ky according to specified # of modes, and
    //  interpolated F(kx,ky)
    T_Spectrum spectrum;
    spectrum = read_spectrum(N_mode); // TO DO:

    // step 1: set eta and h
    foreach() {
      zb[] = -h0;
      eta[] = wave(x, y, N_grid, spectrum);
      double H = wave(x, y, N_grid, spectrum) - zb[];
      foreach_layer() {
        h[] = H/nl;
      } 
    }
    
    // step 2: remap
    vertical_remapping (h, tracers);
    // step 3: set currents
    foreach() {
      double z = zb[];
      foreach_layer() {
        z += h[]/2.;
        u.x[] = u_x(x, y, z, N_grid, spectrum);
        u.y[] = u_y(x, y, z, N_grid, spectrum);
        w[] = u_z(x, y, z, N_grid, spectrum);
        z += h[]/2.;
      }
    }
  }
  else {
    //fprintf(stderr, "restart = %d\n",restart);
    // Restarting from netcdf ...
    // "test, %s", strcat(strcat(my_name,"/"),file_restart));
    // read_nc({zb, h, u, w}, strcat(strcat(my_name,"/"),file_restart));
    fprintf(stderr, "Restarting : reading from file\n"); 
    fprintf(stderr, "->no yet coded ...\n");
    // read_nc({zb, h, u, w}, file_restart);
    //
  }


  fprintf (stderr,"Done initialization!\n");
  // sprintf(fileout, "%0*d.nc", pad, 0); // add the padding
  create_nc({zb, h, u, w, eta}, file_out);

}


// event film surface


// event images (t += dtout) {
//   scalar img[];
//   foreach()
//     img[] = u.x[0,0,l];
//   static FILE * fp = fopen ("grid.ppm", "w");
//
//   output_ppm (img, fp, min = -2, max = 2);
// }
//
// event dumpini (t=0.){
//   dump();
// }

// This event compute layer average of u.x
//
// NOTE: for faster run, I compute only instant horizontal mean
// but then there is the dt/dt_mean that I have to correct after
//
//event compute_horizontal_avg (i++; t<=tend+1e-10){
// event compute_horizontal_avg (t+=dt_mean; t<=tend+1e-10){
//   foreach(reduction(+:u_profile[:nl]))
//     foreach_layer(){
//       //u_profile[point.l] += u.x[] / (N*N) * dt / dt_mean;
//       u_profile[point.l] += u.x[] / (N*N); // * dt / dt_mean;
//
//     }
// }
//
// // This even writes to a file the layer average
// event write_diag(t=0., t+=dt_mean){
//     // main worker is writing the file
//     if (pid()==0) {
//       fp  = fopen("u_profile.dat","a");
//       if (fp == NULL){
//         fprintf(stderr, "Error opening file u_profile.dat");
//         return 2;
//       }
//       for (int i=0; i<nl; ++i) {
//         fprintf (fp, "%f %d %g\n", t, i, u_profile[i]);
//       }
//       fprintf(fp,"\n");
//       fclose(fp);
//     }
//     // Reset the profile for all workers
//     for (int i=0; i<nl; ++i) {
//       u_profile[i] = 0.0;
//     }
//
//
// }

  // vector gradU[];
  // vector gradV[];
  // dudx = u.x[]
  //
  //
  // gradient(u.x, gradU);
  // gradient(u.y, gradV);
  //
  //
  //
  // foreach(reduction(+:dudz), reduction(+:eps)){
  //   foreach_layer(){
  //
  //   }
  // }
 

// event image (t = end) {
//   clear();
//   static FILE * fp = fopen ("image.ppm", "w");
//   output_ppm(eta, fp, min=-0.1, max=0.1);
//
// }

// event snapshot (t = end)
// {
//   clear();
	
  // 3/4 view
  // view (quat = {0.567, 0.137, 0.196, 0.789},
  //     fov = 30, near = 0.01, far = 1000,
  //     tx = 0.048, ty = -0.001, tz = -2.096,
  //     width = 1398, height = 803);
  //
  //
  // char s[80];
  // sprintf (s, "t = %.2f", t);
  // draw_string (s, size = 80);
  // for (double x = -1; x <= 1; x++)
  //   translate (x) {
  //     squares ("eta", linear = true, z = "eta*10", min = -1.0, max = 1.0 , map=gray);
  //   }
  // box ();
  //
  //
  //
  //
  //
  // colorbar(map=gray, label="eta (m)", min=-1.0,max=1.0, pos={-0.95,-0.5},
  //          levels=10);
  // save ("snap_side.png");
  //
  //
  // // top view
  // clear();
  // view (quat = {0.000, 0.000, 0.000, 1.000},
  //     fov = 30, near = 0.01, far = 1000,
  //     tx = 0.000, ty = 0.000, tz = -2.505,
  //     width = 1398, height = 803);
  //
  //
  // sprintf (s, "t = %.2f", t);
  // draw_string (s, size = 80);
  // for (double x = -1; x <= 1; x++)
  //   translate (x) {
  //     squares ("eta", min = -1.0, max = 1.0 , map=gray);
  //   }
  // box ();
  //
  // colorbar(map=gray, label="eta (m)", min=-1.0,max=1.0, pos={-0.95,-0.5},
  //          levels=10);
  // save ("snap_top.png");


//}

// Writing a 4D netcdf file
event output(t = 0.; t<= tend+1e-10; t+=dtout){
  // netcdf output
  write_nc();
  // regular dump
  char dname[100];
  sprintf (dname, "dump_t%g", t);
  dump(dname);
}

// Clean for my diag of layer avg
event cleanup(t=end){
  free(u_profile);
}

/**
Results: plots
~~~gnuplot 
# Plot the heatmap
set pm3d map
set view map
set xlabel "Time (s)"
set ylabel "Layer"
set cblabel "u.x (m/s)"
set yrange [0:29]
set xrange [0:200]
set terminal pngcairo size 800,600 enhanced font 'Verdana,12'
set output 'u_profile.png'
set size 0.9, 0.9
splot "u_profile.dat" using 1:2:3 with pm3d
unset output
~~~
**/

