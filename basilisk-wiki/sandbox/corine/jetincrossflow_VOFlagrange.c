/**
# Atomisation of a liquid jet in cross-flow

We solve the two-phase Navier--Stokes equations with surface
tension and momentum-conserving advection of velocity.
Droplets not expected to breakup anymore are transfered to the Lagrange formalism.*/

#include <sys/types.h>
#include <sys/stat.h>

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "navier-stokes/conserving.h"
#include "tag.h"
#include "lagrange/lagrangenocoales.h"
#include "view.h"
#include "maxruntime.h"
#include "save_data.h"
/** Antoonvh sandbox [scatter2.h](https://basilisk.fr/sandbox/Antoonvh/scatter2.h) */
#include "lagrange/scatter2.h"

/**
## Input parameters 

Isothermal simulations at temperature 20°C, air pressure 5.8 bar, air density 7.19, air velocity 100 m/s, 
momentum flux ratio 6.0, orifice diameter 0.45 mm 155 mm downstream of the inlet in a lean premixed prevaporised channel.
This the baseline case in the paper [website](https://doi.org/10.1016/S1270-9638(01)01135-X).
The VOF simulations covers only a part of the channel near injection.
 */

const double T_END = 0.001; //1e-1

const double DIAMETER = 4.50e-4;

const double SIZE_FACTOR = 80.0;

const double SIZE_FACTOR_INLET_FROM_ORIFICE = 2.5;
const double SIZE_UNREFINE = 10.0;

const double X_ORIFICE = 0.155;

const double PRESSURE = 5.8e5;

/** Properties for Kerosin Jet A-1 and air. */
const double RHO_G = 7.19;
const double RHO_L = 789.269;

const double MU_G = 1.75e-5;
const double MU_L = 1.422e-3;
const double SIGMA = 25.73e-3;

const double VEL_GAS = 100.0;
const double DELTA_BL = 2.5e-3;
/** Velocity of the liquid calculated from momentum flux ratio
$$ q=\frac{\rho_l U_l^2}{\rho_g U_g^2} $$
*/
const double VEL_LIQ_MEAN = 22.99;

/**
The default maximum level of refinement */
int maxlevel = 10;
int inigrid = 128;

/* Error thresholds for adaptation */
double uemax = 1.0;
double femax = 0.001;

/**
## Functions
 */

/** Compute Reynolds number */
double Re_x(double rho, double velocity, double mu, double x_pos)
{
  double Re;
  Re = rho * velocity * x_pos / mu;
  return Re;
}

/** Compute turbulent boundary layer thickness */
double delta_BL_t(double Re_x, double x_pos)
{
  double delta_t;
  delta_t = 0.37 * x_pos / pow(Re_x, 0.2);
  return delta_t;
}

/** Compute velocity at a specific vertical position  $y_pos$  in the turbulent boundary layer (but the inflow is not yet turbulent). */
double vel_BL_t(double U_infinity, double delta_BL, double y_pos)
{
  double vel_at_y;
  vel_at_y = U_infinity * pow((y_pos / delta_BL), (1. / 7.));
  return vel_at_y;
}

/**
## Boundary conditions
 */

/**To impose boundary conditions on a disk we use an auxilliary volume
fraction field *f0* which is one inside the cylinder and zero
outside. */
scalar f0[];

/** Left boundary (air inlet) */
u.n[left] = dirichlet((y < DELTA_BL) ? (VEL_GAS * pow((y / DELTA_BL), (1. / 7.))) : VEL_GAS);
u.t[left] = dirichlet(0);
#if dimension > 2
u.r[left] = dirichlet(0);
#endif
p[left] = neumann(0);
f[left] = 0.;

/** Right boundary (outlet) */
u.n[right] = neumann(0);
p[right] = dirichlet(PRESSURE - RHO_G * VEL_GAS * VEL_GAS);

/** Bottom boundary (wall and inlet for the liquid). Parabolic profile within the liquid jet through the orifice, hydrophobic wall elsewhere (no-slip and f=0). */
u.n[bottom] = dirichlet(sq(x) + sq(z) < sq(0.5 * DIAMETER) ? 2.0 * VEL_LIQ_MEAN * (1.0 - (sq(x) + sq(z)) / (sq(0.5 * DIAMETER))) : 0.0);
u.t[bottom] = dirichlet(0);
#if dimension > 2
u.r[bottom] = dirichlet(0);
#endif
f[bottom] = f0[];

/** Top boundary (symmetry) */
u.n[top] = neumann(0);
u.t[top] = neumann(0);
#if dimension > 2
u.r[top] = neumann(0);
#endif
f[top] = 0.;

/**
## Main driver routine
 */
int main(int argc, char *argv[])
{
/**
Guarantee for a clean end of simulation before walltime exceeds */
  maxruntime(&argc, argv);

/**  The program can take six optional command-line arguments:

the maximum run time with the -m  option (above) 
  in the format H:M:S (hours, minutes and seconds) - At this
  time minus 5 minutes the state of the simulation is dumped 
  in the “restart” file and the program terminates-,
the maximum level for refinement,
the initial grid,
the error threshold on velocity,
the error threshold on the embedded geometry. 
the error threshold on the VOF variable. */
  if (argc > 1)
    maxlevel = atoi(argv[1]);
  if (argc > 2)
    inigrid = atoi(argv[2]);
  if (argc > 3)
    uemax = atof(argv[3]);
  if (argc > 4)
    femax = atof(argv[4]);

/**
The initial domain is discretised with  inigrid  grid points. We set the origin and domain size. */
  init_grid(inigrid);
  origin(-SIZE_FACTOR_INLET_FROM_ORIFICE * DIAMETER, 0.0, -0.5 * SIZE_FACTOR * DIAMETER);
  size(SIZE_FACTOR * DIAMETER);

/**
Set the density and viscosity of each phase as well as the
surface tension coefficient */
  rho1 = RHO_L;
  rho2 = RHO_G;

  mu1 = MU_L;
  mu2 = MU_G;

  f.sigma = SIGMA;

/**
Set timestep for Lagrange Solver. */
  dtlag = 0.5e-5;
  lagstep = false;

  fcut = 1e-4;
/* Start the simulation. */
  run();
}

/**
## Initial conditions
 */
event init(t = 0)
{
/** Prepare output */
  struct stat st_res = {0};
  struct stat st_mov = {0};

/** If it does not exist, create output folder for results, movies and particles with permissions for all*/
  if (stat("results", &st_res) == -1)
  {
    mkdir("results", S_IRWXU | S_IRWXG | S_IRWXO);
  }

  if (stat("movies", &st_mov) == -1)
  {
    mkdir("movies", S_IRWXU | S_IRWXG | S_IRWXO);
  }

  if (stat("particle", &st_mov) == -1)
  {
    mkdir("particle", S_IRWXU | S_IRWXG | S_IRWXO);
  }

  bool restDNS = restore (file = "restart");
  int restLAG = prestore ("restartp",NULL,true);

  if (!restDNS && !restLAG) {
/** Use a static refinement down to *maxlevel* in a cylinder a bit
longer than the initial jet and with twice the radius. */

    refine(y < 2.0 * DIAMETER && sq(x) + sq(z) < 2. * sq(0.5 * DIAMETER) && level < maxlevel);

/** Initialize the auxilliary volume fraction field for a cylinder with constant radius. */

    fraction(f0, sq(0.5 * DIAMETER) - sq(x) - sq(z));
    f0.refine = f0.prolongation = fraction_refine;
    restriction({f0});

/** Set the intial conditions (field values) */
    foreach ()
    {
/** Use the aux. volume field  to define the initial jet and its velocity. */
      f[] = f0[] * (y < DIAMETER);

/** Initialize the x-velocity field with an analytical boundary layer solution.
Note that the distance from the inlet of the test section in the experiment 
to the location if the orifice has to be added in the computation of the BL */
      double Re_at_x = Re_x(RHO_G, VEL_GAS, MU_G, x + X_ORIFICE);
      double bl_thickness_at_x = delta_BL_t(Re_at_x, x + X_ORIFICE);

      u.x[] = (y < bl_thickness_at_x) ? vel_BL_t(VEL_GAS, bl_thickness_at_x, y) : VEL_GAS;
      u.x[] = (sq(x) + sq(z) < sq(0.5 * DIAMETER))
                  ? ( 1.0 - f[]) * u.x[] : u.x[];      


/** Initialize the jet y-velocity with a parabolic profile elsewhere set 0.0
Note that $u_max = 2.0 * u_mean$ in case of a parabolic profile */
      u.y[] = (sq(x) + sq(z) < sq(0.5 * DIAMETER))
                  ? f[] * (2.0 * VEL_LIQ_MEAN * (1.0 - (sq(x) + sq(z)) / (sq(0.5 * DIAMETER))))
                  : 0.0;
    }

    int cellnumber = 0;
    int interfacecells = 0;
    int maxilevel = 0;
    foreach(reduction(+:cellnumber) reduction(+:interfacecells) reduction(max:maxilevel)) {
      cellnumber++;
      if (interfacial (point, f)) {
        interfacecells++;
      }
    if (point.level > maxilevel)
      maxilevel = point.level;
    }
    fprintf(ferr,"cellnumber %d\n",cellnumber);  
    fprintf(ferr,"interfacecells %d\n",interfacecells);  
    fprintf(ferr,"maxilevel %d\n",maxilevel);  

    if (inertial_particles == NULL){
      Particles p1,p2,p3;
      p1  = new_inertial_particles (0);
      p2  = new_inertial_particles (0);
      p3  = new_inertial_particles (0);
    }

  }
  else if(restDNS && !restLAG) {

    fprintf(ferr,"restart loaded, only DNS!!\n");  

    int cellnumber = 0;
    int interfacecells = 0;
    int maxilevel = 0;
    foreach(reduction(+:cellnumber) reduction(+:interfacecells) reduction(max:maxilevel)) {
      cellnumber++;
      if (interfacial (point, f)) {
        interfacecells++;
      }
    if (point.level > maxilevel)
      maxilevel = point.level;
    }
    fprintf(ferr,"cellnumber %d\n",cellnumber);  
    fprintf(ferr,"interfacecells %d\n",interfacecells);  
    fprintf(ferr,"maxilevel %d\n",maxilevel);  

    if (inertial_particles == NULL){
      Particles p1,p2,p3;
      p1  = new_inertial_particles (0);
      p2  = new_inertial_particles (0);
      p3  = new_inertial_particles (0);
    }
    partexist = 0;
  }
  else if(restDNS && restLAG) {

    fprintf(ferr,"restart loaded, DNS and Lagrange!!\n");  
/** realloc the list of Particles and asign a value to the first list. 
This cannot be done in prestore in particle.h, as inertial_particle is not defined there. */
//TODO: restLAG is the number of lists. For restarts with the old code (without 3. list of droplets) set manually 3,
// in future set restLAG.

    inertial_particles = realloc (inertial_particles, 3*sizeof(Particles));
    inertial_particles[0] = 0;
    inertial_particles[1] = 1;
    inertial_particles[2] = 2;
    pn[1] = 0;
    pna[1] = 0;
    pn[2] = 0;
    pna[2] = 0;
    pn[3] = terminate_int;
    pl = realloc (pl, 3*sizeof(particles));
    particle * ptlist1 = malloc (0.);
    pl[1] = ptlist1;
    particle * ptlist2 = malloc (0.);
    pl[2] = ptlist2;

    if (pn[0] > 0) partexist = 1;
#if _MPI     
    MPI_Allreduce (MPI_IN_PLACE, &partexist, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#endif
    
    int cellnumber = 0;
    int interfacecells = 0;
    int maxilevel = 0;
    foreach(reduction(+:cellnumber) reduction(+:interfacecells) reduction(max:maxilevel)) {
      cellnumber++;
      if (interfacial (point, f)) {
        interfacecells++;
      }
    if (point.level > maxilevel)
      maxilevel = point.level;
    }
    fprintf(ferr,"cellnumber %d\n",cellnumber);  
    fprintf(ferr,"interfacecells %d\n",interfacecells);  
    fprintf(ferr,"maxilevel %d\n",maxilevel);  

  }
  else {
    fprintf(ferr,"restart for particles without DNS restart is not possible!! \n");
    exit (1);
  }

  foreach()
    foreach_dimension()
      acoupling.x[] = 0.;

}

/* 
## Log some statistics on the solver.
 */
event logfile(i++)
{
  fprintf(ferr,
            "i t dt mgp.i mgpf.i mgu.i grid->tn perf.t perf.speed\n");
  fprintf(ferr, "%d %g %g %d %d %d %ld %g %g\n",
          i, t, dt, mgp.i, mgpf.i, mgu.i,
          grid->tn, perf.t, perf.speed);
}

/**
## Video output
 */
event video0_output(t += 1e-6)
{
  view (camera = "front",
        fov = 18.0,
        tx = -0.35,
        ty = -0.125,
        width = 362,
        height = 362,
        samples = 1);
  clear();
  box();
  draw_vof("f");
  draw_vof("f", filled = 1, fc = {0.7, 0.7, 0.7});

  if (partexist == 1) {
    scatter(inertial_particles[0], pc = {0.4, 0.6, 0.2});
  }

  //save("movies/mov.mp4");
  char name[80];
  sprintf (name, "movies/vof_and_particles-%g.png", t);
  save(name);

  fprintf(ferr, "\nOutput VOF + particles video!\n");
}

event video_output(t += 5e-5)
{
  scalar velx[];
  foreach()
    velx[] = u.x[];
 
  view (camera = "front",
      fov = 3.0,
      tx = -0.065,
      ty = -0.075,
      width = 362,
      height = 362,
      samples = 1);
  clear();
  draw_vof("f", filled = 1, color = "velx", map = blue_white_red, min = 0, max = 100);

  //save("movies/movZ.mp4");
  char name[80];
  sprintf (name, "movies/vofZ-%g.png", t);
  save(name);

  view (camera = "top",
      fov = 2.0,
      tx = -0.035,
      //ty = -0.075,
      width = 362,
      height = 362,
      samples = 1);
  clear();
  draw_vof("f", filled = 1, color = "velx", map = blue_white_red, min = 0, max = 100);
  //save("movies/movY.mp4");
  sprintf (name, "movies/vofY-%g.png", t);
  save(name);

  view (camera = "left",
      fov = 3.0,
      //tx = -0.035,
      ty = -0.07,
      width = 362,
      height = 362,
      samples = 1);
  clear();
  draw_vof("f", filled = 1, color = "velx", map = blue_white_red, min = 0, max = 100);
  //save("movies/movX.mp4");
  sprintf (name, "movies/vofX-%g.png", t);
  save(name);

  fprintf(ferr, "\nOutput VOF + ux!\n");
}

/**
## Save Snapshots
 */
event snapshot(t += dtlag; t <= T_END)
{
/** Save snapshots of the simulation at regular intervals to
    restart or to post-process with [bview](/src/bview). */
  char name[80];
  sprintf(name, "snapshot-%g", t);
  dump(name);

  //remove old snapshot
  sprintf(name, "snapshot-%g", t-0.5e-05);
  remove(name);

  if (inertial_particles != NULL){
    sprintf (name, "pdump-%g", t);
    pdump (name,inertial_particles,NULL,true,true);
  }

  int cellnumber = 0;
  int interfacecells = 0;
  int maxilevel = 0;
  foreach(reduction(+:cellnumber) reduction(+:interfacecells) reduction(max:maxilevel)) {
    cellnumber++;
    if (interfacial (point, f)) {
      interfacecells++;
    }
  if (point.level > maxilevel)
    maxilevel = point.level;
  }
  fprintf(ferr,"cellnumber %d\n",cellnumber);  
  fprintf(ferr,"interfacecells %d\n",interfacecells);  
  fprintf(ferr,"maxilevel %d\n",maxilevel);

  fprintf(ferr, "snapshot dumped!\n");
}

/**
## Paraview output to vtk for post-processing
 */
event vtk_output(t += 5e-5; t <= T_END)
{
  scalar m[];
  foreach()
    m[] = f[] > fcut;
  int n = tag (m);

  int proc_has_f[npe()];
  for (int j = 0; j < npe(); j++)
    proc_has_f[j] = 0;

  int proc_has_jet[npe()];
  for (int j = 0; j < npe(); j++)
    proc_has_jet[j] = 0;

  foreach_leaf()
    if (m[] > 0) {
      proc_has_f[pid()] = 1; // all the processors with f>fcut
      if (m[] == 1) proc_has_jet[pid()] = 1; // all the processors with cells describing the coherent jet
    }
#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, proc_has_f, npe(), MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, proc_has_jet, npe(), MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#endif

  char vtkname[80] = "allfields";
  char vtkdir[80] = "results/";
  scalar *list = {f, p, rho, m};
  vector *vlist = {u};
  save_data(list, vlist, i, t, vtkname, vtkdir);

  if (pid() ==0 ) {
    FILE * fvoffield;
    FILE * fprocs;
    FILE * fprocs0;
    char filename[80];

    sprintf(filename, "voffield-%03d.pvtu", i);
    fvoffield = fopen(filename, "w");
    sprintf(filename, "procs-%03d.dat", i);
    fprocs = fopen(filename, "w");
    sprintf(filename, "procs0-%03d.dat", i);
    fprocs0 = fopen(filename, "w");

    fputs ("<?xml version=\"1.0\"?>\n"
    "<VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fvoffield);
    fputs ("\t <PUnstructuredGrid GhostLevel=\"0\">\n", fvoffield);
    fputs ("\t\t\t <PCellData Scalars=\"scalars\">\n", fvoffield);
    for (scalar s in list) {
      fprintf (fvoffield,"\t\t\t\t <PDataArray type=\"Float64\" Name=\"%s\" format=\"appended\">\n", s.name);
      fputs ("\t\t\t\t </PDataArray>\n", fvoffield);
    }
    for (vector v in vlist) {
      fprintf (fvoffield,"\t\t\t\t <PDataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"Vect-%s\" format=\"appended\">\n", v.x.name);
      fputs ("\t\t\t\t </PDataArray>\n", fvoffield);
    }
    fputs ("\t\t\t </PCellData>\n", fvoffield);
    fputs ("\t\t\t <PPoints>\n", fvoffield);
    fputs ("\t\t\t\t <PDataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n", fvoffield);
    fputs ("\t\t\t\t </PDataArray>\n", fvoffield);
    fputs ("\t\t\t </PPoints>\n", fvoffield);

    for (int j = 0; j < npe(); j++) {
      if (proc_has_f[j] == 1) {
        // fill in .pvtu file only for processors including jet vof values
        fprintf (fvoffield, "<Piece Source=\"allfields-%03d_n%03d.vtu\"/>\n",i,j);
        // list of processors, comma separated
        fprintf (fprocs, "%03d,",j);
      }
      if (proc_has_jet[j] == 1) {
        fprintf (fprocs0, "%03d,",j);
      }
    }

    fputs ("\t </PUnstructuredGrid>\n", fvoffield);
    fputs ("</VTKFile>\n", fvoffield);

    fclose(fvoffield);
    fclose(fprocs);
    fclose(fprocs0);
  }

  fprintf(ferr, "vtk-output written!\n");
}
      
/**
## Output the list of resolved droplets identified in DNS with the tag function in event droplet_to_particle.
 */
event droplets_out (t = dtlag; t += dtlag)
{
#if dimension == 2
  for (int j=0; j<pn[2]; j++)
    fprintf (fout, "%d %g %d %g %g %g\n", i, t,
             j, pl[2][j].vol, pl[2][j].x, pl[2][j].y, pl[2][j].u.x, pl[2][j].u.y);
#else
/**
For the statistics, we output the volume, position and velocity of each droplet to
standard output.
 */
    if (pn[2] > 0 && pid() == 0)
      for (int j=0; j<pn[2]; j++) {
        printf ("in event droplets: %d, %g, %d, %g, %g, %g, %g, %g, %g, %g\n", i, t,
          j, pl[2][j].vol, pl[2][j].x, pl[2][j].y, pl[2][j].z, pl[2][j].u.x, pl[2][j].u.y, pl[2][j].u.z);
      }
#if _MPI
    if (pid()==0) printf("pid %d, pn[2]: %lu\n", pid(), pn[2]);
#else
    fprintf(ferr,"pn[2]: %lu\n", pn[2]);
#endif

  if (pid() == 0) {
      change_plist_size (inertial_particles[2], -pn[2]);
  }

#endif
  
  fflush (fout);

  fprintf(ferr, "Droplets counted and written to file! (with list)\n");
}

/**
## Output the list of point-droplets in Lagrange formalism.
 */
event particles_outCSV (t = dtlag; t += dtlag)
{
  coord urel;
  if (partexist == 1) {
#if dimension == 2
    for (int j=0; j<pn[0]; j++)
      fprintf (fout, "%d %g %d %g %g %g\n", i, t,
              j, pl[0][j].vol, pl[0][j].x, pl[0][j].y, pl[0][j].u.x, pl[0][j].u.y);
#else

    int numberparticles = pn[0];
#if _MPI
    MPI_Allreduce (MPI_IN_PLACE, &numberparticles, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    int newin = update_mpi_all (inertial_particles[0]);
#endif
    if (pn[0] > 0 && pid() == 0) {

      static FILE * foutpart;
      char name[80];
      sprintf (name, "particle/particles-%d.csv", i);
      foutpart = fopen (name, "a");
      fprintf (foutpart, "j, vol, x, y, z, ux, uy, uz, We\n");

      for (int j=0; j<pn[0]; j++) {
        foreach_dimension()
          urel.x = pl[0][j].uf.x - pl[0][j].u.x;
	      fprintf (foutpart, "%d, %g, %g, %g, %g, %g, %g, %g, %g\n",
          j, pl[0][j].vol, pl[0][j].x, pl[0][j].y, pl[0][j].z, pl[0][j].u.x, pl[0][j].u.y, pl[0][j].u.z, rho2*pl[0][j].dp*sq(normcoord(urel))/f.sigma);
        printf ("in event particles: %d, %g, %g, %g, %g, %g, %g, %g, %g, %g\n", i, t,
          pl[0][j].vol, pl[0][j].x, pl[0][j].y, pl[0][j].z, pl[0][j].u.x, pl[0][j].u.y, pl[0][j].u.z, rho2*pl[0][j].dp*sq(normcoord(urel))/f.sigma);
      }
      fclose(foutpart);
    }
#if _MPI
    assert(numberparticles == pn[0]);
    if (pid()==0) printf("pid %d, pn[0]: %lu\n", pid(), pn[0]);
#else
    fprintf(ferr,"pn[0]: %lu\n", pn[0]);
#endif

#if _MPI
      change_plist_size (inertial_particles[0], -newin);
#endif

#endif
  }
  fflush (fout);

  fprintf(ferr, "Point-droplets counted and written to file! (with list)\n");
}


event adapt(i++)
{
/* Adapt the mesh according to the error on the volume fraction field
    and the velocity. */
#if dimension > 2
  adapt_wavelet2({f, u.x, u.y, u.z}, (double[]){femax, uemax, uemax, uemax}, (int[]){maxlevel,10,10,10});
#else
  adapt_wavelet({f, u.x, u.y}, (double[]){femax, uemax, uemax}, maxlevel);
#endif

/* Set a coarse mesh close to the outlet to avoid backflow */
  unrefine(x > (SIZE_FACTOR - SIZE_UNREFINE - 2.5) * DIAMETER);

  fprintf(ferr, "Mesh adaptation finished!\n");
}

event adapt_lag(t = dtlag; t += dtlag) {	
  adapt_wavelet2({f, u.x, u.y, u.z}, (double[]){femax, uemax, uemax, uemax}, (int[]){maxlevel,10,10,10});
}

/**
## Final message
 */
event the_end(t = T_END)
{
  fprintf(ferr, "\n-----------------------");
  fprintf(ferr, "\n SIMULATION COMPLETED! ");
  fprintf(ferr, "\n-----------------------\n");
}