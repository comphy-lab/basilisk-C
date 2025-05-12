/** 
 <h1> Elastoviscoplastic  Droplet Spreading </h1>
 This file contains the source code used to perform the simulations for the spreading of elastoviscoplastic droplets. 
 The results were published in the paper "[Elasto-viscoplastic spreading: From plastocapillarity to elastocapillarity](https://doi.org/10.1103/PhysRevResearch.6.013226)".
 <img src="https://lh3.googleusercontent.com/d/1-H47znN32JJefnEuvtxaIFzalXd6Mr2u" alt="EVP Droplet Spreading Sketch" width=100%/>
*/

#include <time.h>
#include "axi.h"
#include "navier-stokes/centered.h"
#define FILTERED
#include "two-phaseEVP.h"
#include "navier-stokes/conserving.h"
#include "log-conform.h"
#include "tension.h"
#include "distance.h"
#include "fluidlab_pack.h"

/** 
Nondimensional numbers that define the Saramito Droplet Spreading problem.
*/
double Bo; // Bond number
double Oh; // Ohnesorge number
double Wi; // Weissenberg number
double beta; // Viscosity ratio
double J; // Plastocappilary number
double regularization; // Yield-stress regularization parameter
double h_inf; // Thin-film thickness
double phases_ratio; // Ratio between density and viscosity of droplet to outer fluid

/** 
 Mesh refinement Parameters
*/
int mesh_level; // Max mesh refinement level
int mesh_level_outside = 6;
int mesh_level_inside = 9;
int mesh_level_film = 9;

/**
 We define scalar variables used to visualize certain properties.
*/
scalar lambdav[], mupv[]; // polymer relaxation time and polymeric viscosity (nondimensional: (Wi) and (1-beta))
scalar yielded[]; // Tracking yielded/plugged regions
scalar log_norm_tau_dev_field[], norm_tau_dev_field[]; // Storing the norm of tau_p_dev
scalar D2[]; // Storing the norm of strain rate
scalar tau_s_norm[]; // Storing the norm of tau_s
scalar tau_p_norm[]; // Storing the norm of tau_p
scalar tau_norm[]; // Storing the norm of tau = tau_s + tau_p
scalar flow_type[]; // Storing the flow type parameter (shear=0 extension=1 rotation=-1)
scalar distance_field[]; 

/**
 Boundary conditions:
 
 1. Left: No-slip
 2. Right and top: Open boundary 
 3. Bottom: Axi-symmetry
 */
u.t[left] = dirichlet(0);
f[left] = 0.0;

u.n[right] = neumann(0);
p[right] = dirichlet(0);
pf[right] = dirichlet(0);

u.n[top] = neumann(0);
p[top] = dirichlet(0);
pf[top] = dirichlet(0);

/**
 Reading command line arguments and creating the domain.
 Here is an example command you can use to run:
 
 ./spreading 0.0 0.1 0.0 0.0 1.0 100.0 1e-06 1e-07 1e-04 0.0175 10

 For a full list of all simulation parameters ran in our paper, see the
 spreadsheet in the parameters sub-folder: [parameters](parameters/).
*/
int main (int argc, char * argv[]) {

  if( argc<12 ) {
    printf("Not enough input arguments..\n");
    exit(0);
  }

  DT = 1e-05;
  Bo = atof(argv[1]);
  Oh = atof(argv[2]);
  J = atof(argv[3]);
  Wi = atof(argv[4]);
  beta = atof(argv[5]);
  phases_ratio = atof(argv[6]);
  TOLERANCE = atof(argv[7]);
  regularization = atof(argv[8]);
  DT = atof(argv[9]);
  h_inf = atof(argv[10]);
  mesh_level = atof(argv[11]);

  size(5.0);
  origin(0.0, 0.0);
  init_grid(256);

  /**
   Below we define:
   
   1. Solvent viscosity of both fluids (mu1 and mu2).
   
   2. Density of both fluids (rho1 and rho2).
   
   3. Surface tension coefficient (f.sigma).
   
   4. Variables where the polymeric viscosity and relaxation time will go (mupv and lambdav).
   
   Fluid 1 is EVP droplet and Fluid 2 is the Newtonian air around.
   */
  rho1 = 1.0;
  mu1 = beta*Oh;
  rho2 = rho1/phases_ratio;
  mu2 = mu1/phases_ratio;

  // Setting the surface tension coefficient
  f.sigma = 1.0; 

  mup = mupv;
  lambda = lambdav;

  /**
   We create a folder where all the output of this simulation will go and then start the simulation.
   */
  OpenSimulationFolder("spreading_film%g_reg%.1e_mesh%d_J%g_Oh%g_Bo%g_Wi%g_beta%g", h_inf, regularization, mesh_level, J, Oh, Bo, Wi, beta);
  run();
  CloseSimulationFolder();
}

/**
 This event sets the properties for the Saramito fluid. 
 
 We treat the fluid in two different ways, depending on the value of $Wi$:
 
  1. If $Wi = 0$: We switch into Bingham mode, which means:
    - The log-conform solver is ignored by setting $\lambda_p=1$ and $\mu_p=0$.
    - The solvent viscosity is set to the Bingham-Papanastasiou Generalized Newtonian viscosity. That is:
    $$
    \mu_s = Oh + \frac{\mathcal{J}}{2||\mathbf{D}||}\left( 1 - e^{-||\bm{D}||/\epsilon} \right).
    $$

  2. If $Wi \neq 0$: We go into Saramito mode, which means we want to solve the following equations:
  $$\frac{D\bm{u}}{Dt} = - \nabla p + \nabla \cdot \left( 2 \ \beta \ Oh \bm{D}\right) + \nabla \cdot \bm{\tau}^p + Bo \cdot \bm{g},$$
  $$ 
  De \stackrel{\raisebox{.3ex}{ \hskip -.1cm$\scriptscriptstyle{\bigtriangledown}$}}{\bm{\tau}^p} + \max{\left( 0, 1 - \frac{\mathcal{J}}{||\bm{\tau}^p||} \right)}\bm{\tau}^p = 2 \ (1 - \beta) Oh \ \bm{D}.
  $$
    - Therefore, the solvent viscosity is a simple newtonian constant viscosity:
    $$\mu_s = \beta \ Oh.$$
    - To obtain the constitutive equation above, the log-conform variables $\lambda_p$ and $\mu_p$ are set as
    $$
    \lambda_p = De\  / \ \max{\left( \epsilon, 1 - \frac{\mathcal{J}}{||\bm{\tau}^p||} \right)},
    $$
    $$
    \mu_p = (1 - \beta) \ Oh \  / \ \max{\left( \epsilon, 1 - \frac{\mathcal{J}}{||\bm{\tau}^p||} \right)}.
    $$

*/
event properties (i++) {
  double one_third = 1.0/3.0;

  // If Wi = 0 exactly, I will just use a pure Bingham model (generalized newtonian) with regularization. (is this cheating?)
  // This section was partially taken from Vatsal's sandbox (I changed the regularization method).
  if( Wi==0.0 ) {
    foreach_face(x) {
      double ff = (sf[] + sf[-1])/2.;
      alphav.x[] = fm.x[]/rho(ff);
      face vector muv = mu;
      double D11 = 0.5*(u.y[0,1] - u.y[0,-1] + u.y[-1,1] - u.y[-1,-1])/(2.0*Delta);
      double D22 = (u.y[] + u.y[-1, 0])/(2.0*max(y, 1e-20));
      double D33 = (u.x[] - u.x[-1,0])/Delta;
      double D13 = 0.5*( (u.y[] - u.y[-1, 0])/Delta + 0.5*(u.x[0,1] - u.x[0,-1] + u.x[-1,1] - u.x[-1,-1])/(2.0*Delta) );
      double norm_D = sqrt( sq(D11) + sq(D22) + sq(D33) + 2.0*sq(D13) );
      double apparent_visc = Oh + J/(2.0*norm_D + 1e-10)*(1.0 - exp(-norm_D/regularization));
      // double apparent_visc = Oh + J/(sqrt(2.0)*norm_D + 1e-10)*(1.0 - exp(-norm_D/regularization));
      muv.x[] = fm.x[]*( clamp(ff, 0.0, 1.0)*(apparent_visc - mu2) + mu2 );
    }

    foreach_face(y) {
      double ff = (sf[0,0] + sf[0,-1])/2.;
      alphav.y[] = fm.y[]/rho(ff);
      face vector muv = mu;
      double D11 = (u.y[0,0] - u.y[0,-1])/Delta;
      double D22 = (u.y[0,0] + u.y[0,-1])/(2*max(y, 1e-20));
      double D33 = 0.5*( (u.x[1,0] - u.x[-1,0] + u.x[1,-1] - u.x[-1,-1])/(2.*Delta) );
      double D13 = 0.5*( (u.x[0,0] - u.x[0,-1])/Delta + 0.5*( (u.y[1,0] - u.y[-1,0] + u.y[1,-1] - u.y[-1,-1])/(2.*Delta) ) );
      double norm_D = sqrt( sq(D11) + sq(D22) + sq(D33) + 2.0*sq(D13) );
      double apparent_visc = Oh + J/(2.0*norm_D + 1e-10)*(1.0 - exp(-norm_D/regularization));
      // double apparent_visc = Oh + J/(sqrt(2.0)*norm_D + 1e-10)*(1.0 - exp(-norm_D/regularization));
      muv.y[] = fm.y[]*( clamp(ff, 0.0, 1.0)*(apparent_visc - mu2) + mu2 );
    }

    // Regarding the oldroyd-B solver
    // We just want to force tau_p = 0, so I set mupv=0 and lambdav=1 (lambda is irrelevant if mupv=0)
    foreach() {
      mupv[] = 0.0;
      lambdav[] = 1.0; // This value doesnt really matter because mupv=0
    }
  }
  else { // If Wi!=0, we really have to solve for the polymeric part of the stress tensor
    
    // This loop will set lambda and mup (relaxation time and polymeric viscosity)
    foreach() {
      /// === Norm of the polymeric stress tensor
      double trace_tau = tau_p.x.x[] + tau_p.y.y[] + tau_qq[];
      double tau_dev_xx = tau_p.x.x[] - one_third*trace_tau; 
      double tau_dev_yy = tau_p.y.y[] - one_third*trace_tau; 
      double tau_dev_qq = tau_qq[] - one_third*trace_tau;
      double tau_dev_xy = tau_p.x.y[];
      double norm_tau_dev = sqrt( sq(tau_dev_xx) + sq(tau_dev_yy) + sq(tau_dev_qq) + 2.0*sq(tau_dev_xy) );
      norm_tau_dev_field[] = norm_tau_dev;
      log_norm_tau_dev_field[] = (norm_tau_dev > 0.0) ? log10(norm_tau_dev) : -10.0;
      double yield_term = 1.0;
      yielded[] = 1.0;
      if( J>0.0 ) {
        yield_term = (norm_tau_dev < (J + regularization)) ? regularization : (norm_tau_dev - J)/norm_tau_dev;
        yielded[] = yield_term >= 1e-03;
      }

      // Only the droplet is EVP, which is why we multiply the coefficients by the volume fraction f[]
      mupv[] = (1.0 - beta)*Oh*clamp(f[],0,1)/(yield_term);
      lambdav[] = (Wi/yield_term)*clamp(f[],0,1);
    }

    // This loop will set the solvent viscosity (copy-paste from two-phase.h)
    foreach_face() {
      double ff = (sf[] + sf[-1])/2.;
      alphav.x[] = fm.x[]/rho(ff);
      if (mu1 || mu2) {
        face vector muv = mu;
        muv.x[] = fm.x[]*mu(ff);
      }
    }
  }

  foreach()
    rhov[] = cm[]*rho(sf[]);
 
  sf.prolongation = fraction_refine;
  sf.dirty = true; // boundary conditions need to be updated

  boundary(all);
}

/**
 Function that defines the initial droplet shape.
 */
double initial_shape(double x, double y) {

  double R0 = 1.0;

  // film height
  double h_film = h_inf;

  // Droplet from mazis paper
  double r = y;
  double h = h_film + R0*max(0.0, 1.0 - (r/R0)*(r/R0));
  return h - x;
}

/**
 This functions calculates a pretty bad distance field of each point to the droplet surface.
 It is used for re-meshing in the adapt event.
 */
void CalculateDistanceFunction(scalar d)
{
  // === Alocating memory for the local vertices
  double *local_vertices_x = (double *)malloc( 5000*sizeof(double) );
  double *local_vertices_y = (double *)malloc( 5000*sizeof(double) );

  // === Calculating all polygons again (bad) and saving the vertices to the local arrays
  int index_vertex = 0;
  foreach(serial) {
    if( cfilter (point, f, 1e-05) ) {
      coord n = interface_normal (point, f);
      double alpha = plane_alpha (f[], n);
      coord v[2];
      int m = facets (n, alpha, v);

      local_vertices_x[index_vertex] = local_vertices_y[index_vertex] = 0.0;
      for( int i=0; i<m; i++ ) {
        local_vertices_x[index_vertex] += x + v[i].x*Delta;
        local_vertices_y[index_vertex] += y + v[i].y*Delta;
      }
      local_vertices_x[index_vertex] /= (double)m;
      local_vertices_y[index_vertex] /= (double)m;

      // Ignore vertices too close to the solid surface
      if( local_vertices_x[index_vertex]>0.001 )
        index_vertex++;
    }
  }

  // === Counting how many polys and vertices we have in total between all processes
  int local_count_vertices = index_vertex;
  int total_count_vertices = 0;
  #if _MPI
    MPI_Allreduce(&local_count_vertices, &total_count_vertices, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  #else
    total_count_vertices = local_count_vertices;
  #endif

  // === Sending all the vertices to the global vertex array in process id=0
  double *vertices_x = NULL, *vertices_y = NULL;
  #if _MPI
    vertices_x = (double *)malloc( total_count_vertices*sizeof(double) );
    vertices_y = (double *)malloc( total_count_vertices*sizeof(double) );

    MPI_Gather_Uneven(local_vertices_x, local_count_vertices, MPI_DOUBLE, vertices_x, 0);
    MPI_Gather_Uneven(local_vertices_y, local_count_vertices, MPI_DOUBLE, vertices_y, 0);

    // Sending the total array to other processors as well
    if( pid()==0 ) {
      for(int i=1; i<npe(); i++) {
        MPI_Send(vertices_x, total_count_vertices, MPI_DOUBLE, i, 3000 + i, MPI_COMM_WORLD);
        MPI_Send(vertices_y, total_count_vertices, MPI_DOUBLE, i, 4000 + i, MPI_COMM_WORLD);
      }
    }
    else {
      MPI_Recv(vertices_x, total_count_vertices, MPI_DOUBLE, 0, 3000 + pid(), MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vertices_y, total_count_vertices, MPI_DOUBLE, 0, 4000 + pid(), MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // === Releasing local memory
    free(local_vertices_x);
    free(local_vertices_y);
  #else
    vertices_x = local_vertices_x;
    vertices_y = local_vertices_y;
  #endif

  foreach() {
    d[] = 1e+10;
    for( int i=0; i<total_count_vertices; i++ ) {
      double local_distance = (x - vertices_x[i])*(x - vertices_x[i]) + (y - vertices_y[i])*(y - vertices_y[i]);
      local_distance = sqrt(local_distance);
      if( local_distance<d[] )
        d[] = local_distance;
    }
  }


}

/**
 This function calculates the droplet radius and height for visualization later.
*/
void CalculateDropletRadius(double *Radius, double *Height)
{
  double droplet_radius1 = 1e+10, droplet_radius2 = 1e+10, droplet_radius3 = 1e+10;
  double droplet_height = 0.0;

  // Calculating all the interface line segments
  double fmin = 1e-5; // do not reconstruct fragments smaller than this
  foreach(reduction(min:droplet_radius1) reduction(min:droplet_radius2) reduction(min:droplet_radius3) reduction(max:droplet_height)) {
    if( cfilter (point, f, fmin) ) {
      coord n = interface_normal (point, f);
      double alpha = plane_alpha (f[], n);
      coord v[2];
      facets (n, alpha, v);

      // The coordinates of the two points that form this line segment
      double x0 = x + v[0].x*Delta;
      double y0 = y + v[0].y*Delta;
      double x1 = x + v[1].x*Delta;
      double y1 = y + v[1].y*Delta;

      if( x0>droplet_height )
        droplet_height = x0;
      if( x1>droplet_height )
        droplet_height = x1;

      // Checking if this line segment crosses the line x = x_radius
      double x_radius;

      x_radius = 1.0*h_inf;
      if( (x0<=x_radius && x1>=x_radius) || (x1<=x_radius && x0>=x_radius) ) {
        // Calculating the intersection point (x_radius, y_intersec)
        double lambda = (x_radius - x0)/(x1 - x0 + 1e-10);
        double y_intersec = y0 + lambda*(y1 - y0);
        if( y_intersec<droplet_radius1 )
          droplet_radius1 = y_intersec;
      }
      
      x_radius = 1.1*h_inf;
      if( (x0<=x_radius && x1>=x_radius) || (x1<=x_radius && x0>=x_radius) ) {
        // Calculating the intersection point (x_radius, y_intersec)
        double lambda = (x_radius - x0)/(x1 - x0 + 1e-10);
        double y_intersec = y0 + lambda*(y1 - y0);
        if( y_intersec<droplet_radius2 )
          droplet_radius2 = y_intersec;
      }

      x_radius = 1.2*h_inf;
      if( (x0<=x_radius && x1>=x_radius) || (x1<=x_radius && x0>=x_radius) ) {
        // Calculating the intersection point (x_radius, y_intersec)
        double lambda = (x_radius - x0)/(x1 - x0 + 1e-10);
        double y_intersec = y0 + lambda*(y1 - y0);
        if( y_intersec<droplet_radius3 )
          droplet_radius3 = y_intersec;
      }
    }
  }
  
  Radius[0] = droplet_radius1;
  Radius[1] = droplet_radius2;
  Radius[2] = droplet_radius3;
  *Height = droplet_height;
}

/**
 This event calculates run time of each timestep for bookkeeping.
 I think there is a proper Basilisk function for this. I did not know that yet when I made this...
*/
clock_t start_time = -1.0, end_time = -1.0;
double cpu_time_used = 0.0;
event step_cputime(i++) {
  if( pid() )
    return 0;

  end_time = clock();
  cpu_time_used  = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;
  start_time = end_time;
}

/**
 This event prints some logging information into a file for visualization later.
 It prints: current time, kinetic energy, droplet radius, droplet height
 */
double droplet_radius = 1.0, droplet_height = 1.0;
event logfile (t+=0.001; t<=200.0) {
  // === Calculating the kinetic energy inside the droplet, strain rate tensor, vorticity tensor and flow type parameter
  double kinetic_energy = 0.0;
  foreach(reduction(+:kinetic_energy)){

    // Kinetic energy
    kinetic_energy += f[]*(2.0*pi*y)*0.5*(sq(u.x[]) + sq(u.y[]))*sq(Delta);

    // Calculating velocity derivatives
    double dudx = (u.x[1, 0] - u.x[-1, 0])/(2.0*Delta);
    double dudy = (u.x[0, 1] - u.x[0, -1])/(2.0*Delta);
    double dvdx = (u.y[1, 0] - u.y[-1, 0])/(2.0*Delta);
    double dvdy = (u.y[0, 1] - u.y[0, -1])/(2.0*Delta);
    
    // Calculating the strain rate tensor and its norm
    double D11 = dudx;
    double D22 = dvdy;
    double D33 = u.y[]/max(y, 1e-12);
    double D12 = 0.5*( dudy + dvdx );
    double norm_strain = sqrt( sq(D11) + sq(D22) + sq(D33) + 2.0*sq(D12) );
    D2[] = (norm_strain > 0.0) ? log10(norm_strain) : -10.0;

    // Calculating the vorticity tensor and its norm
    double O11 = 0.0, O22 = 0.0, O33 = 0.0;
    double O12 = 0.5*( dvdx - dudy );
    double O21 = 0.5*( dudy - dvdx );
    double norm_vorticity_rz = sqrt( sq(O11) + sq(O22) + sq(O12) + sq(O21) );
    double norm_strain_rz = sqrt( sq(D11) + sq(D22) + 2.0*sq(D12) );
    flow_type[] = (norm_strain_rz - norm_vorticity_rz)/( norm_strain_rz + norm_vorticity_rz + 1e-12 );

    // Calculating norm of tau_s, tau_p and tau = tau_s + tau_p
    double eta_s = (Wi==0.0) ? Oh : beta*Oh; // Solvent viscosity
    tau_s_norm[] = 2.0*eta_s*norm_strain;
    tau_p_norm[] = sqrt( sq(tau_p.x.x[]) + sq(tau_p.y.y[]) + 2.0*sq(tau_p.x.y[]) + sq(tau_qq[]) );
    tau_norm[] = tau_s_norm[] + tau_p_norm[];
  }

  // Updating the radius and the height of the droplet
  double temp_droplet_radius[3];
  CalculateDropletRadius(temp_droplet_radius, &droplet_height);
  droplet_radius = temp_droplet_radius[0];

  // Updating the distance function used for re-meshing
  CalculateDistanceFunction(distance_field);

  // Printing into a log file
  PrintLog("%d %lf %e %lf %e %lf %lf %lf %lf\n", i, t, dt, cpu_time_used, kinetic_energy, droplet_radius, droplet_height, temp_droplet_radius[1], temp_droplet_radius[2]);
}

/**
 This event creates VTK files with the relevant scalar fields and droplet interface 
 with a certain time frequency (vtk_frequency = 0.05).
 */
event print_solution(t+=0.05) {
  scalar *list = {u.x, u.y, yielded, log_norm_tau_dev_field, D2, flow_type};
  const char *list_names[] = {"ux", "uy", "yielded", "tau_dev", "log_D2", "flow_type"};
  PrintMeshVTK_Binary_Double(i, t, list, list_names);

  PrintInterfaceVTK(i, t);
}

/**
 This event initializes the droplet shape and mesh.
 */
event init (t = 0) {  

  // Adapting mesh based on discretization error of fields
  adapt_wavelet ({f, u.x, u.y}, (double[]){1e-2, 1e-2, 1e-2}, maxlevel = mesh_level, minlevel = mesh_level_outside);

  // Adding some additional mesh refinement near the wall for the thin-film
  mesh_level_film = max(mesh_level, 10);
  refine( (x<=2.0*h_inf) && (y>=droplet_radius)  && (level<mesh_level_film) );

  // Filling in the volume fraction field
  fraction (f, initial_shape(x, y));

  // Refining mesh a bit more inside the droplet (where f>0)
  refine( (f[]>=0.02) && (level<mesh_level) );

  // Refining mesh even more around the interface of the droplet
  CalculateDistanceFunction(distance_field);
  refine( (distance_field[]<0.05) && (level<mesh_level) );

  // Making sure boundary conditions are updated
  boundary(all);
}

/** 
  This event adds gravity acceleration.
 */
event acceleration (i++) {
  face vector av = a;
  foreach_face(x)
    av.x[] -= Bo;
}

/**
 This event perform the mesh adaptation every timestep.
 */
event adapt (i++) {
  // Adapting based on discretization error of fields
  adapt_wavelet ({f, u.x, u.y, norm_tau_dev_field}, (double[]){1e-3, 1e-3, 1e-3, 1e-5}, maxlevel = mesh_level, minlevel = mesh_level_outside);

  // Adding some additional refinement near the wall for the thin-film
  mesh_level_film = max(mesh_level, 10);
  refine( (x<=2.0*h_inf) && (y>=droplet_radius)  && (level<mesh_level_film) );

  // Refining a bit more inside the droplet (where f>0)
  refine( (f[]>=0.02) && (level<mesh_level) );

  // Refining mesh even more around the interface of the droplet
  refine( (distance_field[]<0.05) && (level<mesh_level) );

  boundary(all);
}