/**
# Complex pipe geometry

Here is provided a process for constructing a complex geometry built with the help of single scalars.
The convergence of the set-up with embeded boundaries is tested in a cross-flow, as well as during the restart.
Special care is taken to prevent assertion errors during the restart.
The code is also tested to configure a pipe with a width of one cell.

![Pipe above a flat bottom in a cross-flow
and the $\lambda_2$ vortex-detection isosurfaces](pipe_geometry/pipe_normal.mp4)
 
*/

//The order is important
#include "grid/octree.h" //1
#include "embed.h" //2
#include "navier-stokes/centered.h" //3
#include "fractions.h"
#include "lambda2.h"
#include "view.h"

//Borrowed from acastillo's sandbox with few additions
#include "output_vtu_foreach.h"

//Simulation specifics
int maxlevel;
int minlevel = 4;
double lenght_domain = 6.;
double cross_flow_velocity = 1.;
double timestep_ini =10;
int record_step = 5;
int dump_step = 20;

//Variables
face vector av[], muc[];
scalar pipe_f0[], pipe_wall[], plane[];
scalar l2[];
vector omega[];
int i_for_restart, end_simu;
double eps = 1e-8;
bool onecell = 0, bool_restore = 0;
char file_restart[350];

/**
## Boundary conditions
*/
//Embed BC
u.t[embed] = dirichlet (0.); //no-flow in the solid
u.n[embed] = pipe_wall[] > eps ? dirichlet (0.) : neumann (0.); //no-flow in the solid, free-slip otherwise
u.r[embed] = pipe_wall[] > eps ? dirichlet (0.) : neumann (0.); //no-flow in the solid, free-slip otherwise

//Inflow BC
u.n[left] = dirichlet (cs[]*cross_flow_velocity);
u.t[left] = neumann (0.);
u.r[left] = neumann (0.);

//Outflow BC
u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
u.r[right] = neumann (0.);
uf.n[right] = neumann (0.); //to help convergence

//Top BC
u.n[top] = neumann (0.);
u.t[top] = neumann (0.);
u.r[top] = neumann (0.);

p[top] = dirichlet (0.);    
pf[top] = dirichlet (0.);
p[embed] = neumann (0.);


/**
## Geometry function

Here we define the geometry. The plane and the inner tube are defined first,
which helps to define the "wall" of the pipe, then the scalar cs is defined.
*/
void define_geometry (scalar pipe_f0, scalar pipe_wall, scalar plane,\
  double xcenter, double zcenter, double pipe_height, double radius_exit, double H_plane)
{
  scalar * list_scalar_geometry = {pipe_f0, pipe_wall, plane};
  
  int iteration = 0;
  
  /** Definition of the pipe and the bottom plane */
  do{
    fraction (pipe_f0, y <= pipe_height? (sq(radius_exit) - sq(z-zcenter) - sq(x-xcenter)) : 0);
    fraction (plane, -y + (Y0 + H_plane));
  }
  while (adapt_wavelet ({pipe_f0, plane}, (double []) {1e-15, 1e-15}, maxlevel = maxlevel, maxlevel).nf &&
   iteration++ <= 10);  
  
  /** The pipe is forced to one wherever it's defined and is copied to pipe_wall scalar.*/
  foreach() {
    pipe_wall[] = 0.;
    if (pipe_f0[] > eps) {
      pipe_f0[] = 1.;
      pipe_wall[] = 1.;
    }
    else {pipe_f0[] = 0.;}
  }
  boundary ({pipe_f0, pipe_wall});
  
  /**The contour of the inner tube is detected. */
  int i, j;
  int pipe_width = 2;
  foreach() {
    for (i = -pipe_width ; i <= pipe_width; i++) {
      for (j = -pipe_width ; j <= pipe_width; j++) {
        if (pipe_f0[] < eps) {
          if (pipe_f0[i,0,j] > eps) { //careful, pipe_f0 should be far from boundary (core dumped error)
            pipe_wall[] = 1.;
          }
        }
      }
    }
  }
  boundary ({pipe_wall});
  
  /**The inner tube is "excavated" to keep only the walls. */
  foreach()
    pipe_wall[] = (pipe_wall[] * (1 - pipe_f0[]));
  
  //For BC on levels  
  for (scalar m in list_scalar_geometry) {
    m.refine = m.prolongation = fraction_refine;
    restriction ({m});  
    boundary ({m});
  } 

  /** When restoring the solution, the phi field needs to be cleaned, but HUGE doesn't seem to work,
  so it's initialized with both a random value and HUGE. */
  vertex scalar phi[];
  foreach_vertex()
    phi[] =  15;
  boundary ({phi});  
  fractions (phi, cs, fs);
  fractions_cleanup (cs, fs);
  boundary ({cs, fs});       
  
  /**The vertex scalar field is built from the scalars defined before, to be used for the "fractions" function.*/
  foreach_vertex() {
  //For the case without a restore, phi is also initialized with HUGE
    phi[] =  HUGE;
    if (pipe_wall[] > eps)
      phi[] = 0;
    if (pipe_wall[0,0,1] > eps)
      phi[] = 0; //phi is a vertex scalar
  }
  boundary ({phi});
  
  foreach_vertex()
    if (plane[] > eps && pipe_f0[] < eps )
      phi[] = 0;
  boundary ({phi});
  
  fractions (phi, cs, fs);
  fractions_cleanup (cs, fs);
  cs.refine = cs.prolongation = fraction_refine;
  boundary ({cs, fs});
  restriction ({cs, fs});
  

  /**Somehow cs is shifted, the construction scalars are shifted and corrected accordingly.*/
  for (scalar m in list_scalar_geometry) {
    foreach()
      m[] = m[0,0,1];
    boundary ({m}); 
    foreach()
      if (m[] != m[0,1,0])
        m[] = 0.;
    boundary ({m}); 
  }
}

int main() {

  //Horizontal BC
  periodic (front);

  mu = muc;
  a = av;

  X0 = Y0 = Z0 = 0.;

  for (scalar s in {u, g})
    s.prolongation = refine_embed_linear; 

  sprintf(file_restart, "none");

/**
## Results
*/

/**
### Normal case
*/
  /**
  The first case that produces the video above.
  The perturbation induced by the pipe in the cross-flow can be observed with the $\lambda_2$ negative isosurfaces.
  This case can be used to output a flow from the pipe.
  */
  maxlevel = 8;
  onecell = 0;
  end_simu = 2000;
  i_for_restart = end_simu;
  run();
  
/**
### Restart case
*/
  /**
  The restart case is also tested, to verify that no assertion or convergence error occurs in MPI.
  Special care is taken below to avoid any malfunction even at high resolution.
  Performing a restart with a chain of "run()" prevents the grid from being properly reset,
  and some cells are triggered as a solid when they shouldn't be.
  Interestingly, this doesn't seem to be an issue for the solver convergence and Basilisk.
  
  ![Simulation from a restart](pipe_geometry/pipe_from_restore.mp4)

  ![Anomaly of solid reconstruction](pipe_geometry/pipe_from_restore_anomaly.png)
  */
  end_simu = end_simu + 1000;
  sprintf (file_restart, "dump-%d",i_for_restart);
  run();

/**
### One cell case
*/
  /**
  A particular case where the pipe hole reaches the cell size at maxlevel.
  A little trick is used for this configuration.
  
  ![Pipe of one cell width](pipe_geometry/pipe_one_cell.png)
  */
  sprintf (file_restart, "none");
  maxlevel = 7; //the resolution is kept at a suitable level for the cpu resources
  onecell = 1;
  end_simu = 20;//
  run();
}

/**
## Set-up
*/

event init (t = 0) {

  if (!restore (file = file_restart)) {

    size (lenght_domain);
    
    refine (level < minlevel);
    
    FILE * fp = fopen ("count_cell", "w");
    fprintf (fp, "n_iteration\tn_cell\tdt\n");
    fflush (fp);
    fclose (fp);
    
    
    /**For an exit on one cell only. The pipe in built at malevel+1 and then unrefined to maxlevel.*/
    if (onecell) {
      maxlevel += 1;

      /**The coordinates are given to fall right on the grid (no half-fluid cells).
      There might be some issue in that case where the "boundary_flux()"
      is not assigned to any cell ([User_forum](https://groups.google.com/g/basilisk-fr/c/SLnKb8U0BXs)).
      */
      double dx_maxlevel = L0/pow(2, maxlevel), dx_minlevel = L0/pow(2, minlevel);
      double xcenter = L0/4. - dx_maxlevel, zcenter = L0/2., pipe_height = L0/2. - dx_maxlevel;
      double radius_exit = 0.9*dx_maxlevel;
      double H_plane = 2*dx_minlevel + dx_maxlevel;
      
      /**Mesh preparation for the geometry, the area is refined around the future solid.*/
      refine (level < maxlevel && y < pipe_height \
          && x < xcenter + 2*dx_minlevel && x > xcenter - 2*dx_minlevel \
          && z < zcenter + 2*dx_minlevel && z > zcenter - 2*dx_minlevel );
      
      define_geometry (pipe_f0, pipe_wall, plane, xcenter, zcenter, pipe_height, radius_exit, H_plane);
      
      maxlevel -= 1;
      
      adapt_wavelet ({u, cs}, (double[]) {0.01, 1e-3, 0.01, 0.1}, maxlevel, minlevel);
    }
    /**Otherwise normal exit on 4 cells.*/
    else{
      double dx_maxlevel = L0/pow(2, maxlevel), dx_minlevel = L0/pow(2, minlevel);
      double xcenter = L0/4. - dx_maxlevel, zcenter = L0/2., pipe_height = L0/2. - dx_maxlevel;
      double radius_exit = 0.9*dx_maxlevel;
      double H_plane = 2*dx_minlevel + dx_maxlevel;
    
      //mesh preparation for geometry
      refine (level < maxlevel && y < pipe_height \
          && x < xcenter + 2*dx_minlevel && x > xcenter - 2*dx_minlevel \
          && z < zcenter + 2*dx_minlevel && z > zcenter - 2*dx_minlevel );
          
      define_geometry (pipe_f0, pipe_wall, plane, xcenter, zcenter, pipe_height, radius_exit, H_plane);
    }
    
    //cross-flow initialization
    foreach()
      u.x[] = cs[] < eps || pipe_f0 [] > eps? 0. : cross_flow_velocity;
    boundary ((scalar*) {u});    
  }
  
  else{
    if (pid() == 0) {printf ("****Restoring solution*****\n");}

    /**Refine_embed_linear() may trigger an assertion
    ("refine_embed_linear: Assertion `((double *) ((((Tree *)grid)->L[point.level-1]->m[(p ..."),
    so we follow the workaround from here :
    [User_forum](https://groups.google.com/g/basilisk-fr/c/l7LpA71SEys/m/FyPcgDa4BQAJ).
    */
    
    for (scalar s in {u, g})
      s.prolongation = refine_injection;

    restore (file_restart);

    for (scalar s in {u, g})
      s.prolongation = refine_embed_linear;   

    if (onecell) {
      maxlevel += 1;

      double dx_maxlevel = L0/pow(2, maxlevel), dx_minlevel = L0/pow(2, minlevel);
      double xcenter = L0/4. - dx_maxlevel, zcenter = L0/2., pipe_height = L0/2. - dx_maxlevel;
      double radius_exit = 0.9*dx_maxlevel;
      double H_plane = 2*dx_minlevel + dx_maxlevel;
      
      /**
      At high resolution this step might throw an assertion "Assertion `coarse(cs,child.x,0,0) && coarse(cs,0,child.y,0)' failed.".
      In this case, a solution that should work is to gradually increment the level in the lines below from minlevel to maxlevel.
      Refining everywhere around the solid/topography can also quickly lead to reach the maximum RAM allowed per CPU,
      some efforts can be made to refine precisely at the boundary between the fluid and the embed solid.
      Disabling the assertion in Refine_embed_linear() is also not a good idea, as the solution is more likely to blow up.
      */
      
      //mesh preparation for geometry
      refine (level < maxlevel && y < pipe_height \
          && x < xcenter + 2*dx_minlevel && x > xcenter - 2*dx_minlevel \
          && z < zcenter + 2*dx_minlevel && z > zcenter - 2*dx_minlevel );
      
      define_geometry (pipe_f0, pipe_wall, plane, xcenter, zcenter, pipe_height, radius_exit, H_plane);
      
      maxlevel -= 1;
      
      adapt_wavelet ({u, cs}, (double[]) {0.01, 1e-3, 0.01, 0.1}, maxlevel, minlevel);
    }
    /**Otherwise normal exit on 4 cells.*/
    else{
      double dx_maxlevel = L0/pow(2, maxlevel), dx_minlevel = L0/pow(2, minlevel);
      double xcenter = L0/4. - dx_maxlevel, zcenter = L0/2., pipe_height = L0/2. - dx_maxlevel;
      double radius_exit = 0.9*dx_maxlevel;
      double H_plane = 2*dx_minlevel + dx_maxlevel;
    
      //mesh preparation for geometry
      //Gradually increment the level from minlevel to maxlevel
      int i_level;
      
      for (i_level = minlevel; i_level <= maxlevel; i_level++) {
        refine (level < i_level && y < pipe_height \
            && x < xcenter + 2*dx_minlevel && x > xcenter - 2*dx_minlevel \
            && z < zcenter + 2*dx_minlevel && z > zcenter - 2*dx_minlevel );
      }
          
      define_geometry (pipe_f0, pipe_wall, plane, xcenter, zcenter, pipe_height, radius_exit, H_plane);
    }
    
    bool_restore = 1;
  }
  
  DT=timestep_ini;
}

event defaults (i++) {
  foreach_face()
    muc.x[] = 1e-6 * fs.x[];
  boundary ((scalar*) {muc});
}

event adapt (i++) {
  adapt_wavelet ({u, cs}, (double[]) {0.1, 0.1, 0.1}, maxlevel, minlevel);
  /**To kill the refinement in cs[] - not basilisk friendly*/
  foreach()
    if (cs[] < eps)
      foreach_dimension()
        u.x[] = 0.;
  boundary ((scalar*) {u});
}

event count_cell(i++) {
  int n = 0;
  foreach(reduction(+:n)) {n++;}
  
  if (pid() == 0) {
    FILE * fp = fopen ("count_cell", "a"); 
    fprintf (fp, "%d\t%d\t%f\n", i, n, dt);
    fflush (fp);
    fclose (fp); 
  }
}

event snapshot (i += dump_step)
{
  char name[80];
  sprintf (name, "dump-%d", i);
  p.nodump = pf.nodump = true;
  dump (file = name);

}

/**
## Output
*/

/**An output in MPI for paraview.*/
event export_to_paraview (i++) {
  if (pid() == 0) {printf ("** iteration %d**\n", i);}

  if (i%record_step == 0 || i==0) {
    char name[350], namebis[80];
    sprintf (name, "wip%d",i);
    sprintf (namebis, "wip%d", i);
    
    lambda2 (u, l2);
    boundary ({l2});
  
    foreach()
      foreach_dimension()
        omega.z[] = (u.y[1] - u.y[-1] - u.x[0,1,0] + u.x[0,-1,0]) / (2.*Delta);
    boundary ((scalar*) {omega});
    
    output_vtu((scalar *) {cs, pipe_f0, pipe_wall, plane, u, l2, omega}, (vector *) {0}, name);
    time_output_test ("timetest", namebis, i, end_simu+1, i);
  }
}

/**A movie is generated (above) drawing the geometry and the lambda2 negative isosurface field for
the coherent structures identification.*/

event movie (i++) {
  if (i%record_step == 0 || i==0) {
    if (onecell) {
      clear();
      view (theta = 0., phi = 1.2,fov=11, tx = -0.3, ty = 0.4);
      draw_vof ("cs", "fs", edges = 1);
      draw_vof ("cs", "fs", edges = 0);
      box();

      save ("pipe_one_cell.png");
    }
    else{
      clear();

      view (theta = -0.2, phi = 0.5,fov=29, tx=-0.6, ty=-0.1);
      
      draw_vof ("cs", "fs");
      box();
      translate (y = L0/5.) {
        cells (n = {0, 1, 0});
      }
      translate (x = 0.9*L0) {
        squares ("u.x", n = {1, 0, 0});
      }
      isosurface ("l2", -1.);

      if (bool_restore) {
        save ("pipe_from_restore.mp4");
        if (i == end_simu) {
          clear();
          view (theta = 0., phi = 1.2,fov=11, tx = -0.3, ty = 0.4);
          draw_vof ("cs", "fs", edges = 1);
          draw_vof ("cs", "fs", edges = 0);
          box();

          save ("pipe_from_restore_anomaly.png");
        }
      }
      else
        save ("pipe_normal.mp4");
    }
  }
}

event end (i = end_simu) {}