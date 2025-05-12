#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"

#include "hugo_functions.h"

/**
The density and viscosity ratios are 1000. The Reynolds number based
on the droplet diameter, velocity and viscosity is 5 and the Froude
number is 2.26. */
#define RHO_r 0.001
#define MU_r 0.001
#define LEVEL 7

/**
Dimensionless numbers for the Newtonian simulation. 
*/
double Re = 5.0;
double Fr = 2.26;


scalar lambdav[], mupv[];

/**
The drop comes from the top. We allow the fluid to get through that
boundary. */
u.n[top] = neumann(0);
p[top]   = dirichlet(0);

/**
The wall is at the bottom side. We apply a no-slip boundary condition and a 
non-wetting condition for the VOF tracer.  */
u.t[bottom] = dirichlet(0);
u.r[bottom] = dirichlet(0);
f[bottom]   = 0.0;

int main(int argc, char * argv[]) {

  printf("Starting a simulation with the following parameters: \n");
  printf("Reynolds: %lf\n", Re);
  printf("Froude: %lf\n", Fr);
  

  /** 
    Opening a new folder for the all the output of this simulation
    All the output functions (PrintMesh, PrintInterface, PrintLog) will 
    automatically send stuff to this folder. */
  OpenSimulationFolder ("fall3D_mesh%d_Newt", LEVEL);


  /**
    The domain spans $[-2.5:2.5]\times[0:5]\times[-2.5:2.5]$. 
    The wall is located at y = 0. */
  size (5.0);
  origin(-2.5, 0.0, -2.5);
  init_grid (64);

  
  /**
    Densities and viscosities. */
  rho1 = 1.;
  rho2 = RHO_r;
  mu1 = 1.0/Re;
  mu2 = MU_r/Re;

  /**
    We set a maximum timestep. This is necessary for proper temporal
  resolution of the viscoelastic stresses. */
  DT = 2e-3;
  run();
}

event init (t = 0) {
  /**
    The drop is centered on (0, 2, 0) and has a radius of 0.5. */
  double center_y = 2.0;
  fraction (f, - sq(x - 0.0) - sq(y - center_y) - sq(z - 0.0) + sq(0.5));

  /**
    The initial velocity of the droplet is -1. */
  foreach()
    u.y[] = - f[];
}

/**
We add the acceleration of gravity. */
event acceleration (i++) {
  face vector av = a;
  foreach_face(y)
    av.y[] -= 1./sq(Fr);
}


/**
We adapt the solution at every timestep based on the interface and
velocity errors. */
#if TREE
event adapt (i++) {
  adapt_wavelet ({f, u.x, u.y, u.z}, (double[]){1e-2, 5e-3, 5e-3, 5e-3},
      maxlevel = LEVEL, minlevel = LEVEL - 2);
}
#endif

/**
We track the spreading diameter of the droplet. 
The diameter is printed over time to the screen and to a log file.*/
event logfile (i++; t <= 5.0) {
  scalar pos[];

  position (f, pos, {1, 0, 0});
  double max_x = statsf(pos).max;
  double min_x = statsf(pos).min;
  double spreading_x = max_x - min_x;

  position (f, pos, {0, 0, 1});
  double max_z = statsf(pos).max;
  double min_z = statsf(pos).min;
  double spreading_z = max_z - min_z;

  double avg_spreading = 0.5*( spreading_x + spreading_z );
  
  PrintLog("%g %g %.15lf %lf %lf\n", t, dt, avg_spreading, spreading_x, spreading_z);

  if( !pid() )
    printf("%g %g %.15lf %lf %lf\n", t, dt, avg_spreading, spreading_x, spreading_z);
}


event print_view (t+=0.025) {
  /**
    Printing the mesh to a VTK file. 
    You can use any of the variations below, depending on the type 
    of resulting file you want and which properties you want to print. */

  // === Prints just mesh without any scalar fields in it. Using binary files and float precision
  PrintMeshVTK_Binary_Float(i, t, NULL, NULL);

  // === Prints just mesh without any scalar fields in it. Using binary files and double precision (bigger files!!!)
  // PrintMeshVTK_Binary_Double(i, t, NULL, NULL);

  // // === Prints mesh with two cell-centered scalar fields in it. Using float precision
  // // === (you might want/need double precision depending on what you're solving/printing. watch out.)
  // scalar *list = {u.x, u.y};
  // const char *list_names[] = {"vel-x", "vel-y"};
  // PrintMeshVTK_Binary_Float(i, t, list, list_names);

  // === Prints mesh with two cell-centered scalar fields in it. Using double precision
  // scalar *list = {u.x, u.y};
  // const char *list_names[] = {"vel-x", "vel-y"};
  // PrintMeshVTK_Binary_Double(i, t, list, list_names);
  
  /**
    Printing the droplet interface to a VTK file. 
    Currently printing in ASCII format. BINARY version will be implemented soon (let me know if you desperately need it).*/
  PrintInterfaceVTK(i, t);
}

/**
## Results

The results can be seen in [this video](https://www.youtube.com/watch?v=KVpRKajIoQY). 
(TO DO: I will try to embed the video into the page later)

The video above was created using the python visualization package pyvista. 
You can recreate it by running [this script](fall3D_video_script.py).

Alternatively, the files can be visualized using paraview or 
any other visualization software capable of reading VTK files.
*/
