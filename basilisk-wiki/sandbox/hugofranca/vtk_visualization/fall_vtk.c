/**
# Impact of a viscoelastic drop on a solid
This file is a variation of the original Basilisk example [fall.c](/src/test/fall.c). We use this as a simple example for the VTK output functions and visualization.

We solve the axisymmetric, incompressible, variable-density,
Navier--Stokes equations with two phases and use the log-conformation
method to include viscoelastic stresses. The curvature module is used
to compute interface properties. */

#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "log-conform.h"
#include "curvature.h"

/**
Including my personal functions to output Legacy VTK files */
#include "hugo_functions.h"

/**
The density and viscosity ratios are 1000. The Reynolds number based
on the droplet diameter, velocity and viscosity is 5 and the Froude
number is 2.26. */

#define RHO_r 0.001
#define MU_r 0.001
#define LEVEL 7

/**
Below we set the 4 nondimensional parameters for a simulation (Reynolds,
Froude, Wi, beta). And we solve 6 simulations in total, with varying
Weissenberg numbers. */
#define NUMBER_SIMULATIONS 6
double list_Re[NUMBER_SIMULATIONS] = {5.0, 5.0, 5.0, 5.0, 5.0, 5.0};
double list_Fr[NUMBER_SIMULATIONS] = {2.26, 2.26, 2.26, 2.26, 2.26, 2.26};
double list_Wi[NUMBER_SIMULATIONS] = {0.0, 0.5, 1.0, 5.0, 10.0, 100.0};
double list_beta[NUMBER_SIMULATIONS] = {1.0, 0.1, 0.1, 0.1, 0.1, 0.1};
double Re, Fr, Wi, beta;
scalar lambdav[], mupv[];

/**
The drop comes from the right. We allow the fluid to get through that
boundary. */
u.n[right] = neumann(0);
p[right]   = dirichlet(0);

/** 
No penetration on the volume fraction at the left wall and no slip on the velocity. */
u.t[left] = dirichlet(0);
tau_qq[left] = dirichlet(0);
f[left] = 0.0;

int main(int argc, char * argv[]) {

  size (5.0);

  /**
  The viscoelastic fields will be set below. */
  mup = mupv;
  lambda = lambdav;

  /**
  We set a maximum timestep. This is necessary for proper temporal
  resolution of the viscoelastic stresses. */
  DT = 2e-3;

  /**
  Looping over each simulation and running it. */
  for(int i=0; i<NUMBER_SIMULATIONS; i++) {
    Re = list_Re[i];
    Wi = list_Wi[i];
    beta = list_beta[i];
    Fr = list_Fr[i];
    
    /** 
    Opening a new folder for the all the output of this simulation. <br>
    All the output functions (PrintMesh, PrintInterface, PrintLog) will 
    automatically send stuff to this folder. */
    OpenSimulationFolder ("fall_mesh%d_Wi%g", LEVEL, Wi);

    if( !pid() ) {
      printf ("Starting new simulation with following parameters: \n");
      printf ("Re = %lf  Fr = %lf  Wi = %lf  beta = %lf\n", Re, Fr, Wi, beta);
    }

    init_grid (1 << LEVEL);

    /**
    The densities and viscosities are defined by the parameters above. */
    rho1 = 1.;
    rho2 = RHO_r;
    mu1 = beta/Re;
    mu2 = MU_r/Re;

    run();
  }

  return 0;
}

event init (t = 0) {
  /**
  At a wall of normal $\mathbf{n}$ the component of the viscoelastic
  stress tensor $tau_p_{nn}$ is zero. Since the left boundary is a wall, we
  set $tau_p_{xx}$ equal to zero at that boundary. */ 
  scalar s = tau_p.x.x;
  s[left] = dirichlet(0.);

  /**
  The drop is centered on (2,0) and has a radius of 0.5. */
  fraction (f, - sq(x - 2.0) - sq(y - 0.0) + sq(0.5));

  /**
  The initial velocity of the droplet is -1. */
  foreach()
    u.x[] = - f[];
}

/**
We add the acceleration of gravity. */
event acceleration (i++) {
  face vector av = a;
  foreach_face(x)
    av.x[] -= 1./sq(Fr);
}

/**
We update the viscoelastic properties. Only the droplet is
viscoelastic. */
event properties (i++) {
  foreach() {
    mupv[] = (1. - beta)*clamp(f[],0,1)/Re;
    lambdav[] = Wi*clamp(f[],0,1);
  }
}

/**
We adapt the solution at every timestep based on the interface and
velocity errors. */
#if TREE
event adapt (i++) {
  adapt_wavelet ({f, u.x, u.y}, (double[]){1e-2, 5e-3, 5e-3},
	    maxlevel = LEVEL, minlevel = LEVEL - 2);
}
#endif

/**
We track the spreading diameter of the droplet. 
The diameter is printed over time to the screen and to a log file.*/
event logfile (i++; t <= 5.0) {
  scalar pos[];
  position (f, pos, {0,1});
  double spreading = 2.0*statsf(pos).max;

  PrintLog("%g %g %.15lf\n", t, dt, spreading);

  if( !pid() )
    printf("%g %g %.15lf\n", t, dt, spreading);
}

/**
We print the quadtree mesh with some properties (velocity) to a VTK file. 
We also print the interface to a separate VTK file. */
event print_view (t+=0.025) {
  PrintMeshVTK(i, t, {u.x, u.y}, 
      (const char *[]){"Axial-Velocity", "Radial-Velocity"},
      vtk_type=VTK_TYPE_BINARY, vtk_precision=VTK_PRECISION_FLOAT);
  PrintInterfaceVTK(i, t);
}

/**
## Results

The results can be seen in [this video](https://www.youtube.com/watch?v=9n3WoucL-wc). (TO_DO: I will try to embed the video into the page later)

The video above was created using the python visualization package pyvista. You can recreate it by running [this script](fall_video_script.py).

Alternatively, the files can be visualized using paraview or any other visualization software capable of reading VTK files.
*/