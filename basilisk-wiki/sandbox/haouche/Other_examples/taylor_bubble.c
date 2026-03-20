/**
# Taylor bubble

First, I would like to thank Edoardo Cipriano for his valuable advice and insightful corrections.
In this section, I present the code used to simulate the Taylor bubble. 
For a step-by-step explanation, feel free to watch my video tutorial: [Taylor bubble](https://www.youtube.com/watch?v=yTQWNp18yuY).

We start by including the necessary libraries:
*/

#include "embed.h"
/**
This above library allow us to include the solid and must be above the metric and the main solver.
*/
#include "axi.h"
/**
Note that the metric must be above the Navier-Stokes solver.
*/
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "vtknew_cell.h"

/**
Here, we work in CGS units. Thus, we define all the parameters:

- **Lower density** $\rho_1 = 1cm^3/g$;
- **Upper density** $\rho_2 = 0.0012cm^3/g$;
- **Lower viscosity** $\mu_1 = 0.01P$;
- **Upper viscosity** $\mu_2 = 1.8\cdot10^{-4}P$;
- **Surface tension** $\sigma = 0.5dyn/cm$;
- **Initial length of the Taylor bubble** $L_0 = 0.4cm$;
- **Initial radius of the Taylor bubble** $R_0 = 0.1cm$;
- **Initial velocity** $U_0 = 1cm/s$;
*/

int MAXLEVEL = 8;
int MINLEVEL = 2;

double SIZE = 1. [*];
double RADIUS = 0.1;                      
double DENSITY_1 = 1.;                     
double DENSITY_2 = 0.0012;                    
double DYNAMIC_VISCOSITY_1 = 0.01;
double DYNAMIC_VISCOSITY_2 = 1.8e-4;
double SURFACE_TENSION = 0.5;
double VELOCITY = 1.;

/**
To close the problem, appropriate boundary conditions are prescribed to reproduce the motion of a Taylor bubble rising in a vertical channel.

A unidirectional flow is imposed at the inlet (left boundary), corresponding to a constant axial velocity $U_0$. This mimics a reference frame where the liquid flows upward and the bubble rises due to buoyancy effects. At the outlet (right boundary), a zero-gradient (Neumann) condition is imposed on the velocity field, allowing the flow to leave the domain without artificial reflections. The pressure is fixed to zero at the outlet to ensure a well-posed problem.

On the solid wall (tube boundary), no-slip conditions are enforced, i.e. both normal and tangential components of the velocity vanish. This models the viscous interaction between the liquid and the confining tube.

The boundary conditions are implemented as follows:
*/

u.n[left] = dirichlet(VELOCITY);
u.t[left] = dirichlet(0.);
p[left] = neumann(0.);

u.n[right] = neumann(0.);
u.t[right] = neumann(0.);
p[right] = dirichlet(0.);

/**
The tube is modeled using an embedded boundary, on which a no-slip condition is applied:
*/
u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);
/**
This ensures that the liquid adheres to the tube wall, leading to the formation of a thin lubrication film around the Taylor bubble.
*/

/**
In the `main()` function, the computational domain is initialized as a rectangular box of length `SIZE` discretized using an initial uniform grid of resolution $2^{\text{MAXLEVEL}}$. The physical properties of the two phases (densities, viscosities) as well as the surface tension coefficient are then specified.

The simulation is subsequently launched using `run()`.
*/

int main() {
  rho1 = DENSITY_1;
  rho2 = DENSITY_2;
  mu1 = DYNAMIC_VISCOSITY_1;
  mu2 = DYNAMIC_VISCOSITY_2;
  f.sigma = SURFACE_TENSION;
  
  size(SIZE);
  init_grid(1 << MAXLEVEL);

  system("mkdir -p vtk");
  run();
}

/**
At $t=0$, a Taylor bubble is initialized as an elongated structure aligned with the tube axis. The interface is defined using a combination of cylindrical and spherical cap geometries to approximate the characteristic bullet shape of a Taylor bubble.

The tube confinement is introduced using an embedded boundary of radius $R_{\text{tube}}$.
*/

event init(i = 0) {
  double Rtube = 0.13;
  solid(cs, fs, Rtube - y);
  
  double L = 0.4;               
  double xc = SIZE/3.;          
  fraction(f, -(sq(RADIUS) - sq(y)) * (fabs(x - xc) <= L/2.) + (-(sq(RADIUS) - sq(y)) + sq(fabs(x - xc) - L/2.)) * (fabs(x - xc) > L/2.));
    
#if defined(AXI) && defined(EMBED)
  cm_update(cm, cs, fs);
  fm_update(fm, cs, fs);
  restriction({cs, fs, cm, fm});
#endif

  foreach()
    u.x[] = cs[]*u0;
}

/**
To efficiently capture the interface dynamics and velocity gradients while
limiting the overall computational cost, an adaptive mesh refinement strategy
is employed. The grid adaptation is performed using the wavelet-based
refinement algorithm provided by Basilisk.

The refinement criterion is applied to both the volume fraction field $f$,
which tracks the interface, the embedded boundary field $c_s$ and the velocity field $\mathbf{u}$. Cells are refined or
coarsened according to the estimated discretization error compared to user-
defined thresholds.

The maximum level of refinement is set to `MAXLEVEL`, while the minimum level is
restricted to `MINLEVEL` in order to maintain a sufficiently coarse background
grid.

The adaptation procedure is implemented as follows:
*/

#if TREE
event adapt(i++) {
  adapt_wavelet ({f, cs, u}, (double[]){1e-3, 1e-2, 1e-2, 1e-2}, MAXLEVEL, MINLEVEL);
# if defined(AXI) && defined(EMBED)
  cm_update(cm, cs, fs);
  fm_update(fm, cs, fs);
  restriction({cs, fs, cm, fm});
# endif
}
#endif

/**
Finally, this event handles data output using `output_vtk()`,  
which saves the volume fraction $f$, the embedded boundary field $c_s$, the pressure $p$, and the velocity field $\mathbf{u}$ in `.vtk` format, every 0.001 seconds.
*/

int count = 0;
event vtk(t=0.; t+=0.001; t<=0.2) {
    char filename[100];
    sprintf(filename, "vtk/snap-%d.vtk", count);
    FILE * fileptr = fopen(filename,"w");

    scalar * list = {f, u.x, u.y, p, cs};
    output_vtk(list, fileptr);
    fclose(fileptr);
    
    if (pid() == 0) {
        fprintf(stderr, "t: %-10g\t dt: %-10g\n", t, dt);
    }

    count++;
}