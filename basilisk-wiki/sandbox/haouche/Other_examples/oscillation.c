/**
# Disruption of an interface

In this section, I present the code used to simulate the disruption of an interface. 
For a step-by-step explanation, feel free to watch my video tutorial: [Disruption of an interface](https://www.youtube.com/watch?v=yTQWNp18yuY).

We start by including the necessary libraries:
*/

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "reduced.h"
#include "tension.h"
#include "vtknew_cell.h"

/**
Here, we work in CGS units. Thus, we define all the parameters:

- **Lower density** $\rho_1 = 1cm^3/g$;
- **Upper density** $\rho_2 = 0.001cm^3/g$;
- **Lower viscosity** $\mu_1 = 0.01P$;
- **Upper viscosity** $\mu_2 = 1.8\cdot10^{-4}P$;
- **Surface tension** $\sigma = 72dyn/cm$;
- **Gravity** $g = 981cm/s^2$.

We also set the initial amplitude $a_0$, and the reference length $L$.
*/

double L = 10.;
double a0 = 0.1;
double density_upper = 0.001;
double density_lower = 1.; 
double viscosity_upper = 0.3;
double viscosity_lower = 0.01;
double surface_tension = 72.;
double gravity = 981.;

/**
In the `main()` function, we initialize the domain size (a square domain of side `SIZE_BOX`), the uniform grid resolution $2^8 = 256$, and the physical parameters: densities, viscosities, and surface tension.  
We then call `run()` to launch the simulation.
*/

int main() {
    size(L);
    N = 256;

    rho1 = density_lower;
    rho2 = density_upper;

    mu1 = viscosity_lower;
    mu2 = viscosity_upper;

    f.sigma = surface_tension;
    G.y = -gravity;
    
    system("mkdir -p vtk interface");
    run();
}
/**
We now define the initial condition. The interface is initialized according to: $$\eta\big(x,\,t=0\big)=\frac{L}{2} + a_0\cos\big(2\pi x/L\big)$$
and placed in the middle of the box:
$$F(x, y, t=0) = -y + \bigg(\frac{L}{2} + a_0 \cos(2\pi x/L)\bigg)$$
*/

event init(t=0) {
    fraction(f, -y + L/2. + a0*cos(2*pi*x/L));
}

/**
Finally, this event handles data output using `output_vtk()`,  
which saves the volume fraction $f$, the pressure $p$, and the velocity field $\mathbf{u}$ in `.vtk` format, every 0.01 time units. We also save the shape of the interface $(f=0.5)$ using `output_facets(f, ptr)`.
*/

int count = 0;
event vtk(t=0.; t+=0.01; t<=1.) {
    char filename[100];
    sprintf(filename, "vtk/snap-%d.vtk", count);
    FILE * fileptr = fopen(filename,"w");

    scalar * list = {f, u.x, u.y, p};
    output_vtk(list, fileptr);
    fclose(fileptr);

    char name[80];
    sprintf(name, "interface/snap-%d.dat",count);
    FILE *ptr = fopen(name,"w");
    output_facets(f, ptr);
    fclose(ptr);
    
    if (pid() == 0) {
        fprintf(stderr, "t: %-10g\t dt: %-10g\n", t, dt);
    }

    count++;
}

/**
![Evolution of the volume fraction field](rayleigh-taylor/fields-initial.png)
*/