/**
# Rayleigh–Taylor Simulation Code

In this section, I present the code used to simulate the Rayleigh–Taylor instability.  
For a step-by-step explanation, feel free to watch my video tutorial: [link here].

We start by including the necessary libraries:
*/

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "output_vtu.h"

/**
We define the dimensionless numbers of the problem:

- **Bond number** $Bo = \frac{\rho_l g L^2}{\gamma}$, which compares gravitational forces to surface tension forces;
- **Ohnesorge number** $Oh = \frac{\mu_l}{\sqrt{\rho_l \gamma L}}$, which quantifies the relative importance of viscous forces compared to inertial and capillary effects;
- **Density ratio** $\rho_r = \frac{\rho_1}{\rho_2}$;
- **Viscosity ratio** $\mu_r = \frac{\mu_1}{\mu_2}$.

We also set the wavenumber $k$, the initial amplitude $a_0$, and the reference length $L$.
*/

#define RHOR 0.1383
#define MUR 1.
#define Oh 0.05
#define Bo 100.
#define kk (2.*pi)
#define a0 0.05
#define L 1.

/**
We specify the maximum level of mesh refinement, the size of the domain, and the final simulation time.
*/

int LEVEL = 9;
double SIZE_BOX = 4.;
double FINAL_TIME = 5.;

/**
In the `main()` function, we initialize the domain size (a square domain of side `SIZE_BOX`), the initial grid resolution $2^6 = 64$, and the physical parameters: densities, viscosities, and surface tension.  
We then call `run()` to launch the simulation.
*/

int main() {
    size(SIZE_BOX);
    init_grid(1 << 6);

    rho1 = 1.;
    rho2 = RHOR;
    mu1 = Oh * sqrt(rho1 * L / Bo);
    mu2 = mu1 * MUR;
    f.sigma = 1./Bo;
    
    system("mkdir -p dumpfile");
    run();
}
/**
We now define the initial condition. The interface is initialized according to: $$\eta\big(x,\,t\big)=a_0\cos\big(kx\big)$$
and placed in the middle of the box:
$$F(x, y, t=0) = y - \bigg(\frac{\mathrm{SIZE\_BOX}}{2} + a_0 \cos(kx)\bigg)$$
We also use a mask to change the geometry: instead of a square, we simulate a rectangular domain.
*/

event init (i = 0) {
    mask(x > L ? right : none);
    if (!restore (file = "restart")) {
        refine (y - (1.5*(SIZE_BOX/2. + a0*cos(kk*(x)))) < 0 && level < LEVEL);
        fraction (f, (y)  - (SIZE_BOX/2. + a0*cos(kk*(x))));
    }
}

/**
Gravity is applied using an `acceleration` event.  
Alternatively, you could use `#include "reduced.h"` to implicitly include gravity via the pressure gradient.
*/

event acceleration(i++) {
    face vector av = a;
    foreach_face(y)
        av.y[] += -1.;
    boundary ((scalar *){av});
}

/**
We use `adapt_wavelet()` to dynamically adapt the mesh during the simulation.  
Here, we define error tolerances for the volume fraction and the velocity components to improve numerical accuracy.
*/

event adapt(i++) {
    double femax = 0.001, uxemax = 0.01, uyemax = 0.01;
    adapt_wavelet({f, u}, (double[]){femax, uxemax, uyemax}, LEVEL); 
}

/**
This event is optional. It simply prints runtime information to track the simulation progress.
*/

event running_status(i++) {
    if (pid() == 0)
        fprintf(stderr, "t: %-10g\t dt: %-10g\n", t, dt);
}

/**
Finally, this event handles data output using `output_vtu_bin_foreach()`,  
which saves the volume fraction $f$, the pressure $p$, and the velocity field $\mathbf{u}$ in `.vtu` format, every 0.01 time units.
*/

int pt = 0;
event interface(t += 0.01; t <= FINAL_TIME) {
    char name2[80];
    sprintf(name2, "dumpfile/snap-%d.vtu", pt);
    FILE *ptr2 = fopen(name2,"w");
    output_vtu_bin_foreach((scalar *) {f, p}, (vector *) {u}, ptr2, false);
    fclose(ptr2);
    pt++;
}

/**
![Evolution of the volume fraction field](rayleigh-taylor/fields-initial.png)
*/