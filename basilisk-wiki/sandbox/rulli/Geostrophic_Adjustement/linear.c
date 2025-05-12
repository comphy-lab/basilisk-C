/** 

# Linearised shallow water equations

This code is solving linear shallow water equation for a step as initial condition for the height $\eta$.

~~~gnuplot Comparision between intial and final (steady state) height
set xlabel "X-axis"
set ylabel "Z-axis"
set grid
set yrange [-1.2:1.2]

plot "h0.dat" using 1:2 with lines title 'Initial surface height', \
    "hf.dat" using 1:2 with lines title 'Final surface height (t=200)', \
    "h0.dat" using 1:3 with lines title 'Fond'
~~~

![Height Animation](linear/eta.mp4)(width="200" height="200")

### Boundary conditions

In order to get the steady state as soon as possible, we avoid reflection on the right boundary by using an [open boundary condition](http://basilisk.fr/src/layered/hydro.h#radiation-boundary-conditions) while keeping a dirichlet boundary condition on the left side to prevent the height to decrease.

## Basilisk Code
*/



#include "grid/multigrid.h"
#include "layered/hydro.h"


/**
### Initialisation

We set the initial interface to a step ($\eta(x < 0) = 1$ and $\eta(x > 0 = -1$), and the velocity to 0.
*/


event init(i=0) {

    double h0 = 1;
    foreach() {
        zb[] = -1;
        eta[] = (x < 0)*h0 + (x > 0)*(-h0);
        h[] = eta[] - zb[];
        //eta will be initialised by the hydro.h to zb + h
        u.x[] = 0;
        u.y[] = 0;
    }

    u.n[left] = dirichlet(0);
    u.n[right] = radiation(0);

    // Save the initial interface
    FILE * file = fopen("h0.dat", "w");
    foreach()
    fprintf(file, "%g %g %g %g\n", x, eta[], zb[], t);
}



/**
### Main function
*/
int main() {

    size(8);
    origin(-4, -1);
    init_grid(64);
    
    linearised = true;
    G = 1;
    
    run();

}


/**
## Events
*/
event logs(i++) {
    fprintf(stderr, "%d %g\n", i, t);
}

event animation(t += 0.5) {
    output_ppm(eta, file="eta.mp4", linear=true, min=-1, max=1);
}

event end(t = 100) {
    FILE * file = fopen("hf.dat", "w");
    foreach()
    fprintf(file, "%g %g %g %g\n", x, eta[], zb[], t);
}


