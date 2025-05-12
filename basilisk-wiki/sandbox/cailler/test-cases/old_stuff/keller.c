

/**
~~~gnuplot Interfaces in Physical and Self-Similar Spaces  
array times[1]
times[1] = 10

array thetas[5]
thetas[1] = 40
thetas[2] = 60
thetas[3] = 80
thetas[4] = 100
thetas[5] = 120

t0 = 0


set key top left

set multiplot

    # main plot in self-similar-space [main plot]
    set size 1, 1
    set origin 0, 0
    set xlabel 'xi = x * (rho/sigma.t^2)^(^1^/^3^)'
    set ylabel 'eta = y * (rho/sigma.t^2)^(^1^/^3^)'
    set xrange [0:4]
    set yrange [0:4]

    plot for [i=1:5] sprintf('fig3_%i',thetas[i]) u ($1)/((times[1]-t0)**(2./3.)):($2)/((times[1]-t0)**(2./3.))\
    w l t sprintf('angle = %i °', 180 - thetas[i])


    # Small plot in physical space [insert]
    set size 0.4, 0.4
    set origin 0.58, 0.15

    set xlabel 'x'
    set ylabel 'y'
    set xrange[0:15]
    set yrange [0:15]
    unset key

    plot for [i=1:5] sprintf('fig3_%i',thetas[i]) w l, 0 lt -1

unset multiplot
~~~
*/


#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "contact.h"
#include "vof.h"
#include "tension.h"

#define LEVEL 8

// The interface is represented by the volume fraction field 'c'
scalar c[], * interfaces = {c} ;

vector h[] ;
double theta0 = 45 ;
h.t[bottom] = contact_angle (theta0 * pi / 180.) ;

int main() {
    size(24) ;
    init_grid(1 << LEVEL) ;

    c.height = h ;
    c.sigma = 1. ;

    // Basilisk contact angles are not defined as Kellers
    // with theta = 180 - beta
    for (theta0 = 40; theta0 <= 120; theta0 += 20) {
        h.t[bottom] = contact_angle (theta0 * pi / 180.) ;
        run() ;
    }
}

event init (t = 0){
    fraction(c, y - x * tan(45 * pi / 180.)) ;
}



event profiles (t = 10){
    char filename[20] ;
    sprintf(filename, "fig3_%g", theta0) ;
    FILE * fp = fopen(filename, "w") ;
    output_facets (c, fp) ;
    fclose (fp) ;
}
