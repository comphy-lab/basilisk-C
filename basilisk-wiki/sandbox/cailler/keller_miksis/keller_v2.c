/**
## Visualization of the self-similarity

~~~gnuplot Study in Physical and Self-Similar Spaces for $\beta = 140$°
reset 

array times[8]
times[1] = 5*10**(-5)
times[2] = 1*10**(-4)
times[3] = 2*10**(-4)
times[4] = 5*10**(-4)
times[5] = 1*10**(-3)
times[6] = 2*10**(-3)
times[7] = 5*10**(-3)
times[8] = 1*10**(-2)



array counter[8]
counter[1] = 5
counter[2] = 6
counter[3] = 7
counter[4] = 8
counter[5] = 9
counter[6] = 10
counter[7] = 11
counter[8] = 12

beta = 140 # chosen angle of contact (fixed value)

t0 = 0
xi_0 = 0.4
Y0_keller = 0.03767882187938287 

f(x) = x # bisector function

# rotation angle for Keller & Miksis data:
rot = 359.5 * pi /180.

set size ratio 1
set key top left

set multiplot

    # main plot in self-similar-space [main plot]
    set size 1, 1
    set origin 0, 0
    set xlabel 'xi = (x - x_0) * (rho/sigma.t^2)^(^1^/^3^)'
    set ylabel 'eta = y * (rho/sigma.t^2)^(^1^/^3^)'
    set xrange [0:4]
    set yrange [0:4]
        
    plot for [i=4:8] sprintf('fig3_%d_%i', beta, counter[i])\
        u ($1 - xi_0 * (times[8])**(2./3.))/((times[i]-t0)**(2./3.)):($2)/((times[i]-t0)**(2./3.))\
        w l t sprintf('%i°, t = %g', beta, times[i]) , \
        f(x) linecolor rgb '#dd181f' dashtype 2 lw 0.5 notitle, \
        '../keller_data/keller_digitized_140_raw' u ($1 * cos(rot) - ($2 - Y0_keller) * sin(rot)):($1 * sin(rot) + ($2 - Y0_keller) * cos(rot) )\
         w l dashtype 2 lw 1 lt -1 t "Keller and Miksis"


    # Small plot in physical space [insert]
    set size 0.4, 0.4
    set origin 0.48, 0.15

    set xlabel 'x'
    set ylabel 'y'
    set xrange[0:1*10**(-1)]
    set yrange [0:1*10**(-1)]
    unset key

    plot for [i=4:8] sprintf('fig3_%d_%i', beta, counter[i]) w l

unset multiplot
~~~
*/

#include "navier-stokes/centered.h"
#include "contact.h"
#include "two-phase.h"
#include "tension.h"

#define LEVEL 11 // grid resolution level
#define XI_0 0.4 // offset in the self-similar space
#define T_F 1.e-2 // final time of the simulation
#define SIZE 1 // size of the box
#define THETA_0 45 // initial angle of contact

// The interface is represented by the volume fraction field 'f'
// scalar f[], * interfaces = {f} ; => already defined in "two-phase.h".

// Boundary condition defined for the contact angle:
vector h[] ;
double theta0 = THETA_0 ;
h.t[bottom] = contact_angle (theta0 * pi / 180.) ;

// Power scale for the self-similar shape:
double scale = 2./3. ;

int k = 0 ; // counter

int main() {
    size(SIZE) ;
    init_grid(1 << LEVEL) ;

    rho1 = 1., rho2 = 0.001 ;
    f.height = h ;
    f.sigma = 1. ;

    // Basilisk contact angles are not defined as Kellers
    // with theta = 180 - beta => (theta0 = 40; theta0 <= 120; theta0 += 20)
    for (theta0 = 140; theta0 <= 140; theta0 += 20) {
        k = 0 ;
        h.t[bottom] = contact_angle (theta0 * pi / 180.) ;
        run() ;
    }
}


event init (t = 0){
    /* Offset for the tapered sheet of liquid for avoiding slip-wall condition 
    in the bottom-left corner (constant offset of value XI_0 in self-similar space) */
    double x_0 = XI_0 * pow(T_F, scale) ;

    // Initial free-surface definition:
    fraction(f, - (y - (x - x_0) * tan(THETA_0 * pi / 180.)) ) ;
}


event profiles (t = {1.e-7, 5.e-7, 1.e-6, 5.e-6, 1.e-5, 5.e-5, 
                    1.e-4, 2.e-4, 5.e-4, 
                    1.e-3, 2.e-3, 5.e-3, 
                    1.e-2
                    } 
                ){
                       
    char filename[30] ;
    sprintf(filename, "fig3_%g_%d", theta0, k) ;
    FILE * fp = fopen(filename, "w") ;
    output_facets (f, fp) ;
    fclose (fp) ;
    k += 1 ; 
}


event adapt (i += 10) {
  adapt_wavelet ({f,u}, (double[]){5e-4,1e-3,1e-3}, LEVEL);
}




