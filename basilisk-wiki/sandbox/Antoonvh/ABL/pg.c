/**
# Perlin's surface generator

This (example) script can use command line arguments. E.g., for 7
layers, starting from 15m amplitude, evaluated on a grid of 512 points
in each direction, use:

~~~literatec
$./a.out 7 15 512
~~~

The output is an image `h.png`,

![Height map](pg/h.png)

and an ascii matrix data `h.dat`. This can also be plotted:

~~~gnuplot Mind the aspect ratio
set zr [-20:20]
set pm3d
splot 'h.dat' w l lw 0.001
~~~

We overload the `SMOOTHSTP` function to be a bit more steep. It can be
controlled by `T` > 1. The limit `T` = 1 gives Perlin's smooth step.
 */
double T = 3;
#define SMOOTHSTP(x) (pow((6*cube(x)*sq(x) - 15*sq(sq(x)) + 10*cube(x)), T))
#include "perlin.h"
#include "utils.h"
int ml = 5;     // Pelin layers
double A0 = 10; // Amplitude of first layer
int base = 6;   // First layer frequency 

int main(int argc, char ** argv) { // Command line arguments
  N = 128;
  if (argc > 1)
    ml = atoi (argv[1]);
  if (argc > 2)
    A0 = atof (argv[2]);
  if (argc > 3)
    N = atof (argv[3]);
  foreach_dimension()
    periodic(left);
  L0 = 300;
  init_grid(N);
  srand (time(NULL)); // Non-reproducible output

  scalar h[];
  foreach()
    h[] = 0;
  
  for (int layers = 0; layers < ml; layers++) {
    double nx = pow(2, layers) + base;   // Each new layer has a double frequency above the `base` value.
    double ny = pow(2, layers) + base;
    double A = A0/(double)(layers + 1.); // This function controls the decay of the spectum
    // double A = A0*exp(-layers/3.);
    init_perlin (nx, ny);
    foreach() 
      h[] += A*perlin(x, y, nx, ny);
    free (gradp); //goes after `init_perlin`
  }
  boundary ({h});
  output_ppm (h, file = "h.png", n = 300);
  output_field ({h}, fopen ("h.dat", "w"), N, true);
  
}
