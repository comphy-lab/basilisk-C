/**
This code computes the dynamics of a traveling population $n$ which diffuses, grows, and advects chemotactically toward an attractant $\phi$ which is distributed heterogeneously throughout a fractal terrain. The boundary conditions may be adjusted to simulate traveling wave behavior, and the landscape may be modified to include any arbitrary terrain.

![This movie shows the traveling wave of chemotactic population $n$.](fisherkpp/n.mp4)
![This movie shows the attractant $\phi$ as it is resolved in time as the traveling population senses it.](fisherkpp/phi.mp4)
![This movie shows the terrain as it is resolved in time for the population to pass.](fisherkpp/terrain.mp4)
*/

#include "view.h"
#include "diffusion.h" 
#include "advection.h"
#include "tracer.h"
#include "run.h"

/** 
Define the simulation parameters, including timespan, domain size, species diffusivities, growth rates, carrying capacities, and chemotactic mobility.
*/

#define tfinal 10
#define time_step 1
#define length 1000

#define population_diffusivity 1
#define attractant_diffusivity 10
#define population_growth_rate 0.1
#define attractant_growth_rate 0.1
#define population_carrying_capacity 2
#define attractant_carrying_capacity 2
#define chemotactic_mobility 100 // 1

/**
The midpoint-displacement algorithm for generating a fractal terrain requires size and roughness parameters, which we arbitrarily choose here. The size must be $2^n+1$ for any integer $n$.
*/

#define SIZE 513 
#define ROUGHNESS 2.

/**
Initialize the population $n$, its chemotactic velocity $v_\text{chemo}$, the attractant $\phi$, and the terrain.
*/

scalar n[], phi[], * tracers = {n, phi};
face vector v_chemo[];
scalar terrain[];

/** 
Set up the simulation and initialize the helper variable landscape.
*/

int main() {
  DT = time_step;
  L0 = length; origin (0, 0);
  run();
}
double landscape[SIZE][SIZE];

/**
Implement the midpoint displacement algorithm which generate a fractal terrain and stores it in "landscape".
*/

void midpoint_displacement(int size, double roughness) {
    int step_size = size - 1;
    int x, y;

    landscape[0][0] = fabs(noise());
    landscape[0][size - 1] = fabs(noise());
    landscape[size - 1][0] = fabs(noise());
    landscape[size - 1][size - 1] = fabs(noise());

    while (step_size > 1) {
        int half_step = step_size/2;

        for (y = half_step; y < size - 1; y += step_size) {
            for (x = half_step; x < size - 1; x += step_size) {
                double avg = (landscape[y - half_step][x - half_step] +
                              landscape[y - half_step][x + half_step] +
                              landscape[y + half_step][x - half_step] +
                              landscape[y + half_step][x + half_step])/4;
                double displacement = noise() * roughness * step_size;
                landscape[y][x] = avg + displacement;
            }
        }

        for (y = 0; y < size - 1; y += half_step) {
            for (x = (y + half_step) % step_size; x < size - 1; x += step_size) {
                double avg = (landscape[(y - half_step + size - 1) % (size - 1)][x] +
                              landscape[(y + half_step) % (size - 1)][x] +
                              landscape[y][(x - half_step + size - 1) % (size - 1)] +
                              landscape[y][(x + half_step) % (size - 1)])/4;
                double displacement = noise() * roughness * step_size;
                landscape[y][x] = avg + displacement;
            }
        }

        step_size /= 2;
        roughness *= 1.1;
    }
}

/** 
Because the grid will be adapting, we must interpolate the landscape at each time step.
*/

double interpolate_landscape(double x, double y) {
    int x0 = (int)(x*(SIZE - 1)/length);
    int y0 = (int)(y*(SIZE - 1)/length);
    int x1 = x0 + 1;
    int y1 = y0 + 1;

    double q11 = landscape[x0][y0];
    double q12 = landscape[x0][y1];
    double q21 = landscape[x1][y0];
    double q22 = landscape[x1][y1];

    double x_frac = (x*(SIZE - 1)/length) - x0;
    double y_frac = (y*(SIZE - 1)/length) - y0;

    double r1 = (1 - x_frac)*q11 + x_frac*q21;
    double r2 = (1 - x_frac)*q12 + x_frac*q22;

    return (1 - y_frac)*r1 + y_frac*r2;
}

/**
Initialize the simulation for $t=0$.
*/

event init (t = 0) {

  srand (101);
  midpoint_displacement(SIZE, ROUGHNESS);

  foreach()
    terrain[] = interpolate_landscape(x, y);
  adapt_wavelet ({terrain}, (double[]){10}, 9, 5);

  n[left] = neumann(0); n[right] = neumann(0); 
  phi[left] = neumann(0); phi[right] = neumann(0); 

  foreach() {
    n[] = 1;
    phi[] = 1;
  }
}

/**
Implement the growth and diffusion terms.
*/

event tracer_diffusion (i++) {

  foreach()
    terrain[] = interpolate_landscape(x, y);

  scalar beta[]; // this will hold the reaction terms

  // diffusion of the population
  face vector D_n[];
  foreach_face()
  D_n.x[] = population_diffusivity;
  foreach()
    beta[] = population_growth_rate*(1-n[]/population_carrying_capacity)*(
      (exp(-pow((terrain[]-0)/300 * 0, 2))));
  diffusion (n, dt, D_n, beta=beta);

  // diffusion of the attractant
  face vector D_phi[];
  foreach_face()
    D_phi.x[] = attractant_diffusivity;
  foreach()
    beta[] = attractant_growth_rate*(1 - phi[]/attractant_carrying_capacity)*(
      (exp(-pow((terrain[]-200)/300, 2))));
  diffusion (phi, dt, D_phi, beta=beta);
}

/**
Implement chemotactic advection of $n$ along the gradient of $\phi$.
*/

event tracer_advection (i++) {
  foreach_face() {
    v_chemo.x[] = chemotactic_mobility*(phi[] - phi[-1,0])/Delta;
  }
  advection({n}, v_chemo, DT);
}

/**
Adapt the grid at each time step.
*/

event adapt (i++) {
  adapt_wavelet ({n}, (double[]){1.e-3}, 11, 5);
}

/**
Output data into movies at the specifed time intervals, and end the simulation when desired.
*/

event output (t += time_step) {
  output_ppm (n, file = "n.mp4", n=400);
  output_ppm (phi, file = "phi.mp4", n=400);
  output_ppm (terrain, file = "terrain.mp4", n=400);
}

event stop (t = tfinal) {}

