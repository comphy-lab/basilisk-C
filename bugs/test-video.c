/**
 * Simple test to verify that video generation tools work in Docker
 * This creates a basic Basilisk simulation that outputs video
 */

#include "grid/cartesian2D.h"
#include "run.h"
#include "output.h"

int main (int argc, char * argv[])
{
  // Simple domain
  L0 = 1.;
  N = 32;
  T = 0.1;
  
  run();
  return 0;
}

event init (t = 0) {
  // Initialize a simple scalar field
  scalar s[];
  foreach()
    s[] = x + y;
}

event movie (t += 0.01) {
  // Test video output - this is where the ppm2mp4 issue occurs
  output_ppm (s, file = "test.mp4", n = 32, linear = true);
}