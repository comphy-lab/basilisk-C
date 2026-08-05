/**
# Reductions on GPUs */

#include "utils.h"

double global_sum = 1., global_max = -1.;

int main (int argc, char * argv[])
{
  init_grid (argc > 1 ? atoi(argv[1]) : 1024);

#if !_CUDA // fixme: does not work with CUDA yet
  periodic (right);
  periodic (top);
#endif
  
  size (2.*pi);

  scalar s[];
  
  foreach (serial)
    s[] = sq (cos(2.*x)*cos(2.*y));
  
  timer t = timer_start();

  double sum = 0.;
  int iter;
  for (iter = 0; iter < 400*1024/N; iter++) {
    sum = 0.;
    foreach(reduction(max:sum))
      sum = max (sum, s[]);
  }
  
  double elapsed = timer_elapsed (t);
  printf ("N: %d elapsed: %g speed: %g\n",
	   N, elapsed, grid->tn*iter/elapsed);

  fprintf (stderr, "result: %g\n", sum);
  output_ppm (s, file = "s.png", n = 512, spread = -1);

  /**
  Check that "inout" fields work. */
  
  foreach()
    s[] += sq (0.5*sin(8.*x)*sin(8.*y));
  foreach()
    s[] = s[] + sq (0.5*sin(8.*x)*sin(8.*y));
    
  stats stat = statsf (s);
  fprintf (stderr, "min: %g max: %g\n", stat.min, stat.max);
  output_ppm (s, file = "s1.png", n = 512, spread = -1);

  /**
  Check that reduction on levels works. */

  restriction ({s});
  sum = 0.;
  foreach_level(0, reduction(+:sum))
    sum += s[]*dv();
  fprintf (stderr, "sum: %g %g\n", stat.sum, sum);

  /**
  Check that sum reductions work also when the initial value is not zero. */

  foreach()
    s[] = 1.;
  sum = 1.;
  foreach (reduction(+:sum))
    sum += s[];
  fprintf (stderr, "sum: %f\n", sum);

  /**
  Check reductions on global variables. */

  foreach (reduction(+:global_sum) reduction(max:global_max)) {
    global_sum += s[];
    if (s[] > global_max)
      global_max = s[];
  }
  fprintf (stderr, "global_sum/max: %f %f\n", global_sum, global_max);
}
