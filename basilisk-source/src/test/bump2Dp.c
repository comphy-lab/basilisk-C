// very close to bump2D.c but with additional checks for MPI parallelism
// also used to test parallel dump/restore

#include "saint-venant.h"
#include "check_restriction.h"

#define BGHOSTS 2

int LEVEL = 6;

int main (int argc, char * argv[])
{
  if (argc > 1)
    LEVEL = atoi (argv[1]);
  origin (-0.5, -0.5);
  init_grid (1 << LEVEL);
  mpi.min = 1; // 1 element per process minimum
  mpi.leaves = true; // balance leaves only
  run();
}

event init (i = 0)
{
  double a = 1., b = 200.;
  if (!restore (file = "bump2Dp-restore.dump"))
    foreach()
      h[] = 0.1 + a*exp(- b*(x*x + y*y));
}

event logfile (i++) {
  stats s = statsf (h);
  fprintf (stderr, "%g %d %g %g %.8f\n", t, i, s.min, s.max, s.sum);
}

event outputfile (t <= 2.5; t += 2.5/8);

event image(i++)
{
  scalar pid[];
  foreach()
    pid[] = pid();
  static FILE * fp = fopen ("pid", "w");
  output_ppm (pid, fp, min = 0, max = npe() - 1, map = randomap);
  foreach()
    pid[] = level;
  static FILE * fp1 = fopen ("level", "w");
  output_ppm (pid, fp1, min = 0, max = LEVEL);
}

event adapt (i++) {

#if BGHOSTS == 2
  scalar s[];
  face vector u[];
  reset ({s,u}, undefined);
  foreach()
    s[] = 1;
  foreach_face()
    u.x[] = 1;
#endif
  
  astats st = adapt_wavelet ({h}, (double[]){1e-2}, LEVEL);
  fprintf (stderr, "# refined %d cells, coarsened %d cells\n", st.nf, st.nc);
  restriction ({zb}); // fixme: why is it necessary with MPI?

#if BGHOSTS == 2
  foreach()
    foreach_neighbor()
    assert ((s[] == 1));
  check_restriction (s);
  foreach_face()
    for (int i = -2; i <= 2; i++)
      assert ((u.x[0,i] == 1));
#endif
}

event snapshot (t = 2.5/2.) {
  dump (file = "dump");
}
