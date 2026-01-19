/**
# Sum reductions can be incorrect in MPI

The code below works in serial and OpenMP but fails in MPI. This is because 
the initial value of `s` is incorrectly summed over all processes in MPI. */

int main()
{
  init_grid (4);
  double s = 1;
  foreach(reduction(+:s));
  exit (s - 1);
}