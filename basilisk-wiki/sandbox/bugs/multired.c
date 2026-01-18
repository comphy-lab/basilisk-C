/**
# Multiple successive reductions are incorrect in MPI

The code below works in serial and OpenMP but fails in MPI. */

int main()
{
  init_grid (4);
  double s = 0.;
  foreach(reduction(+:s))
    s += 1.;
#if 1 // this fails in MPI, but is OK in serial and OpenMP
  foreach(reduction(+:s))
    s += 2.;
#else // but this workaround is OK for all versions
  double s1 = 0.;
  foreach(reduction(+:s1))
    s1 += 2.;
  s += s1;
#endif
  assert (s == 48);
}
