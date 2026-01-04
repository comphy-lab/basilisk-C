/**
# Reductions for user-defined macros

Check that (MPI) reductions are correctly applied for user macros. */

macro iterator (int start, int end, int index, Reduce reductions = None) {
  for (int index = start; index <= end; index++)
    {...}
}

int main() {
  int total = 0, end = 16;
  iterator (1, end, i, reduction(+:total))
    total++;
  if (pid() == 0)
    printf ("%d %d\n", total, npe()*end);
  int total2 = 0;
  init_grid (16);
  foreach(reduction(+:total2)) 
    total2++;
  if (pid() == 0)
    printf ("%d %d %d\n", total2, sq(N), npe());
  assert (total2 == sq(N));
  assert (total == npe()*end);
}
