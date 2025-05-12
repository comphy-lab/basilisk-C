/**
# Reductions for user-defined macros

When the code below attempts an automatic reduction on the `iterator`
loop via the macro syntax, it does not work. I noticed that generated
`-source` file does not add the relevant MPI lines to the user-defined
macro, which it does for `foreach()`. See also the [output](reductions/out).
*/
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
