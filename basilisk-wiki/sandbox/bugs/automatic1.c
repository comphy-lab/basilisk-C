/**
# Checks that 'automatic()' works properly */

int main() {
  init_grid (1);
  scalar a = {0};
  for (int i = 0; i < 10; i++) {
    scalar s = automatic (a);
#if 1 // setting this to zero "fixes" the problem
    scalar p[];
#endif
  }
  fprintf (stderr, "%ld\n", datasize/sizeof (double));
  assert (datasize/sizeof (double) <= 2);
}
