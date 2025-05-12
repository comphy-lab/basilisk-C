/**
# `locate() + foreach_dimension()` $\rightarrow$ `input in flex-scanner failed` 

`qcc` gives an error if `foreach_dimension(){...}` is called.

Workarounds include;

1. Use a (dummy) code block
2. Do not use a block with `foreach_dimension()` (for one-line blocks)
 */
#ifndef DUMMY_CODE_BLOCK
#define DUMMY_CODE_BLOCK (0)
#endif

int main() {
  init_grid (N);
  Point point = locate (0.5, 0.6);

#if DUMMY_CODE_BLOCK
  {}           
#endif

  coord a = {x, y, z};
  foreach_dimension() {
      printf ("%g\n", a.x);
      a.x = 0;
  }
}
