/**
# Incompatible arguments in _stencil_function

This test reports a compilation issue for functions with pre-defined inputs.

The error disappears when:

1. We provide all inputs to the function call: `foo (point, f, true);`
1. We do not use the scalar `f` in a function: commenting `test()` it works.
*/

double test (double val) {
  return val;
}

double foo (Point point, scalar f, bool vof = true) {
  test (f[]);
  return 1.;
}

scalar f[];

int main (void) {
  init_grid (1 << 5);
  foreach()
    foo (point, f);
}
