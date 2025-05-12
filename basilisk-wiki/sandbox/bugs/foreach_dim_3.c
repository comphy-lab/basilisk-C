/**
# `foreach_dimension(3)` bug

I expect `foreach_dimension(3)` to rotate twice and do `x->y->z->x`, but it doesn't.
*/
int main() {
  coord n = {1, 2, 3};
  double last = 0, last2 = 0;
  foreach_dimension(3) {
    last = n.x;
    last2 = n.y;
  }
  
  /**
  The assertion below shoud fail: we are 2D so the third dimension is just
  ignored. This is systematic in Basilisk and allows for 3D solvers to
  automatically "degenerate" to 2D. */
  
  assert (last == n.z || last2 == n.z);
}