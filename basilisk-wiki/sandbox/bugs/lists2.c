/**
   # Qcc and list iterator may fail

   This is another bug as [here](lists.c)

   Note that the bug is triggered by the usage of the `x` and `z`
   point variables.
*/

scalar a[],b[],c[], *abc = {a,b,c};

int main() {
  init_grid (N);
  foreach() {
    for (scalar s in abc) {
      //s[]; // <- Uncomment to fix
      double dummy = interpolate_linear (point, s, x, y, z);
    }
  }
}
