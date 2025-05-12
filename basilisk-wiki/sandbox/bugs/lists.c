/**
# Qcc and list iterator may fail

*/
scalar a[],b[],c[], *abc = {a,b,c};

int main() {
  init_grid (N);
  foreach() {
    for (scalar s in abc) {
      //s[]; //<-- Uncomment to fix
      printf ("%s", s.name);
    }
  }
}
