/**

# Qcc and list iterator may fail

*/
scalar a[], b[], c[], * abc = {a, b, c};;

int main() {
  init_grid (N);
  foreach() { //adding `noauto` is ok. 
    scalar s;
    for (s in abc) {
      printf ("%s", s.name);
    }
  }
}
/**
## See similar old bugs

* [v1](lists.c)
* [v2](lists.c)

 */
