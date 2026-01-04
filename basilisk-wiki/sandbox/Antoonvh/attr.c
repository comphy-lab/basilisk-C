/**
# A function as a field attribute

These figures should look like flipped versions

![With attribute function](attr/h.png)

![Without attribute function](attr/f.png)
 */
#include "utils.h"

attribute {
  double (*fun)(double h, double k);
}

double my_fun (double a, double b) {
  return a + b;
}

scalar h[], f[];

int main() {
  h.fun = my_fun;
  init_grid (N);
  foreach() {
    for (scalar s in all) {
      if (s.fun) //check if attribute is set
	s[] = s.fun(x,y);
      else
	s[] = x - y;
    }
  }
  // inspect fields
  output_ppm (h, file = "h.png", min = 0, max = 2, n = 300);
  output_ppm (f, file = "f.png", min = -1, max = 1, n = 300);
}
