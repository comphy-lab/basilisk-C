/**
# output_gfs() swap the component of tensors and vectors.

To plot the components shown in the left side you have to set in
gfsview the right side:

* S_x_x = -S_x_y
* S_x_y =  S_x_x
* S_y_x = -S_y_y
* S_y_y =  S_y_x
* u_x = -u_y
* u_y = u_x

*/

#include "run.h"

tensor S[];
vector u[];

int main() {
  init_grid(2);
  run();
}

event initialize (i = 0) {
  foreach() {
    S.x.x[] = 1;
    S.x.y[] = 2;
    S.y.x[] = 3;
    S.y.y[] = 4;
    u.x[] = 5;
    u.y[] = 6;
  }
}

event gfsview (i = 0) {
  FILE * fp = stdout;
  output_gfs (fp, t=t);
}
