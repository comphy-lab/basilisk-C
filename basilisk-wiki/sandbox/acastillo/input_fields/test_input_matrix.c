/**
# Testing input_matrix()

In this example, we use a matlab routine [test_input_matrix.m](test_input_matrix.m)
to generate some initial conditions 

~~~ {.matlab}
a=0.1
n=256
x = linspace(-0.5,0.5,n);
y = linspace(-0.5,0.5,n);
[X,Y] = meshgrid(x,y);
f=((X.^2+Y.^2).^2 + 4*a*X.*(X.^2+Y.^2) - 4*a^2*Y.^2<=0);
output_matrix(f,n,x,y,'test_input_matrix.bin')
~~~
and this routine [output_matrix.m](output_matrix.m)
to write the results in a format compatible to output_matrix()/input_matrix().


In this case, we take the example from the tutorial (see
[bump.c](http://basilisk.fr/Tutorial)) but instead of the Gaussian bump, we use
an arbitrary condition (a Cardioid) produced by the matlab routine.

*/
#include "saint-venant.h"
#include "auxiliar_input.h"

event init (t = 0) {
  FILE * fp = fopen("test_input_matrix.bin", "r");
  if (!fp) printf("Binary file not found");
  input_matrix(h,fp,N,X0,Y0,L0);
  fclose (fp);
}

event end (i = 10) {
  printf ("i = %d t = %g\n", i, t);
}

int main() {
  init_grid (256);
  origin (0, 0);
  run();
}

event images (i++) {
  output_ppm (h);
}

event end (i = 300){
}

/**

## Output

If we compile and execute this code

~~~bash
qcc -O2 -Wall test_input_matrix.c -o test_input_matrix -lm
./test_input_matrix > test_input_matrix.ppm
animate test_input_matrix.ppm
~~~

we should get something like this

![Bump test case using initial conditions from a file](test_input_matrix.gif)
*/
