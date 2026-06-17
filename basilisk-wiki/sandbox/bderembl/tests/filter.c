/**
Filter function
*/

#include "run.h"


/**
   Main loop */

int main() {
  N = 1024;
  init_grid (N);
  origin(0., 0.);
  L0 = 1;

  scalar s[], w[], s2[];

  s[top] = 0;
  s[bottom] = 0;
  s[right] = 0;
  s[left] = 0;

  w[top] = 0;
  w[bottom] = 0;
  w[right] = 0;
  w[left] = 0;

  s2[top] = 0;
  s2[bottom] = 0;
  s2[right] = 0;
  s2[left] = 0;


/**
   Initial field: random noise
*/
  foreach()
    s[] = noise();
//    s[] = 0.1*sin(2*pi*x/L0)*sin(4*pi*y/L0) + noise();

/**
   Define cutoff legngth scale at all levels

   Lf1 is the cutoff length scale at the bottom of the domain
   Lf2 is the cutoff length scale at the top of the domain

   sig_len is the cutoff length scale (field)
   sig_lev is the weight of the iWT at each level

*/

  double Lf1 = 0.2;
  double Lf2 = 0.02;

  scalar sig_lev[];
  scalar sig_len[];

  foreach()
    sig_len[] = Lf1 + y/L0*(Lf2 - Lf1);

  restriction({sig_len});

  // low pass filter
  for (int l = depth(); l >= 0; l--) {
    foreach_level (l) {
      double ref_flag = 0;
      if (l < depth())
        foreach_child()
          ref_flag += sig_lev[];
      if (ref_flag > 0)
        sig_lev[] = 1;
      else
        if (sig_len[] > 2*Delta)
          sig_lev[] = 0;
        else if (sig_len[] <= 2*Delta && sig_len[] > Delta)
          sig_lev[] = 1-(sig_len[]-Delta)/Delta;
        else
          sig_lev[] = 1;
    }
    boundary_level ({sig_lev}, l);
}

/**
   3 steps:

  1. Wavelet transform
  2. Filer wavelet coefs
  3. inverse wavelet transform


*/

  wavelet (s, w);

  for (int l = 0; l <= depth(); l++) {
    foreach_level (l){
      w[] *= sig_lev[];
    }
    boundary_level ({w}, l);
  }
  inverse_wavelet (s2, w);


  char name[80];
  sprintf (name,"out.dat");
  FILE * fp = fopen (name, "w");
  output_field ({s, s2}, fp, linear=1);
  fclose(fp);
}

/**
We plot the original field

~~~pythonplot

import numpy as np
import matplotlib.pyplot as plt

file1 = 'out.dat'

data = np.loadtxt(file1)

n1, n2 = data.shape

N = int(np.sqrt(n1))
x = data[:,0].reshape((N,N))
y = data[:,1].reshape((N,N))
p0 = data[:,2].reshape((N,N))
p1 = data[:,3].reshape((N,N))

plt.figure()
plt.pcolormesh(x,y,p0)
plt.savefig('field.png')
~~~

and the filtered field
~~~pythonplot

plt.figure()
plt.pcolormesh(x,y,p1)
plt.savefig('filter.png')
~~~
*/
