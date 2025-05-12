/**
Filter function
*/

#include "grid/multigrid.h"
#include "run.h"

/**
The function below will replace the [default one](http://basilisk.fr/src/grid/multigrid-common.h#72).
*/

void waveletb (scalar s, scalar w)
{
  restriction ({s});
  for (int l = depth() - 1; l >= 0; l--) {
    foreach_coarse_level (l) {
      foreach_child()
        w[] = s[];
      s.prolongation (point, s);
      foreach_child() {
        double sp = s[];
        s[] = w[];
        /* difference between fine value and its prolongation */
        w[] -= sp;
      }
    }
    boundary_level ({w}, l + 1);
  }
  /* root cell */
  foreach_level(0) 
    w[] = s[];
  boundary_level ({w}, 0);
}

/**
Inverse wavelet transform */

void inverse_waveletb (scalar s, scalar w)
{
  foreach_level(0) 
    s[] = w[];
  boundary_level ({s}, 0);
  for (int l = 0; l <= depth() - 1; l++) {
    foreach_coarse_level (l) {
      s.prolongation (point, s);
      foreach_child()
        s[] += w[];
    }
    boundary_level ({s}, l + 1);
  }
}

/**
   Main loop */

int main() {
  N = 1 << 8;
  init_grid (N);

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

  foreach()
    s[] = sin(2*pi*x)*sin(10*pi*y);
  boundary({s});

  /**
  Wavelet transform */
  
  waveletb (s, w);
  inverse_waveletb (s2, w);

  /**
  Check that we recover the original signal (almost) exactly. */
  
  foreach()
    assert (fabs(s[] - s2[]) < 1e-10);
  
  char name[80];
  sprintf (name,"out.dat");
  FILE * fp = fopen (name, "w");
  output_field ({s, s2}, fp, linear=1);
  fclose(fp);
}

/**
We plot the filtered field and the filter
   
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
plt.contourf(x,y,p0)
plt.savefig('field.png')

plt.figure()
plt.contourf(x,y,p1)
plt.savefig('filter.png')
~~~
*/
