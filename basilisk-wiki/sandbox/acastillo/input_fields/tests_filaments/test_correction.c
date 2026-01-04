/**
# Dipolar correction for the Batchelor vortex

In this example, we compute the dipolar correction for the Batchelor vortex 
using the approach of [Blanco-Rodriguez et al. (2015)](#blanco2015internal).

~~~pythonplot correction
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('results.txt', delimiter='\t', usecols=[0,2])
fig, ax1 = plt.subplots(ncols=1, sharex=False, figsize=(3.5,3.5))

ax1.plot(data[:,0], data[:,1]/data[:,0], label='$U_c=1.0$');
ax1.legend()
ax1.set_xlabel('$r$')
ax1.set_ylabel('$\\hat{\\psi}/r$')
ax1.set_xlim([0,10])
ax1.set_ylim([0,1.4])
plt.tight_layout()
plt.savefig('test_correction.svg')
~~~

*/

#include "grid/octree.h"
#include "view.h"
#include "acastillo/input_fields/filaments.h"
#include "acastillo/input_fields/draw_filaments.h"


int minlevel = 3;
int maxlevel = 5;

int main() {
  double rho;
  double a = 1.0;
  double U_c = 1.0;
  
  if (pid() == 0){
    FILE *file = fopen("results.txt", "w");
    for (rho = 0.0; rho <= 10.0; rho += 0.05) {
      integration_results results ;
      compute_A_with_derivatives(rho, a, U_c, &results);
      fprintf(file, "%.2f\t\t%.4f\t\t%.4f\t\t%.4f\n", 
          rho, 
          results.A, 
          results.A_p, 
          results.A_pp);
    }
    fclose(file);
  }
  return 0;
}

/**
# References

~~~bib

@article{blanco2015internal,
  title={Internal structure of vortex rings and helical vortices},
  author={Blanco-Rodr{\'\i}guez, Francisco J and Le Diz{\`e}s, St{\'e}phane and Sel{\c{c}}uk, Can and Delbende, Ivan and Rossi, Maurice},
  journal={Journal of Fluid Mechanics},
  volume={785},
  pages={219--247},
  year={2015},
  publisher={Cambridge University Press}
}

~~~
*/