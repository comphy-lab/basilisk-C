/**

# Test case for Vertex Poisson solver

We want to solve the poisson equation

$$
\nabla^2 \psi = \omega
$$

with $\omega$ a known 2d field and $\psi$ the unknown. Both $\psi$ and $\omega$
are vertex fields. For this test case, we first set $\psi$ to a random field,
compute $\omega$ and see if we can recover $\psi$ with the poisson solver.

This test case should work in MPI.

We use a dedicated NetCDF routine to output vertex fields.


~~~pythonplot Numerical solution
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf_file

#plt.ion()

f = netcdf_file('vars.nc')
psi = f.variables['psi'][:,:].copy().squeeze()
f.close()

plt.figure()
plt.imshow(psi,origin='lower')
plt.colorbar()
plt.savefig('psi.svg')
~~~

~~~pythonplot Analytical solution

f = netcdf_file('vars.nc')
exact = f.variables['exact'][:,:].copy().squeeze()
f.close()

plt.figure()
plt.imshow(exact,origin='lower')
plt.colorbar()
plt.savefig('exact.svg')
~~~

~~~pythonplot Difference Analytical - Numerical

plt.figure()
plt.imshow(exact-psi,origin='lower')
plt.colorbar()
plt.savefig('diff.svg')
~~~


*/



#include "grid/multigrid.h"
#define LAYERS 1
#include "run.h"
#include "bderembl/libs/nodal-poisson.h"
#include "bderembl/libs/netcdf_vertex_bas.h"

char* fileout = "vars.nc";

vertex scalar psi[]; 
vertex scalar omega[]; 
vertex scalar diff[]; 
vertex scalar exact[];

#define laplacian(p) (p[1] + p[-1] + p[0,1] + p[0,-1] - 4*p[])/(sq(Delta))

int main() {
  init_grid (64);

  psi[left]   = 0;
  psi[right]  = 0;
  psi[top]    = 0;
  psi[bottom] = 0;
  psi.restriction = restriction_vert;
  psi.prolongation = refine_vert;


  omega[left]   = 0;
  omega[right]  = 0;
  omega[top]    = 0;
  omega[bottom] = 0;
  omega.restriction = restriction_vert;
  omega.prolongation = refine_vert;

  exact[left]   = 0;
  exact[right]  = 0;
  exact[top]    = 0;
  exact[bottom] = 0;
  exact.restriction = restriction_vert;
  exact.prolongation = refine_vert;


  foreach_vertex(){
    psi[] = 0.;
//    omega[] = sin(pi*x)*sin(pi*y);
    exact[] = noise();
  }

  boundary ({exact});

  foreach_vertex(){
    omega[] = laplacian(exact);
  }

  boundary ({psi,omega}); // automatic!

  TOLERANCE=1e-6;
  mgstats mg = vpoisson (psi,omega);

  printf("Cell poisson stats:\n"
	 "mg.i: %d mg.relax: %d mg.resa: %g mg.resb: %g\n",
	 mg.i, mg.nrelax, mg.resa, mg.resb);

  foreach_vertex(){
//    exact[] = -sin(pi*x)*sin(pi*y)/(2*sq(pi));
    diff[] = exact[] - psi[];
  }

/**
   Write netcdf file because: 
   
- ultimately, will need netcdf for the qg model

- no obvious way to write vertex field)

 */

  sprintf (file_nc,"%s", fileout);
  scalar_list_nc = list_copy({psi,exact, diff});
  create_nc();
  write_nc();
  free(scalar_list_nc);
}
