
/**
# Test for vertex boundary conditions


In the present implementation of Basilisk (08/2022), vertex fields have no boundary conditions. Antoon's proposal is to add the same boundary conditions as face vectors.

Here are the stencil for face and vertex fields


![Stencils in Basilisk](/src/figures/stencil_face.svg)

Assuming the vertex cell (right) is the bottom left cell of the domain, we want $\omega[0,0]$ $\omega[1,0]$, and $\omega[0,1]$ to be BC points and $\omega[1,1]$ to be an interior point.

In the exemple below, we set all interior points to 1 and boundary points to zero.

*/


#include "grid/multigrid.h"
#include "run.h"


scalar omega[];

int main() {
  N = 4;
  init_grid (N);

  omega[left] = 0.;
  omega[right] = 0.;
  omega[top] = 0.;
  omega[bottom] = 0.;

  foreach_vertex() {
    omega[] = 1.0;
  }

  boundary({omega});
  foreach_vertex() {
    printf("%g \t %g \t %g\n",x,y,omega[]);
  }

}



/**
With the present implementation of basilisk, this code returns
(first column is x, second is y, third is $\omega$)


	0 	 0 	 1
	0 	 0.25 	 1
	0 	 0.5 	 1
	0 	 0.75 	 1
	0 	 1 	 1
	0.25	 0 	 1
	0.25	 0.25 	 1
	0.25	 0.5 	 1
	0.25	 0.75 	 1
	0.25	 1 	 1
	0.5 	 0 	 1
	0.5 	 0.25 	 1
	0.5 	 0.5 	 1
	0.5 	 0.75 	 1
	0.5 	 1 	 1
	0.75	 0 	 1
	0.75	 0.25 	 1
	0.75	 0.5 	 1
	0.75	 0.75 	 1
	0.75	 1 	 1
	1 	 0 	 1
	1 	 0.25 	 1
	1 	 0.5 	 1
	1 	 0.75 	 1
	1 	 1 	 1

~~~gnuplot
unset key
set size ratio -1
set xrange [-0.1:1.1]
set yrange [-0.1:1.1]
set xtics 0,0.25,1
set ytics 0,0.25,1
set grid
plot 'out' u 1:2:3 w labels
~~~

*/

/**

In a preliminary proposal:

[/sandbox/bderembl/qg-node/basilisk_bc_vertex.patch]()

we were able to change the boundary conditions for NON-MPI codes: here is the output of the same code above when we apply the patch:

	0 	 0 	 0
	0 	 0.25 	 0
	0 	 0.5 	 0
	0 	 0.75 	 0
	0 	 1 	 0
	0.25	 0 	 0
	0.25	 0.25 	 1
	0.25	 0.5 	 1
	0.25	 0.75 	 1
	0.25	 1 	 0
	0.5 	 0 	 0
	0.5 	 0.25 	 1
	0.5 	 0.5 	 1
	0.5 	 0.75 	 1
	0.5 	 1 	 0
	0.75	 0 	 0
	0.75	 0.25 	 1
	0.75	 0.5 	 1
	0.75	 0.75 	 1
	0.75	 1 	 0
	1 	 0 	 0
	1 	 0.25 	 0
	1 	 0.5 	 0
	1 	 0.75 	 0
	1 	 1 	 0


We see that now all boundary points are zero and all interior points are 1.

*/

/**

Now if we compile the code with

CC99='mpicc -std=c99' qcc -D_MPI=1 -grid=multigrid -lm -O3 vertex.c

and run it on 4 processors:
mpirun -np 4 a.out

we get

	0 	 0 	 0
	0 	 0.25 	 0
	0 	 0.5 	 0
	0.25	 0 	 0
	0.25	 0.25 	 1
	0.25	 0.5 	 0
	0.5 	 0 	 0
	0.5 	 0.25 	 0
	0.5 	 0.5 	 0
	0 	 0.5 	 0
	0 	 0.75 	 0
	0 	 1 	 0
	0.25	 0.5 	 0
	0.25	 0.75 	 1
	0.25	 1 	 0
	0.5 	 0.5 	 0
	0.5 	 0.75 	 0
	0.5 	 1 	 0
	0.5 	 0 	 0
	0.5 	 0.25 	 0
	0.5 	 0.5 	 0
	0.75	 0 	 0
	0.75	 0.25 	 1
	0.75	 0.5 	 0
	1 	 0 	 0
	1 	 0.25 	 0
	1 	 0.5 	 0
	0.5 	 0.5 	 0
	0.5 	 0.75 	 0
	0.5 	 1 	 0
	0.75	 0.5 	 0
	0.75	 0.75 	 1
	0.75	 1 	 0
	1 	 0.5 	 0
	1 	 0.75 	 0
	1 	 1 	 0


i.e. all the MPI boundary points are also set to zero.
*/
