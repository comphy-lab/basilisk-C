/**
# Upsample a $\mu$-HH dump file
For $\mu$-HH users it may be usefull to initialize a simulation with a solution from a previous simulation. However, if you want to increase the spatial resolution  between your old and new run you would run into troubles. This file uses the Basilisk octree-grid toolbox in order to upsample a solution to a higher level of refinement. For the moment it only works with grids that have $N^3$ grid cells where $N$ is a power of two. In the future it may be updated to work with $N\times N \times n_z$ grids as well. Special care is taken to make it somewhat convienent to use (see Usage chapter).

## The method
The method relies on the $\mu$-HH to Basilisk data structure header file which contains functions to read and write $\mu$-HH dump file to and from the Basilisk grid. 

*/

#include "grid/octree.h"
#include "utils.h"
#include "muhhtob.h"
/**
We allocate a scalar field and values for the top and bottom boundary condition.
*/
scalar a[];
double ab, at ;
/**
For maximum portability the function relies on command line arguments (see usage).
*/

int main(int argc, char **argv){
  /**
  ## Check the user provided input
 Most lines on this page are dedicated to set things up properly before we take the consecutive read upsample and dump actions.
 
 First we check if the user has provided a sensible number of command line arguments. The function aborts and returns an error and tip if this is not the case.
  */
  if (argc!=4 && argc!=6){
    printf("provide a plausible list of command-line arguments\n");
    printf("e.g.\t$./a.out 8 u.0032400 ud.0000000 0 6.5\n");
    return 1;
  } 
  /**
  Else, we reprint the CLA input for user statisfaction.
  */
  else{
    printf("Level of refinement in the original file data: %d\noriginal file name: %s\nnew file name: %s\n",atoi(argv[1]),argv[2],argv[3]);
  }
   /**
 ## Boundary conditions
 For the upsampling near the domain boundaries we need to fill the values of the ghost cells accoring to boundary conditions.We apply periodic conditions for the lateral directions. Furthermore, the user *may* set dirichlet-condition based values for the top and bottom boundary. Notice that for the boundaries without explictly defined boundary conditions the default Neumann condition will be applied (i.e. $\partial a/ \partial \overrightarrow{n}=0$). 
  */
  periodic(left);
  periodic(front);
  if (argc==6){
    ab = atof(argv[4]);
    at = atof(argv[5]);
    a[bottom]=dirichlet(ab);
    a[top]=dirichlet(at);
    printf("bottom bc: dirichlet(%g)\ntop bc: dirichlet(%g)\n",atof(argv[4]),atof(argv[5]));
  }
  int ilev=atoi(argv[1]);
  /**
  ## The actual reading and dumping
  The main code consists of four steps.  
  
  1. Read a $\mu$-HH file specified by *argv[2]*
  2. Apply the boundary conditions to set the ghost cell values 
  3. Refine the grid and solution using the defailt belinear interpolation
  4. Write the result to a $\mu$-HH compatible file
  
  */
  init_grid(1<<ilev);
  if(readmuhh_with_any_nr_of_threads(a,argv[2],ilev)!=0){ //read
    printf("Something went wrong on thread with pid()=%d\n",pid());
    return 1;
  }
  boundary({a}); //set boundary-ghost-cell values
  refine(level<=ilev); //upsample
  writemuhh_with_any_nr_of_threads(a,argv[3],ilev+1); //write
  return 0;
}
/**
## Usage
If you wish not to install basilisk, the easiest way to run the script is as follows; Download the portable version of the code [here](upsample.c) and for MPI enabled systems [here](MPIupsample.c). The content of this file nicely shows that 99% of the machinery behind this script is taken from Basilisk/src, and coded by Stephane Popinet. 

On a file with gcc installed, copy to file to the folder that contains your data and do something like:

~~~terminal
~/datapath$ gcc upsample.c -lm
~/datapath$ ./a.out 6 u.0032400 u.0000000 0 8
Level of refinement in the original file data: 6
original file name: u.0032400
new file name: u.0000000
bottom bc: dirichlet(0)
top bc: dirichlet(8)
Created a file with name u.0000000
~/datapath$
~~~

I expect you need to do the upsampling equivalent to the fftwplan.000000 and grid.000000 files aswell.
*/
  

