/**
# Call to `boundary()` with 8 threads causes a crash when running Multigrid MPI 
The following code crashes for some levels of refinement when using 8 threads and 9 levels of refinement.
*/

#include "grid/multigrid3D.h"

scalar m[];
int maxlevel = 5;

int main(int argc,char ** argv){
  if (argc>1)
    maxlevel=atoi(argv[1]);
  init_grid(1<<maxlevel);
  printf("thread #%d was here\n",pid());
  boundary({m});
  printf("thread #%d is done\n",pid());
}
/**
## How to crash?
Compile and run the code following like so;

~~~bash
:~$ CC99='mpicc -std=c99' qcc -O2 -Wall -D_MPI=1 MGmpiboundary.c -o maybebug -lm
:~$ mpirun -np 8 ./maybebug 
~~~
There are no issues for this case. However when running with nine levels of refinement;

~~~bash
:~$ mpirun -np 8 ./maybebug 9
~~~
An error is raised: `segmentation fault` that seems to correspond with the call the `boundary()` function as one may get this:

~~~terminal
thread #1 was here
thread #6 was here
thread #3 was here
thread #0 was here
thread #4 was here
thread #2 was here
thread #5 was here
thread #7 was here
--------------------------------------------------------------------------
mpirun noticed that process rank 6 with PID 14254 on node antoon-XPS-15-9550 exited on signal 11 (Segmentation fault).
--------------------------------------------------------------------------
~~~

The bug appears on a test laptop and also on a bigger system (a workstation) that allowed testing with 10 levels of refinement. Here some manupulation of the data associated with the `m[]` field worked fine, untill the call to `boundary()` was made.  

Note that running the executable using a single or 64 threads does not display this issue, i.e. for 9 levels of refinement.  
*/