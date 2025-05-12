/**
We restore on the same level on which the [dump](dumper.c) was
produced. The code is run on the same number of processes using

~~~bash
CC='mpicc -D_MPI=27' make restorer.tst
~~~
*/

#include "grid/multigrid3D.h"
#include "utils.h"

int main()
{
  restore ("../dumper/test");
}

/**
The error messages are:

~~~bash
[restorer.tst]
restorer: /home/popinet/basilisk_1_0/src/output.h:1104: restore: Assertion `0' failed.
[aldebaran:08664] *** Process received signal ***
[aldebaran:08664] Signal: Aborted (6)
...
grid depths do not match. Aborting.
...
~~~
*/
