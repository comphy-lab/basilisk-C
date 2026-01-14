/**
The bug seems to be related to multigrid.h + axi.h.

I tried to run example cases using this combination without success. 

It always returns some MPI error message of the form:
```
    An error occurred in MPI_Recv
    on communicator MPI_COMM_WORLD
    MPI_ERR_TRUNCATE: message truncated
    MPI_ERRORS_ARE_FATAL (processes in this communicator will now abort,
       and potentially your MPI job)
```

Here is a minimal example to show the problem.

If you run it on a single processor, it works:
CC='mpicc -D_MPI=1' make axi_mg_test.tst

If you try with 4, it crashes:
CC='mpicc -D_MPI=4' make axi_mg_test.tst */

#include "grid/multigrid.h"
#include "axi.h"
#include "navier-stokes/centered.h"

int main()
{
  N = 32;
  DT = 0.1 [0,1];
  run();
}

event stop (t = 1.0) {
  return 1;
}