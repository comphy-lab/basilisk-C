/** 
# Default boundary conditions do not work for face vector fields in MPI

You need to run this code using MPI and see the output file for
understanding the bug CC='mpicc -D_MPI=4' make bc_face_vector.tst*/

// uf.x[right] = 0; // including this statement changes the results

face vector uf[];

int main()
{
  init_grid (2);
  reset ((scalar *){uf}, -1);
  foreach_face()
    uf.x[] = pid();

  foreach_face (x) {
    // this correctly triggers the automatic BCs (see the output below)
    double a = uf.x[1] + uf.x[-1] - 2*uf.x[];
    
    /** 
    We output the values of face vector on the boudary between blocks.
    Note that the values in the ghost cells, such as uf[1] for rank 0
    and 1, uf[-1] for rank 2 and 3, are incorrect. */
    
    fprintf (qerr, "pid: %d\tx: %g\ty: %g\tuf: %g\tuf[-1]: %g\tuf[1]: %g\ta: %g\n",
	     pid(), x, y, uf.x[], uf.x[-1], uf.x[1], a);
  }
  fflush (qerr);
  system ("cat log log-* > log.txt");
}

/**
This is the output:
![](bc_face_vector/log.txt){width="100%" height="250px"}
*/
