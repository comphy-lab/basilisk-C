/** You need to run this code using MPI and see the output file for
  understanding the bug CC='mpicc -D_MPI=4' make my_adapt_bug.tst*/

// uf.x[right] = 0; // The automatic BC works again with this statement

face vector uf[];

int main(){
  init_grid(4);

  foreach_face(){
    uf.x[] = -1.;
  }
  boundary ((scalar *) {uf});
  foreach_face(){
    uf.x[] = pid();
  }

  foreach_face(x){
    // this should trigger the automatic boundary condition, but actually not
    double a = uf.x[1] + uf.x[-1] - 2*uf.x[];
    /* We output the values of face vector on the boudary between blocks. 
      Note that the values in the ghost cells, 
      such as uf[1] for rank 0 and 1, uf[-1] for rank 2 and 3, are incorrect. */ 
    if (fabs(x - 0.5) < 1.e-6)
    printf("Face vertor test: pid:%d, x:%g y:%g uf:%g uf[-1]:%g uf[1]:%g\n",
      pid(), x, y, uf.x[], uf.x[-1], uf.x[1]);
  }
  
  return 0;
}

/**
Output:

Face vertor test: pid:0, x:0.5 y:0.125 uf:0 uf[-1]:0 uf[1]:-1

Face vertor test: pid:0, x:0.5 y:0.375 uf:0 uf[-1]:0 uf[1]:-1

Face vertor test: pid:1, x:0.5 y:0.625 uf:1 uf[-1]:1 uf[1]:-1

Face vertor test: pid:1, x:0.5 y:0.875 uf:1 uf[-1]:1 uf[1]:-1

Face vertor test: pid:2, x:0.5 y:0.125 uf:2 uf[-1]:-1 uf[1]:2

Face vertor test: pid:2, x:0.5 y:0.375 uf:2 uf[-1]:-1 uf[1]:2

Face vertor test: pid:3, x:0.5 y:0.625 uf:3 uf[-1]:-1 uf[1]:3

Face vertor test: pid:3, x:0.5 y:0.875 uf:3 uf[-1]:-1 uf[1]:3


*/