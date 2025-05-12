/* You need to run this code using MPI and see the out file for understanding the bug CC='mpicc -D_MPI=4' make my_adapt_bug.tst**/
#define LEVEL 6

int main() {
  init_grid (1 << LEVEL);
  vertex scalar col[];

  unrefine((y > 0.4)  && (level >= LEVEL - 2));
  
  foreach_vertex()
    col[] = 1.;
  boundary({col});

  // refinement takes place inside the cells near block boundary inside
  // the computational domain
  refine((y < 0.5) && (y > 0.45) && (level < LEVEL - 1));

  // some values (on the block boundary) have been change after refinement, 
  // but all of them should equal to 1
  printf("\n");
  foreach_vertex(noauto, serial) {
    if (col[] < 1.) {
      printf("After adapt: x:%8g, y:%8g, col:%3g, pid:%3d, level:%2d\n", x, y, col[], pid(), level);
    }
  }
}