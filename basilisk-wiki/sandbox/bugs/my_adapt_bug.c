/* You need to run this code using MPI and see the output file for understanding the bug CC='mpicc -D_MPI=4' make my_adapt_bug.tst**/
int main() {
  init_grid (1 << 4);
  vertex scalar col[];
  
  foreach_vertex()
    col[] = 1.;

  /**
  Refinement occurs in cells adjacent to the block boundary, even though
  they are entirely within the computational domain. */
  refine((y < 0.5) && (y > 0.45) && (level < 5));

  /**
  After refinement, some values on the block boundary are incorrectly
  reset to zero, although they should all remain equal to one. */
  foreach_vertex(serial)
    if (col[] < 1.)
      fprintf(qerr, "After refine: x:%8g, y:%8g, col:%3g, pid:%3d, level:%2d\n",
       x, y, col[], pid(), level);

  fflush (qerr);
  system ("cat log log-* > log.txt");
}

/**
This is the output (the exact output may differ from run to run):
![](my_adapt_bug/log.txt){width="100%" height="250px"}
*/