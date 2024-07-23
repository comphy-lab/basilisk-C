#define BGHOSTS 2

scalar s[];

int main (int argc, char * argv[])
{
  X0 = Y0 = -0.5;
  init_grid (argc > 1 ? atoi(argv[1]) : 32);

  unrefine (sq(x - 0.1) + sq(y - 0.1) > sq(0.1));
  
  output_cells (stdout);
  
  foreach()
    s[] = 1.;

  foreach()
    foreach_neighbor()
      assert (s[] == 1.);

  // rebalancing
  int nf = 0;
  foreach(reduction(+:nf))
    nf++;
  int npe;
  MPI_Comm_size (MPI_COMM_WORLD, &npe);
  nf = max(1, nf/npe);
  scalar index[];
  z_indexing (index, true);
  foreach()
    fprintf (stderr, "%g %g %d\n", x, y, min(npe - 1, (int)(index[]/nf)));
}
