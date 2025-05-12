void check_restriction (scalar a)
{
  double val = 123;
  int maxlevel = depth();
  mpi_all_reduce (maxlevel, MPI_INT, MPI_MAX);
  for (int l = 0; l <= maxlevel; l++) {
    if (l == 0)
      foreach_level_or_leaf (l)
	a[] = val;
    else
      foreach_level (l) {
	for (int i = -2; i <= 2; i++)
	  for (int j = -2; j <= 2; j++)
	    // fixme: boundary conditions should work too
	    assert ((aparent(i,j).pid < 0) || coarse(a,i,j) == val);
	a[] = val;
      }
    boundary_level ({a}, l);
  }
}
