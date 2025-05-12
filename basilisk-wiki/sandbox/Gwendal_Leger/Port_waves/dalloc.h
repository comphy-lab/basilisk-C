double * dalloc (int N) {
  double * p = malloc (N * sizeof(double));
  int ii;
  OMP_PARALLEL()
    {
      OMP(omp for schedule(static))
      for (ii = 0; ii < N; ii++)
	p[ii] = 0.;
    }// end parallel
  return p;
}
