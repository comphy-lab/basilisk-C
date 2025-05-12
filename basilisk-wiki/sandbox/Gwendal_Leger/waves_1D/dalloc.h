double * dalloc (int N) {
  double * p = malloc (N * sizeof(double));
  int ii;
  for (ii = 0; ii < N; ii++)
    p[ii] = 0.;
  return p;
}
