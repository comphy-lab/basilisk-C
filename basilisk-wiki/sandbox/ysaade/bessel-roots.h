void zeros (int n, double *z) {
  double x, y;
  for (int i = 0; i < n; i++) {
    x = (i + 3./4.)*pi;
    y = j0(x);
    while (fabs(y) > 1e-10) {
      double a = -j1(x);
      double b = y - a*x;
      x = -b/a;
      y =  j0(x);
    }
    *(z + i) = x;
  }
  return;
}
