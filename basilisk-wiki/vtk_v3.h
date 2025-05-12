/**
# VTK Structured grid version 3.0

Generates output to VTK ASCII version 3.0 file format. */
void output_vtk_v3 (scalar * list, int n, FILE * fp, bool linear)
{
  fputs ("# vtk DataFile Version 3.0\n"
         "Basilisk\n"
         "ASCII\n"
         "DATASET STRUCTURED_GRID\n", fp);
  fprintf (fp, "DIMENSIONS %d %d 1\n", n, n);
  fprintf (fp, "POINTS %d double\n", n*n);

  double fn = n;
  double Delta = L0/fn;
  for (int i = 0; i < n; i++) {
    double x = Delta*i + X0 + Delta/2.;
    for (int j = 0; j < n; j++) {
      double y = Delta*j + Y0 + Delta/2.;
      fprintf (fp, "%g %g 0\n", x, y);
    }
  }
  fprintf (fp, "POINT_DATA %d\n", n*n);

  /*
   * todo - Should be more clever way to get the size of the scalar list
   */
  int list_sz=0;
  for(scalar s in list){list_sz++;}

  fprintf (fp, "FIELD FielData %d\n", list_sz);
  for (scalar s in list) {
    fprintf (fp, "%s 1 %d double\n", s.name, n*n );
    double fn = n;
    double Delta = L0/fn;
    for (int i = 0; i < n; i++) {
      double x = Delta*i + X0 + Delta/2.;
      for (int j = 0; j < n; j++) {
        double y = Delta*j + Y0 + Delta/2., v;
        if (linear)
          v = interpolate (s, x, y);
        else {
          Point point = locate (x, y);
          v = point.level >= 0 ? val(s,0,0) : nodata;
        }
        fprintf (fp, "%g\n", v);
      }
    }
  }
  fflush (fp);
}
