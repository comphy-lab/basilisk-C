/*
These function functions are a complement of "vtk.h", which allow to write
results from rectilinear and unstructured grids on legacy VTK format.
Works in multigrid with a square domain.
*/
void output_vtk_rectilinear (scalar * list, int n, FILE * fp, bool linear)
{
  fputs ("# vtk DataFile Version 2.0\n"
  "Basilisk\n"
  "ASCII\n"
  "DATASET RECTILINEAR_GRID\n", fp);
  fprintf (fp, "DIMENSIONS %d %d 1\n", n, n);
  fprintf (fp, "X_COORDINATES %d double\n", n);
  double fn = n;
  double Delta = L0/fn;
  for (int i = 0; i < n; i++) {
    double x = Delta*i + X0 + Delta/2.;
    fprintf (fp, "%g \n", x);
  }
  fprintf (fp, "Y_COORDINATES %d double\n", n);
  for (int j = 0; j < n; j++) {
    double y = Delta*j + Y0 + Delta/2.;
    fprintf (fp, "%g \n", y);
  }
  fprintf (fp, "Z_COORDINATES %d double\n", 1);
  fputs ("0. \n", fp);


  fprintf (fp, "POINT_DATA %d\n", n*n);
  for (scalar s in list) {
    fprintf (fp, "SCALARS %s double\n", s.name);
    fputs ("LOOKUP_TABLE default\n", fp);
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
          v = point.level >= 0 ? val(s) : nodata;
        }
        fprintf (fp, "%g\n", v);
      }
    }
  }
  fflush (fp);
}

void output_vtk_unstructured (scalar * list, int n, FILE * fp, bool linear)
{
  int count = 0;
  double fn = n;
  double Delta = L0/fn;
  double control[n+1][n+1];
  for (int i = 0; i < n+1; i++) {
    for (int j = 0; j < n+1; j++) {
      control[i][j] = count ;
      count+=1;
    }
  }

  fputs ("# vtk DataFile Version 2.0\n"
  "Basilisk\n"
  "ASCII\n"
  "DATASET UNSTRUCTURED_GRID\n", fp);
  fprintf (fp, "POINTS %d double\n", (n+1)*(n+1));
  for (int i = 0; i < n+1; i++) {
    double x = Delta*i + X0;
    for (int j = 0; j < n+1; j++) {
      double y = Delta*j + Y0;
      fprintf (fp, "%g %g 0\n", x, y);
    }
  }
  fputs ("\n", fp);
  fprintf (fp, "CELLS %d %d \n", n*n,n*n*5);
  for (int i = 0; i < n; i++) {
    double x = Delta*i + X0 + Delta/2.;
    for (int j = 0; j < n; j++) {
      double y = Delta*j + Y0 + Delta/2.;
      fprintf (fp, "4 %g %g %g %g\n", control[i][j], control[i+1][j],control[i][j+1],control[i+1][j+1]);
    }
  }
  fputs ("\n", fp);
  fprintf (fp, "CELL_TYPES %d \n", n*n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      fputs ("8 \n", fp);
  fputs ("\n", fp);
  fprintf (fp, "CELL_DATA %d\n", n*n);
  for (scalar s in list) {
    fprintf (fp, "SCALARS %s double\n", s.name);
    fputs ("LOOKUP_TABLE default\n", fp);
    for (int i = 0; i < n; i++) {
      double x = Delta*i + X0 + Delta/2.;
      for (int j = 0; j < n; j++) {
        double y = Delta*j + Y0 + Delta/2., v;
        if (linear)
        v = interpolate (s, x, y);
        else {
          Point point = locate (x, y);
          v = point.level >= 0 ? val(s) : nodata;
        }
        fprintf (fp, "%g\n", v);
      }
    }
  }
  fflush (fp);
}
