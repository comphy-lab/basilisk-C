/**
This is a 3D modification of the orginal output_field() function which outputs multiple fields to uniform grid, written by Shine Sun.
Different from the orginal one, the grid dimension is changed to $N*N*N$.

The arguments and their default values are:

*list*
: list of fields to output. Default is *all*.

*fp*
: file pointer. Default is *stdout*.

*n*
: number of points along each dimension. Default is *N*.

*linear*
: use first-order (default) or bilinear interpolation. 

*box*
: the lower-left and upper-right coordinates of the domain to consider.
 Default is the entire domain. */

void * matrix_new_3d (int nx, int ny, int nz, size_t size) {
  void ** m = qmalloc(nx, void *);
  char * a = qmalloc(nx*ny*nz*size, char);
  for (int i = 0; i < nx; i++)
    m[i] = a + i*ny*nz*size;
  return m;
}

trace
void output_field_3d (scalar * list = all,
		   FILE * fp = stdout,
		   int n = N,
		   bool linear = false,
		   coord box[2] = {{X0, Y0, Z0},{X0 + L0, Y0 + L0, Z0 + L0}})
{
  n++;
  int len = list_len (list);
  double Delta = 0.999999*(box[1].x - box[0].x)/n;
  int ny = (box[1].y - box[0].y)/Delta;
  int nz = (box[1].z - box[0].z)/Delta;
  // fixme: upgrade to foreach_region
  double ** field = (double **) matrix_new_3d (n, ny, nz, len*sizeof(double));
  for (int i = 0; i < n; i++) {
    double xp = Delta*(i + 0.5) + box[0].x;  //the approximate cell centers
    for (int j = 0; j < ny; j++) {
      double yp = Delta*(j + 0.5) + box[0].y;
      for (int k=0; k < nz; k++) {
        double zp = Delta*(k + 0.5)+ box[0].z;
        if (linear) {
          int ind = 0;
          for (scalar s in list)
            field[i][len*j*nz + len*k + ind++] = interpolate (s, xp, yp, zp);
        }
        else {
          Point point = locate (xp, yp, zp);
          int ind = 0;
          for (scalar s in list)
            field[i][len*j*nz + len*k + ind++] = point.level >= 0 ? s[] : nodata;
        }
      }
    }
  }
  
  if (pid() == 0) {
    fprintf (fp, "# 1:x 2:y 3:z");
    int i = 4;
    for (scalar s in list)
      fprintf (fp, " %d:%s", i++, s.name);
    fputc('\n', fp);
    for (int i = 0; i < n; i++) {
      double x = Delta*i + box[0].x;
      for (int j = 0; j < ny; j++) {
	      double y = Delta*j + box[0].y;
        for (int k = 0; k < nz; k++) {
          double z = Delta*k + box[0].z;
          fprintf (fp, "%g %g %g", x, y, z);
          int ind = 0;
          for (scalar s in list)
            fprintf (fp, " %g", field[i][len*j*nz + len*k + ind++]);
          fputc ('\n', fp);
          }
        }
      }
    fflush (fp);
  }
  matrix_free (field);
}