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

struct OutputField_3d {
  scalar * list;
  FILE * fp;
  int n;
  bool linear;
  double box[2][3];  //The lower-left-back and upper-right-front points of the box to be outputed.
};

void * matrix_new_3d (int nx, int ny, int nz, size_t size) {
  void ** m=qmalloc(nx, void *);  //Define a pointer that points to every x coordinate.
  char * a=qmalloc(nx*ny*nz*size, char);
  for (int i=0; i<nx; i++) {
    m[i]=a+i*ny*nz*size;
  }
  return m;
}

trace
void output_field_3d (struct OutputField_3d p)
{
  if (!p.list) p.list = all;  //The variables to be outputed.
  if (p.n == 0) p.n = N;  //changed from N+1.
  if (!p.fp) p.fp = stdout;  //The default output stream is the screen.
  if (p.box[0][0] == 0. && p.box[0][1] == 0. && p.box[0][2]==0.0 &&
     p.box[1][0] == 0. && p.box[1][1] == 0. && p.box[1][2]==0.0) {
    p.box[0][0] = X0;
    p.box[0][1] = Y0;
    p.box[0][2]=Z0;
    p.box[1][0] = X0 + L0;
    p.box[1][1] = Y0 + L0;
    p.box[1][2]=Z0+L0;
  }  //The output box domain
  
  int len = list_len(p.list);  //number of variables to be outputed
  double Delta = 0.999999*(p.box[1][0] - p.box[0][0])/(p.n);  //The output grid spacing. The "0.999999" is used to avoid the interpolated point coinciding with fine grid boundary, which may cause point.level<0.
  int nx = p.n;  //The number of interpolation points in each dimension.
  int ny = p.n;
  int nz = p.n;
  double ** field = (double **) matrix_new_3d (nx, ny, nz, len*sizeof(double));  //create the output matrix.
  for (int i = 0; i < p.n; i++) {
    double xp = Delta*(i+0.5) + p.box[0][0];  //the approximate cell centers
    for (int j = 0; j < ny; j++) {
      double yp = Delta*(j+0.5) + p.box[0][1];
      for (int k=0; k<nz; k++) {
        double zp=Delta*(k+0.5)+p.box[0][2];
        if (p.linear) {
          int ind = 0;
          for (scalar s in p.list)
            field[i][len*j*nz + len*k + ind++] = interpolate (s, xp, yp, zp);
        }  //store the variables in order, every k jumps len variables and every j jumps len*nz variables.
        else {
          Point point = locate (xp, yp, zp);
          int ind = 0;
          for (scalar s in p.list)
            field[i][len*j*nz + len*k + ind++] = point.level >= 0 ? s[] : nodata;
        }
      }
    }
  }

  if (pid() == 0) { // output in the master thread
@if _MPI
    MPI_Reduce (MPI_IN_PLACE, field[0], len*p.n*ny*nz, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);  //Every thread is involved in the interpolation computation of every point, however only one thread contains the meanful data because of the domain decomposition of MPI. The data on other thread are 1e+30, so we get the minmum value from all threads.
@endif
    fprintf (p.fp, "# 1:x 2:y 3:z");  //header for gnuplot or so.
    int i = 4;
    for (scalar s in p.list)
      fprintf (p.fp, " %d:%s", i++, s.name);
    fputc('\n', p.fp);
    for (int i = 0; i < p.n; i++) {
      double xp = Delta*(i+0.5) + p.box[0][0];
      for (int j = 0; j < ny; j++) {
        double yp = Delta*(j+0.5) + p.box[0][1];
        for(int k=0; k<nz; k++) {
          double zp=Delta*(k+0.5)+p.box[0][2];
          fprintf(p.fp, "%g %g %g", xp, yp, zp);
          int ind=0;
          for(scalar s in p.list)
            fprintf(p.fp, " %g", field[i][len*j*nz+len*k+ind++]);
          fputc('\n', p.fp);  //not double quotation""
        }
      }
    }
    fflush (p.fp);
  }
@if _MPI
  else // slave
    MPI_Reduce (field[0], NULL, len*p.n*ny*nz, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
@endif

  matrix_free (field);
}
