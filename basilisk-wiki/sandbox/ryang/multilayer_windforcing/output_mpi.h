/*
This function is an extension of output_matrix() inside "output.h" which allows
to write 2D fields in a gnuplot-complatible format when running in MPI by
performing a MPI_Reduce.
*/

struct OutputMatrix_mpi {
  scalar f;
  FILE * fp;
  int n;
  bool linear;
};

trace
void output_matrix_mpi (struct OutputMatrix_mpi p)
{
  if (p.n == 0) p.n = N;
  if (!p.fp) p.fp = stdout;
  //boundary({p.f});
  if (p.linear) {
    scalar f = p.f;
    boundary ({f});
  }
  
  float fn = p.n, Delta = L0/fn;
  float ** field = matrix_new (p.n, p.n, sizeof(float));

  for (int i = 0; i < p.n; i++) {
    float xp = Delta*i + X0 + Delta/2.;
    for (int j = 0; j < p.n; j++) {
      float yp = Delta*j + Y0 + Delta/2.;
      if (p.linear) {
        field[i][j] = interpolate (p.f, xp, yp);
      }
      else {
        Point point = locate (xp, yp);
        field[i][j] = point.level >= 0 ? val(p.f) : nodata;
      }
    }
  }

  if (pid() == 0) { // master
@if _MPI
    MPI_Reduce (MPI_IN_PLACE, field[0], p.n*p.n, MPI_FLOAT, MPI_MIN, 0,MPI_COMM_WORLD);
@endif

    /*
    fwrite (&fn, sizeof(float), 1, p.fp);
    for (int j = 0; j < p.n; j++) {
      //float yp = Delta*j + Y0 + Delta/2.;
      fwrite (&yp, sizeof(float), 1, p.fp);
    }
    */

    for (int i = 0; i < p.n; i++){
      //float xp = Delta*i + X0 + Delta/2.;
      //fwrite (&xp, sizeof(float), 1, p.fp);
      for (int j = 0; j < p.n; j++) {
        fwrite (&field[i][j], sizeof(float), 1, p.fp);
      }
    }
    fflush (p.fp);
  }
@if _MPI
  else // slave
  MPI_Reduce (field[0], NULL, p.n*p.n, MPI_FLOAT, MPI_MIN, 0,MPI_COMM_WORLD);
@endif

  matrix_free (field);
}

/*
struct OutputMatrix_part {
  scalar f;
  FILE * fp;
  int n;
  bool linear;
  double box[2][2];
};

void output_matrix_part_mpi (struct OutputMatrix_part p)
{
  if (p.n == 0) p.n = N;
  if (!p.fp) p.fp = stdout;
  if (p.box[0][0] == 0. && p.box[0][1] == 0. && 
      p.box[1][0] == 0. && p.box[1][1] == 0.) {
    p.box[0][0] = X0;      p.box[0][1] = Y0;
    p.box[1][0] = X0 + L0; p.box[1][1] = Y0 + L0;
  }
  float fn = p.n, //Delta = L0/fn;
  Delta = (p.box[1][0] - p.box[0][0])/fn;
  float ny = (p.box[1][1] - p.box[0][1])/Delta;
  float nx = (p.box[1][0] - p.box[0][0])/Delta;
  
  float ** field = matrix_new ((int)nx,(int)ny, sizeof(float));

  for (int i = 0; i < nx; i++) {
    float xp = Delta*i + p.box[0][0] + Delta/2.;
    for (int j = 0; j < ny; j++) {
      float yp = Delta*j + p.box[0][1] + Delta/2.;
      if (p.linear) {
        field[i][j] = interpolate (p.f, xp, yp);
      }
      else {
        Point point = locate (xp, yp);
        field[i][j] = point.level >= 0 ? val(p.f) : nodata;
      }
    }
  }

  if (pid() == 0) { // master
@if _MPI
    MPI_Reduce (MPI_IN_PLACE, field[0], nx*ny, MPI_FLOAT, MPI_MIN, 0,MPI_COMM_WORLD);
@endif


    fwrite (&nx, sizeof(float), 1, p.fp);
	fwrite (&ny, sizeof(float), 1, p.fp);
    for (int j = 0; j < ny; j++) {
      float yp = Delta*j + p.box[0][1] + Delta/2.;
      fwrite (&yp, sizeof(float), 1, p.fp);
    }

    for (int i = 0; i < nx; i++){
      float xp = Delta*i + p.box[0][0] + Delta/2.;
      fwrite (&xp, sizeof(float), 1, p.fp);
      for (int j = 0; j < ny; j++) {
        fwrite (&field[i][j], sizeof(float), 1, p.fp);
      }
    }
    fflush (p.fp);
  }
@if _MPI
  else // slave
  MPI_Reduce (field[0], NULL, nx*ny, MPI_FLOAT, MPI_MIN, 0,MPI_COMM_WORLD);
@endif

  matrix_free (field);
}
*/