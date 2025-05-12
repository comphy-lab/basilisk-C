/**

# Obtain a slice of 3D data

Three functions to obtain a square XZ,YZ,ZY slice, respectively, of scalar *s*. The data is outputted in a file with *fname* in a regular grid fashion  (maxlevel$\ ^2$).  

It should work for all grids and parallization styles. 

*/
void sliceXZ(char * fname,scalar s,double yp, int maxlevel){
  FILE *fpver =fopen (fname,"w"); 
  int nn = (1<<maxlevel);
  double ** field = matrix_new (nn, nn, sizeof(double));
  double stp = L0/(double)nn;
  for (int i = 0; i < nn; i++)
    {
      double xp = stp*i + X0 + stp/2.;
      for (int j = 0; j < nn; j++) 
	{
	  double zp = stp*j + Z0 + stp/2.;
	  Point point = locate (xp, yp,zp);
	  field[i][j] = point.level >= 0 ? s[] : nodata;
	}
    }
  if (pid() == 0){ // master
#if _MPI
    MPI_Reduce (MPI_IN_PLACE, field[0], sq(nn), MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
#endif
    for (int i = 0; i < nn; i++) {
      for (int j = 0; j < nn; j++) {
	fprintf (fpver, "%g\t", field[i][j]);
      }
      fputc ('\n', fpver);
    }
    fflush (fpver);
  }
#if _MPI
  else // slave
    MPI_Reduce (field[0], NULL, nn*nn, MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
#endif
  matrix_free (field);
}

void sliceYZ(char * fname,scalar s,double xp, int maxlevel){
  FILE *fpver =fopen (fname,"w"); 
  int nn = (1<<maxlevel);
  double ** field = matrix_new (nn, nn, sizeof(double));
  double stp = L0/(double)nn;
  for (int i = 0; i < nn; i++)
    {
      double yp = stp*i + Y0 + stp/2.;
      for (int j = 0; j < nn; j++) 
	{
	  double zp = stp*j + Z0 + stp/2.;
	  Point point = locate (xp, yp,zp);
	  field[i][j] = point.level >= 0 ? s[] : nodata;
	}
    }
  if (pid() == 0){ // master
#if _MPI
    MPI_Reduce (MPI_IN_PLACE, field[0], sq(nn), MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
#endif
    for (int i = 0; i < nn; i++) {
      for (int j = 0; j < nn; j++) {
	fprintf (fpver, "%g\t", field[i][j]);
      }
      fputc ('\n', fpver);
    }
    fflush (fpver);
  }
#if _MPI
  else // slave
    MPI_Reduce (field[0], NULL, nn*nn, MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
#endif
  matrix_free (field);
}

void sliceXY(char * fname,scalar s,double zp, int maxlevel){
  FILE *fpver =fopen (fname,"w"); 
  int nn = (1<<maxlevel);
  double ** field = matrix_new (nn, nn, sizeof(double));
  double stp = L0/(double)nn;
  for (int i = 0; i < nn; i++)
    {
      double xp = stp*i + X0 + stp/2.;
      for (int j = 0; j < nn; j++) 
	{
	  double yp = stp*j + Y0 + stp/2.;
	  Point point = locate (xp, yp,zp);
	  field[i][j] = point.level >= 0 ? s[] : nodata;
	}
    }
  if (pid() == 0){ // master
#if _MPI
    MPI_Reduce (MPI_IN_PLACE, field[0], sq(nn), MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
#endif
    for (int i = 0; i < nn; i++) {
      for (int j = 0; j < nn; j++) {
	fprintf (fpver, "%g\t", field[i][j]);
      }
      fputc ('\n', fpver);
    }
    fflush (fpver);
  }
#if _MPI
  else // slave
    MPI_Reduce (field[0], NULL, nn*nn, MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
#endif
  matrix_free (field);
}
