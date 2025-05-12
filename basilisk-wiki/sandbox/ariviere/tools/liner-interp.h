void lineX(char * fname,scalar s,double yp, double zp, int maxlevel){
  FILE *fpver =fopen (fname,"a"); 
  int nn = (1<<maxlevel);
  double ** field = matrix_new (1, nn, sizeof(double));
  double stp = L0/(double)nn;
  for (int i = 0; i < nn; i++){
      double xp = stp*i + X0 + stp/2.;
      Point point = locate (xp, yp,zp);
      field[0][i] = point.level >= 0 ? interpolate(s,xp,yp,zp) : nodata;
    }
  if (pid() == 0){ // master
#if _MPI
    MPI_Reduce (MPI_IN_PLACE, field[0], nn, MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
#endif
    for (int i = 0; i < nn; i++) {
			fprintf (fpver, "%g ", field[0][i]);
		}
		fputc ('\n', fpver);
    fflush (fpver);
  }
#if _MPI
    else // slave
    MPI_Reduce (field[0], NULL, nn, MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
#endif
  matrix_free (field);
}

void lineY(char * fname,scalar s,double xp, double zp, int maxlevel){
  FILE *fpver =fopen (fname,"a"); 
  int nn = (1<<maxlevel);
  double ** field = matrix_new (1, nn, sizeof(double));
  double stp = L0/(double)nn;
  for (int i = 0; i < nn; i++)
    {
      double yp = stp*i + Y0 + stp/2.;
      Point point = locate (xp, yp,zp);
      field[0][i] = point.level >= 0 ? interpolate(s,xp,yp,zp) : nodata;
    }
  if (pid() == 0){ // master
#if _MPI
    MPI_Reduce (MPI_IN_PLACE, field[0], nn, MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
#endif
    for (int i = 0; i < nn; i++) {
			fprintf (fpver, "%g ", field[0][i]);
    }
		fputc ('\n', fpver);
    fflush (fpver);
  }
#if _MPI
  else // slave
    MPI_Reduce (field[0], NULL, nn, MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
#endif
  matrix_free (field);
}

void lineZ(char * fname,scalar s,double xp, double yp, int maxlevel){
  FILE *fpver =fopen (fname,"a"); 
  int nn = (1<<maxlevel);
  double ** field = matrix_new (1, nn, sizeof(double));
  double stp = L0/(double)nn;
  for (int i = 0; i < nn; i++)
    {
      double zp = stp*i + Z0 + stp/2.;
      Point point = locate (xp, yp,zp);
      field[0][i] = point.level >= 0 ? interpolate(s,xp,yp,zp) : nodata;
    }
  if (pid() == 0){ // master
#if _MPI
    MPI_Reduce (MPI_IN_PLACE, field[0], nn, MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
#endif
    for (int i = 0; i < nn; i++) {
			fprintf (fpver, "%g ", field[0][i]);
    }
    fputc ('\n', fpver);
    fflush (fpver);
  }
#if _MPI
  else // slave
    MPI_Reduce (field[0], NULL, nn, MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
#endif
  matrix_free (field);
}
