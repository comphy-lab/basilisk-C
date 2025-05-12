void sphere(char * fname,scalar s, coord center, double radius, int maxlevel){
  FILE *fpver = fopen (fname,"a"); 
  int nn = (1<<maxlevel);
  double ** field = matrix_new (nn, 2*nn, sizeof(double));
  double stp = M_PI/(double)nn;
  for (int i = 0; i < nn; i++){//theta
      double thetap = stp*i + stp/2.;
      for (int j = 0; j < 2*nn; j++) {//phi
    	  double phip = stp*j + stp/2.;
	  double xp = center.x + radius*sin(thetap)*cos(phip);
	  double yp = center.y + radius*sin(thetap)*sin(phip);
	  double zp = center.z + radius*cos(thetap);
    	  Point point = locate (xp, yp, zp);
    	  field[i][j] = point.level >= 0 ? interpolate(s,xp,yp,zp) : nodata;
    	}
    }
  if (pid() == 0){ // master
#if _MPI
    MPI_Reduce (MPI_IN_PLACE, field[0], 2*sq(nn), MPI_DOUBLE, MPI_MIN, 0,
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
    MPI_Reduce (field[0], NULL, 2*nn*nn, MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
#endif
  matrix_free (field);
}

