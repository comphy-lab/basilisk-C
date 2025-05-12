struct OutputMatrixNormal {
  scalar f;
  FILE * fp;
  int n;
  double L;
  coord Pos;
  coord nvec;
  coord bvec;
};



trace
void output_matrix_normal (struct OutputMatrixNormal p)
{
  if (p.n == 0) p.n = 32;
  if (!p.fp) p.fp = stdout;
  float fn = p.n;
  coord Pos = p.Pos;
  coord nvec = p.nvec;
  coord bvec = p.bvec;

  float Delta = (float) p.L/fn;
  float ** field = matrix_new (p.n, p.n, sizeof(float));

  for (int i = 0; i < p.n; i++) {
    float bp = (float) (Delta*i - p.L/2);
    for (int j = 0; j < p.n; j++) {
      float np = (float) (Delta*j - p.L/2);

      float xp = (np * nvec.x) + (bp * bvec.x) + Pos.x;
      float yp = (np * nvec.y) + (bp * bvec.y) + Pos.y;
      float zp = (np * nvec.z) + (bp * bvec.z) + Pos.z;

      field[i][j] = interpolate (p.f, xp, yp, zp);
    }
  }

  if (pid() == 0) { // master
@if _MPI
    MPI_Reduce (MPI_IN_PLACE, field[0], p.n*p.n, MPI_FLOAT, MPI_MIN, 0,MPI_COMM_WORLD);
@endif

    fwrite (&fn, sizeof(float), 1, p.fp);
    for (int j = 0; j < p.n; j++) {
      float np = (float) (Delta*j - p.L/2);
      fwrite (&np, sizeof(float), 1, p.fp);
    }

    for (int i = 0; i < p.n; i++){
      float bp = (float) (Delta*i - p.L/2);
      fwrite (&bp, sizeof(float), 1, p.fp);
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
















struct OutputVtsNormal {
  vector * vlist;
  char * name;
  int n;
  double L;
  coord Pos;
  coord tvec;
  coord nvec;
  coord bvec;
};

trace
void output_vts_normal (struct OutputVtsNormal p)
{
  FILE * fp = fopen(p.name, "w");

  if (p.n == 0) p.n = 32;
  double fn = p.n;
  coord Pos = p.Pos;
  coord nvec = p.nvec;
  coord bvec = p.bvec;

  double Delta = (double) p.L/fn;
  double ** field_x = matrix_new (p.n, p.n, sizeof(double));
  double ** field_y = matrix_new (p.n, p.n, sizeof(double));
  double ** field_z = matrix_new (p.n, p.n, sizeof(double));

  if (pid() == 0) { // master

    fputs ("<?xml version=\"1.0\"?>\n"
    "<VTKFile type=\"StructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);
    fprintf (fp,"\t <StructuredGrid WholeExtent=\"0 %d 0 %d 0 %d\">\n", p.n-1, p.n-1, 0);
    fprintf (fp,"\t\t <Piece Extent=\"0 %d 0 %d 0 %d\">\n", p.n-1, p.n-1, 0);
    fputs ("\t\t\t <PointData Scalars=\"scalars\">\n", fp);

    int count = 0;
    for (vector v in p.vlist) {
      fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"3\"  format=\"appended\" offset=\"%d\">\n", v.x.name, count);
      count += ((p.n*p.n*3)*sizeof (double)) + sizeof(uint64_t);
      fputs ("\t\t\t\t </DataArray>\n", fp);
    }
    fputs ("\t\t\t </PointData>\n", fp);
    fputs ("\t\t\t <Points>\n", fp);
    fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\"  format=\"appended\" offset=\"%d\">\n", count);
    fputs ("\t\t\t\t </DataArray>\n", fp);
    fputs ("\t\t\t </Points>\n", fp);
    fputs ("\t\t </Piece>\n", fp);
    fputs ("\t </StructuredGrid>\n", fp);
    fputs ("\t <AppendedData encoding=\"raw\">\n", fp);
    fputs ("_", fp);
  }

  uint64_t block_len=p.n*p.n*sizeof (double)*3;
  for (vector v in p.vlist){
    for (int i = 0; i < p.n; i++) {
      double bp = (double) (Delta*i - p.L/2);
      for (int j = 0; j < p.n; j++) {
        double np = (double) (Delta*j - p.L/2);
        double xp = (np * nvec.x) + (bp * bvec.x) + Pos.x;
        double yp = (np * nvec.y) + (bp * bvec.y) + Pos.y;
        double zp = (np * nvec.z) + (bp * bvec.z) + Pos.z;
        field_x[i][j] = interpolate (v.x, xp, yp, zp);
        field_y[i][j] = interpolate (v.y, xp, yp, zp);
        field_z[i][j] = interpolate (v.z, xp, yp, zp);
      }
    }

@if _MPI
    if (pid() == 0) { // master
      MPI_Reduce (MPI_IN_PLACE, field_x[0], p.n*p.n, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
      MPI_Reduce (MPI_IN_PLACE, field_y[0], p.n*p.n, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
      MPI_Reduce (MPI_IN_PLACE, field_z[0], p.n*p.n, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    }
    else{ // slave
      MPI_Reduce (field_x[0], NULL, p.n*p.n, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
      MPI_Reduce (field_y[0], NULL, p.n*p.n, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
      MPI_Reduce (field_z[0], NULL, p.n*p.n, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    }
@endif

    if (pid() == 0) { // master
      fwrite (&block_len, sizeof (uint64_t), 1, fp);
      for (int i = 0; i < p.n; i++) {
        for (int j = 0; j < p.n; j++) {
          fwrite (&field_x[i][j], sizeof (double), 1, fp);
          fwrite (&field_y[i][j], sizeof (double), 1, fp);
          fwrite (&field_z[i][j], sizeof (double), 1, fp);
      	}
      }
    }
  }


  if (pid() == 0) { // master
    fwrite (&block_len, sizeof (uint64_t), 1, fp);
    for (int i = 0; i < p.n; i++) {
      double bp = (double) (Delta*i - p.L/2);
      for (int j = 0; j < p.n; j++) {
        double np = (double) (Delta*j - p.L/2);

        double xp = (np * nvec.x) + (bp * bvec.x) + Pos.x;
        double yp = (np * nvec.y) + (bp * bvec.y) + Pos.y;
        double zp = (np * nvec.z) + (bp * bvec.z) + Pos.z;

        fwrite (&xp, sizeof (double), 1, fp);
        fwrite (&yp, sizeof (double), 1, fp);
        fwrite (&zp, sizeof (double), 1, fp);
      }
    }
    fputs ("\t\n", fp);
    fputs ("\t </AppendedData>\n", fp);
    fputs ("</VTKFile>\n", fp);
    fflush (fp);
  }

  fclose(fp);

  matrix_free (field_x);
  matrix_free (field_y);
  matrix_free (field_z);
}



trace
void output_vts_normal2 (struct OutputVtsNormal p)
{
  FILE * fp = fopen(p.name, "w");

  if (p.n == 0) p.n = 32;
  double fn = p.n;
  coord Pos = p.Pos, tvec = p.tvec, nvec = p.nvec, bvec = p.bvec;

  double Delta = (double) p.L/fn;
  double ** field_x = matrix_new (p.n, p.n, sizeof(double));
  double ** field_y = matrix_new (p.n, p.n, sizeof(double));
  double ** field_z = matrix_new (p.n, p.n, sizeof(double));

  if (pid() == 0) { // master

    fputs ("<?xml version=\"1.0\"?>\n"
    "<VTKFile type=\"StructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);
    fprintf (fp,"\t <StructuredGrid WholeExtent=\"0 %d 0 %d 0 %d\">\n", p.n-1, p.n-1, 0);
    fprintf (fp,"\t\t <Piece Extent=\"0 %d 0 %d 0 %d\">\n", p.n-1, p.n-1, 0);
    fputs ("\t\t\t <PointData Scalars=\"scalars\">\n", fp);

    int count = 0;
    for (vector v in p.vlist) {
      fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"3\"  format=\"appended\" offset=\"%d\">\n", v.x.name, count);
      count += ((p.n*p.n*3)*sizeof (double)) + sizeof(uint64_t);
      fputs ("\t\t\t\t </DataArray>\n", fp);
    }
    fputs ("\t\t\t </PointData>\n", fp);
    fputs ("\t\t\t <Points>\n", fp);
    fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\"  format=\"appended\" offset=\"%d\">\n", count);
    fputs ("\t\t\t\t </DataArray>\n", fp);
    fputs ("\t\t\t </Points>\n", fp);
    fputs ("\t\t </Piece>\n", fp);
    fputs ("\t </StructuredGrid>\n", fp);
    fputs ("\t <AppendedData encoding=\"raw\">\n", fp);
    fputs ("_", fp);
  }

  uint64_t block_len=p.n*p.n*sizeof (double)*3;
  for (vector v in p.vlist){
    for (int i = 0; i < p.n; i++) {
      double bp = (double) (Delta*i - p.L/2);
      for (int j = 0; j < p.n; j++) {
        double np = (double) (Delta*j - p.L/2);
        double xp = (np * nvec.x) + (bp * bvec.x) + Pos.x;
        double yp = (np * nvec.y) + (bp * bvec.y) + Pos.y;
        double zp = (np * nvec.z) + (bp * bvec.z) + Pos.z;
        field_x[i][j] = interpolate (v.x, xp, yp, zp);
        field_y[i][j] = interpolate (v.y, xp, yp, zp);
        field_z[i][j] = interpolate (v.z, xp, yp, zp);
      }
    }

@if _MPI
    if (pid() == 0) { // master
      MPI_Reduce (MPI_IN_PLACE, field_x[0], p.n*p.n, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
      MPI_Reduce (MPI_IN_PLACE, field_y[0], p.n*p.n, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
      MPI_Reduce (MPI_IN_PLACE, field_z[0], p.n*p.n, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    }
    else{ // slave
      MPI_Reduce (field_x[0], NULL, p.n*p.n, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
      MPI_Reduce (field_y[0], NULL, p.n*p.n, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
      MPI_Reduce (field_z[0], NULL, p.n*p.n, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    }
@endif

    if (pid() == 0) { // master
      fwrite (&block_len, sizeof (uint64_t), 1, fp);
      for (int i = 0; i < p.n; i++) {
        for (int j = 0; j < p.n; j++) {

          double field_t = field_x[i][j]*tvec.x + field_y[i][j]*tvec.y + field_z[i][j]*tvec.z;
          double field_b = field_x[i][j]*bvec.x + field_y[i][j]*bvec.y + field_z[i][j]*bvec.z;
          double field_n = field_x[i][j]*nvec.x + field_y[i][j]*nvec.y + field_z[i][j]*nvec.z;

          fwrite (&field_b, sizeof (double), 1, fp);
          fwrite (&field_n, sizeof (double), 1, fp);
          fwrite (&field_t, sizeof (double), 1, fp);
      	}
      }
    }
  }


  if (pid() == 0) { // master
    double tp = 0.;
    fwrite (&block_len, sizeof (uint64_t), 1, fp);
    for (int i = 0; i < p.n; i++) {
      double bp = (double) (Delta*i - p.L/2);
      for (int j = 0; j < p.n; j++) {
        double np = (double) (Delta*j - p.L/2);
        fwrite (&bp, sizeof (double), 1, fp);
        fwrite (&np, sizeof (double), 1, fp);
        fwrite (&tp, sizeof (double), 1, fp);
      }
    }
    fputs ("\t\n", fp);
    fputs ("\t </AppendedData>\n", fp);
    fputs ("</VTKFile>\n", fp);
    fflush (fp);
  }

  fclose(fp);

  matrix_free (field_x);
  matrix_free (field_y);
  matrix_free (field_z);
}
