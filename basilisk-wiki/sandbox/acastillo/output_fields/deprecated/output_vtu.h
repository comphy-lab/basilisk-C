/*
This function allows to write one XML VTK file per PID process of type
unstructured grid (*.vtu) which can be read using Paraview. File stores scalar
and vector fields defined at the center points. Results are recorded on ASCII
format. Works in multigrid with a square domain. 
*/
void output_vtu_ascii (scalar * list, vector * vlist, int n, FILE * fp, bool linear)
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

  fputs ("<?xml version=\"1.0\"?>\n"
  "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);
  fputs ("\t <UnstructuredGrid>\n", fp);
  fprintf (fp,"\t\t <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", (N+1)*(N+1), N*N);
  fputs ("\t\t\t <CellData Scalars=\"scalars\">\n", fp);

  for (scalar s in list) {
    fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", s.name);
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
    fputs ("\t\t\t\t </DataArray>\n", fp);
  }
  for (vector v in vlist) {
    fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"%s\" format=\"ascii\">\n", v.x.name);
    for (int i = 0; i < n; i++) {
      double x = Delta*i + X0 + Delta/2.;
      for (int j = 0; j < n; j++) {
        double y = Delta*j + Y0 + Delta/2., vx, vy;
        if (linear){
          vx = interpolate (v.x, x, y);
          vy = interpolate (v.y, x, y);
        }
        else {
          Point point = locate (x, y);
          vx = point.level >= 0 ? val(v.x) : nodata;
          vy = point.level >= 0 ? val(v.y) : nodata;
        }
        fprintf (fp, "%g %g 0.\n", vx, vy);
      }
    }
    fputs ("\t\t\t\t </DataArray>\n", fp);
  }
  fputs ("\t\t\t </CellData>\n", fp);
  fputs ("\t\t\t <Points>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n", fp);
  for (int i = 0; i < n+1; i++) {
    double x = Delta*i + X0;
    for (int j = 0; j < n+1; j++) {
      double y = Delta*j + Y0;
      fprintf (fp, "%g %g 0\n", x, y);
    }
  }
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t </Points>\n", fp);
  fputs ("\t\t\t <Cells>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n", fp);
  for (int i = 0; i < n; i++) {
    double x = Delta*i + X0 + Delta/2.;
    for (int j = 0; j < n; j++) {
      double y = Delta*j + Y0 + Delta/2.;
      fprintf (fp, "%g %g %g %g\n", control[i][j], control[i+1][j],control[i][j+1],control[i+1][j+1]);
    }
  }
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n", fp);
  count=0;
  for (int i = 0; i < n; i++) {
    double x = Delta*i + X0 + Delta/2.;
    for (int j = 0; j < n; j++) {
      double y = Delta*j + Y0 + Delta/2.;
      count+=4;
      fprintf (fp, "%d \n", count);
    }
  }
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n", fp);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      fputs ("8 \n", fp);
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t </Cells>\n", fp);
  fputs ("\t\t </Piece>\n", fp);
  fputs ("\t </UnstructuredGrid>\n", fp);
  fputs ("</VTKFile>\n", fp);
  fflush (fp);
}

/*
This function allows to write one XML VTK file per PID process of type
unstructured grid (*.vtu) which can be read using Paraview. File stores scalar
and vector fields defined at the center points. Results are recorded on binary
RAW format.
*/
void output_vtu_bin (scalar * list, vector * vlist, int n, FILE * fp, bool linear)
{
  int count = 0;
  double fn = n;
  double Delta = L0/fn;
  int control[n+1][n+1];
  for (int i = 0; i < n+1; i++) {
    for (int j = 0; j < n+1; j++) {
      control[i][j] = count ;
      count+=1;
    }
  }

  fputs ("<?xml version=\"1.0\"?>\n"
  "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);
  fputs ("\t <UnstructuredGrid>\n", fp);
  fprintf (fp,"\t\t <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", (N+1)*(N+1), N*N);
  fputs ("\t\t\t <CellData Scalars=\"scalars\">\n", fp);

  count = 0;
  for (scalar s in list) {
    fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" Name=\"%s\" format=\"appended\" offset=\"%d\">\n", s.name,count);
    count += ((N*N)+1)*8;
    fputs ("\t\t\t\t </DataArray>\n", fp);
  }
  for (vector v in vlist) {
    fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"3\"  format=\"appended\" offset=\"%d\">\n", v.x.name,count);
    count += ((N*N*3)+1)*8;
    fputs ("\t\t\t\t </DataArray>\n", fp);
  }
  fputs ("\t\t\t </CellData>\n", fp);
  fputs ("\t\t\t <Points>\n", fp);
  fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\"  format=\"appended\" offset=\"%d\">\n",count);
  count += (((N+1)*(N+1)*3)+1)*8;
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t </Points>\n", fp);
  fputs ("\t\t\t <Cells>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n", fp);
  for (int i = 0; i < n; i++) {
    double x = Delta*i + X0 + Delta/2.;
    for (int j = 0; j < n; j++) {
      double y = Delta*j + Y0 + Delta/2.;
      fprintf (fp, "%d %d %d %d\n", control[i][j], control[i+1][j],control[i][j+1],control[i+1][j+1]);
    }
  }
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n", fp);
  count=0;
  for (int i = 0; i < n; i++) {
    double x = Delta*i + X0 + Delta/2.;
    for (int j = 0; j < n; j++) {
      double y = Delta*j + Y0 + Delta/2.;
      count+=4;
      fprintf (fp, "%d \n", count);
    }
  }
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n", fp);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      fputs ("8 \n", fp);
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t </Cells>\n", fp);
  fputs ("\t\t </Piece>\n", fp);
  fputs ("\t </UnstructuredGrid>\n", fp);
  fputs ("\t <AppendedData encoding=\"raw\">\n", fp);
  fputs ("_", fp);
  uint64_t block_len=N*N*8;
  double z=0;
  for (scalar s in list) {
    fwrite (&block_len, sizeof (uint64_t), 1, fp);
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
        fwrite (&v, sizeof (double), 1, fp);
      }
    }
  }
  block_len=N*N*8*3;
  for (vector v in vlist) {
    fwrite (&block_len, sizeof (uint64_t), 1, fp);
    double vz=0;
    for (int i = 0; i < n; i++) {
      double x = Delta*i + X0 + Delta/2.;
      for (int j = 0; j < n; j++) {
        double y = Delta*j + Y0 + Delta/2., vx, vy;
        if (linear){
          vx = interpolate (v.x, x, y);
          vy = interpolate (v.y, x, y);
        }
        else {
          Point point = locate (x, y);
          vx = point.level >= 0 ? val(v.x) : nodata;
          vy = point.level >= 0 ? val(v.y) : nodata;
        }
        fwrite (&vx, sizeof (double), 1, fp);
        fwrite (&vy, sizeof (double), 1, fp);
        fwrite (&vz, sizeof (double), 1, fp);
      }
    }
  }
  block_len=(N+1)*(N+1)*8*3;
  fwrite (&block_len, sizeof (uint64_t), 1, fp);
  for (int i = 0; i < n+1; i++) {
    double x = Delta*i + X0;
    for (int j = 0; j < n+1; j++) {
      double y = Delta*j + Y0;
      fwrite (&x, sizeof (double), 1, fp);
      fwrite (&y, sizeof (double), 1, fp);
      fwrite (&z, sizeof (double), 1, fp);
    }
  }
  fputs ("\t\n", fp);
  fputs ("\t </AppendedData>\n", fp);
  fputs ("</VTKFile>\n", fp);
  fflush (fp);
}
