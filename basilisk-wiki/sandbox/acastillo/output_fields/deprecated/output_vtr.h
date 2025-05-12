/*
This function writes one XML VTK file per PID process of type rectilinear grid
(*.vtr) which can be read using Paraview. File stores scalar
and vector fields defined at the center points. Results are recorded on ASCII
format. Works in multigrid with a square domain.
*/
void output_vtr_ascii (scalar * list, vector * vlist, int n, FILE * fp, bool linear)
{
  fputs ("<?xml version=\"1.0\"?>\n"
  "<VTKFile type=\"RectilinearGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);
  fprintf (fp,"\t <RectilinearGrid WholeExtent=\"0 %d 0 %d 0 %d\">\n", N-1, N-1, 0);
  fprintf (fp,"\t\t <Piece Extent=\"0 %d 0 %d 0 %d\">\n", N-1, N-1, 0);
  fputs ("\t\t\t <PointData Scalars=\"scalars\">\n", fp);

  for (scalar s in list) {
    fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", s.name);
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
    fputs ("\t\t\t\t </DataArray>\n", fp);
  }
  for (vector v in vlist) {
    fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"%s\" format=\"ascii\">\n", v.x.name);
    double fn = n;
    double Delta = L0/fn;
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
  fputs ("\t\t\t </PointData>\n", fp);
  fputs ("\t\t\t <Coordinates>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"Float64\" Name=\"Array 0x47a5770\" format=\"ascii\">\n", fp);
  double fn = n;
  double Delta = L0/fn;
  for (int i = 0; i < n; i++) {
    double x = Delta*i + X0 + Delta/2.;
    fprintf (fp, "%g \n", x);
  }
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"Float64\" Name=\"Array 0x5228a30\" format=\"ascii\">\n", fp);
  for (int j = 0; j < n; j++) {
    double y = Delta*j + Y0 + Delta/2.;
    fprintf (fp, "%g \n", y);
  }
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"Float64\" Name=\"Array 0x52724b0\" format=\"ascii\">\n", fp);
  fputs ("0 \n", fp);
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t </Coordinates>\n", fp);
  fputs ("\t\t </Piece>\n", fp);
  fputs ("\t </RectilinearGrid>\n", fp);
  fputs ("</VTKFile>\n", fp);
  fflush (fp);
}

/*
This function writes one XML VTK file per PID process of type rectilinear grid
(*.vtr) which can be read using Paraview. File stores scalar
and vector fields defined at the center points. Results are recorded on RAW
binary format. Works in multigrid with a square domain.
*/
void output_vtr_bin (scalar * list, vector * vlist, int n, FILE * fp, bool linear)
{
  fputs ("<?xml version=\"1.0\"?>\n"
  "<VTKFile type=\"RectilinearGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);
  fprintf (fp,"\t <RectilinearGrid WholeExtent=\"0 %d 0 %d 0 %d\">\n", N-1, N-1, 0);
  fprintf (fp,"\t\t <Piece Extent=\"0 %d 0 %d 0 %d\">\n", N-1, N-1, 0);
  fputs ("\t\t\t <PointData Scalars=\"scalars\">\n", fp);

  int count = 0;
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
  fputs ("\t\t\t </PointData>\n", fp);
  fputs ("\t\t\t <Coordinates>\n", fp);
  fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" Name=\"Array 0x47a5770\" format=\"appended\" offset=\"%d\">\n",count);
  fputs ("\t\t\t\t </DataArray>\n", fp);
  count += (N+1)*8;
  fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" Name=\"Array 0x5228a30\" format=\"appended\" offset=\"%d\">\n",count);
  fputs ("\t\t\t\t </DataArray>\n", fp);
  count += (N+1)*8;
  fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" Name=\"Array 0x52724b0\" format=\"appended\" offset=\"%d\">\n",count);
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t </Coordinates>\n", fp);
  fputs ("\t\t </Piece>\n", fp);
  fputs ("\t </RectilinearGrid>\n", fp);
  fputs ("\t <AppendedData encoding=\"raw\">\n", fp);
  fputs ("_", fp);
  uint64_t block_len=N*N*8;
  double fn = n, z=0;
  double Delta = L0/fn;
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
  block_len=N*8;
  fwrite (&block_len, sizeof (uint64_t), 1, fp);
  for (int i = 0; i < n; i++) {
    double x = Delta*i + X0 + Delta/2.;
    fwrite (&x, sizeof (double), 1, fp);
  }
  fwrite (&block_len, sizeof (uint64_t), 1, fp);
  for (int j = 0; j < n; j++) {
    double y = Delta*j + Y0 + Delta/2.;
    fwrite (&y, sizeof (double), 1, fp);
  }
  block_len=1*8;
  fwrite (&block_len, sizeof (uint64_t), 1, fp);
  fwrite (&z, sizeof (double), 1, fp);
  fputs ("\t\n", fp);
  fputs ("\t </AppendedData>\n", fp);
  fputs ("</VTKFile>\n", fp);
  fflush (fp);
}
