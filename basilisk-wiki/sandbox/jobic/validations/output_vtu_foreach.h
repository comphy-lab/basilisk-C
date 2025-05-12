#include <stdint.h>
/*
This function writes one XML VTK file per PID process of type unstructured grid
(*.vtu) which can be read using Paraview. File stores scalar and vector fields
defined at the center points. Results are recorded on binary format. If one writes
one *.vtu file per PID process this function may be combined with
output_pvtu() above to read in parallel. Tested in (quad- and oct-)trees
using MPI. Also works with solids (when not using MPI).
*/

void output_vtu_pid (scalar * list, vector * vlist, char * subname)
{
  char name[80];
  sprintf(name, "%s.vtu", subname);
  FILE * fp = fopen(name, "w");

#if defined(_OPENMP)
  int num_omp = omp_get_max_threads();
  omp_set_num_threads(1);
#endif

  // This is to match the number of points and vertices when using periodic conditions
  scalar per_mask[];
  foreach(serial, noauto){
    per_mask[] = 1.;
    if (((Period.x) & (x + Delta > X0 + L0)))
      per_mask[] = 0.;
#if dimension > 1
    if (((Period.y) & (y + Delta > Y0 + L0)))
      per_mask[] = 0.;
#endif
#if dimension > 2
    if (((Period.z) & (z + Delta > Z0 + L0)))
      per_mask[] = 0.;
#endif
  }

  vertex scalar marker[];
  uint64_t no_points = 0, no_cells=0 ;
  foreach_vertex(serial, noauto){
    marker[] = _k;
    if is_vertex(cell)
      no_points += 1;
  }
  foreach(serial, noauto)
    if (per_mask[])
      no_cells += 1;


  fputs ("<?xml version=\"1.0\"?>\n"
  "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);
  fputs ("\t <UnstructuredGrid>\n", fp);
  fprintf (fp,"\t\t <Piece NumberOfPoints=\"%ld\" NumberOfCells=\"%ld\">\n", no_points, no_cells);
  fputs ("\t\t\t <CellData Scalars=\"scalars\">\n", fp);

  uint64_t count = 0;
  for (scalar s in list) {
    fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" Name=\"%s\" format=\"appended\" offset=\"%ld\"/>\n", s.name, count);
    count += (no_cells * sizeof (double)) +  sizeof (uint64_t);
  }
  for (vector v in vlist) {
    fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"3\"  format=\"appended\" offset=\"%ld\"/>\n", v.x.name, count);
    count += (3 * no_cells * sizeof (double)) +  sizeof (uint64_t);
  }
  fputs ("\t\t\t </CellData>\n", fp);
  fputs ("\t\t\t <Points>\n", fp);
  fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\"  format=\"appended\" offset=\"%ld\"/>\n",count);
  fputs ("\t\t\t </Points>\n", fp);
  fputs ("\t\t\t <Cells>\n", fp);
  count += (3 * no_points * sizeof (double)) +  sizeof (uint64_t);

  #if dimension == 2
    //typedef unsigned char uint8_t type=9, noffset=4;
    uint8_t type=9, noffset=4;
  #elif dimension == 3
    uint8_t type=12, noffset=8;
  #endif
  uint64_t offset, connectivity[noffset];

  fprintf (fp,"\t\t\t\t <DataArray type=\"UInt64\" Name=\"offsets\" format=\"appended\" offset=\"%ld\"/>\n", count);
  count +=  (no_cells * sizeof (uint64_t)) + sizeof (uint64_t);

  fprintf (fp,"\t\t\t\t <DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\"%ld\"/>\n", count);
  count +=  (no_cells * sizeof (uint8_t)) + sizeof (uint64_t);

  fprintf (fp,"\t\t\t\t <DataArray type=\"UInt64\" Name=\"connectivity\" format=\"appended\" offset=\"%ld\"/>\n", count);
  count +=  (no_cells * noffset * sizeof (uint64_t)) + sizeof (uint64_t);

  fputs ("\t\t\t </Cells>\n", fp);
  fputs ("\t\t </Piece>\n", fp);
  fputs ("\t </UnstructuredGrid>\n", fp);
  fputs ("\t <AppendedData encoding=\"raw\">\n", fp);
  fputs ("_", fp);

  uint64_t block_len=no_cells*sizeof (double);
#if dimension == 2
  double vz=0;
#endif
  for (scalar s in list) {
    fwrite (&block_len, sizeof (uint64_t), 1, fp);
    foreach(serial, noauto)
      if (per_mask[])
        fwrite (&val(s), sizeof (double), 1, fp);
  }
  block_len=no_cells*3*sizeof (double);
  for (vector v in vlist) {
    fwrite (&block_len, sizeof (uint64_t), 1, fp);
    foreach(serial, noauto){
      if (per_mask[]){
        fwrite (&val(v.x), sizeof (double), 1, fp);
        fwrite (&val(v.y), sizeof (double), 1, fp);
#if dimension == 2
        fwrite (&vz, sizeof (double), 1, fp);
#elif dimension == 3
        fwrite (&val(v.z), sizeof (double), 1, fp);
#endif
      }
    }
  }
  block_len=no_points*3*sizeof (double);
  fwrite (&block_len, sizeof (uint64_t), 1, fp);
  foreach_vertex(serial, noauto){
    fwrite (&x, sizeof (double), 1, fp);
    fwrite (&y, sizeof (double), 1, fp);
    fwrite (&z, sizeof (double), 1, fp);
  }

  block_len=no_cells*1*sizeof(uint64_t);
  fwrite (&block_len, sizeof (uint64_t), 1, fp);
  for (int i = 0; i < no_cells; i++){
    offset = (i+1)*noffset;
    fwrite (&offset, sizeof (int64_t), 1, fp);
  }

  block_len=no_cells*1*sizeof(uint8_t);
  fwrite (&block_len, sizeof (uint64_t), 1, fp);
  for (int i = 0; i < no_cells; i++)
    fwrite (&type, sizeof (int8_t), 1, fp);

  block_len=no_cells*noffset*sizeof(uint64_t);
  fwrite (&block_len, sizeof (uint64_t), 1, fp);
  foreach(serial, noauto){
    if (per_mask[]){
      connectivity[0] = (uint64_t)marker[];
      connectivity[1] = (uint64_t)marker[1];
      connectivity[2] = (uint64_t)marker[1,1];
      connectivity[3] = (uint64_t)marker[0,1];
  #if dimension == 3
      connectivity[4] = (uint64_t)marker[0,0,1];
      connectivity[5] = (uint64_t)marker[1,0,1];
      connectivity[6] = (uint64_t)marker[1,1,1];
      connectivity[7] = (uint64_t)marker[0,1,1];
  #endif
      fwrite (&connectivity, sizeof (uint64_t), noffset, fp);
    }
  }
  fputs ("\t\n", fp);
  fputs ("\t </AppendedData>\n", fp);
  fputs ("</VTKFile>\n", fp);
  fflush (fp);
#if defined(_OPENMP)
  omp_set_num_threads(num_omp);
#endif
  fclose (fp);
}


/*
This function writes one XML file which allows to read the *.vtu files generated
by output_vtu() when used in MPI. Tested in (quad- and oct-)trees
using MPI.
*/

@if _MPI
void output_pvtu (scalar * list, vector * vlist, char * subname)
{
  char name[80];
  FILE * fp;
  if (pid() == 0) {
    sprintf(name, "%s.pvtu", subname);
    fp = fopen(name, "w");

    fputs ("<?xml version=\"1.0\"?>\n"
    "<VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);
    fputs ("\t <PUnstructuredGrid GhostLevel=\"0\">\n", fp);
    fputs ("\t\t\t <PCellData Scalars=\"scalars\">\n", fp);
    for (scalar s in list) {
      fprintf (fp,"\t\t\t\t <PDataArray type=\"Float64\" Name=\"%s\" format=\"appended\">\n", s.name);
      fputs ("\t\t\t\t </PDataArray>\n", fp);
    }
    for (vector v in vlist) {
      fprintf (fp,"\t\t\t\t <PDataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"%s\" format=\"appended\">\n", v.x.name);
      fputs ("\t\t\t\t </PDataArray>\n", fp);
    }
    fputs ("\t\t\t </PCellData>\n", fp);
    fputs ("\t\t\t <PPoints>\n", fp);
    fputs ("\t\t\t\t <PDataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"appended\">\n", fp);
    fputs ("\t\t\t\t </PDataArray>\n", fp);
    fputs ("\t\t\t </PPoints>\n", fp);

    for (int i = 0; i < npe(); i++)
      fprintf (fp, "<Piece Source=\"%s_n%3.3d.vtu\"/> \n", subname, i);

    fputs ("\t </PUnstructuredGrid>\n", fp);
    fputs ("</VTKFile>\n", fp);
    fflush (fp);
    fclose (fp);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  sprintf(name, "%s_n%3.3d", subname, pid());
  output_vtu_pid (list, vlist, name);
}
@endif

// ********************************************************************/
@if _MPI
void output_pvd (scalar * list, vector * vlist, char * subname, double time)
{
  char name[80];
  //char *filename = "collectionParaview.pvd";
  FILE * fp;
  if (pid() == 0) {
    sprintf(name, "%s.pvd", subname);
    //if(i==0.){
      fp = fopen(name, "w");

      fputs ("<?xml version=\"1.0\"?>\n"
    "<VTKFile type='Collection' version='0.1' byte_order='LittleEndian' compressor='vtkZLibDataCompressor'>\n", fp);
      fputs ("\t<Collection>\n", fp);
    //}
    
    for (int i = 0; i < npe(); i++)
      fprintf (fp, "\t\t<DataSet timestep='%g' part='%i' file='%s_n%i.vtu'/>\n", time, i, subname, i);

    fputs ("\t</Collection>\n", fp);
    fputs ("</VTKFile>\n", fp);
    fflush (fp);
    fclose (fp);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  sprintf(name, "%s_n%i", subname, pid());
  output_vtu_pid (list, vlist, name);
}
@endif
// ********************************************************************/
trace
void output_vtu (scalar * list, vector * vlist, char * subname)
{
  @if _MPI
  output_pvtu (list, vlist, subname);
  @else
  output_vtu_pid (list, vlist, subname);
  @endif
}
// ********************************************************************/
trace
void output_vtu2 (scalar * list, vector * vlist, char * subname, double time)
{
  @if _MPI
  output_pvd (list, vlist, subname, time);
  output_pvtu (list, vlist, subname);
  @else
  output_vtu_pid (list, vlist, subname);
  @endif
}
// ********************************************************************/

/**
This function writes a slice defined by $n={n_x, n_y, n_z}$ and $\alpha$
using a similar approach to `view.h` into one XML VTK file per PID process of
type unstructured grid (*.vtu) in binary format which can be read using Paraview.
Only works for x, y, and/or z planes. Here, _alpha is assumed to intersect a cell
face and must be a multiple of L0/1<<MINLEVEL. This also means that scalar and
vector fields are written at the corresponding face values.
*/
#if dimension > 2
void output_vtu_plane_pid (scalar * list, vector * vlist, char * subname, coord n, double _alpha)
{
  char name[80];
  sprintf(name, "%s.vtu", subname);
  FILE * fp = fopen(name, "w");

#if defined(_OPENMP)
  int num_omp = omp_get_max_threads();
  omp_set_num_threads(1);
#endif

  scalar per_mask[];
  foreach(serial, noauto){
    per_mask[] = 0.;
    double alpha = (_alpha - n.x*x - n.y*y - n.z*z)/Delta;
    if (fabs(alpha) > 0.87)
      continue;

    if (alpha > 0.)
      per_mask[] = 1.;

    if ((((Period.x) & (x + Delta > X0 + L0)) || ((Period.y) & (y + Delta > Y0 + L0))) || ((Period.z) & (z + Delta > Z0 + L0)))
      per_mask[] = 0.;
  }

  vertex scalar marker[];
  uint64_t no_points = 0, no_cells=0 ;
  foreach_vertex(serial, noauto){
    marker[] = 0.;
    double alpha = (_alpha - n.x*x - n.y*y - n.z*z)/Delta;

    if (fabs(alpha) > 0.87)
      continue;

    marker[] = no_points;
    no_points += 1;
  }

  foreach(serial, noauto)
    if (per_mask[])
      no_cells += 1;


  fputs ("<?xml version=\"1.0\"?>\n"
  "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);
  fputs ("\t <UnstructuredGrid>\n", fp);
  fprintf (fp,"\t\t <Piece NumberOfPoints=\"%ld\" NumberOfCells=\"%ld\">\n", no_points, no_cells);
  fputs ("\t\t\t <CellData Scalars=\"scalars\">\n", fp);

  uint64_t count = 0;
  for (scalar s in list) {
    fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" Name=\"%s\" format=\"appended\" offset=\"%ld\"/>\n", s.name, count);
    count += (no_cells * sizeof (double)) + sizeof (uint64_t);
  }
  for (vector v in vlist) {
    fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"3\"  format=\"appended\" offset=\"%ld\"/>\n", v.x.name, count);
    count += (3 * no_cells * sizeof (double)) + sizeof (uint64_t);
  }
  fputs ("\t\t\t </CellData>\n", fp);
  fputs ("\t\t\t <Points>\n", fp);
  fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\"  format=\"appended\" offset=\"%ld\"/>\n", count);
  fputs ("\t\t\t </Points>\n", fp);
  fputs ("\t\t\t <Cells>\n", fp);
  count += (3 * no_points * sizeof (double)) + sizeof (uint64_t);

  uint8_t type=9, noffset=4;
  uint64_t offset, connectivity[noffset];

  fprintf (fp,"\t\t\t\t <DataArray type=\"UInt64\" Name=\"offsets\" format=\"appended\" offset=\"%ld\"/>\n", count);
  count +=  (no_cells * sizeof (uint64_t)) + sizeof (uint64_t);

  fprintf (fp,"\t\t\t\t <DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\"%ld\"/>\n", count);
  count +=  (no_cells * sizeof (uint8_t)) + sizeof (uint64_t);

  fprintf (fp,"\t\t\t\t <DataArray type=\"UInt64\" Name=\"connectivity\" format=\"appended\" offset=\"%ld\"/>\n", count);
  count +=  (no_cells * noffset * sizeof (uint64_t)) + sizeof (uint64_t);

  fputs ("\t\t\t </Cells>\n", fp);
  fputs ("\t\t </Piece>\n", fp);
  fputs ("\t </UnstructuredGrid>\n", fp);
  fputs ("\t <AppendedData encoding=\"raw\">\n", fp);
  fputs ("_", fp);

  uint64_t block_len=no_cells*sizeof (double);
  for (scalar s in list) {
    fwrite (&block_len, sizeof (uint64_t), 1, fp);
    foreach(serial, noauto){
      if (per_mask[]){
        double sval;
        if (n.x == 1)
          sval = 0.5*(val(s) + val(s,1,0,0));
        else if (n.y == 1)
          sval = 0.5*(val(s) + val(s,0,1,0));
        else
          sval = 0.5*(val(s) + val(s,0,0,1));
        fwrite (&sval, sizeof (double), 1, fp);
      }
    }
  }
  block_len=no_cells*3*sizeof (double);
  for (vector v in vlist) {
    fwrite (&block_len, sizeof (uint64_t), 1, fp);
    foreach(serial, noauto){
      if (per_mask[]){
        double xval, yval, zval;
        if (n.x == 1) {
          xval = 0.5*(val(v.x) + val(v.x,1,0,0));
          yval = 0.5*(val(v.y) + val(v.y,1,0,0));
          zval = 0.5*(val(v.z) + val(v.z,1,0,0));
        }
        else if (n.y == 1) {
          xval = 0.5*(val(v.x) + val(v.x,0,1,0));
          yval = 0.5*(val(v.y) + val(v.y,0,1,0));
          zval = 0.5*(val(v.z) + val(v.z,0,1,0));
        }
        else {
          xval = 0.5*(val(v.x) + val(v.x,0,0,1));
          yval = 0.5*(val(v.y) + val(v.y,0,0,1));
          zval = 0.5*(val(v.z) + val(v.z,0,0,1));
        }
        fwrite (&xval, sizeof (double), 1, fp);
        fwrite (&yval, sizeof (double), 1, fp);
        fwrite (&zval, sizeof (double), 1, fp);
      }
    }
  }
  block_len=no_points*3*sizeof (double);
  fwrite (&block_len, sizeof (uint64_t), 1, fp);
  foreach_vertex(serial, noauto){

    double alpha = (_alpha - n.x*x - n.y*y - n.z*z)/Delta;

    if (fabs(alpha) > 0.87)
      continue;

    fwrite (&x, sizeof (double), 1, fp);
    fwrite (&y, sizeof (double), 1, fp);
    fwrite (&z, sizeof (double), 1, fp);
  }

  block_len=no_cells*1*sizeof(uint64_t);
  fwrite (&block_len, sizeof (uint64_t), 1, fp);
  for (int i = 0; i < no_cells; i++){
    offset = (i+1)*noffset;
    fwrite (&offset, sizeof (int64_t), 1, fp);
  }

  block_len=no_cells*1*sizeof(uint8_t);
  fwrite (&block_len, sizeof (uint64_t), 1, fp);
  for (int i = 0; i < no_cells; i++)
    fwrite (&type, sizeof (int8_t), 1, fp);

  block_len=no_cells*noffset*sizeof(uint64_t);
  fwrite (&block_len, sizeof (uint64_t), 1, fp);
  foreach(serial, noauto){
    if (per_mask[]){
      if (n.x == 1) {
        connectivity[0] = (uint64_t)marker[1,0,0];
        connectivity[1] = (uint64_t)marker[1,1,0];
        connectivity[2] = (uint64_t)marker[1,1,1];
        connectivity[3] = (uint64_t)marker[1,0,1];
      }
      else if (n.y == 1) {
        connectivity[0] = (uint64_t)marker[0,1,0];
        connectivity[1] = (uint64_t)marker[1,1,0];
        connectivity[2] = (uint64_t)marker[1,1,1];
        connectivity[3] = (uint64_t)marker[0,1,1];
      }
      else {
        connectivity[0] = (uint64_t)marker[0,0,1];
        connectivity[1] = (uint64_t)marker[1,0,1];
        connectivity[2] = (uint64_t)marker[1,1,1];
        connectivity[3] = (uint64_t)marker[0,1,1];
      }
      fwrite (&connectivity, sizeof (uint64_t), noffset, fp);
    }
  }
  fputs ("\t\n", fp);
  fputs ("\t </AppendedData>\n", fp);
  fputs ("</VTKFile>\n", fp);
  fflush (fp);
#if defined(_OPENMP)
  omp_set_num_threads(num_omp);
#endif
}

@if _MPI
void output_pvtu_plane (scalar * list, vector * vlist, char * subname, coord n, double _alpha)
{
  char name[80];
  FILE * fp;
  if (pid() == 0) {
    sprintf(name, "%s.pvtu", subname);
    fp = fopen(name, "w");

    fputs ("<?xml version=\"1.0\"?>\n"
    "<VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);
    fputs ("\t <PUnstructuredGrid GhostLevel=\"0\">\n", fp);
    fputs ("\t\t\t <PCellData Scalars=\"scalars\">\n", fp);
    for (scalar s in list) {
      fprintf (fp,"\t\t\t\t <PDataArray type=\"Float64\" Name=\"%s\" format=\"appended\">\n", s.name);
      fputs ("\t\t\t\t </PDataArray>\n", fp);
    }
    for (vector v in vlist) {
      fprintf (fp,"\t\t\t\t <PDataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"%s\" format=\"appended\">\n", v.x.name);
      fputs ("\t\t\t\t </PDataArray>\n", fp);
    }
    fputs ("\t\t\t </PCellData>\n", fp);
    fputs ("\t\t\t <PPoints>\n", fp);
    fputs ("\t\t\t\t <PDataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"appended\">\n", fp);
    fputs ("\t\t\t\t </PDataArray>\n", fp);
    fputs ("\t\t\t </PPoints>\n", fp);

    for (int i = 0; i < npe(); i++)
      fprintf (fp, "<Piece Source=\"%s_n%3.3d.vtu\"/> \n", subname, i);

    fputs ("\t </PUnstructuredGrid>\n", fp);
    fputs ("</VTKFile>\n", fp);
    fflush (fp);
    fclose (fp);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  sprintf(name, "%s_n%3.3d", subname, pid());
  output_vtu_plane_pid (list, vlist, name, n, _alpha);
}
@endif

trace
void output_vtu_plane (scalar * list, vector * vlist, char * subname, coord n, double _alpha)
{
  @if _MPI
  output_pvtu_plane (list, vlist, subname, n, _alpha);
  @else
  output_vtu_plane_pid (list, vlist, subname, n, _alpha);
  @endif
}

#endif
