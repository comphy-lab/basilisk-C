#include <stdint.h>

/*
This function writes one XML file which allows to read the *.vtu files generated
by output_vtu_ascii() when used in MPI. Tested in (quad- and oct-)trees
using MPI.
*/
void output_pvtu_ascii (scalar * list, vector * vlist, int n, FILE * fp, char * subname)
{
    fputs ("<?xml version=\"1.0\"?>\n"
    "<VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);
    fputs ("\t <PUnstructuredGrid GhostLevel=\"0\">\n", fp);
    fputs ("\t\t\t <PCellData Scalars=\"scalars\">\n", fp);
    for (scalar s in list) {
      fprintf (fp,"\t\t\t\t <PDataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", s.name);
      fputs ("\t\t\t\t </PDataArray>\n", fp);
    }
    for (vector v in vlist) {
      fprintf (fp,"\t\t\t\t <PDataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"%s\" format=\"ascii\">\n", v.x.name);
      fputs ("\t\t\t\t </PDataArray>\n", fp);
    }
    fputs ("\t\t\t </PCellData>\n", fp);
    fputs ("\t\t\t <PPoints>\n", fp);
    fputs ("\t\t\t\t <PDataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n", fp);
    fputs ("\t\t\t\t </PDataArray>\n", fp);
    fputs ("\t\t\t </PPoints>\n", fp);

    for (int i = 0; i < npe(); i++)
      fprintf (fp, "<Piece Source=\"%s_n%3.3d.vtu\"/> \n", subname, i);

    fputs ("\t </PUnstructuredGrid>\n", fp);
    fputs ("</VTKFile>\n", fp);
}

/*
This function writes one XML VTK file per PID process of type unstructured grid
(*.vtu) which can be read using Paraview. File stores scalar and vector fields
defined at the center points. Results are recorded on ASCII format. If one writes
one *.vtu file per PID process this function may be combined with
output_pvtu_ascii() above to read in parallel. Tested in (quad- and oct-)trees
using MPI. Also works with solids (when not using MPI).
*/
void output_vtu_ascii(scalar * list, vector * vlist, int n, FILE * fp)
{
#if defined(_OPENMP)
  int num_omp = omp_get_max_threads();
  omp_set_num_threads(1);
#endif

  //* This is to match the number of points and vertices when using periodic conditions
  scalar per_mask[];
  foreach(){
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
  int no_points = 0, no_cells=0 ;
  foreach_vertex(){
    marker[] = _k;
    no_points += 1;
  }
  foreach()
    if (per_mask[])
      no_cells += 1;

  fputs ("<?xml version=\"1.0\"?>\n"
  "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);
  fputs ("\t <UnstructuredGrid>\n", fp);
  fprintf (fp,"\t\t <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", no_points, no_cells);
  fputs ("\t\t\t <CellData Scalars=\"scalars\">\n", fp);
  for (scalar s in list) {
    fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", s.name);
    foreach()
      if (per_mask[])
        fprintf (fp, "\t\t\t\t\t %g\n", val(s));

    fputs ("\t\t\t\t </DataArray>\n", fp);
  }
  for (vector v in vlist) {
    fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"%s\" format=\"ascii\">\n", v.x.name);
    foreach(){
      if (per_mask[]){
#if dimension == 2
        fprintf (fp, "\t\t\t\t\t %g %g 0.\n", val(v.x), val(v.y));
#elif dimension == 3
        fprintf (fp, "\t\t\t\t\t %g %g %g\n", val(v.x), val(v.y), val(v.z));
#endif
      }
    }
    fputs ("\t\t\t\t </DataArray>\n", fp);
  }
  fputs ("\t\t\t </CellData>\n", fp);
  fputs ("\t\t\t <Points>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n", fp);
  foreach_vertex(){
#if dimension == 2
      fprintf (fp, "\t\t\t\t\t %g %g 0\n", x, y);
#elif dimension == 3
      fprintf (fp, "\t\t\t\t\t %g %g %g\n", x, y, z);
#endif
  }
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t </Points>\n", fp);
  fputs ("\t\t\t <Cells>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n", fp);
  foreach(){
    if (per_mask[]){
#if dimension == 2
      fprintf (fp, "\t\t\t\t\t %u %u %u %u \n", (unsigned int )marker[], (unsigned int )marker[1,0], (unsigned int )marker[1,1], (unsigned int )marker[0,1]);
#elif dimension == 3
      fprintf (fp, "\t\t\t\t\t %u %u %u %u %u %u %u %u\n", (unsigned int )marker[], (unsigned int )marker[1,0,0], (unsigned int )marker[1,1,0], (unsigned int )marker[0,1,0],(unsigned int )marker[0,0,1], (unsigned int )marker[1,0,1], (unsigned int )marker[1,1,1], (unsigned int )marker[0,1,1]);
#endif
    }
  }
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n", fp);

  for (int i = 1; i < no_cells+1; i++){
#if dimension == 2
    fprintf (fp, "\t\t\t\t\t %d \n", i*4);
#elif dimension == 3
    fprintf (fp, "\t\t\t\t\t %d \n", i*8);
#endif
  }
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n", fp);
  foreach(){
    if (per_mask[]){
#if dimension == 2
      fputs ("\t\t\t\t\t 9 \n", fp);
#elif dimension == 3
      fputs ("\t\t\t\t\t 12 \n", fp);
#endif
    }
  }
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t </Cells>\n", fp);
  fputs ("\t\t </Piece>\n", fp);
  fputs ("\t </UnstructuredGrid>\n", fp);
  fputs ("</VTKFile>\n", fp);
  fflush (fp);
#if defined(_OPENMP)
  omp_set_num_threads(num_omp);
#endif
}






/*
This function writes a slice defined by $n={n_x, n_y, n_z}$ and $\alpha$
using a similar approach to `view.h` into one XML VTK file per PID process of
type unstructured grid (*.vtu) in ascii format which can be read using Paraview.
Only works for x, y, and/or z planes. Here, _alpha is assumed to intersect a cell
face and must be a multiple of L0/1<<MINLEVEL. This also means that scalar and
vector fields are written at the corresponding face values.
*/
#if dimension > 2
void output_vtu_plane_ascii (scalar * list, vector * vlist, FILE * fp, coord n, double _alpha)
{
#if defined(_OPENMP)
  int num_omp = omp_get_max_threads();
  omp_set_num_threads(1);
#endif

  scalar per_mask[];
  foreach(){
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
  int no_points = 0, no_cells=0 ;
  foreach_vertex(){
    marker[] = 0.;
    double alpha = (_alpha - n.x*x - n.y*y - n.z*z)/Delta;

    if (fabs(alpha) > 0.87)
      continue;

    marker[] = no_points;
    no_points += 1;
  }

  foreach()
    if (per_mask[])
      no_cells += 1;

  fputs ("<?xml version=\"1.0\"?>\n"
  "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);
  fputs ("\t <UnstructuredGrid>\n", fp);
  fprintf (fp,"\t\t <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", no_points, no_cells);
  fputs ("\t\t\t <CellData Scalars=\"scalars\">\n", fp);
  for (scalar s in list) {
    fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", s.name);

    foreach(){
      if (per_mask[]){
        if (n.x == 1)
          fprintf (fp, "\t\t\t\t\t %g\n", 0.5*val(s) + 0.5*val(s,1,0,0));
        else if (n.y == 1)
          fprintf (fp, "\t\t\t\t\t %g\n", 0.5*val(s) + 0.5*val(s,0,1,0));
        else
          fprintf (fp, "\t\t\t\t\t %g\n", 0.5*val(s) + 0.5*val(s,0,0,1));
      }
    }

    fputs ("\t\t\t\t </DataArray>\n", fp);
  }
  for (vector v in vlist) {
    fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"%s\" format=\"ascii\">\n", v.x.name);

    foreach(){
      if (per_mask[]){
        if (n.x == 1)
          fprintf (fp, "\t\t\t\t\t %g %g %g\n", 0.5*(val(v.x) + val(v.x,1,0,0)), 0.5*(val(v.y) + val(v.y,1,0,0)), 0.5*(val(v.z) + val(v.z,1,0,0)));
        else if (n.y == 1)
          fprintf (fp, "\t\t\t\t\t %g %g %g\n", 0.5*(val(v.x) + val(v.x,0,1,0)), 0.5*(val(v.y) + val(v.y,0,1,0)), 0.5*(val(v.z) + val(v.z,0,1,0)));
        else
          fprintf (fp, "\t\t\t\t\t %g %g %g\n", 0.5*(val(v.x) + val(v.x,0,0,1)), 0.5*(val(v.y) + val(v.y,0,0,1)), 0.5*(val(v.z) + val(v.z,0,0,1)));
      }
    }
    fputs ("\t\t\t\t </DataArray>\n", fp);
  }
  fputs ("\t\t\t </CellData>\n", fp);
  fputs ("\t\t\t <Points>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n", fp);

  foreach_vertex() {
    double alpha = (_alpha - n.x*x - n.y*y - n.z*z)/Delta;
    if (fabs(alpha) > 0.87){
      continue;
    }
    fprintf (fp, "\t\t\t\t\t %g %g %g\n", x, y, z);
  }

  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t </Points>\n", fp);
  fputs ("\t\t\t <Cells>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n", fp);

  foreach(){
    if (per_mask[]){
      if (n.x == 1)
        fprintf (fp, "\t\t\t\t\t %u %u %u %u \n", (unsigned int )marker[1,0,0], (unsigned int )marker[1,1,0], (unsigned int )marker[1,1,1], (unsigned int )marker[1,0,1]);
      else if (n.y == 1)
        fprintf (fp, "\t\t\t\t\t %u %u %u %u \n", (unsigned int )marker[0,1,0], (unsigned int )marker[1,1,0], (unsigned int )marker[1,1,1], (unsigned int )marker[0,1,1]);
      else
        fprintf (fp, "\t\t\t\t\t %u %u %u %u \n", (unsigned int )marker[0,0,1], (unsigned int )marker[1,0,1], (unsigned int )marker[1,1,1], (unsigned int )marker[0,1,1]);
    }
  }

  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n", fp);

  for (int i = 1; i < no_cells+1; i++){
    fprintf (fp, "\t\t\t\t\t %d \n", i*4);
  }

  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n", fp);

  foreach(){
    if (per_mask[]){
      fputs ("\t\t\t\t\t 9 \n", fp);
    }
  }

  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t </Cells>\n", fp);
  fputs ("\t\t </Piece>\n", fp);
  fputs ("\t </UnstructuredGrid>\n", fp);
  fputs ("</VTKFile>\n", fp);
  fflush (fp);
#if defined(_OPENMP)
  omp_set_num_threads(num_omp);
#endif
}
