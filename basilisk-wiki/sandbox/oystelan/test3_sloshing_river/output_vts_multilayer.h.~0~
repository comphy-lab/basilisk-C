#ifndef BASILISK_HEADER_1
#define BASILISK_HEADER_1
#line 1 "./../output_vts_multilayer.h"
/*
This function writes one XML VTK file per PID process of type unstructured grid
(*.vtu) which can be read using Paraview. File stores scalar and vector fields
defined at the center points. Results are recorded on ASCII format. If one writes
one *.vts file per PID process this function may be combined with
output_pvtu_ascii() above to read in parallel. Tested in (quad- and oct-)trees
using MPI. Also works with solids (when not using MPI).
*/
double dt_start = 0;
scalar id_boundary;
extern scalar w;
extern scalar phi;
extern double dry;
//extern int num_omp;
//extern scalar hpres;
// Define alternative traversal order to match with the x ordering of vts (basilisk uses z ordering by default)


macro2 foreach_xorder(char flags = 0, Reduce reductions = None) {
  //OMP_PARALLEL (reductions) {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = depth();
    SET_DIMENSIONS();
    int _k;
    //OMP(omp for schedule(static))
#if dimension == 1
      for (_k = GHOSTS; _k < point.n.x + GHOSTS; _k++) {
	point.i = _k;
#else
    for (_k = GHOSTS; _k < point.n.y + GHOSTS; _k++) {
    point.j = _k;
    for (point.i = GHOSTS; point.i < point.n.x + GHOSTS; point.i++){
#endif
	    {...}
#if dimension > 1
    }
#endif
      }
  //}
}


macro2 foreach_face_xorder (char flags = 0, Reduce reductions = None,
			     const char * order = "xyz")
{
  //OMP_PARALLEL (reductions) {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = depth();
    SET_DIMENSIONS();
    int _k;
    //OMP(omp for schedule(static))
#if dimension == 1
      for (_k = GHOSTS; _k <= point.n.x + GHOSTS; _k++) {
	point.i = _k;
#else
    for (_k = GHOSTS; _k <= point.n.y + GHOSTS; _k++) {
	point.j = _k;
	for (point.i = GHOSTS; point.i <= point.n.x + GHOSTS; point.i++){
#endif
        //POINT_VARIABLES();
	    {...}
#if dimension > 1
    }
#endif
      }
  //}
}


macro2 foreach_vertex_xorder (char flags = 0, Reduce reductions = None,
			     const char * order = "xyz")
{
  //OMP_PARALLEL (reductions) {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = depth();
    SET_DIMENSIONS();
    int _k;
    //OMP(omp for schedule(static))
#if dimension == 1
      for (_k = GHOSTS; _k <= point.n.x + GHOSTS; _k++) {
	point.i = _k;
#else
    for (_k = GHOSTS; _k <= point.n.y + GHOSTS; _k++) {
	point.j = _k;
	for (point.i = GHOSTS; point.i <= point.n.x + GHOSTS; point.i++){
#endif  
        POINT_VARIABLES();
        x -= Delta / 2.;
#if dimension > 1  
    y -= Delta / 2.;
#endif
	    {...}
#if dimension > 1 
    }
#endif        
      }
  //}
}


macro2 foreach_xorder_boundary(char flags = 0, Reduce reductions = None) {
  //OMP_PARALLEL (reductions) {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = depth();
    SET_DIMENSIONS();
    int _k;
    //OMP(omp for schedule(static))
#if dimension == 1
      for (_k = GHOSTS-BGHOSTS; _k < point.n.x + GHOSTS + BGHOSTS; _k++) {
	point.i = _k;
#else
    for (_k = GHOSTS-BGHOSTS; _k < point.n.y + GHOSTS + BGHOSTS; _k++) {
    point.j = _k;
    for (point.i = GHOSTS-BGHOSTS; point.i < point.n.x + GHOSTS + BGHOSTS; point.i++){
#endif
	    {...}
#if dimension > 1
    }
#endif
      }
  //}
}

macro2 foreach_vertex_xorder_boundary (char flags = 0, Reduce reductions = None,
			     const char * order = "xyz")
{
  //OMP_PARALLEL (reductions) {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = depth();
    SET_DIMENSIONS();
    int _k;
    //OMP(omp for schedule(static))
#if dimension == 1
      for (_k = GHOSTS-BGHOSTS; _k <= point.n.x + GHOSTS + BGHOSTS; _k++) {
	point.i = _k;
#else
    for (_k = GHOSTS-BGHOSTS; _k <= point.n.y + GHOSTS + BGHOSTS; _k++) {
	point.j = _k;
	for (point.i = GHOSTS-BGHOSTS; point.i <= point.n.x + GHOSTS + BGHOSTS; point.i++){
#endif  
        POINT_VARIABLES();
        x -= Delta / 2.;
#if dimension > 1  
    y -= Delta / 2.;
#endif
	    {...}
#if dimension > 1 
    }
#endif        
      }
  //}
}


macro2 foreach_xonly(char flags = 0, Reduce reductions = None) {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = depth();
    SET_DIMENSIONS();
    int _k;
    //OMP(omp for schedule(static))
    for (_k = GHOSTS; _k < point.n.x + GHOSTS; _k++) {
	point.i = _k;
	    {...}
    }
}
macro2 foreach_vertex_xonly(char flags = 0, Reduce reductions = None) {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = depth();
    SET_DIMENSIONS();
    int _k;
    //OMP(omp for schedule(static))
    for (_k = GHOSTS; _k <= point.n.x + GHOSTS; _k++) {
	point.i = _k;
    POINT_VARIABLES();
        x -= Delta / 2.;
	    {...}
    }
}



/**
output is a vtk structure file (.vts) in ascii formatting. Cell data
is velocity (and any extra fields you wire in below). Set
`include_boundaries=true` to also emit the BGHOSTS-wide boundary ghost
cells, which is useful for debugging boundary conditions. */
void output_vts_ascii_all_layers (FILE * fp, bool include_boundaries = false)
{
  // The grid extent depends on whether we include the boundary ghost cells.
  // foreach_xorder         iterates interior cells only       (N cells/dim)
  // foreach_xorder_boundary iterates with BGHOSTS extra cells (N+2*BGHOSTS cells/dim)
  int Nx = include_boundaries ? N + 2*BGHOSTS : N;

  fputs("<?xml version=\"1.0\"?>\n"
        "<VTKFile type=\"StructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);

#if dimension == 1
  fprintf(fp, "\t <StructuredGrid WholeExtent=\"%d %d %d %d 0 0\">\n", 0, Nx, 0, nl);
  fprintf(fp, "\t\t <Piece Extent=\"%d %d %d %d  0 0\">\n", 0, Nx, 0, nl);
#endif
#if dimension == 2
  fprintf(fp, "\t <StructuredGrid WholeExtent=\"%d %d %d %d %d %d\">\n", 0, Nx, 0, Nx, 0, nl);
  fprintf(fp, "\t\t <Piece Extent=\"%d %d %d %d %d %d\">\n", 0, Nx, 0, Nx, 0, nl);
#endif

  // Velocity vector field (cell-centered)
  fputs("\t\t\t <CellData Scalars=\"scalars\">\n", fp);
  fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"ascii\">\n");

  if (include_boundaries) {
    foreach_layer() {
      foreach_xorder_boundary(serial, noauto) {
#if dimension == 1
        fprintf(fp, "%g %g 0.\n", u.x[], w[]);
#endif
#if dimension == 2
        fprintf(fp, "%g %g %g\n", u.x[], u.y[], w[]);
#endif
      }
    }
  }
  else {
    foreach_layer() {
      foreach_xorder(serial) {
#if dimension == 1
        fprintf(fp, "%g %g 0.\n", u.x[], w[]);
#endif
#if dimension == 2
        fprintf(fp, "%g %g %g\n", u.x[], u.y[], w[]);
#endif
      }
    }
  }
  fputs("\t\t\t\t </DataArray>\n", fp);
  fputs("\t\t\t </CellData>\n", fp);

  // Vertex coordinates. The interior path uses simple averaging; the
  // boundary path needs to handle the outermost vertices, which only
  // have one neighbour cell on the inside.
  fputs("\t\t\t <Points>\n", fp);
  fputs("\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n", fp);

  scalar hh[];

  if (include_boundaries) {
#if dimension == 1
    foreach_vertex_xorder_boundary(serial, noauto) {
      if (point.i == (GHOSTS - BGHOSTS))
        hh[] = zb[];
      else if (point.i == (point.n.x + GHOSTS + BGHOSTS))
        hh[] = zb[-1];
      else
        hh[] = (zb[] + zb[-1]) / 2.;
      fprintf(fp, "%12.4f %12.4f 0.\n", x, hh[]);
    }
    foreach_layer() {
      foreach_vertex_xorder_boundary(serial, noauto) {
        if (point.i == (GHOSTS - BGHOSTS))
          hh[] += h[];
        else if (point.i == (point.n.x + GHOSTS + BGHOSTS))
          hh[] += h[-1];
        else
          hh[] += (h[] + h[-1]) / 2.;
        fprintf(fp, "%12.4f %12.4f 0.\n", x, hh[]);
      }
    }
#endif
#if dimension == 2
    foreach_vertex_xorder_boundary(serial, noauto) {
      hh[] = (zb[] + zb[-1] + zb[0,-1] + zb[-1,-1]) / 4.;
      fprintf(fp, "%12.4f %12.4f %12.4f\n", x, y, hh[]);
    }
    foreach_layer() {
      foreach_vertex_xorder_boundary(serial, noauto) {
        hh[] += (h[] + h[-1] + h[0,-1] + h[-1,-1]) / 4.;
        fprintf(fp, "%12.4f %12.4f %12.4f\n", x, y, hh[]);
      }
    }
#endif
  }
  else {
#if dimension == 1
    foreach_vertex_xorder(serial) {
      hh[] = (zb[] + zb[-1]) / 2.;
      fprintf(fp, "%12.4f %12.4f 0.\n", x, hh[]);
    }
    foreach_layer() {
      foreach_vertex_xorder(serial) {
        hh[] += (h[] + h[-1]) / 2.;
        fprintf(fp, "%12.4f %12.4f 0.\n", x, hh[]);
      }
    }
#endif
#if dimension == 2
    foreach_vertex_xorder(serial) {
      hh[] = (zb[] + zb[-1] + zb[0,-1] + zb[-1,-1]) / 4.;
      fprintf(fp, "%12.4f %12.4f %12.4f\n", x, y, hh[]);
    }
    foreach_layer() {
      foreach_vertex_xorder(serial) {
        hh[] += (h[] + h[-1] + h[0,-1] + h[-1,-1]) / 4.;
        fprintf(fp, "%12.4f %12.4f %12.4f\n", x, y, hh[]);
      }
    }
#endif
  }

  fputs("\t\t\t\t </DataArray>\n", fp);
  fputs("\t\t\t </Points>\n", fp);
  fputs("\t\t </Piece>\n", fp);

  // Time value as field data
  fprintf(fp, "\t\t <FieldData>\n");
  fprintf(fp, "\t\t\t <DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"ascii\" RangeMin=\"%.3f\" RangeMax=\"%.3f\">\n", t + dt_start, t + dt_start);
  fprintf(fp, "\t\t\t %.3f\n", t + dt_start);
  fprintf(fp, "\t\t\t </DataArray>\n");
  fprintf(fp, "\t\t </FieldData>\n");

  fputs("\t </StructuredGrid>\n", fp);
  fputs("</VTKFile>\n", fp);
  fflush(fp);
}


/*
    does the same as function above, only for a limited area, limited by the bounding box 
    
    bounds[4] = {xmin, xmax, ymin, ymax}

    The function writes only cell where cell centers are within the box.
    (checks if xmin <= x <= xmax and ymin <= y <= ymax.)


*/

/**output is a vtk structure file (.vts) in binary (appended) formating, which is faster to read and process 
 * and stores much more compact . Field data is velocity only */
// void output_vts_bin_all_layers(FILE* fp)
// {

        
//         int count = 0;
//         // MULTIGRID

//         fputs("<?xml version=\"1.0\"?>\n"
//             "<VTKFile type=\"StructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);

// #if dimension == 1
//     fprintf(fp, "\t <StructuredGrid WholeExtent=\"%d %d %d %d 0 0\">\n", 0, N, 0, nl);
// #endif

// #if dimension == 2
//     fprintf(fp, "\t <StructuredGrid WholeExtent=\"%d %d %d %d %d %d\">\n", 0, N, 0, N, 0, nl);
// #endif

    
//     // write time value
//     fprintf(fp, "\t\t <FieldData> \n");
//     fprintf(fp, "\t\t\t <DataArray type = \"Float64\" Name = \"%s\" NumberOfTuples = \"1\" format = \"ascii\" RangeMin = \"%.3f\" RangeMax = \"%.3f\"> \n", "TimeValue", t + dt_start, t + dt_start);
//     fprintf(fp, "\t\t\t %.3f \n", t + dt_start);
//     fprintf(fp, "\t\t\t </DataArray > \n");
//     fprintf(fp, "\t\t </FieldData> \n");
    
// #if dimension == 1
//     fprintf(fp, "\t\t <Piece Extent=\"%d %d %d %d  0 0\">\n", 0, N, 0, nl);
// #endif

// #if dimension == 2
//     fprintf(fp, "\t\t <Piece Extent=\"%d %d %d %d %d %d\">\n", 0, N, 0, N, 0, nl);
// #endif

//     // Loop over velocity data and store kinematics in cell vector stucture
//     fputs("\t\t\t <CellData Scalars=\"scalars\">\n", fp);
//     fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"appended\" offset=\"%d\">\n", count);
//     fputs("\t\t\t\t </DataArray>\n", fp);
//     fputs("\t\t\t </CellData>\n", fp);
    

//     // Coordinates 
//     fputs("\t\t\t <Points>\n", fp);
// #if dimension == 1
//     count += ((N * nl * 3) + 1) * 8;
//     fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\">\n", count);
// #endif
// #if dimension == 2
//     count += ((N * N * nl * 3) + 1) * 8;
//     fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\">\n", count);
// #endif

//     fputs("\t\t\t\t </DataArray>\n", fp);
//     fputs("\t\t\t </Points>\n", fp);
//     fputs("\t\t </Piece>\n", fp);

//     fputs("\t </StructuredGrid>\n", fp);
//     fputs("\t <AppendedData encoding=\"raw\">\n", fp);
//     fputs("_", fp);
    
//     unsigned long long block_len = 8;
//     // append kinematics
// #if dimension == 1
//     block_len = 8 * 3 * N * nl;
//     double dummy = 0.;
// #endif
// #if dimension == 2
//     block_len = 8 * 3 * N * N * nl;
// #endif
//     fwrite(&block_len, sizeof(unsigned long long), 1, fp);
//     foreach_layer() {
//         foreach_xorder(serial, noauto) {
// #if dimension == 1
//             fwrite(&u.x[], sizeof(double), 1, fp);
//             fwrite(&w[], sizeof(double), 1, fp);
//             fwrite(&dummy, sizeof(double), 1, fp);
// #endif
// #if dimension == 2
//             fwrite(&u.x[], sizeof(double), 1, fp);
//             fwrite(&u.y[], sizeof(double), 1, fp);
//             fwrite(&w[], sizeof(double), 1, fp);
// #endif
//         }

//     }
//     // append points
//     // Coordinates 
// #if dimension == 1
//     block_len = 8 * 3 * (N + 1) * (nl + 1);
//     fwrite(&block_len, sizeof(unsigned long long), 1, fp);
//     // first do the bottom coordinates stored in zb
//     scalar hh[];
//     foreach_vertex_xorder(serial, noauto) {
//         hh[] = (zb[] + zb[-1]) / 2.;
//         fwrite(&x, sizeof(double), 1, fp);
//         fwrite(&hh[], sizeof(double), 1, fp);
//         fwrite(&dummy, sizeof(double), 1, fp);
//     }

//     foreach_layer() {
//         //for (scalar h_layer in hl){
//         foreach_vertex_xorder(serial, noauto) {
//             hh[] += (h[] + h[-1]) / 2.;
//             fwrite(&x, sizeof(double), 1, fp);
//             fwrite(&hh[], sizeof(double), 1, fp);
//             fwrite(&dummy, sizeof(double), 1, fp);
//         }
//     }
// #endif
// #if dimension == 2
//     block_len = 8 * 3 * (N + 1) * (N + 1) * (nl + 1);
//     fwrite(&block_len, sizeof(unsigned long long), 1, fp);

//     // first do the bottom coordinates stored in zb
//     scalar hh[];
//     foreach_vertex_xorder(serial, noauto) {
//         hh[] = (zb[] + zb[-1] + zb[0, -1] + zb[-1, -1]) / 4.;
//         fwrite(&x, sizeof(double), 1, fp);
//         fwrite(&y, sizeof(double), 1, fp);
//         fwrite(&hh[], sizeof(double), 1, fp);
//     }

//     foreach_layer() {
//         //for (scalar h_layer in hl){
//         foreach_vertex_xorder(serial, noauto) {
//             hh[] += (h[] + h[-1] + h[0, -1] + h[-1, -1]) / 4.;
//             fwrite(&x, sizeof(double), 1, fp);
//             fwrite(&y, sizeof(double), 1, fp);
//             fwrite(&hh[], sizeof(double), 1, fp);
//         }
//     }
// #endif


//     fputs("\t\n", fp);
//     fputs("\t </AppendedData>\n", fp);
//     fputs("</VTKFile>\n", fp);
//     fflush(fp);
// }

/*
    does the same as function above, only for a limited area, limited by the bounding box

    bounds[4] = {xmin, xmax, ymin, ymax}

    The function writes only cell where cell centers are within the box.
    (checks if xmin <= x <= xmax and ymin <= y <= ymax.)


*/


void output_vts_ascii_all_layers_bounded(FILE* fp, double bounds[])
{
    
        // MULTIGRID

        fputs("<?xml version=\"1.0\"?>\n"
            "<VTKFile type=\"StructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);

#if dimension == 1
    // loop over cells and find which are within bounds
    scalar cells_in_bounds[];
    foreach_xorder(serial, noauto) {
        if (x >= bounds[0] && x <= bounds[1])
            cells_in_bounds[] = 1;
        else
            cells_in_bounds[] = 0;
    }

    vertex scalar vertex_in_bounds[];
    int testcount = 0;
    foreach_vertex_xorder(serial, noauto) {
        if (x >= bounds[0] - Delta / 2. && x <= bounds[1] + Delta / 2.) {
            vertex_in_bounds[] = 1;
            testcount++;
        }
        else {
            vertex_in_bounds[] = 0;
        }
    }
    
    // count number of vertex points in x direction which is within bounds
    int nx = 0;
    double dd = L0 / N;
    for (double xx = X0; xx <= (X0 + L0); xx += dd) {
        if (xx >= (bounds[0] - dd / 2.) && xx <= (bounds[1] + dd / 2.)) {
            nx++;
        }
    }

    fprintf(fp, "\t <StructuredGrid WholeExtent=\"%d %d %d %d 0 0\">\n", 0, nx, 0, nl);
    fprintf(fp, "\t\t <Piece Extent=\"%d %d %d %d  0 0\">\n", 0, nx, 0, nl);
#endif

#if dimension == 2
    // loop over cells and find which are within bounds
    scalar cells_in_bounds[];
    foreach_xorder(serial, noauto) {
        if (x >= bounds[0] && x <= bounds[1] && y >= bounds[2] && y <= bounds[3])
            cells_in_bounds[] = 1;
        else
            cells_in_bounds[] = 0;
    }

    vertex scalar vertex_in_bounds[];
    int testcount = 0;
    foreach_vertex_xorder(serial, noauto) {
        if (x >= bounds[0] - Delta / 2. && x <= bounds[1] + Delta / 2. && y >= bounds[2] - Delta / 2. && y <= bounds[3] + Delta / 2.) {
            vertex_in_bounds[] = 1;
            testcount++;
        }
        else {
            vertex_in_bounds[] = 0;
        }
    }

    // count number of vertex points in x direction which is within bounds
    int nx = 0;
    double dd = L0/N;
    for (double xx = X0; xx <= (X0 + L0); xx += dd) {
        if (xx >= (bounds[0] - dd / 2.) && xx <= (bounds[1] + dd / 2.)) {
            nx++;
        }
    }
    nx--;
    // count number of vertex points in x direction which is within bounds
    int ny = 0;
    for (double yy = Y0; yy <= (Y0 + L0); yy += dd) {
        if (yy >= (bounds[2] - dd / 2.) && yy <= (bounds[3] + dd / 2.)) {
            ny++;
}
    }
    ny--;

    fprintf(stdout, "nx: %d, ny: %d, n: %d, %f %f %f %d\n", nx, ny,N, L0, X0, Y0, testcount);

    fprintf(fp, "\t <StructuredGrid WholeExtent=\"%d %d %d %d %d %d\">\n", 0, nx, 0, ny, 0, nl);
    fprintf(fp, "\t\t <Piece Extent=\"%d %d %d %d %d %d\">\n", 0, nx, 0, ny, 0, nl);
#endif

    // Loop over velocity data and store kinematics in cell vector stucture
    fputs("\t\t\t <CellData Scalars=\"scalars\">\n", fp);
    fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"ascii\">\n");

    //for (u_layer, w_layer in ul, wl){  
    foreach_layer() {
        foreach_xorder(serial, noauto) {


#if dimension == 1
            if (cells_in_bounds[]) {
                fprintf(fp, "%g %g 0.\n", u.x[], w[]);
            }
#endif
#if dimension == 2
            if (cells_in_bounds[]) {
                fprintf(fp, "%g %g %g\n", u.x[], u.y[], w[]);
            }
#endif
        }

    }
    fputs("\t\t\t\t </DataArray>\n", fp);

    fputs("\t\t\t </CellData>\n", fp);


    // Coordinates 
    fputs("\t\t\t <Points>\n", fp);
#if dimension == 1
    fputs("\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n", fp);

    // first do the bottom coordinates stored in zb
    scalar hh[];
    foreach_vertex_xorder(serial, noauto) {
        if (vertex_in_bounds[]) {
            hh[] = (zb[] + zb[-1]) / 2.;
            fprintf(fp, "%12.4f %12.4f 0.\n", x, hh[]);
        }

    }

    foreach_layer() {
        //for (scalar h_layer in hl){
        foreach_vertex_xorder(serial, noauto) {
            if (vertex_in_bounds[]) {
                hh[] += (h[] + h[-1]) / 2.;
                fprintf(fp, "%12.4f %12.4f 0.\n", x, hh[]);
            }
        }
    }
#endif
#if dimension == 2
    fputs("\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n", fp);
    // first do the bottom coordinates stored in zb
    scalar hh[];
    foreach_vertex_xorder(serial, noauto) {
        if (vertex_in_bounds[]) {
            hh[] = (zb[] + zb[-1] + zb[0, -1] + zb[-1, -1]) / 4.;
            fprintf(fp, "%12.4f %12.4f %12.4f\n", x, y, hh[]);
        }
    }

    foreach_layer() {
        //for (scalar h_layer in hl){
        foreach_vertex_xorder(serial, noauto) {
            if (vertex_in_bounds[]) {
                hh[] += (h[] + h[-1] + h[0, -1] + h[-1, -1]) / 4.;
                fprintf(fp, "%12.4f %12.4f %12.4f\n", x, y, hh[]);
            }
        }
    }
#endif

    fputs("\t\t\t\t </DataArray>\n", fp);
    fputs("\t\t\t </Points>\n", fp);
    fputs("\t\t </Piece>\n", fp);
    // write time value
    fprintf(fp, "\t\t <FieldData> \n");
    fprintf(fp, "\t\t\t <DataArray type = \"Float64\" Name = \"%s\" NumberOfTuples = \"1\" format = \"ascii\" RangeMin = \"%.3f\" RangeMax = \"%.3f\"> \n", "TimeValue", t + dt_start, t + dt_start);
    fprintf(fp, "\t\t\t %.3f \n", t + dt_start);
    fprintf(fp, "\t\t\t </DataArray > \n");
    fprintf(fp, "\t\t </FieldData> \n");
    fputs("\t </StructuredGrid>\n", fp);
    fputs("</VTKFile>\n", fp);
    fflush(fp);

}

void output_vts_bin_all_layers_bounded(FILE* fp, double bounds[])
{
        int count = 0;

        // MULTIGRID

        fputs("<?xml version=\"1.0\"?>\n"
            "<VTKFile type=\"StructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);

#if dimension == 1
    // loop over cells and find which are within bounds
    scalar cells_in_bounds[];
    foreach_xorder(serial, noauto) {
        if (x >= bounds[0] && x <= bounds[1])
            cells_in_bounds[] = 1;
        else
            cells_in_bounds[] = 0;
    }

    vertex scalar vertex_in_bounds[];
    int testcount = 0;
    foreach_vertex_xorder(serial, noauto) {
        if (x >= bounds[0] - Delta / 2. && x <= bounds[1] + Delta / 2.) {
            vertex_in_bounds[] = 1;
            testcount++;
        }
        else {
            vertex_in_bounds[] = 0;
        }
    }

    // count number of vertex points in x direction which is within bounds
    int nx = 0;
    double dd = L0 / N;
    for (double xx = X0; xx <= (X0 + L0); xx += dd) {
        if (xx >= (bounds[0] - dd / 2.) && xx <= (bounds[1] + dd / 2.)) {
            nx++;
        }
    }

    fprintf(fp, "\t <StructuredGrid WholeExtent=\"%d %d %d %d 0 0\">\n", 0, nx, 0, nl);
#endif

#if dimension == 2
    // loop over cells and find which are within bounds
    scalar cells_in_bounds[];
    foreach_xorder(serial, noauto) {
        if (x >= bounds[0] && x <= bounds[1] && y >= bounds[2] && y <= bounds[3])
            cells_in_bounds[] = 1;
        else
            cells_in_bounds[] = 0;
    }

    vertex scalar vertex_in_bounds[];
    int testcount = 0;
    foreach_vertex_xorder(serial, noauto) {
        if (x >= bounds[0] - Delta / 2. && x <= bounds[1] + Delta / 2. && y >= bounds[2] - Delta / 2. && y <= bounds[3] + Delta / 2.) {
            vertex_in_bounds[] = 1;
            testcount++;
        }
        else {
            vertex_in_bounds[] = 0;
        }
    }

    // count number of vertex points in x direction which is within bounds
    int nx = 0;
    double dd = L0 / N;
    for (double xx = X0; xx <= (X0 + L0); xx += dd) {
        if (xx >= (bounds[0] - dd / 2.) && xx <= (bounds[1] + dd / 2.)) {
            nx++;
        }
    }
    nx--;
    // count number of vertex points in x direction which is within bounds
    int ny = 0;
    for (double yy = Y0; yy <= (Y0 + L0); yy += dd) {
        if (yy >= (bounds[2] - dd / 2.) && yy <= (bounds[3] + dd / 2.)) {
            ny++;
        }
    }
    ny--;

    fprintf(stdout, "nx: %d, ny: %d, n: %d, %f %f %f %d\n", nx, ny, N, L0, X0, Y0, testcount);

    fprintf(fp, "\t <StructuredGrid WholeExtent=\"%d %d %d %d %d %d\">\n", 0, nx, 0, ny, 0, nl); 
#endif

    fprintf(fp, "\t\t <FieldData> \n");
    fprintf(fp, "\t\t\t <DataArray type = \"Float64\" Name = \"%s\" NumberOfTuples = \"1\" format = \"ascii\" RangeMin = \"%.3f\" RangeMax = \"%.3f\"> \n", "TimeValue", t + dt_start, t + dt_start);
    fprintf(fp, "\t\t\t %.3f \n", t + dt_start);
    fprintf(fp, "\t\t\t </DataArray > \n");
    fprintf(fp, "\t\t </FieldData> \n");

#if dimension == 1
    fprintf(fp, "\t\t <Piece Extent=\"%d %d %d %d  0 0\">\n", 0, nx, 0, nl);
#endif
#if dimension == 2
    fprintf(fp, "\t\t <Piece Extent=\"%d %d %d %d %d %d\">\n", 0, nx, 0, ny, 0, nl);
#endif

    // Loop over velocity data and store kinematics in cell vector stucture
    fputs("\t\t\t <CellData Scalars=\"scalars\">\n", fp);
    fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"appended\" offset=\"%d\">\n", count);
    fputs("\t\t\t\t </DataArray>\n", fp);

    fputs("\t\t\t </CellData>\n", fp);


    // Coordinates 
#if dimension == 1
    count += ((nx * nl * 3) + 1) * 8;
#endif
#if dimension == 2
    count += ((nx * ny * nl * 3) + 1) * 8;
#endif
    fputs("\t\t\t <Points>\n", fp);
    fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\">\n", count);
    fputs("\t\t\t\t </DataArray>\n", fp);
    fputs("\t\t\t </Points>\n", fp);
    fputs("\t\t </Piece>\n", fp);
    // write time value
    fputs("\t </StructuredGrid>\n", fp);
    fputs("\t <AppendedData encoding=\"raw\">\n", fp);
    fputs("_", fp);

    unsigned long long block_len = 8;
    // append kinematics
#if dimension == 1
    block_len = 8 * 3 * nx * nl;
    double dummy = 0.;
#endif
#if dimension == 2
    block_len = 8 * 3 * nx * ny * nl;
#endif
    fwrite(&block_len, sizeof(unsigned long long), 1, fp);
    foreach_layer() {
        foreach_xorder(serial, noauto) {
#if dimension == 1
            if (cells_in_bounds[]) {
                fwrite(&u.x[], sizeof(double), 1, fp);
                fwrite(&w[], sizeof(double), 1, fp);
                fwrite(&dummy, sizeof(double), 1, fp);
            }           
#endif
#if dimension == 2
            if (cells_in_bounds[]) {
                fwrite(&u.x[], sizeof(double), 1, fp);
                fwrite(&u.y[], sizeof(double), 1, fp);
                fwrite(&w[], sizeof(double), 1, fp);
            }          
#endif
        }

    }
    // append points
    // Coordinates 
#if dimension == 1
    block_len = 8 * 3 * (nx + 1) * (nl + 1);
    fwrite(&block_len, sizeof(unsigned long long), 1, fp);
    
    // first do the bottom coordinates stored in zb
    scalar hh[];
    foreach_vertex_xorder(serial, noauto) {
        if (vertex_in_bounds[]) {
            hh[] = (zb[] + zb[-1]) / 2.;
            fwrite(&x, sizeof(double), 1, fp);
            fwrite(&hh[], sizeof(double), 1, fp);
            fwrite(&dummy, sizeof(double), 1, fp);
        }
    }

    foreach_layer() {
        //for (scalar h_layer in hl){
        foreach_vertex_xorder(serial, noauto) {
            hh[] += (h[] + h[-1]) / 2.;
            fwrite(&x, sizeof(double), 1, fp);
            fwrite(&hh[], sizeof(double), 1, fp);
            fwrite(&dummy, sizeof(double), 1, fp);
        }
    }
#endif
#if dimension == 2
    block_len = 8 * 3 * (nx + 1) * (ny + 1) * (nl + 1);
    fwrite(&block_len, sizeof(unsigned long long), 1, fp);

    // first do the bottom coordinates stored in zb
    scalar hh[];
    foreach_vertex_xorder(serial, noauto) {
        if (vertex_in_bounds[]) {
            hh[] = (zb[] + zb[-1] + zb[0, -1] + zb[-1, -1]) / 4.;
            fwrite(&x, sizeof(double), 1, fp);
            fwrite(&y, sizeof(double), 1, fp);
            fwrite(&hh[], sizeof(double), 1, fp);
        }
    }

    foreach_layer() {
        //for (scalar h_layer in hl){
        foreach_vertex_xorder(serial, noauto) {
            if (vertex_in_bounds[]) {
                hh[] += (h[] + h[-1] + h[0, -1] + h[-1, -1]) / 4.;
                fwrite(&x, sizeof(double), 1, fp);
                fwrite(&y, sizeof(double), 1, fp);
                fwrite(&hh[], sizeof(double), 1, fp);
            }
        }
    }
#endif
    fputs("\t\n", fp);
    fputs("\t </AppendedData>\n", fp);
    fputs("</VTKFile>\n", fp);
    fflush(fp);
}


/*
Writes a structured vtk file (.vts) of a single layer.
*/
/*
Writes a structured vtk file (.vts) of a single layer.
*/
/**
Binary (appended raw) VTK writer for a single layer. Drop-in
replacement for output_vts_ascii_single_layer that uses raw fwrite
of doubles instead of fprintf, typically 5-10x faster on large grids.
The output extent and field layout are otherwise identical, so any
ParaView pipeline that reads the ASCII files will read these too. */

void output_vts_bin_single_layer(FILE* fp, int layerid, scalar* list = NULL, vector * vlist = NULL, bool displace_z=false)
{
    int count = 0;
    int layertemp = _layer;
    _layer = layerid;

    fputs("<?xml version=\"1.0\"?>\n<VTKFile type=\"StructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);

#if dimension == 1
    fprintf(fp, "\t <StructuredGrid WholeExtent=\"%d %d %d %d 0 0\">\n", 0, N, 0, 0);
#endif
#if dimension == 2
    fprintf(fp, "\t <StructuredGrid WholeExtent=\"%d %d %d %d %d %d\">\n", 0, N, 0, N, 0, 0);
#endif

    // Time as field data (small, written ascii)
    fprintf(fp, "\t\t <FieldData>\n");
    fprintf(fp, "\t\t\t <DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"ascii\" RangeMin=\"%.3f\" RangeMax=\"%.3f\">\n", t + dt_start, t + dt_start);
    fprintf(fp, "\t\t\t %.3f\n", t + dt_start);
    fprintf(fp, "\t\t\t </DataArray>\n");
    fprintf(fp, "\t\t </FieldData>\n");

#if dimension == 1
    fprintf(fp, "\t\t <Piece Extent=\"%d %d %d %d 0 0\">\n", 0, N, 0, 0);
    int ncells  = N;
    int nverts  = N + 1;
#endif
#if dimension == 2
    fprintf(fp, "\t\t <Piece Extent=\"%d %d %d %d %d %d\">\n", 0, N, 0, N, 0, 0);
    int ncells  = N * N;
    int nverts  = (N + 1) * (N + 1);
#endif

    fputs("\t\t\t <CellData Scalars=\"scalars\">\n", fp);

    // Velocity vector (3 components per cell)
    fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"appended\" offset=\"%d\"/>\n", count);
    count += sizeof(unsigned long long) + ncells * 3 * sizeof(double);

    // Vector fields from vlist
    for (vector v in vlist) {
        fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"%s\" format=\"appended\" offset=\"%d\"/>\n", v.x.name, count);
        count += sizeof(unsigned long long) + ncells * 3 * sizeof(double);
    }

    // Scalar fields from list
    for (scalar s in list) {
        fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"1\" Name=\"%s\" format=\"appended\" offset=\"%d\"/>\n", s.name, count);
        count += sizeof(unsigned long long) + ncells * sizeof(double);
    }
    fputs("\t\t\t </CellData>\n", fp);

    // Points
    fputs("\t\t\t <Points>\n", fp);
    fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"/>\n", count);
    count += sizeof(unsigned long long) + nverts * 3 * sizeof(double);
    fputs("\t\t\t </Points>\n", fp);

    fputs("\t\t </Piece>\n", fp);
    fputs("\t </StructuredGrid>\n", fp);

    // Appended raw data
    fputs("\t <AppendedData encoding=\"raw\">\n", fp);
    fputs("_", fp);

    unsigned long long block_len;
    double dummy = 0.;

    // Velocity
    block_len = ncells * 3 * sizeof(double);
    fwrite(&block_len, sizeof(unsigned long long), 1, fp);
    foreach(serial, noauto) {
#if dimension == 1
        fwrite(&u.x[], sizeof(double), 1, fp);
        fwrite(&w[],   sizeof(double), 1, fp);
        fwrite(&dummy, sizeof(double), 1, fp);
#endif
#if dimension == 2
        fwrite(&u.x[], sizeof(double), 1, fp);
        fwrite(&u.y[], sizeof(double), 1, fp);
        fwrite(&w[],   sizeof(double), 1, fp);
#endif
    }

    // Vector fields
    for (vector v in vlist) {
        block_len = ncells * 3 * sizeof(double);
        fwrite(&block_len, sizeof(unsigned long long), 1, fp);
        foreach(serial, noauto) {
            double vx = val(v.x);
#if dimension == 2
            double vy = val(v.y);
#else
            double vy = 0.;
#endif
            fwrite(&vx,    sizeof(double), 1, fp);
            fwrite(&vy,    sizeof(double), 1, fp);
            fwrite(&dummy, sizeof(double), 1, fp);
        }
    }

    // Scalar fields. For "eta", emit nan in dry columns to match the
    // ascii writer's behaviour.
    for (scalar s in list) {
        block_len = ncells * sizeof(double);
        fwrite(&block_len, sizeof(unsigned long long), 1, fp);
        bool is_eta = (strcmp(s.name, "eta") == 0);
        if (is_eta) {
            double nanval = nan("");
            foreach(serial, noauto) {
                double Htot = 0;
                foreach_layer() Htot += h[];
                double v = (Htot > dry) ? val(s) : nanval;
                fwrite(&v, sizeof(double), 1, fp);
            }
        }
        else {
            foreach(serial, noauto) {
                double v = val(s);
                fwrite(&v, sizeof(double), 1, fp);
            }
        }
    }

    // Points (vertices)
    block_len = nverts * 3 * sizeof(double);
    fwrite(&block_len, sizeof(unsigned long long), 1, fp);
#if dimension == 1
    if (displace_z) {
        scalar z_trans = list[0];
        foreach_vertex(serial, noauto) {
            double xv = x;
            double zv = (z_trans[] + z_trans[-1]) / 2.;
            fwrite(&xv,    sizeof(double), 1, fp);
            fwrite(&zv,    sizeof(double), 1, fp);
            fwrite(&dummy, sizeof(double), 1, fp);
        }
    }
    else {
        foreach_vertex(serial, noauto) {
            double xv = x;
            fwrite(&xv,    sizeof(double), 1, fp);
            fwrite(&dummy, sizeof(double), 1, fp);
            fwrite(&dummy, sizeof(double), 1, fp);
        }
    }
#endif
#if dimension == 2
    if (displace_z) {
        scalar z_trans = list[0];
        foreach_vertex(serial, noauto) {
            double xv = x, yv = y;
            double zv = (z_trans[] + z_trans[-1] + z_trans[0,-1] + z_trans[-1,-1]) / 4.;
            fwrite(&xv, sizeof(double), 1, fp);
            fwrite(&yv, sizeof(double), 1, fp);
            fwrite(&zv, sizeof(double), 1, fp);
        }
    }
    else {
        foreach_vertex(serial, noauto) {
            double xv = x, yv = y;
            fwrite(&xv,    sizeof(double), 1, fp);
            fwrite(&yv,    sizeof(double), 1, fp);
            fwrite(&dummy, sizeof(double), 1, fp);
        }
    }
#endif

    fputs("\n\t </AppendedData>\n", fp);
    fputs("</VTKFile>\n", fp);
    fflush(fp);

    _layer = layertemp;
}


void output_vts_ascii_single_layer(FILE* fp, int layerid, scalar* list = NULL, vector * vlist = NULL, bool displace_z=false)
{

        // MULTIGRID

    fputs("<?xml version=\"1.0\"?>\n<VTKFile type=\"StructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);

#if dimension == 1
    fprintf(fp, "\t <StructuredGrid WholeExtent=\"%d %d %d %d 0 0\">\n", 0, N, 0, 0);
    fprintf(fp, "\t\t <Piece Extent=\"%d %d %d %d  0 0\">\n", 0, N, 0, nl);
#endif

#if dimension == 2
    fprintf(fp, "\t <StructuredGrid WholeExtent=\"%d %d %d %d %d %d\">\n", 0, N, 0, N, 0, 0);
    fprintf(fp, "\t\t <Piece Extent=\"%d %d %d %d %d %d\">\n", 0, N, 0, N, 0, nl);
#endif

    // Loop over velocity data and store kinematics in cell vector stucture
    int layertemp = _layer;
    _layer = layerid;
    fputs("\t\t\t <CellData Scalars=\"scalars\">\n", fp);
    fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"ascii\">\n");
    foreach(serial, noauto) {
#if dimension == 1
        fprintf(fp, "%g %g 0.\n", u.x[], w[]);
#endif
#if dimension == 2
        fprintf(fp, "%g %g %g\n", u.x[], u.y[], w[]);
#endif
    }
    fputs("\t\t\t\t </DataArray>\n", fp);

    // loop over all scalars in scalarlist and store values as cell data
    for (vector v in vlist) {
        fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"%s\" format=\"ascii\">\n", v.x.name);
        foreach(serial, noauto) {
#if dimension == 1
        fprintf(fp, "%g 0. 0.\n", val(v.x));
#endif
#if dimension == 2
        fprintf(fp, "%g %g 0.\n",  val(v.x),  val(v.y));
#endif
        }
        fputs("\t\t\t\t </DataArray>\n", fp);
    }

    // loop over all scalars in scalarlist and store values as cell data
    for (scalar s in list) {
        if (strcmp(s.name, "eta") == 0) {
            fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", s.name);
            foreach(serial, noauto) {
                double Htot = 0;
                foreach_layer()
                    Htot += h[];

                if (Htot > dry) {
                    fprintf(fp, "%g\n", val(s));
                }
                else {
                    fprintf(fp, "nan\n");
                }
            }
            fputs("\t\t\t\t </DataArray>\n", fp);
        }
        else {
            fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", s.name);
            foreach(serial, noauto) {
                fprintf(fp, "%g\n", val(s));
            }
            fputs("\t\t\t\t </DataArray>\n", fp);
        }

    }



    fputs("\t\t\t </CellData>\n", fp);


    // Coordinates 
    fputs("\t\t\t <Points>\n", fp);
#if dimension == 1
    fputs("\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n", fp);

    if (displace_z) {
        scalar z_trans = list[0];
        // set z coordinate to the value of the first scalar field in the provided list.
        double hh;
        foreach_vertex(serial, noauto) {
            hh = (z_trans[] + z_trans[-1]) / 2.;
            fprintf(fp, "%12.4f %12.4f 0.\n", x, hh);

        }
    }
    else {
        // set z coordinate to 0.
        foreach_vertex(serial, noauto) {
            fprintf(fp, "%12.4f 0. 0.\n", x);
        }
    }
#endif
#if dimension == 2
    fputs("\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n", fp);

    if (displace_z) {
        scalar z_trans = list[0];
        // first do the bottom coordinates stored in zb
        double hh;
        foreach_vertex(serial, noauto) {
            hh = (z_trans[] + z_trans[-1] + z_trans[0, -1] + z_trans[-1, -1]) / 4.;
            fprintf(fp, "%12.4f %12.4f %12.4f\n", x, y, hh);
        }
    }
    else {
        // first do the bottom coordinates stored in zb
        foreach_vertex(serial, noauto) {
            fprintf(fp, "%12.4f %12.4f 0.0\n", x, y);
        }
    }
#endif

    fputs("\t\t\t\t </DataArray>\n", fp);
    fputs("\t\t\t </Points>\n", fp);
    fputs("\t\t </Piece>\n", fp);

    // write time value
    fprintf(fp, "\t\t <FieldData> \n");
    fprintf(fp, "\t\t\t <DataArray type = \"Float64\" Name = \"%s\" NumberOfTuples = \"1\" format = \"ascii\" RangeMin = \"%.3f\" RangeMax = \"%.3f\"> \n", "TimeValue", t + dt_start, t + dt_start);
    fprintf(fp, "\t\t\t %.3f \n", t + dt_start);
    fprintf(fp, "\t\t\t </DataArray > \n");
    fprintf(fp, "\t\t </FieldData> \n");
    fputs("\t </StructuredGrid>\n", fp);
    fputs("</VTKFile>\n", fp);
    fflush(fp);

    _layer = layertemp;
    //#endif
}

/**
Binary (appended raw) VTK structured grid writer for all layers.
Drop-in faster replacement for output_vts_ascii_all_layers — same
field semantics but uses raw fwrite of doubles instead of fprintf,
typically 5–10× faster for moderate-to-large grids.

Parameters:
  list                 — extra cell-centered scalar fields (NULL = none)
  vlist                — extra cell-centered vector fields (NULL = none, currently unused)
  flatten_layers       — render the column as a unit-spaced stack (debug)
  include_boundaries   — also emit BGHOSTS-wide boundary ghost cells

The file must be opened in binary mode ("wb"). */
void output_vts_bin_all_layers(FILE* fp, scalar * list = NULL, vector * vlist = NULL, bool flatten_layers = false, bool include_boundaries = false)
{
    int Nx = include_boundaries ? N + 2*BGHOSTS : N;

#if dimension == 1
    int ncells = Nx * nl;
    int nverts = (Nx + 1) * (nl + 1);
#endif
#if dimension == 2
    int ncells = Nx * Nx * nl;
    int nverts = (Nx + 1) * (Nx + 1) * (nl + 1);
#endif

    fputs("<?xml version=\"1.0\"?>\n"
          "<VTKFile type=\"StructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);

#if dimension == 1
    fprintf(fp, "\t <StructuredGrid WholeExtent=\"%d %d %d %d 0 0\">\n", 0, Nx, 0, nl);
#endif
#if dimension == 2
    fprintf(fp, "\t <StructuredGrid WholeExtent=\"%d %d %d %d %d %d\">\n", 0, Nx, 0, Nx, 0, nl);
#endif

    // Time field data (small, ascii)
    fprintf(fp, "\t\t <FieldData>\n");
    fprintf(fp, "\t\t\t <DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"ascii\" RangeMin=\"%.3f\" RangeMax=\"%.3f\">\n", t + dt_start, t + dt_start);
    fprintf(fp, "\t\t\t %.3f\n", t + dt_start);
    fprintf(fp, "\t\t\t </DataArray>\n");
    fprintf(fp, "\t\t </FieldData>\n");

#if dimension == 1
    fprintf(fp, "\t\t <Piece Extent=\"%d %d %d %d 0 0\">\n", 0, Nx, 0, nl);
#endif
#if dimension == 2
    fprintf(fp, "\t\t <Piece Extent=\"%d %d %d %d %d %d\">\n", 0, Nx, 0, Nx, 0, nl);
#endif

    // Compute appended-data offsets. Each block is preceded by an
    // unsigned long long header containing the byte length.
    int count = 0;
    fputs("\t\t\t <CellData Scalars=\"scalars\">\n", fp);

    // Scalar fields from list
    for (scalar s in list) {
        fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"1\" Name=\"%s\" format=\"appended\" offset=\"%d\"/>\n", s.name, count);
        count += sizeof(unsigned long long) + ncells * sizeof(double);
    }
    // Velocity (3 components per cell)
    fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"appended\" offset=\"%d\"/>\n", count);
    count += sizeof(unsigned long long) + ncells * 3 * sizeof(double);
    fputs("\t\t\t </CellData>\n", fp);

    fputs("\t\t\t <Points>\n", fp);
    fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"/>\n", count);
    count += sizeof(unsigned long long) + nverts * 3 * sizeof(double);
    fputs("\t\t\t </Points>\n", fp);

    fputs("\t\t </Piece>\n", fp);
    fputs("\t </StructuredGrid>\n", fp);

    // Appended raw data section
    fputs("\t <AppendedData encoding=\"raw\">\n", fp);
    fputs("_", fp);

    unsigned long long block_len;
#if dimension == 1
    double dummy = 0.;
#endif

    // Cell data: scalar list, then velocity
    for (scalar s in list) {
        block_len = ncells * sizeof(double);
        fwrite(&block_len, sizeof(unsigned long long), 1, fp);
        if (include_boundaries) {
            foreach_layer() {
                foreach_xorder_boundary(serial, noauto) {
                    double v = val(s);
                    fwrite(&v, sizeof(double), 1, fp);
                }
            }
        }
        else {
            foreach_layer() {
                foreach_xorder(serial) {
                    double v = val(s);
                    fwrite(&v, sizeof(double), 1, fp);
                }
            }
        }
    }

    // Velocity
    block_len = ncells * 3 * sizeof(double);
    fwrite(&block_len, sizeof(unsigned long long), 1, fp);
    if (include_boundaries) {
        foreach_layer() {
            foreach_xorder_boundary(serial, noauto) {
#if dimension == 1
                fwrite(&u.x[], sizeof(double), 1, fp);
                fwrite(&w[],   sizeof(double), 1, fp);
                fwrite(&dummy, sizeof(double), 1, fp);
#endif
#if dimension == 2
                fwrite(&u.x[], sizeof(double), 1, fp);
                fwrite(&u.y[], sizeof(double), 1, fp);
                fwrite(&w[],   sizeof(double), 1, fp);
#endif
            }
        }
    }
    else {
        foreach_layer() {
            foreach_xorder(serial) {
#if dimension == 1
                fwrite(&u.x[], sizeof(double), 1, fp);
                fwrite(&w[],   sizeof(double), 1, fp);
                fwrite(&dummy, sizeof(double), 1, fp);
#endif
#if dimension == 2
                fwrite(&u.x[], sizeof(double), 1, fp);
                fwrite(&u.y[], sizeof(double), 1, fp);
                fwrite(&w[],   sizeof(double), 1, fp);
#endif
            }
        }
    }

    // Vertex coordinates
    block_len = nverts * 3 * sizeof(double);
    fwrite(&block_len, sizeof(unsigned long long), 1, fp);

    vertex scalar hh[];

    if (include_boundaries) {
#if dimension == 1
        foreach_vertex_xorder_boundary(serial, noauto) {
            if (point.i == (GHOSTS - BGHOSTS))
                hh[] = zb[];
            else if (point.i == (point.n.x + GHOSTS + BGHOSTS))
                hh[] = zb[-1];
            else
                hh[] = (zb[] + zb[-1]) / 2.;
            double xv = x, zv = hh[];
            fwrite(&xv,    sizeof(double), 1, fp);
            fwrite(&zv,    sizeof(double), 1, fp);
            fwrite(&dummy, sizeof(double), 1, fp);
        }
        foreach_layer() {
            foreach_vertex_xorder_boundary(serial, noauto) {
                if (point.i == (GHOSTS - BGHOSTS))
                    hh[] += h[];
                else if (point.i == (point.n.x + GHOSTS + BGHOSTS))
                    hh[] += h[-1];
                else
                    hh[] += (h[] + h[-1]) / 2.;
                double xv = x, zv = hh[];
                fwrite(&xv,    sizeof(double), 1, fp);
                fwrite(&zv,    sizeof(double), 1, fp);
                fwrite(&dummy, sizeof(double), 1, fp);
            }
        }
#endif
#if dimension == 2
        foreach_vertex_xorder_boundary(serial, noauto) {
            hh[] = (zb[] + zb[-1] + zb[0,-1] + zb[-1,-1]) / 4.;
            double xv = x, yv = y, zv = hh[];
            fwrite(&xv, sizeof(double), 1, fp);
            fwrite(&yv, sizeof(double), 1, fp);
            fwrite(&zv, sizeof(double), 1, fp);
        }
        foreach_layer() {
            foreach_vertex_xorder_boundary(serial, noauto) {
                hh[] += (h[] + h[-1] + h[0,-1] + h[-1,-1]) / 4.;
                double xv = x, yv = y, zv = hh[];
                fwrite(&xv, sizeof(double), 1, fp);
                fwrite(&yv, sizeof(double), 1, fp);
                fwrite(&zv, sizeof(double), 1, fp);
            }
        }
#endif
    }
    else {
#if dimension == 1
        foreach_vertex_xorder(serial) {
            hh[] = (zb[] + zb[-1]) / 2.;
            double xv = x, zv = hh[];
            fwrite(&xv,    sizeof(double), 1, fp);
            fwrite(&zv,    sizeof(double), 1, fp);
            fwrite(&dummy, sizeof(double), 1, fp);
        }
        foreach_layer() {
            foreach_vertex_xorder(serial) {
                hh[] += (h[] + h[-1]) / 2.;
                double xv = x, zv = hh[];
                fwrite(&xv,    sizeof(double), 1, fp);
                fwrite(&zv,    sizeof(double), 1, fp);
                fwrite(&dummy, sizeof(double), 1, fp);
            }
        }
#endif
#if dimension == 2
        foreach_vertex_xorder(serial) {
            if (flatten_layers)
                hh[] = -0.5;
            else
                hh[] = (zb[] + zb[-1] + zb[0,-1] + zb[-1,-1]) / 4.;
            double xv = x, yv = y, zv = hh[];
            fwrite(&xv, sizeof(double), 1, fp);
            fwrite(&yv, sizeof(double), 1, fp);
            fwrite(&zv, sizeof(double), 1, fp);
        }
        foreach_layer() {
            foreach_vertex_xorder(serial) {
                if (flatten_layers)
                    hh[] = 0.5 + 1.*_layer;
                else
                    hh[] += (h[] + h[-1] + h[0,-1] + h[-1,-1]) / 4.;
                double xv = x, yv = y, zv = hh[];
                fwrite(&xv, sizeof(double), 1, fp);
                fwrite(&yv, sizeof(double), 1, fp);
                fwrite(&zv, sizeof(double), 1, fp);
            }
        }
#endif
    }

    fputs("\n\t </AppendedData>\n", fp);
    fputs("</VTKFile>\n", fp);
    fflush(fp);
}


# if dimension == 2

void output_vts_ascii_sliceXZ(scalar * list, FILE* fp, double y_target)
{

    // #if defined(_OPENMP)
    //     //int num_omp = omp_get_max_threads();
    // omp_set_num_threads(1);
    // #endif

    // we compute which cell index point.j that contains our y coordinate
    int jp = floor(N*(y_target-Y0)/L0); NOT_UNUSED(jp);

    int count = 0;
    // MULTIGRID

    fputs("<?xml version=\"1.0\"?>\n"
        "<VTKFile type=\"StructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);

    //fprintf(fp, "\t <StructuredGrid WholeExtent=\"%d %d %d %d 0 0\">\n", 0, N, 0, nl);
    fprintf(fp, "\t <StructuredGrid WholeExtent=\"%d %d %d %d %d %d\">\n", 0, N, 0, 0, 0, nl);


    // write time value
    fprintf(fp, "\t\t <FieldData> \n");
    fprintf(fp, "\t\t\t <DataArray type = \"Float64\" Name = \"%s\" NumberOfTuples = \"1\" format = \"ascii\" RangeMin = \"%.3f\" RangeMax = \"%.3f\"> \n", "TimeValue", t + dt_start, t + dt_start);
    fprintf(fp, "\t\t\t %.3f \n", t + dt_start);
    fprintf(fp, "\t\t\t </DataArray > \n");
    fprintf(fp, "\t\t </FieldData> \n");
    

    fprintf(fp, "\t\t <Piece Extent=\"%d %d %d %d %d %d\">\n", 0, N, 0, 0, 0, nl);


    // Cell data layout (offsets into the appended raw block).
    // Each block is preceded by an unsigned long long byte-length header.
    // Slice has N cells in x, 1 cell in y, nl in z, so ncells = N*nl
    // and the corresponding vertex grid has (N+1)*1*(nl+1) = (N+1)*(nl+1) verts.
    int ncells = N * nl;
    int nverts = (N + 1) * (nl + 1);

    fputs("\t\t\t <CellData Scalars=\"scalars\">\n", fp);
    for (scalar s in list) {
        fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"1\" Name=\"%s\" format=\"appended\" offset=\"%d\"/>\n", s.name, count);
        count += sizeof(unsigned long long) + ncells * sizeof(double);
    }
    // Velocity (3 components per cell)
    fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"appended\" offset=\"%d\"/>\n", count);
    count += sizeof(unsigned long long) + ncells * 3 * sizeof(double);
    fputs("\t\t\t </CellData>\n", fp);

    // Points (vertices)
    fputs("\t\t\t <Points>\n", fp);
    fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"/>\n", count);
    count += sizeof(unsigned long long) + nverts * 3 * sizeof(double);
    fputs("\t\t\t </Points>\n", fp);
    fputs("\t\t </Piece>\n", fp);

    fputs("\t </StructuredGrid>\n", fp);
    fputs("\t <AppendedData encoding=\"raw\">\n", fp);
    fputs("_", fp);
    
    unsigned long long block_len = 8;
    
    // Append scalars
    for (scalar s in list) {
    block_len = 8 * N * nl;

        fwrite(&block_len, sizeof(unsigned long long), 1, fp);
        foreach_layer() {
            foreach_xonly(jp) {
                fwrite(&val(s), sizeof(double), 1, fp);
            }
        }
    }
    // append velocity
    block_len = 8 * 3 * N * nl;
    fwrite(&block_len, sizeof(unsigned long long), 1, fp);
    foreach_layer() {
        foreach_xonly(jp) {
            fwrite(&u.x[], sizeof(double), 1, fp);
            fwrite(&u.y[], sizeof(double), 1, fp);
            fwrite(&w[], sizeof(double), 1, fp);
        }

    }
    // append points
    // Coordinates 
    block_len = 8 * 3 * (N + 1) * (nl + 1);
    fwrite(&block_len, sizeof(unsigned long long), 1, fp);

    // first do the bottom coordinates stored in zb
    scalar hh[];
    foreach_vertex_xonly(jp) {
        hh[] = (zb[] + zb[-1] + zb[0, -1] + zb[-1, -1]) / 4.;
        fwrite(&x, sizeof(double), 1, fp);
        fwrite(&y_target, sizeof(double), 1, fp);
        fwrite(&hh[], sizeof(double), 1, fp);
    }

    foreach_layer() {
        //for (scalar h_layer in hl){
        foreach_vertex_xonly(jp) {
            hh[] += (h[] + h[-1] + h[0, -1] + h[-1, -1]) / 4.;
            fwrite(&x, sizeof(double), 1, fp);
            fwrite(&y_target, sizeof(double), 1, fp);
            fwrite(&hh[], sizeof(double), 1, fp);
        }
    }


    fputs("\t\n", fp);
    fputs("\t </AppendedData>\n", fp);
    fputs("</VTKFile>\n", fp);
    fflush(fp);

    // #if defined(_OPENMP)
    //     omp_set_num_threads(num_omp);
    // #endif
}

#endif



void output_vts_ascii_all_layers_incl_boundary_more(FILE* fp)
{


    // MULTIGRID

    fputs("<?xml version=\"1.0\"?>\n"
        "<VTKFile type=\"StructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);

#if dimension == 1
    fprintf(fp, "\t <StructuredGrid WholeExtent=\"%d %d %d %d 0 0\">\n", 0, N + 2*GHOSTS, 0, nl);
    fprintf(fp, "\t\t <Piece Extent=\"%d %d %d %d  0 0\">\n", 0, N + 2*GHOSTS, 0, nl);
#endif

#if dimension == 2
    fprintf(fp, "\t <StructuredGrid WholeExtent=\"%d %d %d %d %d %d\">\n", 0, N + 2*GHOSTS, 0, N + 2*GHOSTS, 0, nl);
    fprintf(fp, "\t\t <Piece Extent=\"%d %d %d %d %d %d\">\n", 0, N + 2*GHOSTS, 0, N + 2*GHOSTS, 0, nl);
#endif

    // Loop over velocity data and store kinematics in cell vector stucture
    fputs("\t\t\t <CellData Scalars=\"scalars\">\n", fp);
    fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"ascii\">\n");

    //for (u_layer, w_layer in ul, wl){  
    foreach_layer() {
        foreach_xorder_boundary(serial, noauto) {

#if dimension == 1
            fprintf(fp, "%g %g 0.\n", u.x[], w[]);
#endif
#if dimension == 2
            fprintf(fp, "%g %g %g\n", u.x[], u.y[], w[]);
#endif
        }

    }
    fputs("\t\t\t\t </DataArray>\n", fp);



    /*    
    // Output hydrodynamic pressure
    fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"1\" Name=\"%s\" format=\"ascii\">\n", "hpres");
    foreach_layer() {
        foreach_xorder_boundary(serial, noauto) {
            fprintf(fp, "%g\n", hpres[]);
        }

    }
    fputs("\t\t\t\t </DataArray>\n", fp);
    */
    fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"1\" Name=\"%s\" format=\"ascii\">\n", "phi");
    foreach_layer() {
        foreach_xorder_boundary(serial, noauto) {
            fprintf(fp, "%g\n", phi[]);
        }

    }
    fputs("\t\t\t\t </DataArray>\n", fp);

    fputs("\t\t\t </CellData>\n", fp);



    // Coordinates 
    fputs("\t\t\t <Points>\n", fp);
#if dimension == 1
    fputs("\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n", fp);

    // first do the bottom coordinates stored in zb
    // we cannot us the default vertex or face arrays here, because we get trouble with the boundaries.\OLAND
    double hh[N+2*GHOSTS+1];
    //array = (double*) malloc( * sizeof(int));
    foreach_vertex_xorder_boundary(serial, noauto) {

        
        if (point.i == 0){
            //fprintf(stdout,"xiter1: %d, %d, %g\n", point.i, point.n, zb[]);
            hh[point.i] = zb[];
        }
        else if(point.i == point.n.x+2*GHOSTS){
            //fprintf(stdout,"xiter2: %d, %d, %g\n", point.i, point.n, zb[-1]);
            hh[point.i] = zb[-1];
            
        }
        else{
            //fprintf(stdout,"xiter3: %d, %d, %g %g\n", point.i, point.n, zb[], zb[-1]);
            hh[point.i] = (zb[] + zb[-1]) / 2.;
        }
        //fprintf(stdout,"xiter4: %d, %d\n", point.i, point.n);
        fprintf(fp, "%12.4f %12.4f 0.\n", x, hh[point.i]);

    }
    //fprintf(stdout, "ape\n");

    foreach_layer() {
        //for (scalar h_layer in hl){
        foreach_vertex_xorder_boundary(serial, noauto) {
            if (point.i == 0){
                hh[point.i] += h[];
            }
            else if(point.i == point.n.x+2*GHOSTS){
                hh[point.i] += h[-1];
            }
            else{
                hh[point.i] += (h[] + h[-1]) / 2.;
            }

            fprintf(fp, "%12.4f %12.4f 0.\n", x, hh[point.i]);
        }
    }
#endif
#if dimension == 2
// Todo: fix to reflect changes that is done for 1D. current 2d implementation will occationally result in seg faults
    fputs("\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n", fp);
    // first do the bottom coordinates stored in zb
    int NN2 = (N+2*GHOSTS+1);
    double hh[NN2*NN2];
    foreach_vertex_xorder_boundary(serial, noauto) {
        if (point.i == 0){
            //fprintf(stdout,"xiter1: %d, %d, %g\n", point.i, point.n, zb[]);
            hh[point.j*NN2+point.i] = zb[];
        }
        else if(point.i == point.n.x+2*GHOSTS){
            //fprintf(stdout,"xiter2: %d, %d, %g\n", point.i, point.n, zb[-1]);
            hh[point.j*NN2+point.i] = zb[-1];
            
        }
        if (point.j == 0){
            //fprintf(stdout,"xiter1: %d, %d, %g\n", point.i, point.n, zb[]);
            hh[point.j*NN2+point.i] = zb[];
        }
        else if(point.j == point.n.y+2*GHOSTS){
            //fprintf(stdout,"xiter2: %d, %d, %g\n", point.i, point.n, zb[-1]);
            hh[point.j*NN2+point.i] = zb[0,-1];
            
        }
        else{
            //fprintf(stdout,"xiter3: %d, %d, %g %g\n", point.i, point.n, zb[], zb[-1]);
            hh[point.j*NN2+point.i] = (zb[] + zb[-1] + zb[0, -1] + zb[-1, -1]) / 4.;
        }
        //fprintf(stdout,"xiter4: %d, %d\n", point.i, point.n);
        //fprintf(fp, "%12.4f %12.4f 0.\n", x, hh[point.i]);
        fprintf(fp, "%12.4f %12.4f %12.4f\n", x, y, hh[point.j*NN2+point.i]);
    }

    foreach_layer() {
        //for (scalar h_layer in hl){
        foreach_vertex_xorder_boundary(serial, noauto) {
            if (point.i == 0){
                //fprintf(stdout,"xiter1: %d, %d, %g\n", point.i, point.n, zb[]);
                hh[point.j*NN2+point.i] = h[];
            }
            else if(point.i == point.n.x+2*GHOSTS){
                //fprintf(stdout,"xiter2: %d, %d, %g\n", point.i, point.n, zb[-1]);
                hh[point.j*NN2+point.i] = h[-1];
                
            }
            if (point.j == 0){
                //fprintf(stdout,"xiter1: %d, %d, %g\n", point.i, point.n, zb[]);
                hh[point.j*NN2+point.i] = h[];
            }
            else if(point.j == point.n.y+2*GHOSTS){
                //fprintf(stdout,"xiter2: %d, %d, %g\n", point.i, point.n, zb[-1]);
                hh[point.j*NN2+point.i] = h[0,-1];
                
            }
            else{
                //fprintf(stdout,"xiter3: %d, %d, %g %g\n", point.i, point.n, zb[], zb[-1]);
                hh[point.j*NN2+point.i] = (h[] + h[-1] + h[0, -1] + h[-1, -1]) / 4.;
                fprintf(fp, "%12.4f %12.4f %12.4f\n", x, y, hh[point.j*NN2+point.i]);
            }
        }
    }
#endif

    fputs("\t\t\t\t </DataArray>\n", fp);
    fputs("\t\t\t </Points>\n", fp);
    fputs("\t\t </Piece>\n", fp);
    // write time value
    fprintf(fp, "\t\t <FieldData> \n");
    fprintf(fp, "\t\t\t <DataArray type = \"Float64\" Name = \"%s\" NumberOfTuples = \"1\" format = \"ascii\" RangeMin = \"%.3f\" RangeMax = \"%.3f\"> \n", "TimeValue", t + dt_start, t + dt_start);
    fprintf(fp, "\t\t\t %.3f \n", t + dt_start);
    fprintf(fp, "\t\t\t </DataArray > \n");
    fprintf(fp, "\t\t </FieldData> \n");
    fputs("\t </StructuredGrid>\n", fp);
    fputs("</VTKFile>\n", fp);
    fflush(fp);
}
#endif
