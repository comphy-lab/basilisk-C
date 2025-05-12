/**

# XML File Formats

## preamble: define some useful macros

To make the routines (a little bit) easier to follow, we'll define the following 
macros to deal with the reconstruction of the VOF interface

*/


#include "geometry.h"
#include "fractions.h"

#if dimension == 1
coord mycs(Point point, scalar c)
{
  coord n = {1.};
  return n;
}
#elif dimension == 2
#include "myc2d.h"
#define mfacets int m = facets(n, alpha, v);
#else // dimension == 3
#include "myc.h"
#define mfacets int m = facets(n, alpha, v, 1.);
#endif

/** Macro to simplify facets calculation */ 
#define shortcut_facets               \
  coord n;                            \
  n = mycs(point, c);                 \
  double alpha = plane_alpha(c[], n); \
  coord v[12];                        \
  mfacets;

/** Macro to simplify dealing with slices */ 
#define shortcut_slice(n, _alpha)                                \
  double alpha = (_alpha - n.x * x - n.y * y - n.z * z) / Delta; \
  if (fabs(alpha) > 0.87)                                        \
    continue;

/** and some useful functions */

#include "output_vtu_helpers.h"

/**
## output_vtu(): Exports full 2D (or 3D) fields.

This function writes one [VTK XML
file](https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html) which
can be read using Paraview. The file
stores scalar and vector fields defined at the center points are stored 
at the cell center of an unstructured grid with the following structure: 

<center>
<table>
<tr>
<td>![](vtk_quad.png){ width="80%" }</td>
<td>![](vtk_hexahedron.png){ width="100%" }</td>
</tr>
<tr>
<td><center>if 2D</center></td> 
<td><center>if 3D</center></td> 
</tr>
</table>
</center>

In MPI environments, each task writes its own `.vtu` file, linked together by a
`.pvtu` file. Results are recorded in binary format. 

The arguments and their default values are:

*list*
: pointer to the list of scalar fields to be exported.

*vlist*
: pointer to the list of vector fields to be exported.

*subname*
: subname to be used for the output file.

### Example Usage

```c
scalar * slist = {a,b};
vector * vlist = {c,d};
char *subname = "domain";
output_vtu(slist, vlist, subname);
```

see, also [example](test_output0_vtu.c).


*/

// Function prototypes
void output_vtu_pid(scalar *list, vector *vlist, char *subname);
#ifdef _MPI
void output_pvtu(scalar *list, vector *vlist, char *subname);
#endif

trace    
void output_vtu(scalar *list, vector *vlist, char *subname){
// Check if MPI is defined
@if _MPI
  // If MPI is defined, call output_pvtu to handle parallel VTU output
  output_pvtu(list, vlist, subname);
@else
  // If MPI is not defined, call output_vtu_pid to handle serial VTU output
  output_vtu_pid(list, vlist, subname);
@endif
}

/** ### *output_vtu_pid()*: writes one `.vtu` file for the current process */
void output_vtu_pid(scalar *list, vector *vlist, char *subname) {

#if defined(_OPENMP)
  int num_omp = omp_get_max_threads(); // Get the number of OpenMP threads
  omp_set_num_threads(1);              // Set the number of OpenMP threads to 1
#endif

  char name[111];                   // Buffer for file name construction
  sprintf(name, "%s.vtu", subname); // Construct the VTU filename

  FILE *fp = fopen(name, "w"); // Open the VTU file for writing
  if (!fp){
    fprintf(stderr, "Error opening file %s for writing.\n", name);
  }

  /** Define a scalar field for periodic conditions */ 
  scalar per_mask[];
  foreach () {
    per_mask[] = 1.;
  }

  /** Obtain the number of points, cells, and marker used for connectivity */ 
  vertex scalar marker[];
  long no_points = 0, no_cells = 0;
  foreach_vertex(serial, noauto){
    #if TREE
      marker[] = _k;
    #else
      #if dimension == 2
        marker[] = (point.i - 2) * ((1 << point.level) + 1) + (point.j - 2);
      #else
        marker[] = (point.i - 2) * sq((1 << point.level) + 1) + (point.j - 2) * ((1 << point.level) + 1) + (point.k - 2);
      #endif
    #endif
    no_points++; // Increment the number of points
  }

  foreach (serial, noauto){
    if (per_mask[]){
      no_cells++; // Increment the number of cells
    }
  }

  /** VTK cell types: VTK_QUAD (in 2D) or VTK_HEXAHEDRON (in 3D)  */
  #if dimension == 2
    char type = 9, noffset = 4;
  #elif dimension == 3
    char type = 12, noffset = 8;
  #endif

  /** Write the light data of the VTU file with data blocks specified by offset */ 
  long count = 0;
  write_vtu_header(fp, no_points, no_cells);                        // Write the VTU file header
  write_scalar_light_data(fp, list, vlist, &count, no_cells);       // Write scalar data arrays
  write_points_light_data(fp, &count, no_points);                   // Write points data array
  write_cells_light_data(fp, &count, no_cells, no_cells * noffset); // Write cells data arrays
  write_vtu_appended(fp);                                           // Write the VTU appended data section

  /** Write the heavy data blocks */ 
  write_scalar_heavy_data(fp, list, per_mask, no_cells);            // Write scalar field data
  write_vector_heavy_data(fp, vlist, per_mask, no_cells);           // Write vector field data
  write_points_heavy_data(fp, no_points);                           // Write points data
  write_cell_offsets(fp, no_cells, noffset);                        // Write cell offsets
  write_cell_types(fp, no_cells, type);                             // Write cell types
  write_cell_connectivity(fp, marker, per_mask, no_cells, noffset); // Write cell connectivity

  /** and close the file */ 
  fputs("\t\n", fp);
  fputs("\t</AppendedData>\n", fp);
  fputs("</VTKFile>\n", fp);
  fflush(fp); // Flush the file buffer
  fclose(fp); // Close the VTU file

#if defined(_OPENMP)
  omp_set_num_threads(num_omp); // Restore the original number of OpenMP threads
#endif
}

/** ### *output_pvtu()*: if MPI, writes one `.pvtu` and `.vtu` for each process */
@ if _MPI 
void output_pvtu(scalar *list, vector *vlist, char *subname)
{
  char name[112]; // Buffer for file name construction
  FILE *fp;       // File pointer for file operations

  if (pid() == 0){

    sprintf(name, "%s.pvtu", subname); // Construct the PVTU filename
    fp = fopen(name, "w");             // Open the PVTU file for writing
    if (!fp){
      fprintf(stderr, "Error opening file %s for writing.\n", name);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Write sections of the PVTU file
    write_pvtu_header(fp);                     // Write the header
    write_scalar_light_pdata(fp, list, vlist); // Write scalar data arrays
    write_points_light_pdata(fp);              // Write points data array
    write_pieces_light_pdata(fp, subname);     // Write piece references

    // Close the PVTU file
    fflush(fp); // Flush the file buffer
    fclose(fp); // Close the file
  }
  MPI_Barrier(MPI_COMM_WORLD);
  sprintf(name, "%s_n%3.3d", subname, pid());
  output_vtu_pid(list, vlist, name);
}
@endif





/**
## *output_slice_vtu()*: Exports a 2D slice of a 3D field

This function takes a slice defined by `n={n_x, n_y, n_z}` and  `_alpha` as as
in [view.h](view.h). Only works for `x`, `y`, and/or `z` planes. Here, `_alpha`
is assumed to intersect a cell face and **must** be a multiple of `L0/1<<MINLEVEL`.
This also means that scalar and vector fields are written at the corresponding
face values. 

Naturaly, this routine only works in 3D.

The arguments and their default values are:

*list*
: pointer to the list of scalar fields to be exported.

*vlist*
: pointer to the list of vector fields to be exported.

*subname*
: subname to be used for the output file.

*n*
: vector $\vec{n}$ normal to the plane.

*alpha*
: coordinate $\alpha$ intersecting the plane.

### Example Usage

```c
scalar * slist = {a,b};
vector * vlist = {c,d};
output_slice_vtu(slist, vlist, "slice_x", (coord){1,0,0}, 0);
output_slice_vtu(slist, vlist, "slice_y", (coord){0,1,0}, 0);
output_slice_vtu(slist, vlist, "slice_z", (coord){0,0,1}, 0);
```

see, also [example](test_output0_vtu.c).

*/

#if dimension > 2

// Function prototypes
void output_slice_vtu_pid(scalar *list, vector *vlist, char *subname, coord n, double _alpha);
#ifdef _MPI
void output_slice_pvtu(scalar *list, vector *vlist, char *subname, coord n, double _alpha);
#endif

trace    
void output_slice_vtu(scalar *list, vector *vlist, char *subname, coord n = {0, 0, 1}, double _alpha = 0){
// Check if MPI is defined
@ if _MPI
  // If MPI is defined, call output_pvtu to handle parallel VTU output
  output_slice_pvtu(list, vlist, subname, n, _alpha);
@ else
  // If MPI is not defined, call output_vtu_pid to handle serial VTU output
  output_slice_vtu_pid(list, vlist, subname, n, _alpha);
@endif
}



/**
### *output_slice_vtu_pid()*: writes one `.vtu` file for the current process */


void output_slice_vtu_pid(scalar *slist, vector *vlist, char *subname, coord n = {0, 0, 1}, double _alpha = 0){

#if defined(_OPENMP)
  int num_omp = omp_get_max_threads();  // Get the number of OpenMP threads
  omp_set_num_threads(1);               // Set the number of OpenMP threads to 1
#endif

  char name[111];                       // Buffer for file name construction
  sprintf(name, "%s.vtu", subname);     // Construct the VTU filename

  FILE *fp = fopen(name, "w");          // Open the VTU file for writing
  if (!fp){
    fprintf(stderr, "Error opening file %s for writing.\n", name);
  }

  /** Define a scalar field to deal with solids and periodic conditions */ 
  scalar per_mask[];
  foreach (){
    per_mask[] = 0.;
    shortcut_slice(n, _alpha);
    if (alpha > 0.){
#if EMBED
      per_mask[] = cs[];
#else
      per_mask[] = 1.;
#endif
    }
  }

  /** Obtain the number of points, cells, and marker used for connectivity*/ 
  vertex scalar marker[];
  long no_points = 0, no_cells = 0;
  foreach_vertex(serial, noauto){
    marker[] = 0.;
    shortcut_slice(n, _alpha);  // if not in the requested plane, we cycle
    marker[] = no_points;
    no_points++;                // Increment the number of points
  }

  foreach (serial, noauto){
    if (per_mask[]){
      no_cells++;               // Increment the number of cells
    }
  }

  /** Every time we write someting we increase the counter to track where the
  data array starts. Cells are defined as VTK_QUAD(=9) defined by 4 points. 
  */
  char type = 9, noffset = 4;
  long count = 0;
  
  /** Write the light data. */
  write_vtu_header(fp, no_points, no_cells);                        // Write the VTU file header
  write_scalar_light_data(fp, slist, vlist, &count, no_cells);      // Write scalar data arrays
  write_points_light_data(fp, &count, no_points);                   // Write points data array
  write_cells_light_data(fp, &count, no_cells, no_cells * noffset); // Write cells data arrays
  write_vtu_appended(fp);                                           // Write the VTU appended data section

  /**   Write the heavy data blocks */
  write_scalar_heavy_data_slice(fp, slist, per_mask, no_cells, n, _alpha);   // Write scalar field data
  write_vector_heavy_data_slice(fp, vlist, per_mask, no_cells, n, _alpha);   // Write vector field data
  write_points_heavy_data_slice(fp, no_points, n, _alpha);                   // Write points data
  write_cell_offsets(fp, no_cells, noffset);                                 // Write cell offsets
  write_cell_types(fp, no_cells, type);                                      // Write cell types
  write_cell_connectivity_slice(fp, marker, per_mask, no_cells, noffset, n); // Write cell connectivity

  /** and close the file */ 
  fputs("\t\n", fp);
  fputs("\t </AppendedData>\n", fp);
  fputs("</VTKFile>\n", fp);
  fflush(fp); // Flush the file buffer
  fclose(fp); // Close the VTU file

#if defined(_OPENMP)
  omp_set_num_threads(num_omp); // Restore the original number of OpenMP threads
#endif
}



/** ### *output_slice_pvtu()*: if MPI, writes one `.pvtu` and `.vtu` for each process */
@ if _MPI 
void output_slice_pvtu(scalar *list, vector *vlist, char *subname, coord n = {0, 0, 1}, double _alpha = 0)
{
  char name[112]; // Buffer for file name construction
  FILE *fp;       // File pointer for file operations

  if (pid() == 0){

    sprintf(name, "%s.pvtu", subname); // Construct the PVTU filename
    fp = fopen(name, "w");             // Open the PVTU file for writing
    if (!fp){
      fprintf(stderr, "Error opening file %s for writing.\n", name);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Write sections of the PVTU file
    write_pvtu_header(fp);                     // Write the header
    write_scalar_light_pdata(fp, list, vlist); // Write scalar data arrays
    write_points_light_pdata(fp);              // Write points data array
    write_pieces_light_pdata(fp, subname);     // Write piece references

    // Close the PVTU file
    fflush(fp); // Flush the file buffer
    fclose(fp); // Close the file
  }
  MPI_Barrier(MPI_COMM_WORLD);
  sprintf(name, "%s_n%3.3d", subname, pid());
  output_slice_vtu_pid(list, vlist, name, n, _alpha);
}
@endif
#endif



/** 
## output_facets_vtu(): Exports the VOF interface.

This function writes a surface from field data using the built-in VOF PLIC
surface. The file also stores a scalar field at the cell center of an 
unstructured grid with the following structure: 

<center>
<table>
<tr>
<td>![](vtk_poly_line.png){ width="80%" }</td>
<td>![](vtk_polygon.png){ width="100%" }</td>
</tr>
<tr>
<td><center>if 2D</center></td> 
<td><center>if 3D</center></td> 
</tr>
</table>
</center>

This extends the routines in
[`output_surfaces.h`](http://basilisk.fr/sandbox/oystelan/output_surfaces.h). by
Oystein Lande in two points: one, this is extended to write a single .vtu file
when using MPI; and two, heavy data is stored in raw binary format.

The arguments and their default values are:

*c*
: vof scalar field.

*s*
: scalar field to store, for instance, curvature.

*subname*
: subname to be used for the output file.

### Example Usage

```c
scalar kappa[];
curvature (f, kappa);
output_facets_vtu(f, kappa, "Interface");
```

see, also [example](test_output0_vtu2.c).

*/

// Function prototypes
void output_facets_pid(scalar c, scalar s, char *subname);
#ifdef _MPI
void output_facets_pvtu(scalar c, scalar s, char *subname);
#endif

trace
void output_facets_vtu(scalar c, scalar s, char *subname){
// Check if MPI is defined
@ if _MPI
  // If MPI is defined, call output_pvtu to handle parallel VTU output
  output_facets_pvtu(c, s, subname);
@ else
  // If MPI is not defined, call output_vtu_pid to handle serial VTU output
  output_facets_pid(c, s, subname);
@endif
}


/** ### *output_facets_pid()*: writes one `.vtu` file for the current process */
void output_facets_pid(scalar c, scalar kappa, char *subname)
{

#if defined(_OPENMP)
  // Save current number of OpenMP threads and set it to 1 for this function
  int num_omp = omp_get_max_threads();
  omp_set_num_threads(1);
#endif

  /** Construct the filename and open the file for writing */ 
  char name[99];
  sprintf(name, "%s.vtu", subname);
  FILE *fp = fopen(name, "w");

  /** Obtain the number of vertices (points) and assets (cells) */ 
  long nverts = 0, nfacets = 0;
  count_vertices_and_facets(c, &nverts, &nfacets);

  /** Declare arrays to store vertex coordinates and facet offsets */ 
  double *pt_array_x = malloc(nverts * sizeof(double));
  double *pt_array_y = malloc(nverts * sizeof(double));
  #if dimension <= 2
    double *pt_array_z = NULL;
  #else
    double *pt_array_z = malloc(nverts * sizeof(double));
  #endif
  long *offsets = malloc(nfacets * sizeof(long));

  /** Populate the arrays with vertex coordinates and facet offsets */ 
  populate_vertex_and_offset_arrays(c, nverts, nfacets, pt_array_x, pt_array_y, pt_array_z, offsets);

  /** Populate the arrays with the curvature kappa */ 
  double *fc_array_kappa = malloc(nfacets * sizeof(double));
  populate_facet_arrays(c, kappa, nverts, nfacets, fc_array_kappa);

  /**  Set cell type: VTK_POLY_LINE for 2D, VTK_POLYGON for 3D */
  #if dimension <= 2
    char type = 4; // VTK_POLY_LINE (=4)
  #else
    char type = 7; // VTK_POLYGON (=7)
  #endif

  /** Counter to keep track of the byte offsets in the VTK file */ 
  long count = 0.;

  /** Write the light data */ 
  write_vtu_header(fp, nverts, nfacets);                       // Write the VTU file header
  write_points_light_data(fp, &count, nverts);                 // Write points data array
  write_cells_light_data(fp, &count, nfacets, nverts);         // Write cells data arrays
  write_scalar_light_data(fp, {kappa}, NULL, &count, nfacets); // Write scalar data arrays
  write_vtu_appended(fp);                                      // Write the VTU appended data section

  /** Write the points (vertices) data to the VTK file */ 
  write_points_heavy_data_array(fp, nverts, pt_array_x, pt_array_y, pt_array_z);
  write_cell_offsets2(fp, nfacets, offsets);                  // Write cell offsets
  write_cell_types(fp, nfacets, type);                        // Write cell types
  write_cell_offsets3(fp, nverts);                            // Write cell offsets
  write_scalar_heavy_data_array(fp, nfacets, fc_array_kappa); // Write cell curvature

  /** and close the file */
  fputs("\t\n", fp);
  fputs("\t </AppendedData>\n", fp);
  fputs("</VTKFile>\n", fp);

  fflush(fp); // Flush the file buffer
  fclose(fp); // Close the VTU file

  free(offsets);
  free(pt_array_x);
  free(pt_array_y);
#if dimension == 3
  free(pt_array_z);
#endif
  free(fc_array_kappa);

#if defined(_OPENMP)
  omp_set_num_threads(num_omp);
#endif
}

/** ### *output_facets_pvtu()*: if MPI, writes one `.pvtu` and `.vtu` for each process */

@ if _MPI
void output_facets_pvtu(scalar c, scalar kappa, char *subname){

#if defined(_OPENMP)
  // Save current number of OpenMP threads and set it to 1 for this function
  int num_omp = omp_get_max_threads();
  omp_set_num_threads(1);
#endif

  /** Obtain the number of vertices (points) and assets (cells) */ 
  long nverts_long = 0, nfacets_long = 0;
  count_vertices_and_facets(c, &nverts_long, &nfacets_long);
  int nverts = (int)nverts_long, nfacets = (int)nfacets_long;

  /** Declare arrays to store vertex coordinates and facet offsets */ 
  double *pt_array_x = malloc(nverts * sizeof(double));
  double *pt_array_y = malloc(nverts * sizeof(double));
  #if dimension <= 2
    double *pt_array_z = NULL;
  #else
    double *pt_array_z = malloc(nverts * sizeof(double));
  #endif
  long *offsets = malloc(nfacets * sizeof(long));

  /** Populate the arrays with vertex coordinates and facet offsets */ 
  populate_vertex_and_offset_arrays(c, nverts_long, nfacets_long, pt_array_x, pt_array_y, pt_array_z, offsets);

  /** Populate the arrays with the curvature kappa */ 
  double *fc_array_kappa = malloc(nfacets * sizeof(double));
  populate_facet_arrays(c, kappa, nverts, nfacets, fc_array_kappa);

  /** Global variables to store the total number of vertices and facets across all processes */ 
  long glob_nverts = 0, glob_nfacets = 0;
  MPI_Allreduce(&nverts_long, &glob_nverts, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&nfacets_long, &glob_nfacets, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

  /** Arrays to store the number of vertices and facets per process */ 
  int list_nverts[npe()], pos_nverts[npe()];
  int list_nfacets[npe()], pos_nfacets[npe()];

  /** Gather the number of vertices and facets from all processes */ 
  MPI_Allgather(&nverts, 1, MPI_INT, list_nverts, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Allgather(&nfacets, 1, MPI_INT, list_nfacets, 1, MPI_INT, MPI_COMM_WORLD);

  /** Calculate the starting positions for vertices and facets for each process */ 
  pos_nverts[0] = 0;
  pos_nfacets[0] = 0;
  for (int i = 1; i < npe(); i++){
    pos_nverts[i] = pos_nverts[i - 1] + list_nverts[i - 1];
    pos_nfacets[i] = pos_nfacets[i - 1] + list_nfacets[i - 1];
  }

  /** Adjust facet offsets for the current process */ 
  for (int i = 0; i < nfacets; i++)
    offsets[i] = offsets[i] + pos_nverts[pid()];

  // Pointers for global arrays (only used by the root process) 
  long *glob_offsets = NULL;
  double *glob_pt_array_x = NULL;
  double *glob_pt_array_y = NULL;
  double *glob_pt_array_z = NULL;
  double *glob_fc_array_kappa = NULL;

  /** Gather the data from each process */
  const int root = 0;
  if (pid() == root){
    // Allocate memory for global arrays on the root process
    glob_offsets = malloc(glob_nfacets * sizeof(long));
    glob_pt_array_x = malloc(glob_nverts * sizeof(double));
    glob_pt_array_y = malloc(glob_nverts * sizeof(double));
    #if dimension == 3
      glob_pt_array_z = malloc(glob_nverts * sizeof(double));
    #endif
    glob_fc_array_kappa = malloc(glob_nfacets * sizeof(double));
  }

  // Gather vertex and facet data from all processes to the root process
  MPI_Gatherv(offsets, nfacets, MPI_LONG, glob_offsets, list_nfacets, pos_nfacets, MPI_LONG, root, MPI_COMM_WORLD);
  MPI_Gatherv(pt_array_x, nverts, MPI_DOUBLE, glob_pt_array_x, list_nverts, pos_nverts, MPI_DOUBLE, root, MPI_COMM_WORLD);
  MPI_Gatherv(pt_array_y, nverts, MPI_DOUBLE, glob_pt_array_y, list_nverts, pos_nverts, MPI_DOUBLE, root, MPI_COMM_WORLD);
#if dimension == 3
  MPI_Gatherv(pt_array_z, nverts, MPI_DOUBLE, glob_pt_array_z, list_nverts, pos_nverts, MPI_DOUBLE, root, MPI_COMM_WORLD);
#endif
  MPI_Gatherv(fc_array_kappa, nfacets, MPI_DOUBLE, glob_fc_array_kappa, list_nfacets, pos_nfacets, MPI_DOUBLE, root, MPI_COMM_WORLD);

  // Free local arrays
  free(offsets);
  free(pt_array_x);
  free(pt_array_y);
  #if dimension == 3
    free(pt_array_z);
  #endif
  free(fc_array_kappa);

  /** Only the root process writes the data to the VTK file */ 
  if (pid() == root){

    #if dimension <= 2
      char type = 4; // VTK_POLY_LINE (=4)
    #else
      char type = 7; // VTK_POLYGON (=7)
    #endif

    // Construct the filename and open the file for writing
    char name[99];
    sprintf(name, "%s.vtu", subname);
    FILE *fp = fopen(name, "w");

    // Counter to keep track of the byte offsets in the VTK file
    long count = 0;
    write_vtu_header(fp, glob_nverts, glob_nfacets);               // Write the VTU file header
    write_points_light_data(fp, &count, glob_nverts);              // Write points data array
    write_cells_light_data(fp, &count, glob_nfacets, glob_nverts); // Write cells data arrays
    write_scalar_light_data(fp, {kappa}, NULL, &count, nfacets);   // Write scalar data arrays
    write_vtu_appended(fp);                                        // Write the VTU appended data section

    // Write the points (vertices) data to the VTK file
    write_points_heavy_data_array(fp, glob_nverts, glob_pt_array_x, glob_pt_array_y, glob_pt_array_z);
    write_cell_offsets2(fp, glob_nfacets, glob_offsets);                  // Write cell offsets
    write_cell_types(fp, glob_nfacets, type);                             // Write cell types
    write_cell_offsets3(fp, glob_nverts);                                 // Write cell offsets
    write_scalar_heavy_data_array(fp, glob_nfacets, glob_fc_array_kappa); // Write cell curvature

    fputs("\t\n", fp);
    fputs("\t </AppendedData>\n", fp);
    fputs("</VTKFile>\n", fp);
    fflush(fp); // Flush the file buffer
    fclose(fp); // Close the VTU file

    // Free allocated memory
    free(glob_offsets);
    free(glob_pt_array_x);
    free(glob_pt_array_y);
#if dimension == 3
    free(glob_pt_array_z);
#endif
    free(glob_fc_array_kappa);
  }

#if defined(_OPENMP)
  omp_set_num_threads(num_omp);
#endif
}
@endif



/** ## postamble: delete macros */
#undef shortcut_slice
#undef shortcut_facets
#undef mfacets