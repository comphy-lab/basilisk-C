#ifndef OUTPUT_VTU_FACETS_H
#define OUTPUT_VTU_FACETS_H

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
[`output_surfaces.h`](/sandbox/oystelan/output_surfaces.h). by
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

#include "output_vtu_helpers.h"

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
  write_points_heavy_data_array(fp, nverts, pt_array_x, pt_array_y, pt_array_z); // Write points data
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

#endif
