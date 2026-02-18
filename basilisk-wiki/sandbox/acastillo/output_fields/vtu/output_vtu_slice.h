#ifndef OUTPUT_VTU_SLICE_H
#define OUTPUT_VTU_SLICE_H

/**
# Export specific slices 
*/

#include "output_vtu_helpers.h"

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


void output_slice_vtu_pid(scalar *list, vector *vlist, char *subname, coord n = {0, 0, 1}, double _alpha = 0){

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
  write_scalar_light_data(fp, list, vlist, &count, no_cells);      // Write scalar data arrays
  write_points_light_data(fp, &count, no_points);                   // Write points data array
  write_cells_light_data(fp, &count, no_cells, no_cells * noffset); // Write cells data arrays
  write_vtu_appended(fp);                                           // Write the VTU appended data section

  /**   Write the heavy data blocks */
  write_scalar_heavy_data_slice(fp, list, per_mask, no_cells, n, _alpha);   // Write scalar field data
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
