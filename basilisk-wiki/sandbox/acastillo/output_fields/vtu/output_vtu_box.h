/**

# XML File Formats

## preamble: define some useful macros

To make the routines (a little bit) easier to follow, we'll define the following 
macros to deal with the reconstruction of the VOF interface

*/

#include "output_vtu_helpers.h"

/**
## output_vtu_box(): Exports full 2D (or 3D) fields.

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

*box*
: we only save data inside the region defined by box.

### Example Usage

```c
scalar * slist = {a,b};
vector * vlist = {c,d};
char *subname = "domain";
coord box[2] = {{-0.25, -0.25}, {0.25, 0.25}};
output_vtu_box(slist, vlist, subname, box);
```

*/

// Function prototypes
void output_vtu_box_pid(scalar *list, vector *vlist, char *subname, coord box[2]);
#ifdef _MPI
void output_pvtu(scalar *list, vector *vlist, char *subname, coord box[2]);
#endif

trace    
void output_vtu_box(scalar *list, vector *vlist, char *subname, coord box[2]){
// Check if MPI is defined
@if _MPI
  // If MPI is defined, call output_pvtu to handle parallel VTU output
  output_pvtu(list, vlist, subname, box);
@else
  // If MPI is not defined, call output_vtu_box_pid to handle serial VTU output
  output_vtu_box_pid(list, vlist, subname, box);
@endif
}

/** ### *output_vtu_box_pid()*: writes one `.vtu` file for the current process */
void output_vtu_box_pid(scalar *list, vector *vlist, char *subname, coord box[2]) {

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

  /** Define a scalar field for cell selection with consistent boundaries */ 
  scalar cell_mask[];
  foreach () {
    cell_mask[] = 0.;  // Initialize to 0
#if dimension == 2
    if (x >= box[0].x && x < box[1].x &&
        y >= box[0].y && y < box[1].y)   
#elif dimension == 3
    if (x >= box[0].x && x < box[1].x &&
        y >= box[0].y && y < box[1].y &&
        z >= box[0].z && z < box[1].z)
#endif
    {  
      cell_mask[] = 1.;
    }
  }

  /** Obtain the number of points, cells, and marker used for connectivity */ 
  vertex scalar marker[];
  vertex scalar vertex_needed[];
  foreach_vertex(){
    marker[] = 0;
    vertex_needed[] = 0;
  }
  
  long no_cells = 0;
  foreach (serial, noauto){
    if (cell_mask[] > 0.5){
      no_cells++;
      // Mark all vertices of this cell as needed      
      vertex_needed[0] = 1;
      vertex_needed[1] = 1;
      vertex_needed[1,1] = 1;
      vertex_needed[0,1] = 1;
#if dimension == 3
      vertex_needed[0,0,1] = 1;
      vertex_needed[1,0,1] = 1;
      vertex_needed[1,1,1] = 1;
      vertex_needed[0,1,1] = 1;
#endif      
    }
  }
  
  // Now count and number the needed vertices
  long no_points = 0;
  foreach_vertex(serial, noauto){
    if (vertex_needed[] > 0.5) {
      marker[] = no_points;
      no_points++;
    }
  }
  
  // Debug output to verify counts
  fprintf(stderr, "VTU Box Output: %ld points, %ld cells\n", no_points, no_cells);
  
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
  write_scalar_heavy_data(fp, list, cell_mask, no_cells);            // Write scalar field data
  write_vector_heavy_data(fp, vlist, cell_mask, no_cells);           // Write vector field data
  write_points_heavy_data_masked(fp, no_points, vertex_needed);      // Write points data
  write_cell_offsets(fp, no_cells, noffset);                        // Write cell offsets
  write_cell_types(fp, no_cells, type);                             // Write cell types
  write_cell_connectivity(fp, marker, cell_mask, no_cells, noffset); // Write cell connectivity

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
void output_pvtu(scalar *list, vector *vlist, char *subname, coord box[2])
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
  output_vtu_box_pid(list, vlist, name, box);
}
@endif