/**
 This file contains the functions required to pass on data from Basilisk to
 Paraview/Catalyst2 as unstructured grids. This is done in two parts: grid and
 attributes, named here as cgrid and cattributes. */

typedef struct CGrid {
  long no_points;
  long no_cells;
  double* points;
  long* cells;
} CGrid;

typedef struct CAttributes {
  // A structure for generating and storing point and cell fields.
  // TODO: make it work with lists.
  double* Velocity;
  double* Pressure;
  double* Density;
  CGrid* GridPtr;
} CAttributes;

 /** Like in output_vtu() we need
 special treatment for periodic conditions to match the number of points
 inside foreach_vertex() with the number of cells inside foreach(). To this end,
 we define a mask (and a corresponding macro for readability).
*/

#define shortcut_periodic()     \
  if (Period.x)                 \
    foreach()                   \
      if (x + Delta > X0 + L0)  \
        per_mask[] = 0.;        \
  if (Period.y)                 \
    foreach()                   \
      if (y + Delta > Y0 + L0)  \
        per_mask[] = 0.;        \
  if (Period.z)                 \
    foreach()                   \
      if (z + Delta > Z0 + L0)  \
        per_mask[] = 0.;        \

void InitializeGrid(CGrid* cgrid){
  // We have to modify the way we pass on data. This is based on the output_vtu
  cgrid->no_points = 0;
  cgrid->no_cells = 0;
  cgrid->points = 0;
  cgrid->cells = 0;

  // First, we need to obtain the number of points and the topology
  // This a cheap way to match the number of points and vertices when using periodic conditions since vertex are not defined in that case.
  scalar per_mask[];
  foreach ()
    per_mask[] = 1.;
  shortcut_periodic();

  // We obtain the number of points, cells and a marker used to recover the
  // connectivity.
  vertex scalar marker[];
  long no_points = 0, no_cells=0 ;
  foreach_vertex(serial, noauto){
#if TREE
    marker[] = _k;
#else // Ugly cheat if using multigrid
# if dimension == 2
    marker[] = (point.i-2)*((1 << point.level) + 1) + (point.j-2);
# else
    marker[] = (point.i-2)*sq((1 << point.level) + 1) + (point.j-2)*((1 << point.level) + 1) + (point.k-2);
# endif
#endif
    no_points++;
  }

  foreach (serial, noauto)
    if (per_mask[])
      no_cells++;

  // Now, we need to write the data in the buffers.
  // First, we write points
  if (cgrid->points != 0)
    free(cgrid->points);
  cgrid->points = (double*)malloc(3 * sizeof(double) * no_points);

	foreach_vertex(serial, noauto){
#if !TREE
# if dimension == 2
		int _k = (point.i - 2) * ((1 << point.level) + 1) + (point.j - 2);
# else
		int _k = (point.i - 2) * sq((1 << point.level) + 1) + (point.j - 2) * ((1 << point.level) + 1) + (point.k - 2);
# endif
#endif
		int counter = _k * 3;
		cgrid->points[counter + 0] = x;
		cgrid->points[counter + 1] = y;
#if dimension == 2
		cgrid->points[counter + 2] = 0.;
#endif
#if dimension == 3
		cgrid->points[counter + 2] = z;
#endif
	}
  cgrid->no_points = no_points;

  // Now, we create the hex cells
  if (cgrid->cells != 0)
    free(cgrid->cells);

  cgrid->cells = (long*)malloc( pow(2, dimension) * sizeof(long) * no_cells);
	foreach (serial, noauto){
		if (per_mask[]){
#if !TREE
# if dimension == 2
			int _k = (point.i - 2) * ((1 << point.level)) + (point.j - 2);
# else
			int _k = (point.i - 2) * sq((1 << point.level)) + (point.j - 2) * ((1 << point.level)) + (point.k - 2);
# endif
#endif
			int counter = _k * pow(2, dimension);
			cgrid->cells[counter + 0] = (long)marker[];
			cgrid->cells[counter + 1] = (long)marker[1, 0];
			cgrid->cells[counter + 2] = (long)marker[1, 1];
			cgrid->cells[counter + 3] = (long)marker[0, 1];
#if dimension == 3
			cgrid->cells[counter + 4] = (long)marker[0, 0, 1];
			cgrid->cells[counter + 5] = (long)marker[1, 0, 1];
			cgrid->cells[counter + 6] = (long)marker[1, 1, 1];
			cgrid->cells[counter + 7] = (long)marker[0, 1, 1];
#endif
		}
	}
  cgrid->no_cells = no_cells;
}

void FinalizeGrid(CGrid* cgrid){
  if (cgrid->points) {
    free(cgrid->points);
    cgrid->points = 0;
  }
  if (cgrid->cells) {
    free(cgrid->cells);
    cgrid->cells = 0;
  }
  cgrid->no_points = 0;
  cgrid->no_cells = 0;
}

































void InitializeAttributes(CAttributes* cattributes, CGrid* cgrid) {
  cattributes->GridPtr = cgrid;
  cattributes->Velocity = 0;
  cattributes->Pressure = 0;
}

void FinalizeAttributes(CAttributes* cattributes) {
  if (cattributes->Velocity) {
    free(cattributes->Velocity);
    cattributes->Velocity = 0;
  }
  if (cattributes->Pressure) {
    free(cattributes->Pressure);
    cattributes->Pressure = 0;
  }
}

void UpdateFields(CAttributes* cattributes, double time, scalar p, vector v) {
  scalar per_mask[];
  foreach ()
    per_mask[] = 1.;
  shortcut_periodic();


  int no_cells = cattributes->GridPtr->no_cells;
  if (cattributes->Velocity != 0)
    free(cattributes->Velocity);
  cattributes->Velocity = (double*)malloc(sizeof(double) * no_cells * 3);

  if (cattributes->Pressure != 0)
    free(cattributes->Pressure);
  cattributes->Pressure = (double*)malloc(sizeof(double) * no_cells);

	foreach (serial, noauto){
		if (per_mask[]) {
#if !TREE
# if dimension == 2
			int _k = (point.i - 2) * ((1 << point.level)) + (point.j - 2);
# else
			int _k = (point.i - 2) * sq((1 << point.level)) + (point.j - 2) * ((1 << point.level)) + (point.k - 2);
# endif
#endif
      int counter = _k;
      cattributes->Velocity[counter + 0*no_cells] = v.x[];
      cattributes->Velocity[counter + 1*no_cells] = v.y[];
#if dimension == 2
      cattributes->Velocity[counter + 2*no_cells] = 0.;
#else
      cattributes->Velocity[counter + 2*no_cells] = v.z[];
#endif
      cattributes->Pressure[counter] = p[];
		}
	}
}

#undef shortcut_periodic