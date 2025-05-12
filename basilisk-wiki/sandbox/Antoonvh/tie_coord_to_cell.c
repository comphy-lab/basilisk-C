/**
# Is it possible to reference arrays from cell data?

Yes. Here is an example of how array data (a coord) can be assigned to
a cell. See [the output](tie_coord_to_cell/out).
*/

scalar data_pointer[];  

typedef union{
  double* p;
  double d;
} datic;

// Our data we wish to store in some cell
coord a = {0.3123, 0.6123, 35};

int main() {
  // This is an important assertion:
  assert (sizeof(double*) == sizeof (double)); 
  init_grid (N);
  Point point = locate (a.x, a.y, a.z); {        // Find cell
    // datic ap = {.p = (double*)&a};            // option 1 (not used)
    datic ap = {.p = malloc (3*sizeof(double))}; // allocate array
    memcpy (ap.p, &a, sizeof(coord));            // assign values
    data_pointer[] = ap.d;                       // store pointer
  }
  // Search for the data using the `data_pointer` field
  foreach() {
    if (data_pointer[]) {
      datic ap = {.d = data_pointer[]};
      coord * ar = (coord*)ap.p, cc = {x, y, z};
      puts ("A position was found:");
      foreach_dimension(3) 
	printf ("%g is near %g?\n", ar->x, cc.x);
    }
  }
  // cleanup
  foreach() {
    if (data_pointer[]) {
      free((datic){.d = data_pointer[]}.p); // Not needed with option 1.
      data_pointer[] = 0.;
    }
  }
}
