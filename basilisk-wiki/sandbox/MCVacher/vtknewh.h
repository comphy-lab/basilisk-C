/**
Useful to visualize .vtk files. Copied from Daniel Fulster's sandbox. 
*/

struct OutputVTK {
  scalar * list;
  FILE * fp;
  double box[2][2];
};

void output_vtk (struct OutputVTK ov)
{

  scalar * list = ov.list;
  FILE * fp = ov.fp;

  if (ov.box[0][0] == 0. && ov.box[0][1] == 0. && 
      ov.box[1][0] == 0. && ov.box[1][1] == 0.) {
    ov.box[0][0] = X0;      ov.box[0][1] = Y0;
    ov.box[1][0] = X0 + L0; ov.box[1][1] = Y0 + L0;
  }

  fputs ("# vtk DataFile Version 2.0\n"
	 "Basilisk\n"
	 "ASCII\n"
	 "DATASET UNSTRUCTURED_GRID\n \n", fp);

  vertex scalar psi[];
  int np=0;
  foreach_vertex () {
     if (x >= ov.box[0][0] && x <= ov.box[1][0] &&
         y >= ov.box[0][1] && y <= ov.box[1][1] ) {
     	psi[] = np;
     	np++;
     }
  }

  fprintf (fp, "POINTS %d double\n", np);

  foreach_vertex () {
     if (x >= ov.box[0][0] && x <= ov.box[1][0] &&
         y >= ov.box[0][1] && y <= ov.box[1][1] ) 
	      fprintf (fp, "%g %g 0\n", x, y);
  }

  fprintf (fp, "\n");

  int ncells=0;
  foreach () {
     if ( (x - Delta) >= ov.box[0][0] && (x + Delta) <= ov.box[1][0] &&
          (y - Delta) >= ov.box[0][1] && (y + Delta) <= ov.box[1][1] ) 
	ncells++;
  }

  fprintf (fp, "CELLS %i %i \n", ncells, ncells*5);

  foreach () {
     if ( (x - Delta) >= ov.box[0][0] && (x + Delta) <= ov.box[1][0] &&
          (y - Delta) >= ov.box[0][1] && (y + Delta) <= ov.box[1][1] ) 
      fprintf (fp, "4 %i %i %i %i \n", (int) psi[0,0], (int) psi[0,1], (int) psi[1,0], (int) psi[1,1]);
  }
  fprintf (fp, "\n");


  fprintf (fp, "CELL_TYPES %i \n", ncells);
  for (int i = 0; i < ncells; i++)
	fprintf (fp, "8 \n");
  fprintf (fp, "\n");

  fprintf (fp, "POINT_DATA %d\n", np);
  for (scalar s in list) {
    fprintf (fp, "SCALARS %s double\n", s.name);
    fputs ("LOOKUP_TABLE default\n", fp);
    foreach_vertex () 
     if (x >= ov.box[0][0] && x <= ov.box[1][0] &&
         y >= ov.box[0][1] && y <= ov.box[1][1] ) 
		fprintf (fp, "%g\n", s[]);
    fprintf (fp, "\n");
  }

  fflush (fp);
}