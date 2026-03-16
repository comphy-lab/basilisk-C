/**
I would like to thank Pr. Palas Kumar farsoiya for this `file.h` which allow us to visualize our data with a `vtk` format. Here is, his sandbox : [P. K. Farsoiya](https://basilisk.fr/sandbox/farsoiya/).
*/
/** 
.vtk output for AMR: does not interpolate the data on a nxn grid as in the standard outputs of Basilisk, but keep and register the current grid.

This function is almost identical to the one in [D. Fuster sandbox](./../../fuster/Miscellaneous/vtknew.h). The difference lies in the way the data are stored: in D. Fuster code, the data are interpolated on the vertex and register as a vertex field. In my code, the data are registered as cell centered data (this is especially important when plotting discrete quantities such as the level of refinement).

N.B. : the current code doesn't work with periodic() boundary condition, as the vertex shared by the periodic BC doesn't exist twice. (In, other words, when using periodic(right) on a cartesian mesh, the vertex of the right boundaries "doesn't exist anymore", they redirect to the vertex of the left boundary). See the end of this file to a ugly fix (which is for now commented: it works with Basilisk versions prior to the automatic boundary conditions patch and has not been tested with the current version)

[Example](./../AMR_examples/basilisk_logo/AMRlogo.c)

*/


// struct OutputVTK {
//    scalar * list;
//    FILE * fp;
// #if dimension==2
//    double box[2][2];
// #else // dimension=3
//    double box[2][3];
// #endif
// };

// void output_vtk (struct OutputVTK ov)
void output_vtk (scalar * list, FILE * fp)
{
   #if dimension==2
   double box[2][2];
   #else // dimension=3
      double box[2][3];
   #endif

   // scalar * list = ov.list;
   // FILE * fp = ov.fp;

// #if dimension==2
//    if (ov.box[0][0] == 0. && ov.box[0][1] == 0. && 
//        ov.box[1][0] == 0. && ov.box[1][1] == 0.) {
//       ov.box[0][0] = X0;      ov.box[0][1] = Y0;
//       ov.box[1][0] = X0 + L0; ov.box[1][1] = Y0 + L0;
//    }
// #else //dimension=3
// if (ov.box[0][0] == 0. && ov.box[0][1] == 0. && 
//        ov.box[1][0] == 0. && ov.box[1][1] == 0. &&
//        ov.box[0][2] == 0. && ov.box[1][2] == 0.) {
//       ov.box[0][0] = X0;      ov.box[0][1] = Y0;      ov.box[0][2] = Z0;
//       ov.box[1][0] = X0 + L0; ov.box[1][1] = Y0 + L0; ov.box[1][2] = Z0 + L0;
//    }
// #endif

#if dimension==2
   // if (box[0][0] == 0. && box[0][1] == 0. && 
   //     box[1][0] == 0. && box[1][1] == 0.) {
      box[0][0] = X0;      box[0][1] = Y0;
      box[1][0] = X0 + L0; box[1][1] = Y0 + L0;
   // }
   // printf("stop");
#else //dimension=3
// if (box[0][0] == 0. && box[0][1] == 0. && 
//        box[1][0] == 0. && box[1][1] == 0. &&
//        box[0][2] == 0. && box[1][2] == 0.) {
      box[0][0] = X0;      box[0][1] = Y0;      box[0][2] = Z0;
      box[1][0] = X0 + L0; box[1][1] = Y0 + L0; box[1][2] = Z0 + L0;
   // }
#endif
 
   fputs ("# vtk DataFile Version 2.0\n"
	  "Basilisk\n"
	  "ASCII\n"
	  "DATASET UNSTRUCTURED_GRID\n \n", fp);

   vertex scalar psi[];
   int np=0;
   foreach_vertex (serial) {
#if dimension==2
      if (x >= box[0][0] && x <= box[1][0] &&
	  y >= box[0][1] && y <= box[1][1] ) {
#else // dimension=3
        if (x >= box[0][0] && x <= box[1][0] &&
	  y >= box[0][1] && y <= box[1][1] &&
	  z >= box[0][2] && z <= box[1][2] ) {
#endif
	 psi[] = np;
	 np++;
      }
   }
  
   fprintf (fp, "POINTS %d double\n", np);

   foreach_vertex () {
      if (x >= box[0][0] && x <= box[1][0] &&
	  y >= box[0][1] && y <= box[1][1]
#if dimension > 2
	  && z >= box[0][2] && z <= box[1][2]
#endif
	) 
        fprintf (fp, "%g %g %g\n", x, y, z);
   }
     
   fprintf (fp, "\n");

   int ncells=0;
   foreach (serial) {
      if ( (x - Delta/2) >= box[0][0] && (x + Delta/2) <= box[1][0] &&
	   (y - Delta/2) >= box[0][1] && (y + Delta/2) <= box[1][1]
#if dimension > 3
	   && (z - Delta/2) >= box[0][2] && (z + Delta/2) <= box[1][2]
#endif
	 ) 
	 ncells++;
   }

#if dimension==2
   fprintf (fp, "CELLS %i %i\n", ncells, ncells*5);
#else // dimension=3
   fprintf (fp, "CELLS %i %i\n", ncells, ncells*9);
#endif
   
   // for new basilisk version, with automatic boundaries
   psi.dirty = false;   // https://groups.google.com/g/basilisk-fr/c/ASR5Mu83Un0/m/atK7ZjhPCgAJ
   foreach () {
#if dimension==2
      if ( (x - Delta/2) >= box[0][0] && (x + Delta/2) <= box[1][0] &&
	   (y - Delta/2) >= box[0][1] && (y + Delta/2) <= box[1][1] ) 
	 fprintf (fp, "4 %i %i %i %i \n", (int) psi[0,0], (int) psi[0,1], (int) psi[1,0], (int) psi[1,1]);
#else // dimension==3
      if ( (x - Delta/2) >= box[0][0] && (x + Delta/2) <= box[1][0] &&
	   (y - Delta/2) >= box[0][1] && (y + Delta/2) <= box[1][1] &&
	   (z - Delta/2) >= box[0][2] && (z + Delta/2) <= box[1][2]) 
        fprintf (fp, "8 %i %i %i %i %i %i %i %i\n", (int) psi[0,0,0], (int) psi[0,0,1], (int) psi[0,1,0], (int) psi[0,1,1], (int) psi[1,0,0], (int) psi[1,0,1], (int) psi[1,1,0], (int) psi[1,1,1]);
#endif
   }
   fprintf (fp, "\n");

  
   fprintf (fp, "CELL_TYPES %i\n", ncells);
   for (int i = 0; i < ncells; i++)
#if dimension==2
     fprintf (fp, "8\n");  // 8 : VTK_PIXEL
#else // dimension=3
   fprintf (fp, "11\n");  // 11 : VTK_VOXEL
#endif
   
   fprintf (fp, "\n");

   fprintf (fp, "CELL_DATA %d\n", ncells);
   for (scalar s in list) {
      fprintf (fp, "SCALARS %s double\n", s.name);
      fputs ("LOOKUP_TABLE default\n", fp);
      foreach () 
	 if ( (x - Delta/2) >= box[0][0] && (x + Delta/2) <= box[1][0] &&
	      (y - Delta/2) >= box[0][1] && (y + Delta/2) <= box[1][1]
#if dimension>2
	      && (z - Delta/2) >= box[0][2] && (z + Delta/2) <= box[1][2]
#endif
	      ) 
	    fprintf (fp, "%g\n", s[]);
      fprintf (fp, "\n");
   }

   fflush (fp);
}