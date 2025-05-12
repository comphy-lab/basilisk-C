
/** 
.vtk output for AMR: does not interpolate the data on a nxn grid as in the standard outputs of Basilisk, but keep and register the current grid.

This function is almost identical to the one in [D. Fuster sandbox](./../../fuster/Miscellaneous/vtknew.h). The difference lies in the way the data are stored: in D. Fuster code, the data are interpolated on the vertex and register as a vertex field. In my code, the data are registered as cell centered data (this is especially important when plotting discrete quantities such as the level of refinement).

N.B. : the current code doesn't work with periodic() boundary condition, as the vertex shared by the periodic BC doesn't exist twice. (In, other words, when using periodic(right) on a cartesian mesh, the vertex of the right boundaries "doesn't exist anymore", they redirect to the vertex of the left boundary). See the end of this file to a ugly fix (which is for now commented: it works with Basilisk versions prior to the automatic boundary conditions patch and has not been tested with the current version)

[Example](./../AMR_examples/basilisk_logo/AMRlogo.c)

*/


struct OutputVTK {
   scalar * list;
   FILE * fp;
#if dimension==2
   double box[2][2];
#else // dimension=3
   double box[2][3];
#endif
};

void output_vtk (struct OutputVTK ov)
{

   scalar * list = ov.list;
   FILE * fp = ov.fp;

#if dimension==2
   if (ov.box[0][0] == 0. && ov.box[0][1] == 0. && 
       ov.box[1][0] == 0. && ov.box[1][1] == 0.) {
       int maxlev=0;
       foreach(reduction(max:maxlev))
          if (level>maxlev)
             maxlev=level;
      ov.box[0][0] = X0 -1./(2.*pow(2.,maxlev));      ov.box[0][1] = Y0 -1./(2.*pow(2.,maxlev));
      ov.box[1][0] = X0 + L0 +1./(2.*pow(2.,maxlev)); ov.box[1][1] = Y0 + L0 +1./(2.*pow(2.,maxlev));
   }
#else //dimension=3
if (ov.box[0][0] == 0. && ov.box[0][1] == 0. && 
       ov.box[1][0] == 0. && ov.box[1][1] == 0. &&
       ov.box[0][2] == 0. && ov.box[1][2] == 0.) {
       int maxlev=0;
       foreach(reduction(max:maxlev))
          if (level>maxlev)
             maxlev=level;
      ov.box[0][0] = X0 -1./(2.*pow(2.,maxlev));      ov.box[0][1] = Y0 -1./(2.*pow(2.,maxlev));      ov.box[0][2] = Z0 -1./(2.*pow(2.,maxlev));
      ov.box[1][0] = X0 + L0 +1./(2.*pow(2.,maxlev)); ov.box[1][1] = Y0 + L0 +1./(2.*pow(2.,maxlev)); ov.box[1][2] = Z0 + L0 +1./(2.*pow(2.,maxlev));
   }
#endif
 
   fputs ("# vtk DataFile Version 2.0\n"
	  "Basilisk\n"
	  "ASCII\n"
	  "DATASET UNSTRUCTURED_GRID\n \n", fp);

   vertex scalar psi[];
   int np=0;
   foreach_vertex (serial) {
#if dimension==2
      if (x >= ov.box[0][0] && x <= ov.box[1][0] &&
	  y >= ov.box[0][1] && y <= ov.box[1][1] ) {
#else // dimension=3
        if (x >= ov.box[0][0] && x <= ov.box[1][0] &&
	  y >= ov.box[0][1] && y <= ov.box[1][1] &&
	  z >= ov.box[0][2] && z <= ov.box[1][2] ) {
#endif
	 psi[] = np;
	 np++;
      }
   }
  
   fprintf (fp, "POINTS %d double\n", np);

   foreach_vertex () {
      if (x >= ov.box[0][0] && x <= ov.box[1][0] &&
	  y >= ov.box[0][1] && y <= ov.box[1][1]
#if dimension > 2
	  && z >= ov.box[0][2] && z <= ov.box[1][2]
#endif
	) 
        fprintf (fp, "%g %g %g\n", x, y, z);
   }
     
   fprintf (fp, "\n");

   int ncells=0;
   foreach (serial) {
      if ( (x - Delta/2) >= ov.box[0][0] && (x + Delta/2) <= ov.box[1][0] &&
	   (y - Delta/2) >= ov.box[0][1] && (y + Delta/2) <= ov.box[1][1]
#if dimension > 3
	   && (z - Delta/2) >= ov.box[0][2] && (z + Delta/2) <= ov.box[1][2]
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
      if ( (x - Delta/2) >= ov.box[0][0] && (x + Delta/2) <= ov.box[1][0] &&
	   (y - Delta/2) >= ov.box[0][1] && (y + Delta/2) <= ov.box[1][1] ) 
	 fprintf (fp, "4 %i %i %i %i \n", (int) psi[0,0], (int) psi[0,1], (int) psi[1,0], (int) psi[1,1]);
#else // dimension==3
      if ( (x - Delta/2) >= ov.box[0][0] && (x + Delta/2) <= ov.box[1][0] &&
	   (y - Delta/2) >= ov.box[0][1] && (y + Delta/2) <= ov.box[1][1] &&
	   (z - Delta/2) >= ov.box[0][2] && (z + Delta/2) <= ov.box[1][2]) 
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
	 if ( (x - Delta/2) >= ov.box[0][0] && (x + Delta/2) <= ov.box[1][0] &&
	      (y - Delta/2) >= ov.box[0][1] && (y + Delta/2) <= ov.box[1][1]
#if dimension>2
	      && (z - Delta/2) >= ov.box[0][2] && (z + Delta/2) <= ov.box[1][2]
#endif
	      ) 
	    fprintf (fp, "%g\n", s[]);
      fprintf (fp, "\n");
   }

   fflush (fp);
}





/**
Some automatization for MPI run, but we can do better

*/

#if _MPI

struct OutputVTK_MPI {
  scalar * list;
};


event initialize_output_vtk (i=0){
   if (pid()==0){
   FILE * fpvtk = fopen("fields.visit", "w");
   fprintf(fpvtk,"!NBLOCKS %i\n",npe());
   fclose (fpvtk);
   }
}

int ivtk=0;
//event myivtk(i++){
//   ivtk=i;  // remove and add ivtk++ at the end of void output_vtk_MPI ??
//}

void output_vtk_MPI (struct OutputVTK_MPI ov)
{
   scalar * list = ov.list;

   if (pid()==0){
   FILE * fpvtk = fopen("fields.visit", "a");
   for (int k=0;k<=npe()-1;k++)
      fprintf(fpvtk,"fields_%i_%i.vtk\n",ivtk,k);   
   fclose (fpvtk);
   }


   char fname[80];
   sprintf(fname,"fields_%i_%i.vtk",ivtk,pid());  
   FILE * fpvtk = fopen(fname, "w");
   output_vtk ((scalar*){list}, fpvtk);
   fclose (fpvtk);

   ivtk++;   // add if event myivtk removed
}


#endif // _MPI









/**
For periodic(right); simulations. The issue with the above output_vtk function when using periodic(right) is the cell-column lying on the right boundary : this column is not shown when displaying vtk result file with visit.
The problem comes from the (correct) implementation of the periodic(right) boundary condition, which makes "disappear" the vertex on the right boundary. That implies logically that the right boundary is registered as the same coordinates as the left boundary.

The function below simply adds coordinates and scalar values on the right boundary. It may simply be extended to peridodic(top) case (not done here). 

At least the line with the .dirty flag must be added to work with the Basilisk version after the automatic boundary conditions patch.

*/


/* void output_vtk_periodic (struct OutputVTK ov) */
/* { */

/*   scalar * list = ov.list; */
/*   FILE * fp = ov.fp; */

/*   if (ov.box[0][0] == 0. && ov.box[0][1] == 0. &&  */
/*       ov.box[1][0] == 0. && ov.box[1][1] == 0.) { */
/*     ov.box[0][0] = X0;      ov.box[0][1] = Y0; */
/*     ov.box[1][0] = X0 + L0; ov.box[1][1] = Y0 + L0; */
/*   } */

/*   fputs ("# vtk DataFile Version 2.0\n" */
/* 	 "Basilisk\n" */
/* 	 "ASCII\n" */
/* 	 "DATASET UNSTRUCTURED_GRID\n \n", fp); */

/*   vertex scalar psi[]; */
/*   vector psi2[]; */
/*   int np=0; */
/*   foreach_vertex () { */
/*      if (x >= ov.box[0][0] && x <= ov.box[1][0] && */
/*          y >= ov.box[0][1] && y <= ov.box[1][1] ) { */
/*      	psi[] = np; */
/*      	np++; */
/*      } */
/*      if (x>=X0+L0-Delta){   // periodic(right) : we want to add the right boundary */
/* 	psi2.x[] = np; // left vertex */
/* 	np++; */
/* 	psi2.y[] = np;  // right vertex */
/* 	np++; */
/*      } */
/*   } */
  
/*   fprintf (fp, "POINTS %d double\n", np); */

/*   foreach_vertex () { */
/*      if (x >= ov.box[0][0] && x <= ov.box[1][0] && */
/*          y >= ov.box[0][1] && y <= ov.box[1][1] )  */
/* 	      fprintf (fp, "%g %g 0\n", x, y); */
/*      if (x>=X0+L0-Delta){ */
/* 	fprintf (fp, "%g %g 0\n", x, y); */
/* 	fprintf (fp, "%g %g 0\n", x+Delta, y); */
/*      } */
/*   } */
     
/*   fprintf (fp, "\n"); */

/*   int ncells=0; */
/*   foreach () { */
/*      if ( (x - Delta) >= ov.box[0][0] && (x + Delta) <= ov.box[1][0] && */
/*           (y - Delta) >= ov.box[0][1] && (y + Delta) <= ov.box[1][1] )  */
/* 	ncells++; */
/*   } */

/*   fprintf (fp, "CELLS %i %i \n", ncells, ncells*5); */

/*   foreach () { */
/*      if ( (x - Delta) >= ov.box[0][0] && (x + Delta) <= ov.box[1][0] && */
/*           (y - Delta) >= ov.box[0][1] && (y + Delta) <= ov.box[1][1] && */
/* 	  x < X0+L0-Delta )  */
/* 	fprintf (fp, "4 %i %i %i %i \n", (int) psi[0,0], (int) psi[0,1], (int) psi[1,0], (int) psi[1,1]); */
/*      else if ( (x - Delta) >= ov.box[0][0] && (x + Delta) <= ov.box[1][0] && */
/* 	       (y - Delta) >= ov.box[0][1] && (y + Delta) <= ov.box[1][1] ) */
/* 	fprintf (fp, "4 %i %i %i %i \n", (int) psi2.x[0,0], (int) psi2.x[0,1], (int) psi2.y[0,0], (int) psi2.y[0,1]); */
/*   } */
/*   fprintf (fp, "\n"); */


  
  
/*   fprintf (fp, "CELL_TYPES %i \n", ncells); */
/*   for (int i = 0; i < ncells; i++) */
/* 	fprintf (fp, "8 \n"); */
/*   fprintf (fp, "\n"); */

/*   fprintf (fp, "POINT_DATA %d\n", np); */
/*   for (scalar s in list) { */
/*     fprintf (fp, "SCALARS %s double\n", s.name); */
/*     fputs ("LOOKUP_TABLE default\n", fp); */
/*     foreach_vertex () { */
/*      if (x >= ov.box[0][0] && x <= ov.box[1][0] && */
/*          y >= ov.box[0][1] && y <= ov.box[1][1] )  */
/* 		fprintf (fp, "%g\n", s[]); */
/*      if (x>=X0+L0-Delta){ */
/* 	fprintf (fp, "%g\n", s[]); */
/* 	fprintf (fp, "%g\n", s[1,0]); */
/*      } */
/*     } */
/*     fprintf (fp, "\n"); */
/*   } */

/*   fflush (fp); */
/* } */










/**
Below: 2D periodic(top) vtk export. it works with the automatic boundary conditions
*/ 





















/* void output_vtk_periodic_top (struct OutputVTK ov) */
/* { */

/*    scalar * list = ov.list; */
/*    FILE * fp = ov.fp; */

/*    if (ov.box[0][0] == 0. && ov.box[0][1] == 0. &&  */
/*        ov.box[1][0] == 0. && ov.box[1][1] == 0.) { */
/*       ov.box[0][0] = X0;      ov.box[0][1] = Y0; */
/*       ov.box[1][0] = X0 + L0; ov.box[1][1] = Y0 + L0; */
/*    } */

/*    fputs ("# vtk DataFile Version 2.0\n" */
/* 	  "Basilisk\n" */
/* 	  "ASCII\n" */
/* 	  "DATASET UNSTRUCTURED_GRID\n \n", fp); */

/*    vertex scalar psi[]; */
/*    vertex scalar  psi2[];  // we don't have access to psi[..,1] on the top boundary due to periodic. Thus we use a vector to store the corresponding data in psi2[..,0] */
/*    int np=0; */
/*    foreach_vertex () { */
/*       if (x >= ov.box[0][0] && x <= ov.box[1][0] && */
/* 	  y >= ov.box[0][1] && y <= ov.box[1][1]  && */
/* 	  y<Y0+L0-Delta ) { */
/* 	 psi[] = np; */
/* 	 np++; */
/*       } */
/*       if (x >= ov.box[0][0] && x <= ov.box[1][0] && */
/* 	  y >= ov.box[0][1] && y <= ov.box[1][1] && */
/* 	  y>=Y0+L0-Delta){   // periodic(top) : we want to add the top boundary  */
/*  	psi[] = np; // bottom vertex  */
/*  	np++;  */
/*  	psi2[] = np;  // top vertex  */
/*  	np++;  */
/*       } */
/*    } */
  
/*    fprintf (fp, "POINTS %d double\n", np); */

/*    foreach_vertex () { */
/*       if (x >= ov.box[0][0] && x <= ov.box[1][0] && */
/* 	  y >= ov.box[0][1] && y <= ov.box[1][1] && */
/* 	  y<Y0+L0-Delta )  */
/* 	 fprintf (fp, "%g %g 0\n", x, y); */
/*       if (x >= ov.box[0][0] && x <= ov.box[1][0] && */
/* 	  y >= ov.box[0][1] && y <= ov.box[1][1] && */
/* 	  y>=Y0+L0-Delta){  */
/*          fprintf (fp, "%g %g 0\n", x, y);  */
/*          fprintf (fp, "%g %g 0\n", x, y+Delta);  */
/*       }  */
/*    } */
     
/*    fprintf (fp, "\n"); */

/*    int ncells=0; */
/*    foreach () { */
/*       if ( (x - Delta/2) >= ov.box[0][0] && (x + Delta/2) <= ov.box[1][0] && */
/* 	   (y - Delta/2) >= ov.box[0][1] && (y + Delta/2) <= ov.box[1][1] )  */
/* 	 ncells++; */
/*    } */

/*    fprintf (fp, "CELLS %i %i \n", ncells, ncells*5); */

/*    psi.dirty = false;   // for new basilisk version, with automatic boundaries : https://groups.google.com/g/basilisk-fr/c/ASR5Mu83Un0/m/atK7ZjhPCgAJ */
/*    foreach () { */
/*       if ( (x - Delta/2) >= ov.box[0][0] && (x + Delta/2) <= ov.box[1][0] && */
/* 	   (y - Delta/2) >= ov.box[0][1] && (y + Delta/2) <= ov.box[1][1] &&  */
/*  	   y<Y0+L0-Delta)  */
/* 	 fprintf (fp, "4 %i %i %i %i \n", (int) psi[0,0], (int) psi[0,1], (int) psi[1,0], (int) psi[1,1]); */
/*       else if ( (x - Delta/2) >= ov.box[0][0] && (x + Delta/2) <= ov.box[1][0] &&  */
/*  	       (y - Delta/2) >= ov.box[0][1] && (y + Delta/2) <= ov.box[1][1] )  */
/*  	fprintf (fp, "4 %i %i %i %i \n", (int) psi[0,0], (int) psi2[0,0], (int) psi[1,0], (int) psi2[1,0]);  */
/*    } */
/*    fprintf (fp, "\n"); */

  
/*    fprintf (fp, "CELL_TYPES %i \n", ncells); */
/*    for (int i = 0; i < ncells; i++) */
/*       fprintf (fp, "8 \n"); */
/*    fprintf (fp, "\n"); */

/*    fprintf (fp, "CELL_DATA %d\n", ncells); */
/*    for (scalar s in list) { */
/*       fprintf (fp, "SCALARS %s double\n", s.name); */
/*       fputs ("LOOKUP_TABLE default\n", fp); */
/*       foreach ()  */
/* 	 if ( (x - Delta/2) >= ov.box[0][0] && (x + Delta/2) <= ov.box[1][0] && */
/* 	      (y - Delta/2) >= ov.box[0][1] && (y + Delta/2) <= ov.box[1][1] )  */
/* 	    fprintf (fp, "%g\n", s[]); */
/*       fprintf (fp, "\n"); */
/*    } */

/*    fflush (fp); */
/* } */








