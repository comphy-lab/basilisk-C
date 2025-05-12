int findintersectionpoints(double nx, double ny, double alpha, double xc, double yc, double delta, double xi[2], double yi[2], char typei[2])
{
	int ppp = 0;
	double xtmp[2], ytmp[2], underflow = 1.0e-6;
	if (fabs(nx) < underflow)
	{
		ytmp[0] = (alpha - nx * (-0.50)) / ny;
		ytmp[1] = (alpha - nx * (+0.50)) / ny;
		xi[ppp] = xc + (-0.50) * delta;
		yi[ppp] = yc + (ytmp[0]) * delta;
		typei[ppp] = 'l';
		(ppp)++;
		xi[ppp] = xc + (+0.50) * delta;
		yi[ppp] = yc + (ytmp[1]) * delta;
		typei[ppp] = 'r';
		(ppp)++;
	}
	else if (fabs(ny) < underflow)
	{
		xtmp[0] = (alpha - ny * (-0.50)) / nx;
		xtmp[1] = (alpha - ny * (+0.50)) / nx;
		xi[ppp] = xc + (xtmp[0]) * delta;
		yi[ppp] = yc + (-0.50) * delta;
		typei[ppp] = 'b';
		(ppp)++;
		xi[ppp] = xc + (xtmp[1]) * delta;
		yi[ppp] = yc + (+0.50) * delta;
		typei[ppp] = 't';
		(ppp)++;
	}
	else
	{
		xtmp[0] = (alpha - ny * (-0.50)) / nx;
		xtmp[1] = (alpha - ny * (+0.50)) / nx;
		ytmp[0] = (alpha - nx * (-0.50)) / ny;
		ytmp[1] = (alpha - nx * (+0.50)) / ny;

		if (-0.50 <= ytmp[0] && ytmp[0] <= +0.50)
		{
			xi[ppp] = xc + (-0.50) * delta;
			yi[ppp] = yc + (ytmp[0]) * delta;
			typei[ppp] = 'l';
			(ppp)++;
		}
		if (-0.50 <= xtmp[0] && xtmp[0] <= +0.50)
		{
			xi[ppp] = xc + (xtmp[0]) * delta;
			yi[ppp] = yc + (-0.50) * delta;
			typei[ppp] = 'b';
			(ppp)++;
		}
		if (-0.50 <= ytmp[1] && ytmp[1] <= +0.50)
		{
			xi[ppp] = xc + (+0.50) * delta;
			yi[ppp] = yc + (ytmp[1]) * delta;
			typei[ppp] = 'r';
			(ppp)++;
		}
		if (-0.50 <= xtmp[1] && xtmp[1] <= +0.50)
		{
			xi[ppp] = xc + (xtmp[1]) * delta;
			yi[ppp] = yc + (+0.50) * delta;
			typei[ppp] = 't';
			(ppp)++;
		}
	}
	return true;
}


void output_paraview_IF(char *name, double time, double rotationangle, scalar intrfc)
{
	int interfacepoints, cellnumber, i;
	char typei[2], *typeintersect;
	double xi[2], yi[2], *xintersect, *yintersect;
	FILE *fp;
	const double angle = rotationangle * pi / 180.0;
	scalar alpha[];
	vector n[];
	;
	cellnumber = 0;
	foreach ()
	cellnumber++;
	xintersect = (double *)calloc(sizeof(double), cellnumber);
	yintersect = (double *)calloc(sizeof(double), cellnumber);
	typeintersect = (char *)calloc(sizeof(char), cellnumber);
	reconstruction(intrfc, n, alpha);
	interfacepoints = 0;
	foreach_leaf()
	{
		if (intrfc[] > 0.01/*R_VOFLIMIT*/ && intrfc[] < 1.0 - 0.01/*R_VOFLIMIT*/)
		{
			findintersectionpoints(n.x[], n.y[], alpha[], x, y, Delta, xi, yi, typei);
			xintersect[interfacepoints] = xi[0] * cos(angle) - yi[0] * sin(angle);
			yintersect[interfacepoints] = yi[0] * cos(angle) + xi[0] * sin(angle);
			typeintersect[interfacepoints] = typei[0];
			interfacepoints++;
			xintersect[interfacepoints] = xi[1] * cos(angle) - yi[1] * sin(angle);
			yintersect[interfacepoints] = yi[1] * cos(angle) + xi[1] * sin(angle);
			typeintersect[interfacepoints] = typei[1];
			interfacepoints++;
		}
	};
	fp = fopen(name, "w");
	fprintf(fp, "# vtk DataFile Version 2.0\r\n");
	fprintf(fp, "Unstructured Grid\r\n");
	fprintf(fp, "ASCII\r\n");
	fprintf(fp, "DATASET UNSTRUCTURED_GRID\r\n");
	switch (interfacepoints)
	{
	case 0:
	{
		fprintf(fp, "POINTS %d float\r\n", 2);
		foreach_leaf()
		{
			fprintf(fp, "%.9f %.9f 0.0\r\n", x, y);
			fprintf(fp, "%.9f %.9f 0.0\r\n", x, y);
			break;
		}
		fprintf(fp, "CELLS %d %d\r\n", 1, 3);
		fprintf(fp, "%d %d %d\r\n", 2, 0, 1);
		fprintf(fp, "CELL_TYPES %d\r\n", 1);
		fprintf(fp, "3\r\n");
		break;
	}
	default:
	{
		fprintf(fp, "POINTS %d float\r\n", interfacepoints);
		for (i = 0; i < interfacepoints; i++)
			fprintf(fp, "%.9f %.9f 0.0\r\n", xintersect[i], yintersect[i]);
		fprintf(fp, "CELLS %d %d\r\n", interfacepoints / 2, interfacepoints + interfacepoints / 2);
		for (i = 0; i < interfacepoints / 2; i++)
			fprintf(fp, "%d %d %d\r\n", 2, i * 2 + 0, i * 2 + 1);
		fprintf(fp, "CELL_TYPES %d\r\n", interfacepoints / 2);
		for (i = 0; i < interfacepoints / 2; i++)
			fprintf(fp, "3\r\n");
		break;
	}
	};
	free(xintersect);
	free(yintersect);
    free(typeintersect);
	fclose(fp);
	return;
}


void output_xmf_ascii_foreach (scalar * list, vector * vlist, int n, FILE * fp, bool linear)
{
#if defined(_OPENMP)
  int num_omp = omp_get_max_threads();
  omp_set_num_threads(1);
#endif

  vertex scalar marker[];
  int no_points = 0, no_cells=0 ;
  foreach_vertex(){
    marker[] = _k;
    no_points += 1;
  }
  foreach(){
    no_cells += 1;
  }

  fputs ("<?xml version=\"1.0\"?>\n"
	       "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n"
  	     "<Xdmf xmlns:xi=\"http://www.w3.org/2003/XInclude\" Version=\"2.2\">\n", fp);
  fputs ("\t <Domain>\n", fp);
  fputs ("\t\t <Grid GridType=\"Uniform\">\n", fp);
#if dimension == 2
  fprintf (fp,"\t\t\t <Topology TopologyType=\"Quadrilateral\" Dimensions=\"%d\">\n", no_cells);
  fprintf (fp,"\t\t\t\t <DataItem Dimensions=\"%d 4\" NumberType=\"Int\" Precision=\"8\" Format=\"XML\">\n", no_cells);    
  foreach(){
    fprintf (fp, "\t\t\t\t\t %g %g %g %g \n", marker[], marker[1,0], marker[1,1], marker[0,1]);
  }
#endif
#if dimension == 3
  fprintf (fp,"\t\t\t <Topology TopologyType=\"Hexahedron\" Dimensions=\"%d\">\n", no_cells);
  fprintf (fp,"\t\t\t\t <DataItem Dimensions=\"%d 8\" NumberType=\"Int\" Precision=\"8\" Format=\"XML\">\n", no_cells);    
  foreach(){
    fprintf (fp, "\t\t\t\t\t %g %g %g %g %g %g %g %g\n", marker[], marker[1,0,0], marker[1,1,0], marker[0,1,0],marker[0,0,1], marker[1,0,1], marker[1,1,1], marker[0,1,1]);
  }
#endif
  fputs ("\t\t\t\t </DataItem>\n", fp);      
  fputs ("\t\t\t </Topology>\n", fp);      
  fputs ("\t\t\t <Geometry GeometryType=\"XYZ\">\n", fp);      
  fprintf (fp,"\t\t\t\t <DataItem Dimensions=\"%d 3\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\">\n", no_points);
  foreach_vertex(){
#if dimension == 2
    fprintf (fp, "\t\t\t\t\t %g %g 0\n", x, y);
#endif
#if dimension == 3
    fprintf (fp, "\t\t\t\t\t %g %g %g\n", x, y, z);
#endif
  }
  fputs ("\t\t\t\t </DataItem>\n", fp);      
  fputs ("\t\t\t </Geometry>\n", fp);   
  for (scalar s in list) {
    fprintf (fp,"\t\t\t <Attribute Name=\"%s\" Active=\"1\" AttributeType=\"Scalar\" Center=\"Cell\">\n", s.name);
  	fprintf (fp,"\t\t\t\t <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\">\n", no_cells);    
    foreach(){
      fprintf (fp, "\t\t\t\t\t %g\n", val(s));
    }
    fputs ("\t\t\t\t </DataItem>\n", fp);
    fputs ("\t\t\t </Attribute>\n", fp);
  }
  for (vector v in vlist) {
    fprintf (fp,"\t\t\t <Attribute Name=\"%s\" Active=\"1\" AttributeType=\"Vector\" Center=\"Cell\">\n", v.x.name);
  	fprintf (fp,"\t\t\t\t <DataItem Dimensions=\"%d 3\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\">\n", no_cells); 
    foreach(){
#if dimension == 2
      fprintf (fp, "\t\t\t\t\t %g %g 0.\n", val(v.x), val(v.y));
#endif
#if dimension == 3
      fprintf (fp, "\t\t\t\t\t %g %g %g\n", val(v.x), val(v.y), val(v.z));
#endif
    }
    fputs ("\t\t\t\t </DataItem>\n", fp);
    fputs ("\t\t\t </Attribute>\n", fp);
  }
  fputs ("\t\t </Grid>\n", fp);
  fputs ("\t </Domain>\n", fp);
  fputs ("</Xdmf>\n", fp);      

#if defined(_OPENMP)
  omp_set_num_threads(num_omp);
#endif
}




#include <hdf5.h>
#define FAIL -1
#define verbose 0
#define MESG(x) if (verbose) printf("%s\n", x);
void create_output_h5()
{
  hid_t   file_id;
	if (pid() == 0){
		file_id = H5Fcreate("aa.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT); 
		assert(file_id != FAIL); 
		MESG("H5Fcreate succeed");
		H5Fclose(file_id);   
	}
}


void write_xml_h5 (scalar * list, vector * vlist, int n, FILE * fp, char * stamp, int no_cells, int no_points, char * time){
	/* Write the xmf file. */
  fputs ("<?xml version=\"1.0\"?>\n"
	       "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n"
  	     "<Xdmf xmlns:xi=\"http://www.w3.org/2003/XInclude\" Version=\"2.2\">\n", fp);
  fputs ("\t <Domain>\n", fp);
  fputs ("\t\t <Grid GridType=\"Uniform\">\n", fp);
  fprintf(fp, "\t\t<Time Value =\"%s\" />", time);
#if dimension == 2
  fprintf (fp,"\t\t\t <Topology TopologyType=\"Quadrilateral\" Dimensions=\"%d\">\n", no_cells);
  fprintf (fp,"\t\t\t\t <DataItem Dimensions=\"%d 4\" NumberType=\"Int\" Precision=\"8\" Format=\"HDF\">\n", no_cells); 
#endif
#if dimension == 3   
  fprintf (fp,"\t\t\t <Topology TopologyType=\"Hexahedron\" Dimensions=\"%d\">\n", no_cells);
  fprintf (fp,"\t\t\t\t <DataItem Dimensions=\"%d 8\" NumberType=\"Int\" Precision=\"8\" Format=\"HDF\">\n", no_cells); 
#endif
  char h5_filename[1024]; 
  sprintf(h5_filename, "out_%s.h5", stamp);
  
  fprintf (fp,"\t\t\t\t\t %s:/%s/Topology\n", h5_filename,stamp);    
  fputs ("\t\t\t\t </DataItem>\n", fp);
  fputs ("\t\t\t </Topology>\n", fp);
  fputs ("\t\t\t <Geometry GeometryType=\"XYZ\">\n", fp);      
  fprintf (fp,"\t\t\t\t <DataItem Dimensions=\"%d 3\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", no_points);
  fprintf (fp,"\t\t\t\t\t %s:/%s/Geometry/Points\n", h5_filename,stamp);    
  fputs ("\t\t\t\t </DataItem>\n", fp);
  fputs ("\t\t\t </Geometry>\n", fp);
  for (scalar s in list) {
    fprintf (fp,"\t\t\t <Attribute Name=\"%s\" Active=\"1\" AttributeType=\"Scalar\" Center=\"Cell\">\n", s.name);
  	fprintf (fp,"\t\t\t\t <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", no_cells);
  	fprintf (fp,"\t\t\t\t\t %s:/%s/Cell/%s\n", h5_filename,stamp,s.name);     
	  fputs ("\t\t\t\t </DataItem>\n", fp);   
	  fputs ("\t\t\t </Attribute>\n", fp);   
	}
  for (vector v in vlist) {
    fprintf (fp,"\t\t\t <Attribute Name=\"%s\" Active=\"1\" AttributeType=\"Vector\" Center=\"Cell\">\n", v.x.name);
  	fprintf (fp,"\t\t\t\t <DataItem Dimensions=\"%d 3\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", no_cells);
  	fprintf (fp,"\t\t\t\t\t %s:/%s/Cell/%s\n", h5_filename,stamp,v.x.name);     
	  fputs ("\t\t\t\t </DataItem>\n", fp);   
	  fputs ("\t\t\t </Attribute>\n", fp);   
	}
  fputs ("\t\t </Grid>\n", fp);
  fputs ("\t </Domain>\n", fp);
  fputs ("</Xdmf>\n", fp); 
}

// event init (t=0) {
// 	create_output_h5();
// }


@if _MPI
void output_xmf_h5_foreach (scalar * list, vector * vlist, int n, FILE * fp, char * stamp, char * time)
{
  hid_t   acc_tpl1;		            /* File access templates */
  hid_t   file_id;                /* HDF5 file IDs */
  hid_t   group1_id, group2_id;   /* HDF5 group IDs */
  hid_t   dataspace_id;           /* File dataspace ID */
  hid_t   memspace_id;            /* memory dataspace ID */
  hid_t   filespace_id;            /* memory dataspace ID */
  hid_t   dataset_id;             /* Dataset ID */
  hsize_t dims1[1], dims2[2];  		/* dataspace dim sizes */
  herr_t  status;                 /* Generic return value */
  hid_t	  plist_id;               /* property list identifier */
  H5L_info_t info;

	hsize_t	count[2];	          		/* hyperslab selection parameters */
	hsize_t	offset[2];

  vertex scalar marker[];
  int no_points = 0, no_cells=0, no_points_loc=0, no_cells_loc=0;
  foreach_vertex(reduction(+:no_points)){
    no_points += 1;
		no_points_loc += 1;
  }
  foreach(reduction(+:no_cells)){
    no_cells += 1;
		no_cells_loc += 1;
  }

	/* Find the total number of cells/points and per sub-domain */ 
	int no_points_Loc[npe()], no_points_Off[npe()];
	int no_cells_Loc[npe()], no_cells_Off[npe()];
	for (int i = 0; i < npe(); ++i ){
		no_points_Loc[i] = 0.;	
		no_cells_Loc[i] = 0.;	
	}

	no_points_Loc[pid()] = no_points_loc;
	no_cells_Loc[pid()] = no_cells_loc;

	MPI_Allreduce(no_points_Loc, no_points_Off, npe(), MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(no_cells_Loc, no_cells_Off, npe(), MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	
	offset[0] = 0.;
	if (pid()!=0){
		for (int  i = 1; i <= pid(); ++i ){
			offset[0] += no_points_Off[i-1];	
		}
	}

	foreach_vertex(){
    marker[] = _k+offset[0];
	}

	/* Write the xmf file. */
	write_xml_h5 (list, vlist, n, fp, stamp, no_cells, no_points, time);

  /* setup file access template with parallel IO access. */
  acc_tpl1 = H5Pcreate(H5P_FILE_ACCESS); 
	assert(acc_tpl1 != FAIL); 
	MESG("H5Pcreate access succeed");
  /* set Parallel access with communicator */
  status = H5Pset_fapl_mpio(acc_tpl1, MPI_COMM_WORLD, MPI_INFO_NULL); 
	assert(status != FAIL); 
	MESG("H5Pset_fapl_mpio succeed");
  /* Open an existing file collectively. */
  char h5_filename[1024]; 
  sprintf(h5_filename, "out_%s.h5", stamp);
  file_id = H5Fcreate(h5_filename, H5F_ACC_TRUNC, H5P_DEFAULT, acc_tpl1);
	assert(file_id != FAIL); 
	MESG("H5Fopen access succeed");
  /* Release file-access template */
  status = H5Pclose(acc_tpl1); assert(status != FAIL);


	/* Create (or open) a group for this snapshot. */
	status = H5Lexists (file_id, stamp, H5P_DEFAULT);  
	group1_id = (int)status == 0 ? H5Gcreate(file_id, stamp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) : H5Gopen(file_id, stamp, H5P_DEFAULT);
	assert(group1_id != FAIL); MESG("H5Gcreate/H5Gopen succeed"); 

	/****************************************************************************/
	/* Create the dataspace and write the Topology dataset. */
  /* a) Create the dataspace for the dataset. */
  char substamp[1024]; 
  sprintf(substamp, "/%s/Topology", stamp);
#if dimension == 2
	dims2[0] = no_cells; dims2[1] = 4;
#endif
#if dimension == 3
	dims2[0] = no_cells; dims2[1] = 8;
#endif
	dataspace_id = H5Screate_simple(2, dims2, NULL);
	assert(dataspace_id != FAIL); MESG("H5Screate_simple succeed"); 


  /* b) Create the dataset with default properties and close dataspace. */
	status = H5Lexists (group1_id, substamp, H5P_DEFAULT);  
	if ((int)status != 0) H5Ldelete(group1_id, substamp, H5P_DEFAULT);
	assert(status != FAIL); MESG("H5Lexists succeed"); 
	dataset_id = H5Dcreate2(group1_id, substamp, H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	assert(dataset_id != FAIL); MESG("H5Dcreate2 succeed"); 
  H5Sclose(dataspace_id);

	/* c)  Each process defines dataset in memory and writes it to the hyperslab
     in the file. */
#if dimension == 2
	count[0] = no_cells_loc; count[1] = 4;	offset[0] = 0; 	offset[1] = 0;
#endif
#if dimension == 3
	count[0] = no_cells_loc; count[1] = 8;	offset[0] = 0; 	offset[1] = 0;
#endif
	if (pid()!=0){
		for (int  i = 1; i <= pid(); ++i ){
			offset[0] += no_cells_Off[i-1];	
		}
	}
	memspace_id = H5Screate_simple(2, count, NULL);
	assert(memspace_id != FAIL); MESG("H5Screate_simple succeed"); 

	/* d)  Select hyperslab in the file. */
	dataspace_id = H5Dget_space(dataset_id);
	assert(dataspace_id != FAIL); MESG("H5Dget_space succeed"); 

	status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
	assert(status != FAIL); MESG("H5Sselect_hyperslab succeed"); 

#if dimension == 2
  	int *topo_dset = (int *)malloc(no_cells_loc * 4 * sizeof(int));
  foreach(){
		topo_dset[_k*4+0] = marker[];
		topo_dset[_k*4+1] = marker[1,0];
		topo_dset[_k*4+2] = marker[1,1];
		topo_dset[_k*4+3] = marker[0,1];
	}
#endif
#if dimension == 3
  	int *topo_dset = (int *)malloc(no_cells_loc * 8 * sizeof(int));
  foreach(){
		topo_dset[_k*8+0] = marker[];
		topo_dset[_k*8+1] = marker[1,0,0];
		topo_dset[_k*8+2] = marker[1,1,0];
		topo_dset[_k*8+3] = marker[0,1,0];
		topo_dset[_k*8+4] = marker[0,0,1];
		topo_dset[_k*8+5] = marker[1,0,1];
		topo_dset[_k*8+6] = marker[1,1,1];
		topo_dset[_k*8+7] = marker[0,1,1];
	}
#endif

	/* e) Create property list for collective dataset write. */
	acc_tpl1 = H5Pcreate(H5P_DATASET_XFER);
	status = H5Pset_dxpl_mpio(acc_tpl1, H5FD_MPIO_COLLECTIVE);
	assert(status != FAIL); MESG("H5Pset_dxpl_mpio succeed"); 

	status = H5Dwrite(dataset_id, H5T_NATIVE_INT, memspace_id, dataspace_id, acc_tpl1, topo_dset);
	assert(status != FAIL); MESG("H5Dwrite succeed"); 

  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);
  H5Sclose(memspace_id);
  H5Pclose(acc_tpl1);
	/****************************************************************************/

	/****************************************************************************/
  /* Check if subgroup Geometry exists already */
  status = H5Lexists (group1_id, "Geometry", H5P_DEFAULT);  
  if ((int)status == 0)
    group2_id = H5Gcreate(group1_id, "Geometry", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); 
  else {
    group2_id = H5Gopen(group1_id, "Geometry", H5P_DEFAULT); 
  }
	/****************************************************************************/
	/* Create the dataspace and write the Points dataset. */
  /* a) Create the dataspace for the dataset. */
  sprintf(substamp, "/%s/Geometry/Points", stamp);
	dims2[0] = no_points; dims2[1] = 3;

	dataspace_id = H5Screate_simple(2, dims2, NULL); 
	assert(dataspace_id != FAIL); MESG("H5Screate_simple succeed"); 

  /* b) Create the dataset with default properties and close dataspace. */
	status = H5Lexists (group2_id, substamp, H5P_DEFAULT);  
	if ((int)status != 0) H5Ldelete(group2_id, substamp, H5P_DEFAULT);
	assert(status != FAIL); MESG("H5Lexists succeed"); 
	dataset_id = H5Dcreate2(group2_id, substamp, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); 
	assert (dataset_id != FAIL); MESG("H5Dcreate2 succeed");
  H5Sclose(dataspace_id);

	/* c)  Each process defines dataset in memory and writes it to the hyperslab
     in the file. */

	count[0] = no_points_loc; count[1] = 3;	offset[0] = 0; 	offset[1] = 0;
	if (pid()!=0){
		for (int  i = 1; i <= pid(); ++i ){
			offset[0] += no_points_Off[i-1];	
		}
	}

	memspace_id = H5Screate_simple(2, count, NULL);
	assert(memspace_id != FAIL); MESG("H5Screate_simple succeed"); 

	/* d)  Select hyperslab in the file. */
	dataspace_id = H5Dget_space(dataset_id);
	assert(dataspace_id != FAIL); MESG("H5Dget_space succeed"); 

	status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
	assert(status != FAIL); MESG("H5Sselect_hyperslab succeed"); 

  	double *points_dset = (double *)malloc(no_points_loc * 3 * sizeof(double));
  foreach_vertex(){
		points_dset[_k*3+0]=x;    
		points_dset[_k*3+1]=y;    
#if dimension == 2
		points_dset[_k*3+2]=0.; 
#endif
#if dimension == 3  
		points_dset[_k*3+2]=z; 
#endif
	}
	/* e) Create property list for collective dataset write. */
	acc_tpl1 = H5Pcreate(H5P_DATASET_XFER);
	status = H5Pset_dxpl_mpio(acc_tpl1, H5FD_MPIO_COLLECTIVE);
	assert(status != FAIL); MESG("H5Pset_dxpl_mpio succeed"); 

	status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, acc_tpl1, points_dset);
	assert(status != FAIL); MESG("H5Dwrite succeed"); 

  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);
  H5Sclose(memspace_id);
  H5Pclose(acc_tpl1);
	H5Gclose(group2_id); 
/****************************************************************************/
  /* Check if subgroup Cell exists already */
  status = H5Lexists (group1_id, "Cell", H5P_DEFAULT);  
  if ((int)status == 0)
    group2_id = H5Gcreate(group1_id, "Cell", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); 
  else {
    group2_id = H5Gopen(group1_id, "Cell", H5P_DEFAULT); 
  }
/****************************************************************************/
	/* Create the dataspace and write the scalar dataset. */
	/* Write the scalar datasets. */
  double *scalar_dset = (double *)malloc(no_cells_loc * 1 * sizeof(double));
  for (scalar s in list) {
	  /* a) Create the dataspace for the dataset. */
	  sprintf(substamp, "/%s/Cell/%s", stamp,s.name);
		dims2[0] = no_cells; dims2[1] = 1;

		dataspace_id = H5Screate_simple(2, dims2, NULL); 
		assert(dataspace_id != FAIL); MESG("H5Screate_simple succeed");

		/* b) Create the dataset with default properties and close dataspace. */
		status = H5Lexists (group2_id, substamp, H5P_DEFAULT);  
		if ((int)status != 0) H5Ldelete(group2_id, substamp, H5P_DEFAULT);
		assert(status != FAIL); MESG("H5Lexists succeed"); 
		dataset_id = H5Dcreate2(group2_id, substamp, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); 
		assert (dataset_id != FAIL); MESG("H5Dcreate2 succeed");
		H5Sclose(dataspace_id);

		/* c)  Each process defines dataset in memory and writes it to the hyperslab
		   in the file. */

		count[0] = no_cells_loc; count[1] = 1; offset[0] = 0; 	offset[1] = 0;
		if (pid()!=0){
			for (int  i = 1; i <= pid(); ++i ){
				offset[0] += no_cells_Off[i-1];	
			}
		}

		memspace_id = H5Screate_simple(2, count, NULL);
		assert(memspace_id != FAIL); MESG("H5Screate_simple succeed"); 

		/* d)  Select hyperslab in the file. */
		dataspace_id = H5Dget_space(dataset_id);
		assert(dataspace_id != FAIL); MESG("H5Dget_space succeed"); 

		status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
		assert(status != FAIL); MESG("H5Sselect_hyperslab succeed"); 

		foreach()
			scalar_dset[_k]=s[];   

		/* e) Create property list for collective dataset write. */
		acc_tpl1 = H5Pcreate(H5P_DATASET_XFER);
		status = H5Pset_dxpl_mpio(acc_tpl1, H5FD_MPIO_COLLECTIVE);
		assert(status != FAIL); MESG("H5Pset_dxpl_mpio succeed"); 

		status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, acc_tpl1, scalar_dset);
		assert(status != FAIL); MESG("H5Dwrite succeed"); 

		H5Dclose(dataset_id);
		H5Sclose(dataspace_id);
		H5Sclose(memspace_id);
		H5Pclose(acc_tpl1);
/****************************************************************************/		
	}
/****************************************************************************/
	/* Create the dataspace and write the vector datasets. */
	/* Write the scalar datasets. */
  double *vector_dset = (double *)malloc(no_cells_loc * 3 * sizeof(double));
  for (vector v in vlist) {
	  /* a) Create the dataspace for the dataset. */
	  sprintf(substamp, "/%s/Cell/%s", stamp,v.x.name);
		dims2[0] = no_cells; dims2[1] = 3;

		dataspace_id = H5Screate_simple(2, dims2, NULL); 
		assert(dataspace_id != FAIL); MESG("H5Screate_simple succeed");

		/* b) Create the dataset with default properties and close dataspace. */
		status = H5Lexists (group2_id, substamp, H5P_DEFAULT);  
		if ((int)status != 0) H5Ldelete(group2_id, substamp, H5P_DEFAULT);
		assert(status != FAIL); MESG("H5Lexists succeed"); 
		dataset_id = H5Dcreate2(group2_id, substamp, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); 
		assert (dataset_id != FAIL); MESG("H5Dcreate2 succeed");
		H5Sclose(dataspace_id);

		/* c)  Each process defines dataset in memory and writes it to the hyperslab
		   in the file. */

		count[0] = no_cells_loc; count[1] = 3; offset[0] = 0; 	offset[1] = 0;
		if (pid()!=0){
			for (int  i = 1; i <= pid(); ++i ){
				offset[0] += no_cells_Off[i-1];	
			}
		}

		memspace_id = H5Screate_simple(2, count, NULL);
		assert(memspace_id != FAIL); MESG("H5Screate_simple succeed"); 

		/* d)  Select hyperslab in the file. */
		dataspace_id = H5Dget_space(dataset_id);
		assert(dataspace_id != FAIL); MESG("H5Dget_space succeed"); 

		status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
		assert(status != FAIL); MESG("H5Sselect_hyperslab succeed"); 

		foreach(){
			vector_dset[_k*3+0]=v.x[];   
			vector_dset[_k*3+1]=v.y[];   
			vector_dset[_k*3+2]=v.z[];   
		}

		/* e) Create property list for collective dataset write. */
		acc_tpl1 = H5Pcreate(H5P_DATASET_XFER);
		status = H5Pset_dxpl_mpio(acc_tpl1, H5FD_MPIO_COLLECTIVE);
		assert(status != FAIL); MESG("H5Pset_dxpl_mpio succeed"); 

		status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, acc_tpl1, vector_dset);
		assert(status != FAIL); MESG("H5Dwrite succeed"); 

		H5Dclose(dataset_id);
		H5Sclose(dataspace_id);
		H5Sclose(memspace_id);
		H5Pclose(acc_tpl1);
/****************************************************************************/		
	}
	H5Gclose(group2_id); 
	H5Gclose(group1_id); 
  H5Fclose(file_id);
}
@else
void output_xmf_h5_foreach (scalar * list, vector * vlist, int n, FILE * fp, char * stamp, char * time)
{
#if defined(_OPENMP)
  int num_omp = omp_get_max_threads();
  omp_set_num_threads(1);
#endif

  hid_t   file_id;                /* HDF5 file IDs */
  hid_t   group1_id, group2_id;   /* HDF5 group IDs */
  hid_t   dataspace_id;           /* File dataspace ID */
  hid_t   memspace_id;            /* memory dataspace ID */
  hid_t   dataset_id;             /* Dataset ID */
  hsize_t dims1[1], dims2[2];  		/* dataspace dim sizes */
  herr_t  status;                 

  vertex scalar marker[];
  int no_points = 0, no_cells=0 ;
  int i=0;
  foreach_vertex(){
//     marker[] = _k;
//     printf("_k is %d \n",_k);
    marker[]=i;
    i++;
    no_points += 1;
  }
  foreach(){
    no_cells += 1;
  }

	/* Write the xmf file. */
	write_xml_h5 (list, vlist, n, fp, stamp, no_cells, no_points, time);

	/* Write the corresponding H5 file. */
	/* Open existing H5 file. */
    char h5_filename[1024]; 
    sprintf(h5_filename, "out_%s.h5", stamp);
    file_id = H5Fcreate(h5_filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	assert(file_id != FAIL); MESG("H5Fopen succeed"); 

	/* Create (or open) a group for this snapshot. */
	status = H5Lexists (file_id, stamp, H5P_DEFAULT);  
	group1_id = (int)status == 0 ? H5Gcreate(file_id, stamp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) : H5Gopen(file_id, stamp, H5P_DEFAULT);
	assert(group1_id != FAIL); MESG("H5Gcreate/H5Gopen succeed"); 

	/* Create the dataspace and write the Topology dataset. */
  char substamp[1024]; 
  sprintf(substamp, "/%s/Topology", stamp);
#if dimension == 2
	dims2[0] = no_cells; dims2[1] = 4;
#endif
#if dimension == 3
	dims2[0] = no_cells; dims2[1] = 8;
#endif
	dataspace_id = H5Screate_simple(2, dims2, NULL);
	assert(dataspace_id != FAIL); MESG("H5Screate_simple succeed"); 

	status = H5Lexists (group1_id, substamp, H5P_DEFAULT);  
	if ((int)status != 0) H5Ldelete(group1_id, substamp, H5P_DEFAULT);
	assert(status != FAIL); MESG("H5Lexists succeed"); 

	dataset_id = H5Dcreate2(group1_id, substamp, H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	assert(dataset_id != FAIL); MESG("H5Dcreate2 succeed"); 

#if dimension == 2
    	int *topo_dset = (int *)malloc(no_cells * 4 * sizeof(int));
        i=0;
  foreach(){
// 		topo_dset[_k*4+0] = marker[];
// 		topo_dset[_k*4+1] = marker[1,0];
// 		topo_dset[_k*4+2] = marker[1,1];
// 		topo_dset[_k*4+3] = marker[0,1];
        topo_dset[i*4+0] = marker[];
		topo_dset[i*4+1] = marker[1,0];
		topo_dset[i*4+2] = marker[1,1];
		topo_dset[i*4+3] = marker[0,1];
        i++;
	}
#endif
#if dimension == 3
    	int *topo_dset = (int *)malloc(no_cells * 8 * sizeof(int));
  foreach(){
		topo_dset[_k*8+0] = marker[];
		topo_dset[_k*8+1] = marker[1,0,0];
		topo_dset[_k*8+2] = marker[1,1,0];
		topo_dset[_k*8+3] = marker[0,1,0];
		topo_dset[_k*8+4] = marker[0,0,1];
		topo_dset[_k*8+5] = marker[1,0,1];
		topo_dset[_k*8+6] = marker[1,1,1];
		topo_dset[_k*8+7] = marker[0,1,1];
	}
#endif
	H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, topo_dset);
    free(topo_dset);
 	H5Sclose(dataspace_id);
  H5Dclose(dataset_id);

  /* Check if subgroup Geometry exists already */
  status = H5Lexists (group1_id, "Geometry", H5P_DEFAULT);  
  if ((int)status == 0)
    group2_id = H5Gcreate(group1_id, "Geometry", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); 
  else {
    group2_id = H5Gopen(group1_id, "Geometry", H5P_DEFAULT); 
  }

	/* Write the points dataset. */
  sprintf(substamp, "/%s/Geometry/Points", stamp);
	dims2[0] = no_points; dims2[1] = 3;
	dataspace_id = H5Screate_simple(2, dims2, NULL); assert (dataspace_id != FAIL); MESG("H5Screate_simple succeed");

	status = H5Lexists (group2_id, substamp, H5P_DEFAULT);  
	if ((int)status != 0) H5Ldelete(group2_id, substamp, H5P_DEFAULT);
	dataset_id = H5Dcreate2(group2_id, substamp, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); 
	assert (dataset_id != FAIL); MESG("H5Dcreate2 succeed");


    	double *points_dset = (double *)malloc(no_points * 3 * sizeof(double));
        i=0;
  	foreach_vertex(){
		points_dset[i/*_k*/*3+0]=x;    
		points_dset[i/*_k*/*3+1]=y;    
#if dimension == 2
		points_dset[i/*_k*/*3+2]=0.; 
#endif
#if dimension == 3  
		points_dset[i/*_k*/*3+2]=z; 
#endif
        i++;
	}
	H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, points_dset);
    free(points_dset);
 	H5Sclose(dataspace_id);
	H5Dclose(dataset_id);
	H5Gclose(group2_id); 

  	/* Check if subgroup Cell exists already */
  	status = H5Lexists (group1_id, "Cell", H5P_DEFAULT);  
  	if ((int)status == 0)
    		group2_id = H5Gcreate(group1_id, "Cell", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); 
  	else {
		group2_id = H5Gopen(group1_id, "Cell", H5P_DEFAULT); 
	}
	/* Write the scalar datasets. */
    double *scalar_dset = (double *)malloc(no_cells * sizeof(double));
  	for (scalar s in list) {
	  sprintf(substamp, "/%s/Cell/%s", stamp,s.name);
		dims2[0] = no_cells; dims2[1] = 1;
		dataspace_id = H5Screate_simple(2, dims2, NULL); assert (dataspace_id != FAIL); MESG("H5Screate_simple succeed");

		status = H5Lexists (group2_id, substamp, H5P_DEFAULT);  
		if ((int)status != 0) H5Ldelete(group2_id, substamp, H5P_DEFAULT);
		dataset_id = H5Dcreate2(group2_id, substamp, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); 
		assert (dataset_id != FAIL); MESG("H5Dcreate2 succeed");
        int i=0;
		foreach(){
			scalar_dset[i/*_k*/]=s[]; 
            i++;
        }
		
		H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, scalar_dset);
	 	H5Sclose(dataspace_id);
		H5Dclose(dataset_id);
	}
    free(scalar_dset);
	/* Write the vector datasets. */
    	double *vector_dset = (double *)malloc(no_cells * 3 * sizeof(double));
  	for (vector v in vlist) {
	  	sprintf(substamp, "/%s/Cell/%s", stamp,v.x.name);
		dims2[0] = no_cells; dims2[1] = 3;
		dataspace_id = H5Screate_simple(2, dims2, NULL); assert (dataspace_id != FAIL); MESG("H5Screate_simple succeed");

		status = H5Lexists (group2_id, substamp, H5P_DEFAULT);  
		if ((int)status != 0) H5Ldelete(group2_id, substamp, H5P_DEFAULT);
		dataset_id = H5Dcreate2(group2_id, substamp, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); 
		// assert (dataset_id != FAIL); MESG("H5Dcreate2 succeed");
        int i=0;
		foreach(){
			vector_dset[i/*_k*/*3+0]=v.x[];   
			vector_dset[i/*_k*/*3+1]=v.y[];   
			vector_dset[i/*_k*/*3+2]=v.z[]; 
            i++;
		}
		H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vector_dset);
	 	H5Sclose(dataspace_id);
		H5Dclose(dataset_id);
	}
    free(vector_dset);
	status =H5Gclose(group2_id); 
    assert(status==0);
	status =H5Gclose(group1_id);
    assert(status==0);
	status =H5Fclose(file_id);
    assert(status==0);
    
#if defined(_OPENMP)
  omp_set_num_threads(num_omp);
#endif
}
@endif
