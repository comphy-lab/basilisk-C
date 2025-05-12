/**
** Deprecated! use: output_htg.h instead! *

This headerfile enables the export of AMR-Meshes (Quadtree/Octree) with or without MPI. In this current version Paraview v5.9.0 and v5.9.1 are compatible as well as the Paraview-Nighthly-Builds (use *output_htg_data_v2_xml20(...)* in *output_htg(...)* for this case).
The HTG file format is still under development so compatibility with future Paraview versions might be limited.
Advantage of the HTG file format is the ~7x smaller size due to implicit point location. Disatvantage are besided the currently known bugs/errors (see below) the fact that x-z axis are swapped.

# *output_htg(...)*
 * is the functions to use in your basilisk code. It uses following subroutines
 1. *output_htg_data_v2()*
    - This functions writes the actual data
 2. *output pvd()* [Only pid()==0]
    - This function writes the pvd-files (overwrites the tail with the new values)
  		including the correct time step value
                
## Example Usage

    event report(i+=10){
        boundary(all);
        restriction (all);
        char path[]="htg/";
        char prefix[80];
        sprintf(prefix, "data_%06d", i);  
        output_htg((scalar *) {cs, d,p},(vector *){u}, path, prefix, i, t);
    }

## VTK-Format
 * The VTK-Format is a mixed "ascii" and "appended" format 
 * *output_htg_data_v2()* is the default and compatible with HTG file format version 1.0
 * *output_htg_data_v2_xml20()* enables the export in HTG file format version 2.0 (currently only in Paraview nightly)

## Known bugs/errors (this exporter / htg-fileformat)
 * Contour filter does not work
    - Workaround: HyperTreeGridToUnstructuredGrid -> ResampleToImage
 * HyperTreeGridToDualGrid works only for 1 time step
*/


#include "output_pvd.h"
	  
void output_htg(scalar * list, vector * vlist, char* path, char* prefix, int i, double t);
void output_htg_data_v2(scalar * list, vector * vlist, FILE * fp);
void output_htg_data_v2_xml20(scalar * list, vector * vlist, FILE * fp);

void output_htg(scalar * list, vector * vlist, char* path, char* prefix, int i, double t) {
	FILE * fp ;

	char htg_name[80];  	  
	sprintf(htg_name, "%s%s.htg", path, prefix);  

	fp = fopen(htg_name, "w");	  
	if(fp == NULL){
		printf("output_htg_pv_590_mpi.h : %s could not be opened\n Does the Folder exist?\n", htg_name);
		exit(1);
	}  
	output_htg_data_v2((scalar *) list,(vector *)vlist,fp);
	//~ output_htg_data_v2_xml20((scalar *) list,(vector *)vlist,fp);
	fclose(fp);

	if(pid()==0){
		bool firstTimeWritten = false;
		char pvd_name[]="hypertree.pvd";	  
		fp = fopen(pvd_name, "r+");
		if( (i == 0) ||  (fp == NULL) ) {
			fp = fopen(pvd_name,"w");
			firstTimeWritten = true;
		}
		output_pvd(htg_name, t, fp, firstTimeWritten);
		fclose(fp);
	}
	#if _MPI  
	MPI_Barrier(MPI_COMM_WORLD);    
	#endif
}

/**
# *output_htg_data_v2*
 * Uses MPI_Accumulate to calculate the correct byte offset
 * XML-Headers are only written by [pid()==0]
 * The data is written using the *foreach_level(lvl)* iterator
 */
void output_htg_data_v2(scalar * list, vector * vlist, FILE * fp)
{
	#if defined(_OPENMP)
		int num_omp = omp_get_max_threads();
		omp_set_num_threads(1);
	#endif


	/* Caluculate Offsets */
	long verticesPerLevelLokal[grid->maxdepth+1];  
	memset(verticesPerLevelLokal, 0,(grid->maxdepth+1)*sizeof(long));
  
	for (int lvl = 0; lvl <= depth(); lvl++) {
		foreach_level(lvl) 
			if(is_local(cell))
				verticesPerLevelLokal[lvl]++;	  		
	}
  
	long offset = 0;

	/*
	*  myOffsetPerLevel[2*lvl+0] contains the level depending offset to write
	*  myOffsetPerLevel[2*lvl+1] contains the offset for the next level (total number of written data per level)  
	* 		this equals the number of Cells Per Level */
	long * myOffsetPerLevel;      
	int window_elements = 2*(grid->maxdepth+1);  
    #if _MPI
	MPI_Win win;
	MPI_Win_allocate((MPI_Aint)window_elements*sizeof(long),\
	sizeof(long), MPI_INFO_NULL, MPI_COMM_WORLD, &myOffsetPerLevel, &win);				
	MPI_Win_fence(0, win);
	memset(myOffsetPerLevel, 0, window_elements*sizeof(long));
	MPI_Barrier(MPI_COMM_WORLD);	
	MPI_Win_fence(0, win);

	for(int lvl = 0; lvl < grid->maxdepth+1; ++lvl){
		for(int pe = 0; pe < npe();++pe){
			// Offset for next level  in myOffsetPerLevel[1]
			MPI_Accumulate(&verticesPerLevelLokal[lvl], 1, MPI_LONG, pe, 2*lvl+1, 1, MPI_LONG, MPI_SUM, win);    
			if( pe > pid()){
				// Offset in myOffsetPerLevel[0]
				MPI_Accumulate(&verticesPerLevelLokal[lvl], 1, MPI_LONG, pe, 2*lvl+0, 1, MPI_LONG, MPI_SUM, win);
			}
		}
	}
	
	MPI_Win_fence(0, win);	
    #else
    myOffsetPerLevel = (long*)malloc(window_elements*sizeof(long));
    for(int lvl = 0; lvl < grid->maxdepth+1; ++lvl){
		myOffsetPerLevel[2*lvl+1]=verticesPerLevelLokal[lvl];
		myOffsetPerLevel[2*lvl+0]=0;
	}
    #endif
    // Count global vertices
    long vertices = 0;
    for(int lvl=0; lvl<=grid->maxdepth;++lvl){
	vertices +=myOffsetPerLevel[2*lvl+1];
    }

    // Count Descriptor Bits
    long num_descbit = 0;
    for (int lvl = 0; lvl < grid->maxdepth; ++lvl) {
  	num_descbit +=myOffsetPerLevel[2*lvl+1];
    }
	
	/* File Header vtkHyperTreeGrid */
	if (pid() == 0){	  
		fputs ("<VTKFile type=\"HyperTreeGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n",fp);
		#if dimension==2
		fprintf(fp, "\t<HyperTreeGrid BranchFactor=\"2\" TransposedRootIndexing=\"0\" Dimensions=\"%d %d %d\">\n", 2,2,1);
		#elif dimension==3
		fprintf(fp, "\t<HyperTreeGrid BranchFactor=\"2\" TransposedRootIndexing=\"0\" Dimensions=\"%d %d %d\">\n", 2,2,2);
		#endif
		fputs ("\t\t<Grid>\n",fp);
		fprintf(fp, "\t\t\t<DataArray type=\"Float32\" Name=\"XCoordinates\" NumberOfTuples=\"2\" format=\"ascii\" RangeMin=\"%g\" RangeMax=\"%g\">\n",X0, X0+L0);
		fprintf (fp, "\t\t\t\t%g %g",X0, X0+L0);
		fputs ("\n\t\t\t</DataArray>\n",fp);
		fprintf(fp, "\t\t\t<DataArray type=\"Float32\" Name=\"YCoordinates\" NumberOfTuples=\"2\" format=\"ascii\" RangeMin=\"%g\" RangeMax=\"%g\">\n",Y0, Y0+L0);
		fprintf (fp, "\t\t\t\t%g %g",Y0, Y0+L0);
		fputs ("\n\t\t\t</DataArray>\n",fp);	  
		fprintf(fp, "\t\t\t<DataArray type=\"Float32\" Name=\"ZCoordinates\" NumberOfTuples=\"2\" format=\"ascii\" RangeMin=\"%g\" RangeMax=\"%g\">\n",Z0, Z0+L0);
		#if dimension==2
		fprintf (fp, "\t\t\t\t%g %g", 0., 0.);	  
		#elif dimension==3
		fprintf (fp, "\t\t\t\t%g %g",Z0, Z0+L0);	  
		#endif
		fputs ("\n\t\t\t</DataArray>\n",fp);	  
		fputs ("\t\t</Grid>\n",fp);
		fputs ("\t\t<Trees>\n",fp);

		// Tree Begin
		fprintf (fp,"\t\t\t<Tree Index=\"0\" NumberOfLevels=\"%d\" NumberOfVertices=\"%li\">\n", grid->maxdepth +1, vertices);

		// Array Descriptor Bits
		fprintf (fp,"\t\t\t\t<DataArray type=\"Bit\" Name=\"Descriptor\" NumberOfTuples=\"%li\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"1\">\n", num_descbit);
		offset = ftell(fp);	  
	} // pid()==0
	#if _MPI
	MPI_Barrier(MPI_COMM_WORLD);	
	MPI_Bcast(&offset, 1, MPI_LONG, 0, MPI_COMM_WORLD);
	#endif
	/* Print Bitmask per Process per Level */
	int cell_size=2*sizeof(char);    
	for(int lvl=0; lvl<grid->maxdepth;++lvl){ // Bitmask is "0" at grid->maxdepth	and can therefore be neglected
		fseek (fp, offset + myOffsetPerLevel[2*lvl+0]*cell_size, SEEK_SET); 
		foreach_level(lvl) {
			if (is_local(cell)) {
				if(is_leaf(cell)){
					fputs("0 ",fp);
				} else {
					fputs("1 ",fp);
				}			
			} else {
				continue;
			}
		}				
		offset += myOffsetPerLevel[2*lvl+1]*cell_size;
	}
	fseek (fp, offset, SEEK_SET);
    
     /* Print VerticesByLevel */     		 
	if(pid()==0){
		fputs ("\n\t\t\t\t</DataArray>\n",fp);
		// Array NbVerticesByLevel
		fprintf(fp,"\t\t\t\t<DataArray type=\"Int64\" Name=\"NbVerticesByLevel\" NumberOfTuples=\"%d\" format=\"ascii\" RangeMin=\"1\" RangeMax=\"%li\" >\n" ,grid->maxdepth + 1, myOffsetPerLevel[2*grid->maxdepth + 1]);
		fputs("\t\t\t\t\t", fp);

		for (int lvl = 0; lvl <= grid->maxdepth; lvl++) {
			fprintf(fp,"%li ", myOffsetPerLevel[2*lvl+1]);
		}
		
		fputs ("\n\t\t\t\t</DataArray>\n",fp);
		//~ fputs ("\t\t\t\t<PointData>\n",fp);
		//~ fputs ("\t\t\t\t</PointData>\n",fp);
		
		// Cell Data Begin
		fputs ("\t\t\t\t<CellData>\n",fp);	
	  
		long byte_offset = 0;
		fprintf(fp,"\t\t\t\t\t<DataArray type=\"UInt8\" Name=\"Level\" NumberOfTuples=\"%li\" format=\"appended\" RangeMin=\"0\" RangeMax=\"%i\" offset=\"%li\"/>\n", vertices, grid->maxdepth, byte_offset);
		
		byte_offset += vertices * sizeof(u_int8_t) + sizeof(u_int32_t);
		for (scalar s in list) {			
			fprintf(fp,"\t\t\t\t\t<DataArray type=\"Float32\" Name=\"%s\" NumberOfTuples=\"%li\" format=\"appended\" offset=\"%li\"/>\n", s.name, vertices, byte_offset);			
			byte_offset += vertices * sizeof(float_t) + sizeof(u_int32_t);
		}
		
		for (vector v in vlist) {
			char *vname = strtok(v.x.name, ".");
			fprintf (fp, "\t\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"%i\" Name=\"%s\" NumberOfTuples=\"%li\" format=\"appended\" offset=\"%li\"/>\n", dimension,vname, vertices, byte_offset);
			byte_offset += vertices * dimension * sizeof(float_t) + sizeof(u_int32_t);
		} 
		// Cell Data End
		fputs ("\t\t\t\t</CellData>\n",fp);

		// Tree End
		fputs ("\t\t\t</Tree>\n",fp);
		// File Tail
		fputs ("\t\t</Trees>\n",fp);
		fputs ("\t</HyperTreeGrid>\n",fp);
		fputs ("\t<AppendedData encoding=\"raw\">\n",fp);
		fputs ("_",fp);	  	 
	} // end pid() == 0
	
	cell_size=sizeof(u_int8_t);
	if(pid()==0){	  
		u_int32_t prepend_size = vertices * cell_size;
		fwrite(&prepend_size, sizeof(u_int32_t), 1, fp); 
		offset = ftell(fp);	  
	} 
	#if _MPI
	MPI_Barrier(MPI_COMM_WORLD);	
	MPI_Bcast(&offset, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
	#endif
   
	for(int lvl=0; lvl<=grid->maxdepth;++lvl){
		fseek (fp, offset + myOffsetPerLevel[2*lvl+0]*cell_size, SEEK_SET); 		
		foreach_level(lvl) {
			if (is_local(cell)) {
				fwrite(&lvl, cell_size, 1, fp);
			} else {
			}
		}
		offset += myOffsetPerLevel[2*lvl+1]*cell_size;		
	}
	fseek (fp, offset, SEEK_SET);
  
	for (scalar s in list) {
		cell_size=sizeof(float_t);
		if(pid()==0){		  
			u_int32_t prepend_size = vertices * cell_size;
			fwrite(&prepend_size, sizeof(u_int32_t), 1, fp); 		  
	    } 
		offset += sizeof(u_int32_t);
		for(int lvl=0; lvl<=grid->maxdepth;++lvl){
			// offset per level			
			fseek (fp, offset + myOffsetPerLevel[2*lvl+0]*cell_size, SEEK_SET); 
			foreach_level(lvl) {
				if (is_local(cell)) {
					float_t value[1] = {s[]};
					
					fwrite(value, cell_size, 1, fp);
				} else {
				}
			}
			offset += myOffsetPerLevel[2*lvl+1]*cell_size;
		}
		fseek (fp, offset, SEEK_SET);	  
	}
	
	for (vector v in vlist) {
		cell_size = dimension*sizeof(float_t);
		if(pid()==0){		  
			u_int32_t prepend_size = vertices * cell_size;
			fwrite(&prepend_size, sizeof(u_int32_t), 1, fp); 		  
		} 
		offset += sizeof(u_int32_t);
		for(int lvl=0; lvl<=grid->maxdepth;++lvl){
			// offset per level
			fseek (fp, offset + myOffsetPerLevel[2*lvl+0]*cell_size, SEEK_SET); 			
			foreach_level(lvl) {
				if (is_local(cell)) {
					#if dimension == 2
					float_t data[2] = {v.y[],v.x[]};
					#elif dimension == 3
					float_t data[3] = {v.z[],v.y[],v.x[]};
					#endif
					fwrite(data, cell_size, 1, fp);
				} else {
				}
			}
			offset += myOffsetPerLevel[2*lvl+1]*cell_size;
		}
		fseek (fp, offset, SEEK_SET);	  
	}  
	// Write Tail
	if (pid() == 0)	{		
		fprintf(fp, "ENDBINARY\n\t</AppendedData>\n</VTKFile>\n");
	}
  
  #if _MPI  
    MPI_Barrier(MPI_COMM_WORLD);  
    MPI_Win_free(&win); 
  #else 
    free(myOffsetPerLevel);
  #endif
  #if defined(_OPENMP)
    omp_set_num_threads(num_omp);
  #endif
}

/**
# *output_htg_data_v2_xml20(...)*
 * xml-version 2.0, for Paraview-nightly builds
*/
void output_htg_data_v2_xml20(scalar * list, vector * vlist, FILE * fp)
{
	#if defined(_OPENMP)
		int num_omp = omp_get_max_threads();
		omp_set_num_threads(1);
	#endif
	
	long verticesPerLevelLokal[grid->maxdepth+1];  
	memset(verticesPerLevelLokal, 0,(grid->maxdepth+1)*sizeof(long));
  
	for (int lvl = 0; lvl <= depth(); lvl++) {
		foreach_level(lvl) 
			if(is_local(cell))
				verticesPerLevelLokal[lvl]++;	  		
	}
  
	long offset = 0;

	/*
	*  myOffsetPerLevel[2*lvl+0] contains the level depending offset to write
	*  myOffsetPerLevel[2*lvl+1] contains the offset for the next level (total number of written data per level)  
	* 		this equals the number of Cells Per Level*/
	long * myOffsetPerLevel;      
	int window_elements = 2*(grid->maxdepth+1);  
    #if _MPI
	MPI_Win win;
	MPI_Win_allocate((MPI_Aint)window_elements*sizeof(long),\
	sizeof(long), MPI_INFO_NULL, MPI_COMM_WORLD, &myOffsetPerLevel, &win);				
	MPI_Win_fence(0, win);
	memset(myOffsetPerLevel, 0, window_elements*sizeof(long));
	MPI_Win_fence(0, win);

	for(int lvl = 0; lvl < grid->maxdepth+1; ++lvl){
		for(int pe = 0; pe < npe();++pe){
			// Offset for next level  in myOffsetPerLevel[1]
			MPI_Accumulate(&verticesPerLevelLokal[lvl], 1, MPI_LONG, pe, 2*lvl+1, 1, MPI_LONG, MPI_SUM, win);    
			if( pe > pid()){
				// Offset in myOffsetPerLevel[0]
				MPI_Accumulate(&verticesPerLevelLokal[lvl], 1, MPI_LONG, pe, 2*lvl+0, 1, MPI_LONG, MPI_SUM, win);
			}
		}
	}
	
	MPI_Win_fence(0, win);	
    #else
    myOffsetPerLevel = (long*)malloc(window_elements*sizeof(long));
    for(int lvl = 0; lvl < grid->maxdepth+1; ++lvl){
		myOffsetPerLevel[2*lvl+1]=verticesPerLevelLokal[lvl];
		myOffsetPerLevel[2*lvl+0]=0;
	}
    #endif
    
	long vertices = 0;
	for(int lvl=0; lvl<=grid->maxdepth;++lvl){
		vertices +=myOffsetPerLevel[2*lvl+1];
	}

	// Count Descriptor Bits
	long num_descbit = 0;
	for (int lvl = 0; lvl < grid->maxdepth; ++lvl) {
		num_descbit +=myOffsetPerLevel[2*lvl+1];
	}

	// File Header vtkHyperTreeGrid
	if (pid() == 0){	  
		fputs ("<VTKFile type=\"HyperTreeGrid\" version=\"2.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n",fp);
		#if dimension==2
		fprintf(fp, "\t<HyperTreeGrid BranchFactor=\"2\" TransposedRootIndexing=\"0\" Dimensions=\"%d %d %d\">\n", 2,2,1);
		#elif dimension==3
		fprintf(fp, "\t<HyperTreeGrid BranchFactor=\"2\" TransposedRootIndexing=\"0\" Dimensions=\"%d %d %d\">\n", 2,2,2);
		#endif
		fputs ("\t\t<Grid>\n",fp);
		fprintf(fp, "\t\t\t<DataArray type=\"Float32\" Name=\"XCoordinates\" NumberOfTuples=\"2\" format=\"ascii\" RangeMin=\"%g\" RangeMax=\"%g\">\n",X0, X0+L0);
		fprintf (fp, "\t\t\t%g %g",X0, X0+L0);
		fputs ("\n\t\t\t</DataArray>\n",fp);
		fprintf(fp, "\t\t\t<DataArray type=\"Float32\" Name=\"YCoordinates\" NumberOfTuples=\"2\" format=\"ascii\" RangeMin=\"%g\" RangeMax=\"%g\">\n",Y0, Y0+L0);
		fprintf (fp, "\t\t\t%g %g",Y0, Y0+L0);
		fputs ("\n\t\t\t</DataArray>\n",fp);	  
		fprintf(fp, "\t\t\t<DataArray type=\"Float32\" Name=\"ZCoordinates\" NumberOfTuples=\"2\" format=\"ascii\" RangeMin=\"%g\" RangeMax=\"%g\">\n",Z0, Z0+L0);
		#if dimension==2
		fprintf (fp, "\t\t\t%g %g", 0., 0.);	  
		#elif dimension==3
		fprintf (fp, "\t\t\t%g %g",Z0, Z0+L0);	  
		#endif
		fputs ("\n\t\t\t</DataArray>\n",fp);	  
		fputs ("\t\t</Grid>\n",fp);
		// Trees Begin
		fputs ("\t\t<Trees>\n",fp);

		// Array Descriptor Bits
		fprintf (fp,"\t\t\t<DataArray type=\"Bit\" Name=\"Descriptors\" NumberOfTuples=\"%li\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"1\">\n", num_descbit);
		fputs("\t\t\t", fp);
		offset = ftell(fp);	  
	} // pid()==0
	#if _MPI
	MPI_Bcast(&offset, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
	#endif
	   /* Print Bitmask per Process per Level */
	int cell_size=2*sizeof(char);    
	for(int lvl=0; lvl<grid->maxdepth;++lvl){ // Bitmask is "0" at grid->maxdepth. Therefore lvl < grid->maxdepth
		fseek (fp, offset + myOffsetPerLevel[2*lvl+0]*cell_size, SEEK_SET); 
		foreach_level(lvl) {
			if (is_local(cell)) {
				if(is_leaf(cell)){
					fputs("0 ",fp);
				} else {
					fputs("1 ",fp);
				}			
			} else {
				continue;
			}
		}				
		offset += myOffsetPerLevel[2*lvl+1]*cell_size;
	}
	fseek (fp, offset, SEEK_SET);

     /* Print VerticesByLevel */         		 
	if(pid()==0){
		fputs ("\n\t\t\t</DataArray>\n",fp);
		// Array NbVerticesByLevel
		fprintf(fp,"\t\t\t<DataArray type=\"Int64\" Name=\"NumberOfVerticesPerDepth\" NumberOfTuples=\"%d\" format=\"ascii\" RangeMin=\"1\" RangeMax=\"%li\" >\n" ,grid->maxdepth + 1, myOffsetPerLevel[2*grid->maxdepth + 1]);
		fputs("\t\t\t", fp);

		for (int lvl = 0; lvl <= grid->maxdepth; lvl++) {
			fprintf(fp,"%li ", myOffsetPerLevel[2*lvl+1]);
			}
		fputs ("\n\t\t\t</DataArray>\n",fp);

		fprintf(fp,"\t\t\t<DataArray type=\"Int64\" Name=\"TreeIds\" NumberOfTuples=\"1\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\" >\n");
		fputs("\t\t\t0", fp);
		fputs ("\n\t\t\t</DataArray>\n",fp);
		fprintf(fp,"\t\t\t<DataArray type=\"UInt32\" Name=\"DepthPerTree\" NumberOfTuples=\"1\" format=\"ascii\" RangeMin=\"%i\" RangeMax=\"%i\" >\n" ,grid->maxdepth + 1,grid->maxdepth + 1);
		fputs("\t\t\t", fp);
		fprintf(fp, "%i", grid->maxdepth + 1);
		fputs ("\n\t\t\t</DataArray>\n",fp);
		fputs ("\t\t</Trees>\n",fp);
		// Cell Data Begin
		fputs ("\t\t<CellData>\n",fp);
	  
		long byte_offset = 0;
		fprintf(fp,"\t\t\t<DataArray type=\"UInt8\" Name=\"Level\" NumberOfTuples=\"%li\" format=\"appended\" RangeMin=\"0\" RangeMax=\"%i\" offset=\"%li\"/>\n", vertices, grid->maxdepth, byte_offset);
		byte_offset += vertices * sizeof(u_int8_t) + sizeof(u_int32_t);
		for (scalar s in list) {

			fprintf(fp,"\t\t\t<DataArray type=\"Float32\" Name=\"%s\" NumberOfTuples=\"%li\" format=\"appended\" offset=\"%li\"/>\n", s.name, vertices, byte_offset);
			byte_offset += vertices * sizeof(float_t) + sizeof(u_int32_t);
		}
		// Vectors
		for (vector v in vlist) {
			char *vname = strtok(v.x.name, ".");
			fprintf (fp, "\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"%i\" Name=\"%s\" NumberOfTuples=\"%li\" format=\"appended\" RangeMin=\"\" RangeMax=\"\" offset=\"%li\"/>\n", dimension,vname, vertices, byte_offset);
			byte_offset += vertices * dimension * sizeof(float_t) + sizeof(u_int32_t);
		} 
		// Cell Data End	  
		fputs ("\t\t</CellData>\n",fp);

		// File Tail
		
		fputs ("\t</HyperTreeGrid>\n",fp);
		fputs ("\t<AppendedData encoding=\"raw\">\n",fp);
		fputs ("_",fp);	  	 
	} // end pid() == 0
	
	cell_size=sizeof(u_int8_t);
	if(pid()==0){	  
		u_int32_t prepend_size = vertices * cell_size;
		fwrite(&prepend_size, sizeof(u_int32_t), 1, fp); 
		offset = ftell(fp);	  
	} 
	#if _MPI
	MPI_Bcast(&offset, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
	#endif
   
	for(int lvl=0; lvl<=grid->maxdepth;++lvl){
		fseek (fp, offset + myOffsetPerLevel[2*lvl+0]*cell_size, SEEK_SET); 		
		foreach_level(lvl) {
			if (is_local(cell)) {
				fwrite(&lvl, cell_size, 1, fp);
			} else {
			}
		}
		offset += myOffsetPerLevel[2*lvl+1]*cell_size;		
	}
	fseek (fp, offset, SEEK_SET);
  
	for (scalar s in list) {
		cell_size=sizeof(float_t);
		if(pid()==0){		  
			u_int32_t prepend_size = vertices * cell_size;
			fwrite(&prepend_size, sizeof(u_int32_t), 1, fp); 		  
	    } 
		offset += sizeof(u_int32_t);
		for(int lvl=0; lvl<=grid->maxdepth;++lvl){
			// offset per level			
			fseek (fp, offset + myOffsetPerLevel[2*lvl+0]*cell_size, SEEK_SET); 
			foreach_level(lvl) {
				if (is_local(cell)) {
					float_t value[1] = {s[]};
					fwrite(value, cell_size, 1, fp);
				} else {
				}
			}
			offset += myOffsetPerLevel[2*lvl+1]*cell_size;
		}
		fseek (fp, offset, SEEK_SET);	  
	}
	
	for (vector v in vlist) {
		cell_size = dimension*sizeof(float_t);
		if(pid()==0){		  
			u_int32_t prepend_size = vertices * cell_size;
			fwrite(&prepend_size, sizeof(u_int32_t), 1, fp); 		  
		} 
		offset += sizeof(u_int32_t);
		for(int lvl=0; lvl<=grid->maxdepth;++lvl){
			// offset per level
			fseek (fp, offset + myOffsetPerLevel[2*lvl+0]*cell_size, SEEK_SET); 			
			foreach_level(lvl) {
				if (is_local(cell)) {
					#if dimension == 2
					float_t data[2] = {v.y[],v.x[]};
					#elif dimension == 3
					float_t data[3] = {v.z[],v.y[],v.x[]};
					#endif
					fwrite(data, cell_size, 1, fp);
				} else {
				}
			}
			offset += myOffsetPerLevel[2*lvl+1]*cell_size;
		}
		fseek (fp, offset, SEEK_SET);	  
	}  
	// Write Tail
	if (pid() == 0)	{		
		fprintf(fp, "ENDBINARY\n\t</AppendedData>\n</VTKFile>\n");
	}
  
  #if _MPI  
    MPI_Barrier(MPI_COMM_WORLD);  
    MPI_Win_free(&win);    
  #else 
    free(myOffsetPerLevel);
  #endif
  #if defined(_OPENMP)
    omp_set_num_threads(num_omp);
  #endif
}