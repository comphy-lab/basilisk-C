/**
 Export  VTK-(P)-ImageData
 
 * output_vti()  
 * output_pvti() 
 * output_vti_data_v1()
 
## Usage
        event log(i++)
        {
		char path[]="vti/";
		char prefix[80];
		sprintf(prefix, "data_%06d", i);  
		output_vti((scalar *) {cs, d,p,mpi_color},(vector *){u}, path, prefix, i, t);
	}
 
 */
#include "output_pvd.h"
void output_vti(scalar * list, vector * vlist, char* path, char* prefix, int i, double t);
void output_pvti(scalar * list, vector * vlist, FILE * fp, char* subname);
void output_pvd(char* name, double t, FILE* fp,bool firstTimeWritten);
void output_vti_data_v1(scalar * list, vector * vlist, FILE * fp);

void output_vti(scalar * list, vector * vlist, char* path, char* prefix, int i, double t){
	
	FILE * fp ;

	char vti_name[80];  	  
	sprintf(vti_name, "%s%s_n%04d.vti", path, prefix, pid());  

	fp = fopen(vti_name, "w");
	if(fp == NULL){
		printf("output_vti_pv_590_mpi.h : %s could not be opend\n Does the Folder exist?\n", vti_name);
		exit(1);
	}  	  
	output_vti_data_v1((scalar *) list,(vector *)vlist,fp);
	fclose(fp);

	if(0 == pid()){		  
		char pvd_name[80];  
		char pvti_name[80];  

		sprintf(pvti_name, "%s%s.pvti", path, prefix);

		fp = fopen(pvti_name, "w");		  
		output_pvti((scalar *) list,(vector *)vlist,fp, prefix);
		fclose(fp);

		bool firstTimeWritten = false;
		sprintf(pvd_name, "imageData.pvd");	  
		fp = fopen(pvd_name, "r+");
		if( ((int)i == 0) || (fp == NULL) ){
			fp = fopen(pvd_name,"w");
			firstTimeWritten = true;
			printf("firstTimeWritten\n");
		}
		  
		output_pvd(pvti_name, t, fp, firstTimeWritten);
		fclose(fp);
	  }
	#if _MPI  
	MPI_Barrier(MPI_COMM_WORLD);    
	#endif
}

void output_pvti(scalar * list, vector * vlist, FILE * fp, char* subname)
{
	#if defined(_OPENMP)
		int num_omp = omp_get_max_threads();
		omp_set_num_threads(1);
	#endif

	long byte_offset=0;
	long vertices= grid->tn;
	
	fprintf(fp,"<VTKFile type=\"PImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
	#if MULTIGRID_MPI 
		int a3 = (1<<depth())*mpi_dims[0];
		int a2 = (1<<depth())*mpi_dims[1];
		#if dimension > 2
			int a1 = (1<<depth())*mpi_dims[2];
		#else
			int a1 = 0;
		#endif
		double spac3 = (double)L0/((double)(1<<depth()))/(double)mpi_dims[0];		
		double spac2 = spac3;
		#if dimension > 2		
			double spac1 = spac2;
		#else
			double spac1 = 0.;
		#endif
	
	#else // NO_MULTIGRID_MPI				
		int a3 = (1<<depth());
		int a2 = a3;
		#if dimension == 3
			int a1 = a2;
		#else
			int a1 = 0;
		#endif
		double spac3 = (double)L0/((double)(1<<depth()));
		double spac2 = spac3;
		#if dimension == 3
			double spac1 = spac2;
		#else
			double spac1 = 0;
		#endif

	#endif
	fprintf(fp,"\t<PImageData WholeExtent=\"%i %i %i %i %i %i\" GhostLevel=\"0\" Origin=\"%g %g %g\" Spacing=\"%g %g %g\" Direction=\"1 0 0 0 1 0 0 0 1\">\n",\
							0,a1,0,a2,0,a3,X0,Y0,Z0,spac1,spac2,spac3);
	fprintf(fp,"\t\t<PCellData>\n");
	
	for(scalar s in list) {		
		fprintf(fp,"\t\t\t<PDataArray type=\"Float32\" Name=\"%s\"/>\n", s.name);		
		byte_offset += vertices * sizeof(float_t) + sizeof(u_int64_t);
	}
	for (vector v in vlist) {
		char *vname = strtok(v.x.name, ".");
		fprintf (fp, "\t\t\t<PDataArray type=\"Float32\" NumberOfComponents=\"%i\" Name=\"%s\"/>\n", dimension,vname);
		byte_offset += vertices * dimension * sizeof(float_t) + sizeof(u_int64_t);
	} 
	fprintf(fp,"\t\t</PCellData>\n");
	for(int pe=0; pe<npe(); ++pe){
		#if MULTIGRID_MPI
		#if dimension > 2
		int nz = pe % (mpi_dims[2]); 
		int ny = (pe/mpi_dims[2]) % (mpi_dims[1]);
		int nx = pe/(mpi_dims[2]*mpi_dims[1]);
		#else
		int nz = 0;
		int ny = (pe) % (mpi_dims[1]);
		int nx = pe/(mpi_dims[1]);
		#endif
		int nxoffset = ((1 << depth())*nx);
		int nyoffset = ((1 << depth())*ny);
		int nzoffset = ((1 << depth())*nz);
		#else
		int nxoffset = 0;
		int nyoffset = 0;
		int nzoffset = 0;
		#endif
		
		#if dimension > 2
		fprintf(fp,"\t\t<Piece Extent=\"%i %i %i %i %i %i\" Source=\"%s_n%04d.vti\"/>\n",\
					nzoffset,nzoffset+(1 << depth()),nyoffset,nyoffset+(1 << depth()),nxoffset,nxoffset+(1 << depth()),\
					subname, pe); 
		#else 
		fprintf(fp,"\t\t<Piece Extent=\"%i %i %i %i %i %i\" Source=\"%s_n%04d.vti\"/>\n",\
					0,0,nyoffset,nyoffset+(1 << depth()),nxoffset,nxoffset+(1 << depth()),\
					subname, pe); 
		#endif
	}
	fprintf(fp,"\t</PImageData>\n");
	fprintf(fp,"</VTKFile>\n");
	fflush(fp);		
  #if defined(_OPENMP)
    omp_set_num_threads(num_omp);
  #endif	

}

void output_vti_data_v1(scalar * list, vector * vlist, FILE * fp)
{
	#if defined(_OPENMP)
		int num_omp = omp_get_max_threads();
		omp_set_num_threads(1);
	#endif

	long verticesLokal= grid->n;  	 
	//~ long vertices= grid->tn;

	#if MULTIGRID_MPI 	
	#if dimension > 2
	int nz = pid() % (mpi_dims[2]); 
	int ny = (pid()/mpi_dims[2]) % (mpi_dims[1]);
	int nx = pid()/(mpi_dims[2]*mpi_dims[1]);
	#else
	int nz = 0;
	int ny = (pid()) % (mpi_dims[1]);
	int nx = (pid()) / (mpi_dims[1]);
	#endif
	int nxoffset = ((1 << depth() )*nx);
	int nyoffset = ((1 << depth() )*ny);
	int nzoffset = ((1 << depth() )*nz);
		
	double spac3 = (double)L0/((double)(1<<depth()))/(double)mpi_dims[0];
	double spac2 = spac3;
	#if dimension > 2
	double spac1 = spac2;
	#else
	double spac1 = 0.;
	#endif
	#else // NO_MULTIGRID_MPI	
	double spac3 = (double)L0/((double)(1<depth()));
	double spac2 = spac3;
	double spac1 = spac2;
	int nxoffset = 0;
	int nyoffset = 0;
	int nzoffset = 0;
	#endif

	fprintf(fp,"<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");	
	
	#if dimension > 2
	
	fprintf(fp,"\t<ImageData WholeExtent=\"%i %i %i %i %i %i\" Origin=\"%g %g %g\" Spacing=\"%g %g %g\" Direction=\"1 0 0 0 1 0 0 0 1\">\n",\
									nzoffset,nzoffset+(1 << depth()),nyoffset,nyoffset+(1 << depth()),nxoffset,nxoffset+(1 << depth()),\
									X0,Y0,Z0,spac1,spac2,spac3);
	#else
	
	fprintf(fp,"\t<ImageData WholeExtent=\"%i %i %i %i %i %i\" Origin=\"%g %g %g\" Spacing=\"%g %g %g\" Direction=\"1 0 0 0 1 0 0 0 1\">\n",\
									0,0,nyoffset,nyoffset+(1 << depth()),nxoffset,nxoffset+(1 << depth()),\
									X0,Y0,Z0,spac1,spac2,spac3);
	#endif
	
	long byte_offset = 0;
	#if dimension > 2		
	fprintf(fp,"\t\t<Piece Extent=\"%i %i %i %i %i %i\">\n",nzoffset,nzoffset+(1 << depth()),nyoffset,nyoffset+(1 << depth()),nxoffset,nxoffset+(1 << depth())); // Wie richtige indizieren?
	#else
	fprintf(fp,"\t\t<Piece Extent=\"%i %i %i %i %i %i\">\n",0,0,nyoffset,nyoffset+(1 << depth()),nxoffset,nxoffset+(1 << depth())); // Wie richtige indizieren?
	#endif
	fprintf(fp,"\t\t\t<PointData>\n");				
	fprintf(fp,"\t\t\t</PointData>\n");				
	fprintf(fp,"\t\t\t<CellData>\n");				
	for (scalar s in list) {		
		fprintf(fp,"\t\t\t\t<DataArray type=\"Float32\" Name=\"%s\"  format=\"appended\" offset=\"%li\"/>\n", s.name, byte_offset);		
		byte_offset += verticesLokal * sizeof(float_t) + sizeof(u_int64_t);
	}
	for (vector v in vlist) {
		char *vname = strtok(v.x.name, ".");
		fprintf (fp, "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"%i\" Name=\"%s\" format=\"appended\" offset=\"%li\"/>\n", dimension,vname, byte_offset);
		byte_offset += verticesLokal * dimension * sizeof(float_t) + sizeof(u_int64_t);
	} 
	fprintf(fp,"\t\t\t</CellData>\n");
	fprintf(fp,"\t\t</Piece>\n");
	fprintf(fp,"\t</ImageData>\n");
	fprintf(fp,"\t<AppendedData encoding=\"raw\">\n_");
	
	int cell_size;
	for (scalar s in list) {
		cell_size=sizeof(float_t);	
		u_int64_t prepend_size = verticesLokal * cell_size;
		fwrite(&prepend_size, sizeof(u_int64_t), 1, fp); 		  	
		foreach() {
			float_t value[1] = {s[]};
			fwrite(value, cell_size, 1, fp);
		}
	}
	for (vector v in vlist) {
		cell_size = dimension*sizeof(float_t);	
		u_int64_t prepend_size = verticesLokal * cell_size;
		fwrite(&prepend_size, sizeof(u_int64_t), 1, fp); 		  	
		foreach() {
			#if dimension == 2
			float_t data[2] = {v.y[],v.x[]};
			#elif dimension == 3
			float_t data[3] = {v.z[],v.y[],v.x[]};
			#endif
			fwrite(data, cell_size, 1, fp);					
		}	
	}  
	// Write Tail	
	fprintf(fp,"ENDBINARY\n\t</AppendedData>\n</VTKFile>\n");
	
  #if defined(_OPENMP)
    omp_set_num_threads(num_omp);
  #endif
}


