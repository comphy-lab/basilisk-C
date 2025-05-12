/**
  These functions do the same as "sandbox/Antoonvh/readxyz.h", except it 
  is not limited to one node.
   
  These functions allow to read
  
  	(a) binary x-y-z field into Basilisk
  		read_xyz_float_v1()
  		read_xyz_float_v2()
  	
  
  	(b) reading of data,that was created by Basilisk's "foreach()" operator
  		read_foreach_float_v1()
  		read_foreach_float_v2()
  
  
 */ 

void read_xyz_float_v1( FILE * fp, vertex scalar cs, int LEVEL)
{
	unsigned long long nvox = (1 << LEVEL);		
	float* data = (float*) malloc(sizeof(float)*nvox*nvox*nvox);	
	fread(data, sizeof(float), nvox*nvox*nvox, fp);	//< C order!!!	
	
	#if MULTIGRID_MPI 
	int nz = pid() % (mpi_dims[2]); 
	int ny = (pid()/mpi_dims[2]) % (mpi_dims[1]);
	int nx = pid()/(mpi_dims[2]*mpi_dims[1]);
	int nxoffset = ((1 << depth())*nx);
	int nyoffset = ((1 << depth())*ny);
	int nzoffset = ((1 << depth())*nz);
	#else // Non MG-MPI, no offset 
	int nxoffset = 0;
	int nyoffset = 0;
	int nzoffset = 0;
	#endif
	
	long long index_c;
	int o = -2; // alternative: "-BGHOSTS - 1;"	
	
	
	foreach(){		
		index_c = ((nxoffset+point.i+o)*nvox + (nyoffset+point.j+o))*nvox + (nzoffset+point.k+o);
		cs[] = (double)data[index_c];	//< row major, C

		//~ index_fortran = ((point.k+o)*nvox + (point.j+o))*nvox + (point.i+o);		
		//~ cs[] = (double)data[index_fortran];	//< col major, FORTRAN
  }  
  free(data);
}

void read_xyz_float_v2( FILE * fp, vertex scalar cs, int LEVEL)
{
	unsigned long long nvox = (1 << LEVEL);		
	MPI_Win win;
	float* snd_buf;
	if (pid()==0){		
		MPI_Win_allocate((MPI_Aint)nvox*nvox*nvox*sizeof(float),\
		sizeof(float), MPI_INFO_NULL, MPI_COMM_WORLD, &snd_buf, &win);				
	} else {
		MPI_Win_allocate((MPI_Aint)0,\
		sizeof(float), MPI_INFO_NULL, MPI_COMM_WORLD, &snd_buf, &win);
	}
	MPI_Win_fence(0, win);
	MPI_Barrier (MPI_COMM_WORLD);
    if (pid() == 0){ 
		MPI_Win_lock (MPI_LOCK_EXCLUSIVE, 0, MPI_MODE_NOCHECK, win);
		fread(snd_buf, sizeof(float), nvox*nvox*nvox, fp);			
		MPI_Win_unlock (0, win);
    }
    MPI_Barrier (MPI_COMM_WORLD);
    
    #if MULTIGRID_MPI 
	int nz = pid() % (mpi_dims[2]); 
	int ny = (pid()/mpi_dims[2]) % (mpi_dims[1]);
	int nx = pid()/(mpi_dims[2]*mpi_dims[1]);
	int nxoffset = ((1 << depth())*nx);
	int nyoffset = ((1 << depth())*ny);
	int nzoffset = ((1 << depth())*nz);
	#else // Non MG-MPI, no offset 
	int nxoffset = 0;
	int nyoffset = 0;
	int nzoffset = 0;
	#endif
    
	long long index_c;
	int o = -2;
	
	float temp;
	
	MPI_Win_fence(0, win);
	foreach(){
		if (is_local(cell)){
			index_c = ((nxoffset+point.i+o)*nvox + (nyoffset+point.j+o))*nvox + (nzoffset+point.k+o);		//< row major, C
			MPI_Get(&temp, 1, MPI_FLOAT, 0, (MPI_Aint)index_c, 1, MPI_FLOAT, win);
			MPI_Win_fence(0, win);
			cs[] = (double)temp;
		}
	} 

	MPI_Win_fence(0, win);
	MPI_Win_free(&win);
	snd_buf=NULL;  
}

void read_xyz_float_v3( MPI_File * fh, vertex scalar cs, int LEVEL)
{
	MPI_Status status;
	unsigned long long nvox = (1 << LEVEL);		
	        
    #if MULTIGRID_MPI 
	int nz = pid() % (mpi_dims[2]); 
	int ny = (pid()/mpi_dims[2]) % (mpi_dims[1]);
	int nx = pid()/(mpi_dims[2]*mpi_dims[1]);
	int nxoffset = ((1 << depth())*nx);
	int nyoffset = ((1 << depth())*ny);
	int nzoffset = ((1 << depth())*nz);
	#else // Non MG-MPI, no offset 
	int nxoffset = 0;
	int nyoffset = 0;
	int nzoffset = 0;
	#endif
    
	long long index_c;
	int o = -2;
	
	float temp;
	
	foreach(){
		if (is_local(cell)){
			index_c = ((nxoffset+point.i+o)*nvox + (nyoffset+point.j+o))*nvox + (nzoffset+point.k+o);		//< row major, C			
			MPI_File_read_at(*fh, sizeof(float)*index_c, &temp, 1, MPI_FLOAT, &status);									
			cs[] = (double)temp;
		}
	} 

}

void read_foreach_float_v1( FILE * fp, scalar d, int LEVEL)
{
	unsigned long long nvox = (1 << LEVEL);		
	MPI_Win win;
	float* snd_buf;
	if (pid()==0){		
		MPI_Win_allocate((MPI_Aint)nvox*nvox*nvox*sizeof(float),\
		sizeof(float), MPI_INFO_NULL, MPI_COMM_WORLD, &snd_buf, &win);				
	} else {
		MPI_Win_allocate((MPI_Aint)0,\
		sizeof(float), MPI_INFO_NULL, MPI_COMM_WORLD, &snd_buf, &win);
	}
	MPI_Win_fence(0, win);
	MPI_Barrier (MPI_COMM_WORLD);
    if (pid() == 0){ 
		MPI_Win_lock (MPI_LOCK_EXCLUSIVE, 0, MPI_MODE_NOCHECK, win);
		fread(snd_buf, sizeof(float), nvox*nvox*nvox, fp);			
		MPI_Win_unlock (0, win);
    }
    MPI_Barrier (MPI_COMM_WORLD);
	float temp;
		 
	unsigned long n=(unsigned long)grid->n;
	
	unsigned long offsetPerThread[npe()];
	MPI_Allgather(&n, 1, MPI_UNSIGNED_LONG, &offsetPerThread, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    unsigned long offset=0;
    unsigned long ntotal=0;
    for (int i = 0; i < npe(); i++) {
		ntotal += offsetPerThread[i];
		if (i < pid())
			offset += offsetPerThread[i];
    }
	
	MPI_Win_fence(0, win);
	int i = 0, j=0;
	//~ foreach_vertex(){ // Problem bei foreach_offset: Bei Parallelen Programmen gibt es (an den Kanten) mehr vertices! -> nicht mapbar
	foreach(){				
		MPI_Get(&temp, 1, MPI_FLOAT, 0, (MPI_Aint)offset, 1, MPI_FLOAT, win);			
		MPI_Win_fence(0, win);
		d[] = (double)temp;
		offset++;
	}  
;
	MPI_Win_fence(0, win);
	MPI_Win_free(&win);
	snd_buf=NULL;

}

// Using one Sided "MPI_Accumulate" to pass offsets to processes
void read_foreach_float_v2( FILE * fp, scalar d, int LEVEL)
{
	unsigned long long nvox = (1 << LEVEL);		
	
	MPI_Win win1;
	float* snd_buf;
	if (pid()==0){		
		MPI_Win_allocate((MPI_Aint)nvox*nvox*nvox*sizeof(float),\
		sizeof(float), MPI_INFO_NULL, MPI_COMM_WORLD, &snd_buf, &win1);				
	} else {
		MPI_Win_allocate((MPI_Aint)0,\
		sizeof(float), MPI_INFO_NULL, MPI_COMM_WORLD, &snd_buf, &win1);
	}
	MPI_Win_fence(0, win1);
	MPI_Barrier (MPI_COMM_WORLD);
    if (pid() == 0){ 
		MPI_Win_lock (MPI_LOCK_EXCLUSIVE, 0, MPI_MODE_NOCHECK, win1);
		fread(snd_buf, sizeof(float), nvox*nvox*nvox, fp);			
		MPI_Win_unlock (0, win1);
    }
    MPI_Barrier (MPI_COMM_WORLD);
	float temp;
		 
	unsigned long n=(unsigned long)grid->n;

    MPI_Win win2;
    unsigned long * myOffset;
    MPI_Win_allocate((MPI_Aint)2*sizeof(unsigned long),\
		sizeof(unsigned long), MPI_INFO_NULL, MPI_COMM_WORLD, &myOffset, &win2);				
	MPI_Win_fence(0, win2);
	memset(myOffset, 0, 2*sizeof(unsigned long));
	MPI_Win_fence(0, win2);
	
	for(int i = 0; i<npe();++i){
		// Total Number of leaves in myOffset[1]
		MPI_Accumulate(&n, 1, MPI_UNSIGNED_LONG, i, 1, 1, MPI_UNSIGNED_LONG, MPI_SUM, win2);
		if( i > pid()){
			// Offset in myOffset[0]
			MPI_Accumulate(&n, 1, MPI_UNSIGNED_LONG, i, 0, 1, MPI_UNSIGNED_LONG, MPI_SUM, win2);
		}
	}
	MPI_Win_fence(0, win2);	
	printf("pid:%i, Offset:%lu\n",pid(),myOffset[0]);
	MPI_Win_fence(0, win1);	
	int i = 0;	
	foreach(){				
		MPI_Get(&temp, 1, MPI_FLOAT, 0, (MPI_Aint)myOffset[0], 1, MPI_FLOAT, win1);			
		MPI_Win_fence(0, win1);	
		d[] = (double)temp;
		myOffset[0]++;
	}  	
	MPI_Win_fence(0, win1);
	MPI_Win_free(&win1);
	
	MPI_Win_fence(0, win2);	
	MPI_Win_free(&win2);
	snd_buf=NULL;

}

void read_foreach_float_v3(MPI_File * fh, scalar d, int LEVEL)
{
	MPI_Status status;
	unsigned long long nvox = (1 << LEVEL);		

	float temp;
		 
	unsigned long n=(unsigned long)grid->n;

    MPI_Win win;
    unsigned long * myOffset;
    MPI_Win_allocate((MPI_Aint)2*sizeof(unsigned long),\
		sizeof(unsigned long), MPI_INFO_NULL, MPI_COMM_WORLD, &myOffset, &win);				
	MPI_Win_fence(0, win);
	memset(myOffset, 0, 2*sizeof(unsigned long));
	MPI_Win_fence(0, win);
	
	for(int i = 0; i<npe();++i){
		// Total Number of leaves in myOffset[1]
		MPI_Accumulate(&n, 1, MPI_UNSIGNED_LONG, i, 1, 1, MPI_UNSIGNED_LONG, MPI_SUM, win);
		if( i > pid()){
			// Offset in myOffset[0]
			MPI_Accumulate(&n, 1, MPI_UNSIGNED_LONG, i, 0, 1, MPI_UNSIGNED_LONG, MPI_SUM, win);
		}
	}
	MPI_Win_fence(0, win);	
	printf("pid:%i, Offset:%lu\n",pid(),myOffset[0]);
	
	int i = 0;	
	foreach(){			
		MPI_File_read_at(*fh, sizeof(float)*myOffset[0]++, &temp, 1, MPI_FLOAT, &status);
		d[] = (double)temp;
	}  	
	
	MPI_Win_fence(0, win);	
	MPI_Win_free(&win);	
}


