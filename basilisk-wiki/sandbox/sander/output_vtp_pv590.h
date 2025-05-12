/**
This headerfile enables the export of particle data in the Paraview format "*.vtp".

## Example Usage
        ...
        #include "tracer-particles.h"
        #include "stokes-particles.h" 
        #include "output_vtp_pv590.h"
	
        Particles flow, heavy;
	int np = 100; // number particle
	double ti = 0.5; // time injection
        ...
        event init(i=0){
            ...
            new_tracer_particles (0);
	    new_inertial_particles (0);
        }
        
        event init_particle(t=ti){
            ...
            flow  = new_tracer_particles(np);
	    heavy = new_inertial_particles(np);
            ...
        }
        
	if(t>=ti){
	    char path[]="vtp/";
	    char prefix[80];
	    sprintf(prefix, "flow_%06d", i); // underscore in name is necessary for name splitting (see "strtok(...))"! 
	    output_vtp(flow, false, path, prefix, i, t, ti);
	    sprintf(prefix, "heavy_%06d", i);  
	    output_vtp(heavy, true, path, prefix, i, t, ti);
	}

# VTP-Export
*/

#include "output_pvd.h"
void output_vtp_v1(Particles P, bool particleChar,FILE * fp);
void output_vtp (Particles P, bool particleChar,char* path, char* prefix, int i, double t,double ti);

void output_vtp (Particles P, bool particleChar,char* path, char* prefix, int i, double t,double ti){
	FILE * fp ;

	char vtp_name[80];  	  
	sprintf(vtp_name, "%s%s.vtp", path, prefix);  

	fp = fopen(vtp_name, "w");
	if(fp == NULL){
		printf("output_vtp_pv_590.h : %s could not be opend\n Does the Folder exist?\n", vtp_name);
		exit(1);
	}  
	output_vtp_v1(P, particleChar, fp);
	fclose(fp);

	if(pid()==0){
		bool firstTimeWritten = false;
		char pvd_name[80];
		char *partname = strtok(prefix, "_");	  
		sprintf(pvd_name,"%s_particles.pvd",partname);	  
		fp = fopen(pvd_name, "r+");
		if( (i == 0) || (t == ti) ||  fp == NULL ){
			fp = fopen(pvd_name,"w");
			firstTimeWritten = true;
		}
		 
		output_pvd(vtp_name, t, fp, firstTimeWritten);
		fclose(fp);
	}
	#if _MPI  
	MPI_Barrier(MPI_COMM_WORLD);    
	#endif
	
	
}

void output_vtp_v1(Particles P, bool particleChar,FILE * fp)
{
  #if defined(_OPENMP)
    int num_omp = omp_get_max_threads();
    omp_set_num_threads(1);
  #endif

  // number of local particles
  long unsigned int nparts_all = pn[P];
  long unsigned int parts_off = 0;
#if _MPI
   
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Win win;
  unsigned long * myOffset;
  MPI_Win_allocate((MPI_Aint)2*sizeof(unsigned long),\
	sizeof(unsigned long), MPI_INFO_NULL, MPI_COMM_WORLD, &myOffset, &win);				
  MPI_Win_fence(0, win);
  memset(myOffset, 0, 2*sizeof(unsigned long));
  MPI_Win_fence(0, win);
  
  for(int i = 0; i<npe();++i){
  // Total Number  in myOffset[1]
    MPI_Accumulate(&pn[P], 1, MPI_UNSIGNED_LONG, i, 1, 1, MPI_UNSIGNED_LONG, MPI_SUM, win);
    
    if( i > pid()){
	  // Offset in myOffset[0]
	  MPI_Accumulate(&pn[P], 1, MPI_UNSIGNED_LONG, i, 0, 1, MPI_UNSIGNED_LONG, MPI_SUM, win);
    }
    
  }
  MPI_Win_fence(0, win);	
  parts_off= myOffset[0];
  nparts_all = myOffset[1];  
  MPI_Win_free(&win);

#endif

  // Write Header
  if (pid() == 0)
  {
    unsigned long offset = 0;
    // File vtkPolyData
    fprintf(fp, "<?xml version=\"1.0\"?>\n<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n  <PolyData>\n    <Piece NumberOfPoints=\"%lu\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\" >\n      <Verts>\n      </Verts>\n      <Lines>\n      </Lines>\n      <Strips>\n      </Strips>\n      <Polys>\n      </Polys>\n      <CellData>\n      </CellData>\n      <Points>\n        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"appended\" RangeMin=\"\" RangeMax=\"\" offset=\"%lu\" />\n      </Points>\n      <PointData>\n", nparts_all, offset);

    offset += nparts_all*3*sizeof(float) + 8;
    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"3\" format=\"appended\" RangeMin=\"\" RangeMax=\"\" offset=\"%lu\" />\n", "u", offset);
    offset += nparts_all*3*sizeof(float) + 8;
    
    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"1\" format=\"appended\" RangeMin=\"\" RangeMax=\"\" offset=\"%lu\" />\n", "id", offset);    
    offset += nparts_all*1*sizeof(float) + 8;
    
    if (particleChar == true){
		fprintf(fp, "        <DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"1\" format=\"appended\" RangeMin=\"\" RangeMax=\"\" offset=\"%lu\" />\n", "density", offset);
		offset += nparts_all*1*sizeof(float) + 8;
		fprintf(fp, "        <DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"1\" format=\"appended\" RangeMin=\"\" RangeMax=\"\" offset=\"%lu\" />\n", "radius", offset);
		offset += nparts_all*1*sizeof(float) + 8;
		fprintf(fp, "        <DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"1\" format=\"appended\" RangeMin=\"\" RangeMax=\"\" offset=\"%lu\" />\n", "tau", offset);
		offset += nparts_all*1*sizeof(float) + 8;
	}
    // Appended Data
    fprintf(fp, "      </PointData>\n    </Piece>\n  </PolyData>\n  <AppendedData encoding=\"raw\">\n    _");
  }

  // Write Points
  unsigned long header_off = 0;
  u_int64_t data_size = nparts_all*3*sizeof(float);
  if (pid() == 0)
  {
    // Write Particle Positions Size
    fwrite(&data_size, sizeof(u_int64_t), 1, fp);
    header_off = ftell(fp);
  }
#if _MPI
  MPI_Bcast(&header_off, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
#endif
  //~ printf("POINTS %d npart: %lu header %lu partoff %lu fileloc %lu\n", pid(),pn[P], header_off, parts_off, header_off + parts_off*3*sizeof(float));

  // Go to Write Position
  fseek(fp, header_off + parts_off*3*sizeof(float), SEEK_SET);

  // Write Points Data
  foreach_particle_in (P) {
    float pos[3] = {z, y, x};	//< must convert from float to double for VTK Float32, see Header | "zyx" for compatibility with htg
    fwrite (pos, sizeof(float), 3, fp);
  }

  // Move to End of Data
  header_off += data_size;

  // Write u
  if (pid() == 0)
  {
    fseek(fp, header_off, SEEK_SET);
    fwrite(&data_size, sizeof(u_int64_t), 1, fp);
  }

  // Go to Write Position
  header_off += 8;	//< uint 64 size
  fseek(fp, header_off + parts_off*3*sizeof(float), SEEK_SET);

  //~ printf("UUUUUu %d npart: %lu header %lu partoff %lu fileloc %lu\n", pid(),pn[P], header_off, parts_off, header_off + parts_off*3*sizeof(float));

  // Write u Data
  foreach_particle_in (P) {
    float val[3] = {pl[l][j].u.z, pl[l][j].u.y, pl[l][j].u.x};	//< must convert from float to double for VTK Float32, see Header | "zyx" for compatibility with htg
    fwrite (val, sizeof(float), 3, fp);
  }

  // Move to End of Data
  header_off += data_size;

	data_size = nparts_all*1*sizeof(float);
	if (pid() == 0)
	{
		fseek(fp, header_off, SEEK_SET);
		fwrite(&data_size, sizeof(u_int64_t), 1, fp);
	}
	header_off += 8;	//< uint 64 size
	fseek(fp, header_off+ parts_off*1*sizeof(float), SEEK_SET);
	// Write density Data
	foreach_particle_in (P) {
		float val[1] = {pl[l][j].tag};	//< must convert from float to double for VTK Float32, see Header | "zyx" for compatibility with htg
		fwrite (val, sizeof(float), 1, fp);
	}
	
	// Move to End of density Data
	header_off += data_size;

  if (particleChar == true){
	/** Density  */
	data_size = nparts_all*1*sizeof(float);
	if (pid() == 0)
	{
		fseek(fp, header_off, SEEK_SET);
		fwrite(&data_size, sizeof(u_int64_t), 1, fp);
	}
	header_off += 8;	//< uint 64 size
	fseek(fp, header_off+ parts_off*1*sizeof(float), SEEK_SET);
	// Write density Data
	foreach_particle_in (P) {
		float val[1] = {pl[l][j].u2.x};	//< must convert from float to double for VTK Float32, see Header | "zyx" for compatibility with htg
		fwrite (val, sizeof(float), 1, fp);
	}
	
	// Move to End of density Data
	header_off += data_size;
	/** Radius */	
	data_size = nparts_all*1*sizeof(float);
	if (pid() == 0)
	{
		fseek(fp, header_off, SEEK_SET);
		fwrite(&data_size, sizeof(u_int64_t), 1, fp);
		//~ printf("Radius, npartsall:%lu\n",nparts_all);
	}
	header_off += 8;	//< uint 64 size
	fseek(fp, header_off+ parts_off*1*sizeof(float), SEEK_SET);
	// Write radius Data
	foreach_particle_in (P) {
		
		float val[1] = {pl[l][j].u2.y};	//< must convert from float to double for VTK Float32, see Header | "zyx" for compatibility with htg
		//~ printf("%f\n",val[0]);
		fwrite (val, sizeof(float), 1, fp);
	}
	
	// Move to End of radius Data
	header_off += data_size;	
	/** tau */
	data_size = nparts_all*1*sizeof(float);
	if (pid() == 0)
	{
		fseek(fp, header_off, SEEK_SET);
		fwrite(&data_size, sizeof(u_int64_t), 1, fp);
	}
	header_off += 8;	//< uint 64 size
	fseek(fp, header_off+ parts_off*1*sizeof(float), SEEK_SET);
	// Write tau Data
	foreach_particle_in (P) {
		float val[1] = {pl[l][j].u2.z};	//< must convert from float to double for VTK Float32, see Header | "zyx" for compatibility with htg
		fwrite (val, sizeof(float), 1, fp);
	}
	// Move to End of tau Data
	header_off += data_size;
  
  }
	  
  // Write Tail
  if (pid() == 0)
  {
    fseek(fp, header_off, SEEK_SET);
    fprintf(fp, "ENDBINARY\n  </AppendedData>\n</VTKFile>\n");
  }

  #if defined(_OPENMP)
    omp_set_num_threads(num_omp);
  #endif
}
