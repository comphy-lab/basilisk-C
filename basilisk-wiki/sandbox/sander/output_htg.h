#include "output_pvd.h"
#include "utils.h"

/**
This field enables the export of  AMR-meshes  into the HyperTreeGrid Dataformat.
It works on Quad-/Octrees with and without MPI. 
The HTG file format is still under development so compatibility with future Paraview versions 
might be limited. Advantage of the HTG file format is the ~7x smaller size due to implicit point location.
The data is stored under the `path` variable, a time series can be opend
using the corresponding collection file (.pvd) with name "output_<path>.pvd"



`output_htg()` is a wrapper function that additionally creates a collection 
file pointing to different time steps of an htg simulation. If MPI is used 
the HTG is written with collective MPI_File_IO. If no MPI is used standard
POSIX fopen/fclose operations are used.

Be sure to used the updated, MPI_IO compatible `output_pvd.h`!

## Example Usage ##
    event report(i+=10){
      char path[]="htg"; // no slash at the end!!
      char prefix[80];
      sprintf(prefix, "data_%06d", i);
      output_htg((scalar ) {cs, d,p},(vector ){u}, path, prefix, i, t);
    }
  The Collection File Name includes the used path. For the example above
   the Collection File is called "output_htg". This allows to run a parameter
   study based on the same binary if you include the current case configuration 
   in the path name. (e.g. Level, Reynolds number, velocity, acceleration, ...)

## Flags ##
    #define HEADER_MIN_MAX_VAL 1
If this flag is set the corresponding min/max values for each values are included
in the header. This is probably not needed, the adavantage is unknown. 
It probably increases the posprocessing speed.
Defaults to 1 (true).

    #define WRITE_HTG_LEVEL 1
If this flag is set the level of each cell is exported as well and selectable during postprocessing.
Defaults to 0 (false). 

    #define HTG_SPEED_STATS 1
If this flag is set information about the write speed are posted (to stdout).
Defaults to 0 (false).

    #define VTK_FILE_VERSION 20
If this flag is set to 20, the created HyperTreeGrid has the VTKFile version 2.0. It is supported
starting with Paraview Version 5.10.
Defaults to 10 (VTKFile version 1.0)

## Known HTG Problems ##
* HyperTreeGridToDualGrid stops working when advancing time step (Paraview related)
* Contour Filter does not work on HTG with only one tree (like this exporter creates) (Paraview related)
* x-z Axis are swapped (3D), x-y Axis are swapped (2D)


## Details of the parallel implementation ##
The mapping of the Basilisk internal tree structure to the HyperTreeGrid
structure is not trivial. The HTG structure is defined level wise. Multiple processors
can write a single tree, if each processes writes its cells level wise after those from
the preceeding process. The offsets and cells per level have to recalculated every 
output step.
These offsets and cells per level are used to create a derived datatype (MPI_Type_indexed).
Each process copys the values of all its nodes into a continous array. The MPI_Type_indexed 
datatype allows to write this continous array as a correclty spaced array to the file, which 
allows other processes to write their data in the spaces left.
The HyperTreeGrid further needs the size of the following data in bytes in addition to the actual data.

   
     _| size_in_bytes_of_the_following_blob |_actual_data_scalar1 | size...blob |_ac...data |...
   

It is therefore  convenient to define a struct consisting of the an integer with the 
blob size and the array of the data to write. This can be done using an MPI_Type_struct
datatype. While the struct contains the interger and a pointer to a continous space in
memory the data on the file will be the integer followed by the correctly spaced MPI_Type_indexed.

This enables using MPI_File_set_view and MPI_File_write_all (collective operation) for good
IO-performance.

The header and the tail is written only by process 0.

The actual tree structure of the HTG ist defined using a bitmask. If a node has children
its value is 1, else it is 0. This is done level wise as well. On the next level the value
of each (still existing) node is described using 1 (is refined/has childen) or 0 
(not refined/ does not have childen) as well. This is repeated for each level except 
the last level because these values are all 0 due to being the last level and are neglected.

The problem arises when storing this bitmask, as it is only possible the write a byte (8 bits).
If a process on a level does not have a number of nodes that is devisible by 8 (a whole byte), logical gaps
in the bitmask appear which are not allowed. This problem arises immediately at level 0,
as process 0 has one 1 bit (level zero has children), but writing  a `1000 0000` byte is wrong
as it already describes the values of the 7 following nodes as well, whose values are unknown to
process 0.
Therefore the describtor bits (stored as bytes (u_int8_t) in memory)  have to be moved to the next
process (or level) until a number divisible by 8 is resident on the process. A 'wave' of excess 
descriptor bytes is swept through all processes and levels til the last process of the last level
who fills the missing bytes with 0.

To allow this moving of upto 7 bytes a special array structure is chosen:

 |8_empty_byt|descr_byt_lvl0|7_empty_byt|descr_byt_lvl1|...|7_empty_byt|
 0           8              

Each process has to space to accommodate upto 7 bytes from the preceeding process/level
infront of the own descriptor bytes. Each process keeps track of how many bytes to send/
to recv/ and stay on a level and on which index the level begins taking the recvieved 
bytes into account. 
After this every process has a number of descriptor bytes that can be transformed into
descriptor bits. This is actually done in the same array starting a index 0 which is guaranteed
to not contain any data. In byte 0 the information of the next 8 relevant bytes is written as bits.

With information about offset and length again a MPI_Type_indexed and a MPI_Type struct
can be defined. This again is used to write the bitmask in a collective operation.

## Detail of the serial implementation ##
Each process caches the data on each level in a contingues array which is written
to disk sequentially with only one fwrite() funtion call.
The Bitmask is written in "appened" mode as well resulting in an uncluttered file
 */ 

//~ #define HEADER_MIN_MAX_VAL 0
//~ #define WRITE_HTG_LEVEL 1
//~ #define HTG_SPEED_STATS 1
//~ #define VTK_FILE_VERSION 20 // 2.0 

#ifndef HEADER_MIN_MAX_VAL
#define HEADER_MIN_MAX_VAL 1
#endif
#ifndef WRITE_HTG_LEVEL
#define WRITE_HTG_LEVEL 0
#endif
#ifndef HTG_SPEED_STATS
#define HTG_SPEED_STATS 0
#endif
#ifndef VTK_FILE_VERSION
#define VTK_FILE_VERSION 10 // 1.0
#endif

void output_htg(scalar * list, vector * vlist, const char* path, char* prefix, int i, double t);

#if _MPI
void output_htg_data_mpiio(scalar * list, vector * vlist, MPI_File fp);
#else
void output_htg_data(scalar * list, vector * vlist, FILE * fp);
#endif
 
#if _MPI
void output_htg(scalar * list, vector * vlist, const char* path, char* prefix, int i, double t) {
	MPI_File fp;
	
	int ec;
	char htg_name[80];  	  
	sprintf(htg_name, "%s/%s.htg", path, prefix);  

	ec = MPI_File_open(MPI_COMM_WORLD, htg_name, \
		MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fp);
  // Overwrite File 
	if (ec == MPI_ERR_FILE_EXISTS) {
    printf("ERR, htg_name exists!\n");
		// MPI_File_delete(htg_name, MPI_INFO_NULL);
		MPI_File_open(MPI_COMM_WORLD, htg_name, \
		  MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fp);
    MPI_File_set_size(fp, 0);
	}
	
	if(ec != MPI_SUCCESS){
		printf("output_htg.h : %s could not be opened\n Does the Folder exist?\n", htg_name);
		MPI_Abort(MPI_COMM_WORLD, 2);
	}  
	
	output_htg_data_mpiio((scalar *) list,(vector *)vlist,fp);

	MPI_File_close(&fp);
	
  if(pid()==0){
		bool firstTimeWritten = false;
		char pvd_name[80];	  
		sprintf(pvd_name,"output_%s.pvd",path);
	
		ec = MPI_File_open(MPI_COMM_SELF, pvd_name, \
			MPI_MODE_RDWR, MPI_INFO_NULL, &fp);
      
		if( (i == 0) ||  (ec == MPI_ERR_NO_SUCH_FILE) ) {
			// Overwrite File 	
			// Best implementation??		
			MPI_File_open(MPI_COMM_SELF, pvd_name, \
			  MPI_MODE_WRONLY|MPI_MODE_CREATE, MPI_INFO_NULL, &fp);			
			MPI_File_set_size(fp, 0);
			firstTimeWritten = true;
		}
		output_pvd_mpiio(htg_name, t, fp, firstTimeWritten);
		
		MPI_File_close(&fp);
	}	
	MPI_Barrier(MPI_COMM_WORLD);    
}

#else // no MPI

void output_htg(scalar * list, vector * vlist, const char* path, char* prefix, int i, double t) {
	FILE * fp ;

	char htg_name[80];  	  
	sprintf(htg_name, "%s/%s.htg", path, prefix);  

	fp = fopen(htg_name, "w");	  
	if(!fp){
		printf("output_htg.h : %s could not be opened\n Does the Folder exist?\n", htg_name);
		exit(1);
	}
  
  output_htg_data((scalar *) list,(vector *) vlist, fp);

  fclose(fp);

  bool firstTimeWritten = false;
  char pvd_name[80];	  
  sprintf(pvd_name,"output_%s.pvd",path);
  fp = fopen(pvd_name, "r+");
  if( (i == 0) ||  (fp == NULL) ) {
    fp = fopen(pvd_name,"w");
    firstTimeWritten = true;
  }
  output_pvd(htg_name, t, fp, firstTimeWritten);
  fclose(fp);	
}

#endif


#if _MPI

#define Write2File(x) do { \
  (x); \
  MPI_File_write(fp,&buffer, strlen(buffer), MPI_CHAR, MPI_STATUS_IGNORE); \
  } while(0)
  
void output_htg_data_mpiio(scalar * list, vector * vlist, MPI_File fp)
{
	#if defined(_OPENMP)
		int num_omp = omp_get_max_threads();
		omp_set_num_threads(1);
	#endif
    
	unsigned int vertices_local = 0;
	unsigned int descBits_local = 0;  
  unsigned int vertices_local_pL[grid->maxdepth+1];  // pL = per Level 
  
  unsigned int descBits;
  unsigned int vertices;  
  unsigned int vertices_pL[grid->maxdepth+1];  
	    
	for (int lvl = 0; lvl < grid->maxdepth+1; ++lvl){
    vertices_local_pL[lvl] = 0;
		foreach_level(lvl,serial) 
			if(is_local(cell)) 
				vertices_local_pL[lvl]++;
    
    vertices_local += vertices_local_pL[lvl];         
  }
  
  descBits_local = vertices_local - vertices_local_pL[grid->maxdepth]; // 
  
  MPI_Reduce(&vertices_local_pL[0], &vertices_pL[0], grid->maxdepth + 1,\
              MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);  
 
#if HTG_SPEED_STATS
	MPI_Barrier(MPI_COMM_WORLD);
	double start = MPI_Wtime();
#endif	
  MPI_Offset offset = 0;

	int vertices_global_offset[grid->maxdepth+1];
  vertices_global_offset[0] = 0;
  unsigned int carryover = 0;
  for (int lvl = 0; lvl <= grid->maxdepth; ++lvl){
    MPI_Exscan (&vertices_local_pL[lvl], \
                &vertices_global_offset[lvl], \
                1, \
                MPI_UNSIGNED, \
                MPI_SUM, \
                MPI_COMM_WORLD);    
    // Pass Offset to process 0 for Exscan on next level
    if (pid() == (npe() - 1)) {
      unsigned int next_offset;
      next_offset = vertices_global_offset[lvl] + vertices_local_pL[lvl];
      MPI_Ssend(&next_offset, 1, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD);            
    }
    if (pid() == 0) {
      vertices_local_pL[lvl] -= carryover; // temporarily remove offset added in previous level
      
      MPI_Recv(&carryover, 1, MPI_UNSIGNED, npe() - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);      
      
      if (lvl < grid->maxdepth) {
        vertices_local_pL[lvl+1] += carryover; // temporarily use array to carry over offset to next level
        vertices_global_offset[lvl+1] = carryover;
      }      
      else
        vertices = carryover;
      if (lvl == grid->maxdepth - 1)
        descBits = carryover;
    }
  }
  
  MPI_Bcast(&vertices, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  MPI_Bcast(&descBits, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  

 /**
   * Calculate min,max for xml header (prob. not needed?)
   */
#if HTG_SPEED_STATS     
  MPI_Barrier(MPI_COMM_WORLD);
  double t_stats_start = MPI_Wtime();
#endif  
  double min_val[list_len(list)];
  double max_val[list_len(list)];
  {
    int i = 0;
    for (scalar s in list) {
#if HEADER_MIN_MAX_VAL
      stats stat = statsf(s);
      min_val[i]=stat.min;
      max_val[i]=stat.max;
#else
      min_val[i]=0.;
      max_val[i]=0.;
#endif
      i++;
    }
	}
  double min_val_v[vectors_len(vlist)];
  double max_val_v[vectors_len(vlist)];
  { 
    int i = 0;	
    for (vector v in vlist) {
#if HEADER_MIN_MAX_VAL
      min_val_v[i]= HUGE;
      max_val_v[i]= -HUGE;
      foreach_dimension(){
        stats stat = statsf(v.x);	
        min_val_v[i] = min(stat.min,min_val_v[i]);
        max_val_v[i]=  max(stat.max,max_val_v[i]);			
      }
#else
      min_val_v[i]= 0.;
      max_val_v[i]= 0.;
#endif
      i++;
    }
	}
#if HTG_SPEED_STATS    
  MPI_Barrier(MPI_COMM_WORLD);
  double t_stats = MPI_Wtime() - t_stats_start ;
	double t_header_start = MPI_Wtime();
#endif  
	char buffer[256];
	/** File Header vtkHyperTreeGrid */
	if (pid() == 0){
    //~ WriteHeader(&offset, infoStruct);
    int maj_v = VTK_FILE_VERSION/10, min_v = VTK_FILE_VERSION%10;
    
    Write2File(sprintf(buffer,\
      "<VTKFile %s version=\"%i.%i\" %s %s>\n",\
        "type=\"HyperTreeGrid\"", \
        maj_v, min_v , \
        "byte_order=\"LittleEndian\" ",\
        "header_type=\"UInt32\"")); 
  
  
  #if dimension==2
		Write2File(sprintf(buffer,\
      "\t<HyperTreeGrid BranchFactor=\"2\" TransposedRootIndexing=\"0\" Dimensions=\"%d %d %d\">\n",\
      2,2,1));    
  #elif dimension==3
		Write2File(sprintf(buffer,\
      "\t<HyperTreeGrid BranchFactor=\"2\" TransposedRootIndexing=\"0\" Dimensions=\"%d %d %d\">\n", 2,2,2));
  #endif
		Write2File(sprintf(buffer,"\t\t<Grid>\n"));
  #if dimension==2
		Write2File(sprintf(buffer, "\t\t\t<DataArray type=\"Float64\" Name=\"XCoordinates\" NumberOfTuples=\"2\" format=\"ascii\" RangeMin=\"%g\" RangeMax=\"%g\">\n",Y0, Y0+L0));
    Write2File(sprintf(buffer, "\t\t\t\t%g %g",Y0, Y0+L0));
    Write2File(sprintf(buffer,"\n\t\t\t</DataArray>\n"));
    Write2File(sprintf(buffer, "\t\t\t<DataArray type=\"Float64\" Name=\"YCoordinates\" NumberOfTuples=\"2\" format=\"ascii\" RangeMin=\"%g\" RangeMax=\"%g\">\n",X0, X0+L0));
    Write2File(sprintf(buffer, "\t\t\t\t%g %g",X0, X0+L0));
    Write2File(sprintf(buffer,"\n\t\t\t</DataArray>\n"));	  
		Write2File(sprintf(buffer, "\t\t\t<DataArray type=\"Float64\" Name=\"ZCoordinates\" NumberOfTuples=\"2\" format=\"ascii\" RangeMin=\"%g\" RangeMax=\"%g\">\n",Z0, Z0+L0));		
		Write2File(sprintf(buffer, "\t\t\t\t%g %g", 0., 0.));	  
  #elif dimension==3
		Write2File(sprintf(buffer, "\t\t\t<DataArray type=\"Float64\" Name=\"XCoordinates\" NumberOfTuples=\"2\" format=\"ascii\" RangeMin=\"%g\" RangeMax=\"%g\">\n",\
        Z0, Z0+L0));
		Write2File(sprintf(buffer, "\t\t\t\t%g %g",Z0, Z0+L0));
		Write2File(sprintf(buffer,"\n\t\t\t</DataArray>\n"));
		Write2File(sprintf(buffer, "\t\t\t<DataArray type=\"Float64\" Name=\"YCoordinates\" NumberOfTuples=\"2\" format=\"ascii\" RangeMin=\"%g\" RangeMax=\"%g\">\n",\
        Y0, Y0+L0));
		Write2File(sprintf(buffer, "\t\t\t\t%g %g",Y0, Y0+L0));
		Write2File(sprintf(buffer,"\n\t\t\t</DataArray>\n"));	  
		Write2File(sprintf(buffer, "\t\t\t<DataArray type=\"Float64\" Name=\"ZCoordinates\" NumberOfTuples=\"2\" format=\"ascii\" RangeMin=\"%g\" RangeMax=\"%g\">\n",\
        X0, X0+L0));		
		Write2File(sprintf(buffer, "\t\t\t\t%g %g",X0, X0+L0));	  
  #endif
		Write2File(sprintf(buffer,"\n\t\t\t</DataArray>\n"));	  
		Write2File(sprintf(buffer,"\t\t</Grid>\n"));
		Write2File(sprintf(buffer,"\t\t<Trees>\n"));
    
		// Tree Begin
    unsigned int byte_offset = 0;
  #if VTK_FILE_VERSION == 10
		Write2File(sprintf(buffer,"\t\t\t<Tree Index=\"0\" NumberOfLevels=\"%d\" NumberOfVertices=\"%u\">\n",\
        grid->maxdepth + 1, vertices));   	
		
    // Array Descriptor Bits    
    Write2File(sprintf(buffer,"\t\t\t\t<DataArray type=\"Bit\" Name=\"Descriptor\" NumberOfTuples=\"%u\" format=\"appended\" RangeMin=\"0\" RangeMax=\"1\" offset=\"%u\"/>\n",\
      descBits, byte_offset));        
    byte_offset += (descBits/8+1) * sizeof(u_int8_t) + sizeof(u_int32_t);
    
		// Array NbVerticesByLevel
		Write2File(sprintf(buffer,"\t\t\t\t<DataArray type=\"Int64\" Name=\"NbVerticesByLevel\" NumberOfTuples=\"%d\" format=\"ascii\" RangeMin=\"1\" RangeMax=\"%u\" >\n\t\t\t\t\t",\
        grid->maxdepth + 1, vertices_pL[grid->maxdepth]));
    		
    for (int lvl = 0; lvl <= grid->maxdepth; lvl++) {
			Write2File(sprintf(buffer,"%u ", vertices_pL[lvl]));
		}		
		Write2File(sprintf(buffer,"\n\t\t\t\t</DataArray>\n"));
  #endif     
  
  #if VTK_FILE_VERSION == 20
    // Array Descriptor Bits    
    Write2File(sprintf(buffer,"\t\t\t\t<DataArray type=\"Bit\" Name=\"Descriptors\" NumberOfTuples=\"%u\" format=\"appended\" RangeMin=\"0\" RangeMax=\"1\" offset=\"%u\"/>\n",\
      descBits, byte_offset));        
    byte_offset += (descBits/8+1) * sizeof(u_int8_t) + sizeof(u_int32_t);
    
    Write2File(sprintf(buffer,"\t\t\t\t<DataArray type=\"Int64\" Name=\"NumberOfVerticesPerDepth\" NumberOfTuples=\"%d\" format=\"ascii\" RangeMin=\"1\" RangeMax=\"%u\" >\n\t\t\t\t\t",\
        grid->maxdepth + 1, vertices_pL[grid->maxdepth]));
          
    for (int lvl = 0; lvl <= grid->maxdepth; lvl++) {
      Write2File(sprintf(buffer,"%u ", vertices_pL[lvl]));
    }		
    Write2File(sprintf(buffer,"\n\t\t\t\t</DataArray>\n"));
    
    Write2File(sprintf(buffer,"\t\t\t\t<DataArray type=\"Int64\" Name=\"TreeIds\" NumberOfTuples=\"1\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\" >\n"));
    Write2File(sprintf(buffer,"\t\t\t\t\t%i\n",0));
    Write2File(sprintf(buffer,"\t\t\t\t</DataArray>\n"));
    
    Write2File(sprintf(buffer,"\t\t\t\t<DataArray type=\"UInt32\" Name=\"DepthPerTree\" NumberOfTuples=\"1\" format=\"ascii\" RangeMin=\"%i\" RangeMax=\"%i\" >\n",\
    grid->maxdepth+1, grid->maxdepth+1));
    Write2File(sprintf(buffer,"\t\t\t\t\t%i\n",grid->maxdepth+1));
    Write2File(sprintf(buffer,"\t\t\t\t</DataArray>\n"));
            
    Write2File(sprintf(buffer,"\t\t</Trees>\n"));
  #endif  
		// Cell Data Begin
		Write2File(sprintf(buffer,"\t\t\t\t<CellData>\n"));
    
#if WRITE_HTG_LEVEL    
		Write2File(sprintf(buffer,"\t\t\t\t\t<DataArray type=\"UInt8\" Name=\"Level\" NumberOfTuples=\"%u\" format=\"appended\" RangeMin=\"0\" RangeMax=\"%i\" offset=\"%u\"/>\n",\
        vertices, grid->maxdepth, byte_offset));		
		byte_offset += vertices * sizeof(u_int8_t) + sizeof(u_int32_t);
#endif    
		{
      int i = 0;
      for (scalar s in list)
      {			
        Write2File(sprintf(buffer,"\t\t\t\t\t<DataArray type=\"Float32\" Name=\"%s\" NumberOfTuples=\"%u\" format=\"appended\" RangeMin=\"%g\" RangeMax=\"%g\" offset=\"%u\"/>\n",\
            s.name, vertices, min_val[i], max_val[i], byte_offset));			    
        byte_offset += vertices * sizeof(float_t) + sizeof(u_int32_t);
        i++;
      }
    }
    {
      int i = 0;
      for (vector v in vlist)
      {
        char *vname = strtok(v.x.name, ".");
        Write2File(sprintf(buffer, "\t\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"%i\" Name=\"%s\" NumberOfTuples=\"%u\" format=\"appended\"  RangeMin=\"%g\" RangeMax=\"%g\" offset=\"%u\"/>\n",\
            3,vname, vertices, min_val_v[i],max_val_v[i],byte_offset));
        byte_offset += vertices * 3 * sizeof(float_t) + sizeof(u_int32_t);
        i++;
      } 
    }
		// Cell Data End  		
		Write2File(sprintf(buffer,"\t\t\t\t</CellData>\n"));
#if VTK_FILE_VERSION == 10
    Write2File(sprintf(buffer,"\t\t\t</Tree>\n\t\t</Trees>\n"));
#endif

    Write2File(sprintf(buffer,"\t</HyperTreeGrid>\n\t<AppendedData encoding=\"raw\">\n_"));
    MPI_Offset offset_tmp;
    MPI_File_get_position(fp, &offset_tmp);
    offset += offset_tmp;    
	} // end pid() == 0
	
	MPI_Bcast(&offset, 1, MPI_OFFSET, 0, MPI_COMM_WORLD);
#if HTG_SPEED_STATS  
  	double t_header = MPI_Wtime() - t_header_start;
#endif
  /** Write Bitmask */
  int cell_size;
  
  { 
    cell_size=sizeof(u_int8_t);
    int vertices_local_pL_offset[grid->maxdepth+1];
    /** Create an array large enough to hold data + spacing between the levels */
    long length_w_spacing = descBits_local + 7*(grid->maxdepth+1) + 8;    
    u_int8_t* mask = (u_int8_t*)calloc( length_w_spacing, cell_size); //calloc needed?
    
    long index = 8; // start a 8
    for(int lvl=0; lvl < grid->maxdepth;++lvl) { // iterate through array and place lokal mask.
      vertices_local_pL_offset[lvl] = index;
      foreach_level(lvl,serial) {
        if (is_local(cell)) {
          mask[index++] = (u_int8_t)(!is_leaf(cell));
        }
      }
      index += 7; // make space for byte movement
    }
    vertices_local_pL_offset[grid->maxdepth] = index;
    
    assert(length_w_spacing > index);
  
    int vertices_local_pL_corr[grid->maxdepth+1];
    
    for (int lvl = 0; lvl < grid->maxdepth; ++lvl) {
      vertices_local_pL_corr[lvl] = vertices_local_pL[lvl];
    }
    
    
    int vertices_local_pL_offset_corr[grid->maxdepth]; // new offset
    long vertices_local_corr = 0;
    
    int num_from_prev = 0;
    u_int8_t tmp[7];
    for (int lvl = 0; lvl < grid->maxdepth; ++lvl) {
      for (int pe=0; pe < npe();  ++pe) {        
        int send_rank = pe;
        int recv_rank = (pe+1) % npe(); // next higher rank + loop
                 
        if (pid() == send_rank) {
          // work on prev send data          
          vertices_local_pL_corr[lvl] += num_from_prev;          
          int trg_position = vertices_local_pL_offset[lvl]-num_from_prev;
          vertices_local_pL_offset_corr[lvl] = trg_position;            
          for (int tmp_cnt = 0; tmp_cnt < num_from_prev; ++tmp_cnt)
          {
            // printf("pid: %i, lvl: %i trg_pos: %i, tmp_cnt: %i, lws: %li\n", pid(), lvl, trg_position, tmp_cnt,length_w_spacing); fflush(stdout);
            mask[trg_position + tmp_cnt] = tmp[tmp_cnt];
          }
          int num_to_next = (vertices_local_pL_corr[lvl]%8);
         
          vertices_local_pL_corr[lvl] -= num_to_next;
          
          // if last level and last process, append space 
          if ((lvl ==  grid->maxdepth -1 ) && (pid() == npe() -1) ) {
            vertices_local_pL_corr[lvl] += 8;            
          }
          
          vertices_local_corr += vertices_local_pL_corr[lvl];               
          
          // MPI_Ssend(&num_to_next,  1, MPI_INT, recv_rank, 0, MPI_COMM_WORLD); 
          MPI_Send(&num_to_next,  1, MPI_INT, recv_rank, 0, MPI_COMM_WORLD); 
                      
          int src_position = vertices_local_pL_offset[lvl+1] - 7 - num_to_next; // sum of vertices incl thsi level          
          // MPI_Ssend(&mask[src_position],  num_to_next, MPI_UINT8_T, recv_rank, 0, MPI_COMM_WORLD); 
          MPI_Send(&mask[src_position],  num_to_next, MPI_UINT8_T, recv_rank, 1, MPI_COMM_WORLD); 
        }     
        if (pid() == recv_rank) {
          MPI_Recv(&num_from_prev, 1, MPI_INT, send_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
          MPI_Recv(&tmp[0], num_from_prev, MPI_UINT8_T, send_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);           
        }
      }                        
    }
    
        
    if (vertices_local_corr %8 != 0) // sanity check
      MPI_Abort(MPI_COMM_WORLD, 2);
      
    //* TEMP EXPORT FILE STRUCTURE */
    #if 0
    FILE * fpt;
    char exp_name[80];
    sprintf(exp_name, "structure_pid%04i.dat", pid());
    fpt = fopen(exp_name, "w");
    fprintf(fpt, "%li\n",  length_w_spacing);
    for (int lvl = 0; lvl < grid->maxdepth; ++lvl) {       
      fprintf(fpt, "%i\n",  vertices_local_pL_offset_corr[lvl]);
      fprintf(fpt, "%i\n",  vertices_local_pL_corr[lvl]) ;
    }
    fclose(fpt);
    #endif
    //* END TEMP EXPORT FILE STRUCTURE */
    
    // zero copy bitmask creation
    
    int i = 0, cnt = 0;
    for (int lvl = 0; lvl < grid->maxdepth; ++lvl) {      
      int displacement = vertices_local_pL_offset_corr[lvl]; // offset of correct array within my array
      int count = vertices_local_pL_corr[lvl]; // corrected count        
      for (int c = 0; c < count; ++c) {
        mask[i] |= mask[displacement + c] << (7-cnt);        
        if (++cnt%8 == 0) {
           mask[++i] = 0; // incr i and clear next byte
           cnt=0; 
         }
      }
    }
        
    /** We need to calculate the global offset of the bit packets for the
     *  MPI_Type_indexed datatype describing the tree data */ 
    int vertices_global_offset_corr[grid->maxdepth]; // not "+1"
    vertices_global_offset_corr[0] = 0;
    unsigned int carryover = 0;
    for (int lvl = 0; lvl < grid->maxdepth; ++lvl){
    MPI_Exscan (&vertices_local_pL_corr[lvl], \
                &vertices_global_offset_corr[lvl], \
                1, \
                MPI_INT, \
                MPI_SUM, \
                MPI_COMM_WORLD);    
      // Send Data to Process 0 for the Exscan on next Level
      if (pid() == (npe() - 1)) {
        unsigned int next_offset;
        next_offset = vertices_global_offset_corr[lvl] + vertices_local_pL_corr[lvl];
        MPI_Ssend(&next_offset, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);            
      }
      if (pid() == 0) {
        vertices_local_pL_corr[lvl] -= carryover; // remove temporary offset added in previous level
        
        MPI_Recv(&carryover, 1, MPI_INT, npe() - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);            
        if (lvl + 1 < grid->maxdepth ) {
          vertices_local_pL_corr[lvl+1] += carryover; // use array temporarily to carry over offset to next level
          vertices_global_offset_corr[lvl+1] = carryover;
        }      
      } 
    }
  
    /** Create a struct which will be written to RAM */
    struct descBit_t{
      u_int32_t size;
      u_int8_t* data;
    }descBit_struct;
        
    descBit_struct.size = (descBits/8 + 1) * cell_size;
    descBit_struct.data =  &mask[0];
 
    /** Define where the data is located in Memory */
 
    MPI_Aint m_displacements[2];
    MPI_Aint base_address;
    MPI_Get_address(&descBit_struct,         &base_address);
    MPI_Get_address(&descBit_struct.size,    &m_displacements[0]);
    MPI_Get_address(&descBit_struct.data[0], &m_displacements[1]);
    m_displacements[0] = MPI_Aint_diff(m_displacements[0], base_address);
    m_displacements[1] = MPI_Aint_diff(m_displacements[1], base_address);
    
    MPI_Datatype m_types[2] = { MPI_UINT32_T, MPI_BYTE };    
    //~ /** Tree has to have some data, else there is a segmentationfault when setting view */
    int m_lengths[2] = {1, vertices_local_corr/8}; // Bits
        
    MPI_Datatype m_view; // memory view
    MPI_Type_create_struct(2, m_lengths, m_displacements, m_types, &m_view);
    MPI_Type_commit(&m_view);    
    
    int lengths[grid->maxdepth];
    int displacements[grid->maxdepth];
    /** Define the File view for writing the bitmask [1/2]*/
    for(int lvl = 0; lvl < grid->maxdepth; ++lvl) {
      lengths[lvl]       = (int)vertices_local_pL_corr[lvl]/8; // bits, not bytes
      displacements[lvl] = (int)vertices_global_offset_corr[lvl]/8; // bits, not bytes
    }    
    
    MPI_Datatype tree_type_descBit;
    MPI_Type_indexed(grid->maxdepth, lengths, displacements, MPI_BYTE, &tree_type_descBit);
    MPI_Type_commit(&tree_type_descBit);
    
    /** Define the File view for writing the bitmask [2/2]*/
    MPI_Aint f_displacements[2] = {0, sizeof(u_int32_t)};
    int f_lengths[2] = {1, 1};
    MPI_Datatype f_types[2] = { MPI_UINT32_T, tree_type_descBit};    
    
    MPI_Datatype f_view; // file_view
    MPI_Type_create_struct(2, f_lengths, f_displacements, f_types, &f_view);
    MPI_Type_commit(&f_view);

    /** Set the view according to f_fiew (file_view)*/
    MPI_File_set_view(fp, offset, f_view, f_view, "native", MPI_INFO_NULL);
    /** Write from memory according to m_view (memory_view) */
    MPI_File_write_all(fp, &descBit_struct, 1, m_view, MPI_STATUS_IGNORE);     
    
    offset += (descBits/8+1) * sizeof(u_int8_t) + sizeof(u_int32_t);
    
    MPI_Type_free(&m_view);
    MPI_Type_free(&f_view);           
    MPI_Type_free(&tree_type_descBit);
    
    free(mask); mask = NULL;       
    descBit_struct.data = NULL;   
  }

  /** Write LEVEL */
#if WRITE_HTG_LEVEL
  {
    struct level_t{
      u_int32_t size;
      u_int8_t* data;
    }level_struct;
        
    cell_size=sizeof(u_int8_t);
    level_struct.size =  vertices * cell_size;
    level_struct.data = (u_int8_t*) malloc(vertices_local*cell_size);
 
    MPI_Aint m_displacements[2];
    MPI_Aint base_address;
    MPI_Get_address(&level_struct,         &base_address);
    MPI_Get_address(&level_struct.size,    &m_displacements[0]);
    MPI_Get_address(&level_struct.data[0], &m_displacements[1]);
    m_displacements[0] = MPI_Aint_diff(m_displacements[0], base_address);
    m_displacements[1] = MPI_Aint_diff(m_displacements[1], base_address);
    
    MPI_Datatype m_types[2] = { MPI_UINT32_T, MPI_UINT8_T};
    int m_lengths[2] = {1, vertices_local};
    
    MPI_Datatype m_view;
    MPI_Type_create_struct(2, m_lengths, m_displacements, m_types, &m_view);
    MPI_Type_commit(&m_view);    
    
    int lengths[grid->maxdepth+1];
    int displacements[grid->maxdepth+1];
    for(int lvl = 0; lvl < grid->maxdepth + 1; ++lvl){
      lengths[lvl]       = (int) vertices_local_pL[lvl];
      displacements[lvl] = (int) vertices_global_offset[lvl];
    }
    MPI_Datatype tree_type_level;
    MPI_Type_indexed(grid->maxdepth + 1, lengths, displacements, MPI_UINT8_T, &tree_type_level);
    MPI_Type_commit(&tree_type_level);
    
    MPI_Aint f_displacements[2] = {0, sizeof(u_int32_t)};
    int f_lengths[2] = {1, 1};
    MPI_Datatype f_types[2] = { MPI_UINT32_T, tree_type_level};
    
    MPI_Datatype f_view;
    MPI_Type_create_struct(2, f_lengths, f_displacements, f_types, &f_view);
    MPI_Type_commit(&f_view);
    
    long index = 0;   
    for(int lvl=0; lvl<=grid->maxdepth;++lvl)		
      foreach_level(lvl,serial) 
        if (is_local(cell)) 
          level_struct.data[index++] = (u_int8_t)lvl;
                          
    MPI_File_set_view(fp, offset, f_view,f_view, "native", MPI_INFO_NULL);        
    MPI_File_write_all(fp, &level_struct, 1, m_view, MPI_STATUS_IGNORE);

    offset += vertices*cell_size + sizeof(u_int32_t);

    free(level_struct.data); level_struct.data = NULL;
    MPI_Type_free(&m_view);       
    MPI_Type_free(&f_view);       
    MPI_Type_free(&tree_type_level);
  }
#endif // end WRITE_HTG_LEVEL
  {
    struct scalar_t{
      u_int32_t size;
      float* data;
    }scalar_struct;

    cell_size=sizeof(float);
    scalar_struct.size = vertices*cell_size;
    scalar_struct.data = (float*) malloc(vertices_local*cell_size);
    
    MPI_Aint m_displacements[2];
    MPI_Aint base_address;
    MPI_Get_address(&scalar_struct,         &base_address);
    MPI_Get_address(&scalar_struct.size,    &m_displacements[0]);
    MPI_Get_address(&scalar_struct.data[0], &m_displacements[1]);
    m_displacements[0] = MPI_Aint_diff(m_displacements[0], base_address);
    m_displacements[1] = MPI_Aint_diff(m_displacements[1], base_address);
    
    MPI_Datatype m_types[2] = { MPI_UINT32_T, MPI_FLOAT};
    int m_lengths[2] = {1, vertices_local};
    
    MPI_Datatype m_view;
    MPI_Type_create_struct(2, m_lengths, m_displacements, m_types, &m_view);
    MPI_Type_commit(&m_view);    
    
    int lengths[grid->maxdepth+1];
    int displacements[grid->maxdepth+1];
    for(int lvl = 0; lvl < grid->maxdepth + 1; ++lvl){
      lengths[lvl]       = (int) vertices_local_pL[lvl];
      displacements[lvl] = (int) vertices_global_offset[lvl];
    }
    MPI_Datatype tree_type_scalar;
    MPI_Type_indexed(grid->maxdepth + 1, lengths, displacements, MPI_FLOAT, &tree_type_scalar);
    MPI_Type_commit(&tree_type_scalar);
    
    MPI_Aint f_displacements[2] = {0, sizeof(u_int32_t)};
    int f_lengths[2] = {1, 1};
    MPI_Datatype f_types[2] = { MPI_UINT32_T, tree_type_scalar};
    
    MPI_Datatype f_view;
    MPI_Type_create_struct(2, f_lengths, f_displacements, f_types, &f_view);
    MPI_Type_commit(&f_view);
    
    for (scalar s in list) {
      long index = 0;
      for(int lvl=0; lvl<=grid->maxdepth;++lvl)
        foreach_level(lvl,serial)
          if (is_local(cell))
             scalar_struct.data[index++] = (float)val(s);
      
      MPI_File_set_view(fp, offset, f_view,f_view, "native", MPI_INFO_NULL);      
      MPI_File_write_all(fp, &scalar_struct, 1, m_view, MPI_STATUS_IGNORE); 
      offset += vertices*cell_size + sizeof(u_int32_t);
    }
    free(scalar_struct.data); scalar_struct.data = NULL;
    MPI_Type_free(&m_view);       
    MPI_Type_free(&f_view);
    MPI_Type_free(&tree_type_scalar);  
  }   
  {
    struct vector_t{
      u_int32_t size;
      float* data;
    }vector_struct;

    cell_size=3*sizeof(float);
    vector_struct.size = vertices*cell_size;
    vector_struct.data = (float*) malloc(vertices_local*cell_size);
    
    MPI_Aint m_displacements[2];
    MPI_Aint base_address;
    MPI_Get_address(&vector_struct,         &base_address);
    MPI_Get_address(&vector_struct.size,    &m_displacements[0]);
    MPI_Get_address(&vector_struct.data[0], &m_displacements[1]);
    m_displacements[0] = MPI_Aint_diff(m_displacements[0], base_address);
    m_displacements[1] = MPI_Aint_diff(m_displacements[1], base_address);
    
    MPI_Datatype m_types[2] = { MPI_UINT32_T, MPI_FLOAT};
    int m_lengths[2] = {1, 3*vertices_local};
    
    MPI_Datatype m_view;
    MPI_Type_create_struct(2, m_lengths, m_displacements, m_types, &m_view);
    MPI_Type_commit(&m_view);    
    
    int lengths[grid->maxdepth+1];
    int displacements[grid->maxdepth+1];
    for(int lvl = 0; lvl < grid->maxdepth + 1; ++lvl) {
      lengths[lvl]       = 3*(int)vertices_local_pL[lvl];
      displacements[lvl] = 3*(int)vertices_global_offset[lvl];
    }
    MPI_Datatype tree_type_vector;
    MPI_Type_indexed(grid->maxdepth + 1, lengths, displacements, MPI_FLOAT, &tree_type_vector);
    MPI_Type_commit(&tree_type_vector);
    
    MPI_Aint f_displacements[2] = {0, sizeof(u_int32_t)};
    int f_lengths[2] = {1, 1};
    MPI_Datatype f_types[2] = { MPI_UINT32_T, tree_type_vector};
    
    MPI_Datatype f_view;
    MPI_Type_create_struct(2, f_lengths, f_displacements, f_types, &f_view);
    MPI_Type_commit(&f_view);        

    for (vector v in vlist) {
      long index = 0;
        for(int lvl=0; lvl<=grid->maxdepth;++lvl)
          foreach_level(lvl,serial) 
            if (is_local(cell)) {
              #if dimension == 2					
              vector_struct.data[index]   = (float)val(v.y);
              vector_struct.data[index+1] = (float)val(v.x);
              vector_struct.data[index+2] = (float)0.;
              index += 3;
              #elif dimension == 3					
              vector_struct.data[index]   = (float)val(v.z);
              vector_struct.data[index+1] = (float)val(v.y);
              vector_struct.data[index+2] = (float)val(v.x);
              index += 3;
              #endif					
            } 			            
      MPI_File_set_view(fp, offset, f_view,f_view, "native", MPI_INFO_NULL);
      MPI_File_write_all(fp, &vector_struct, 1, m_view, MPI_STATUS_IGNORE);       
      offset += vertices*cell_size + sizeof(u_int32_t);
    }
    free(vector_struct.data); vector_struct.data = NULL;
    MPI_Type_free(&m_view);       
    MPI_Type_free(&f_view);
    MPI_Type_free(&tree_type_vector);
  } 
  
  MPI_File_set_view(fp, offset, MPI_BYTE,MPI_BYTE, "native", MPI_INFO_NULL);
	
	if (pid() == 0)	
		Write2File(sprintf(buffer, "ENDBINARY\n\t</AppendedData>\n</VTKFile>\n"));    
  
#if HTG_SPEED_STATS  
  MPI_Barrier(MPI_COMM_WORLD);  
	double end = MPI_Wtime();  
	fprintf (stdout,"Write Time Taken: %f ms, Throughput: %.2f MB/s, Stats %.2f%%, Header: %.2f%%\n",\
    (end-start)*1000., (double)(offset+strlen(buffer))/(end-start)/(double)sq(1024),\
    t_stats/(end-start)*100.,t_header/(end-start)*100.);    
  fflush(stdout);
#endif    
  MPI_File_sync(fp);
  MPI_Barrier(MPI_COMM_WORLD);  
  #if defined(_OPENMP)
    omp_set_num_threads(num_omp);
  #endif
}
#else
void output_htg_data(scalar * list, vector * vlist, FILE * fp)
{
	#if defined(_OPENMP)
		int num_omp = omp_get_max_threads();
		omp_set_num_threads(1);
	#endif  
	
  unsigned int vertices_local = 0;
	unsigned int descBits_local = 0;  
  unsigned int vertices_local_pL[grid->maxdepth+1];  // pL = per Level 

    
	for (int lvl = 0; lvl < grid->maxdepth+1; ++lvl){
    vertices_local_pL[lvl] = 0;
		foreach_level(lvl,serial) 
			if(is_local(cell)) 
				vertices_local_pL[lvl]++;
    
    vertices_local += vertices_local_pL[lvl];         
  }
  
  descBits_local = vertices_local - vertices_local_pL[grid->maxdepth];
  
#if HTG_SPEED_STATS
  struct timespec start;
  
  if( clock_gettime( CLOCK_REALTIME, &start) == -1 ) {
    perror( "clock gettime" );    
  }
#endif
  
 /**
   * Calculate min,max for xml header (prob. not needed?)
   */
#if HTG_SPEED_STATS       
  struct timespec t_stats_start;  
  if( clock_gettime( CLOCK_REALTIME, &t_stats_start) == -1 ) {
    perror( "clock gettime" );
     
  } 
#endif  
  double min_val[list_len(list)];
  double max_val[list_len(list)];
  {
    int i = 0;
    for (scalar s in list) {
#if HEADER_MIN_MAX_VAL
      stats stat = statsf(s);
      min_val[i]=stat.min;
      max_val[i]=stat.max;
#else
      min_val[i]=0.;
      max_val[i]=0.;
#endif
      i++;
    }
	}
  double min_val_v[vectors_len(vlist)];
  double max_val_v[vectors_len(vlist)];
  { 
    int i = 0;	
    for (vector v in vlist) {
#if HEADER_MIN_MAX_VAL
      min_val_v[i]= HUGE;
      max_val_v[i]= -HUGE;
      foreach_dimension(){
        stats stat = statsf(v.x);	
        min_val_v[i] = min(stat.min,min_val_v[i]);
        max_val_v[i]=  max(stat.max,max_val_v[i]);			
      }
#else
      min_val_v[i]= 0.;
      max_val_v[i]= 0.;
#endif
      i++;
    }
	}
#if HTG_SPEED_STATS    
  struct timespec t_stats_stop;  
  if( clock_gettime( CLOCK_REALTIME, &t_stats_stop) == -1 ) {
    perror( "clock gettime" );
     
  }  
  // to do t_stats_stop    
  double t_stats = ( t_stats_stop.tv_sec - t_stats_start.tv_sec )
                  + (double)( t_stats_stop.tv_nsec - t_stats_start.tv_nsec )
                    / (double)1000000000L;
	
  struct timespec t_header_start;
  if( clock_gettime( CLOCK_REALTIME, &t_header_start) == -1 ) {
    perror( "clock gettime" );
     
  }  
#endif  	
	/** File Header vtkHyperTreeGrid */
	
    //~ WriteHeader(&offset, infoStruct);
  int maj_v = VTK_FILE_VERSION/10, min_v = VTK_FILE_VERSION%10;
  
  fprintf(fp,\
    "<VTKFile %s version=\"%i.%i\" %s %s>\n",\
      "type=\"HyperTreeGrid\"", \
      maj_v, min_v , \
      "byte_order=\"LittleEndian\" ",\
      "header_type=\"UInt32\""); 

#if dimension==2
  fprintf(fp,\
    "\t<HyperTreeGrid BranchFactor=\"2\" TransposedRootIndexing=\"0\" Dimensions=\"%d %d %d\">\n",\
    2,2,1);    
#elif dimension==3
  fprintf(fp,\
    "\t<HyperTreeGrid BranchFactor=\"2\" TransposedRootIndexing=\"0\" Dimensions=\"%d %d %d\">\n", 2,2,2);
#endif
  fprintf(fp,"\t\t<Grid>\n");
#if dimension==2
  fprintf(fp, "\t\t\t<DataArray type=\"Float64\" Name=\"XCoordinates\" NumberOfTuples=\"2\" format=\"ascii\" RangeMin=\"%g\" RangeMax=\"%g\">\n",Y0, Y0+L0);
  fprintf(fp, "\t\t\t\t%g %g",Y0, Y0+L0);
  fprintf(fp,"\n\t\t\t</DataArray>\n");
  fprintf(fp, "\t\t\t<DataArray type=\"Float64\" Name=\"YCoordinates\" NumberOfTuples=\"2\" format=\"ascii\" RangeMin=\"%g\" RangeMax=\"%g\">\n",X0, X0+L0);
  fprintf(fp, "\t\t\t\t%g %g",X0, X0+L0);
  fprintf(fp,"\n\t\t\t</DataArray>\n");	  
  fprintf(fp, "\t\t\t<DataArray type=\"Float64\" Name=\"ZCoordinates\" NumberOfTuples=\"2\" format=\"ascii\" RangeMin=\"%g\" RangeMax=\"%g\">\n",Z0, Z0+L0);		
  fprintf(fp, "\t\t\t\t%g %g", 0., 0.);	  
#elif dimension==3
  fprintf(fp, "\t\t\t<DataArray type=\"Float64\" Name=\"XCoordinates\" NumberOfTuples=\"2\" format=\"ascii\" RangeMin=\"%g\" RangeMax=\"%g\">\n",\
      Z0, Z0+L0);
  fprintf(fp, "\t\t\t\t%g %g",Z0, Z0+L0);
  fprintf(fp,"\n\t\t\t</DataArray>\n");
  fprintf(fp, "\t\t\t<DataArray type=\"Float64\" Name=\"YCoordinates\" NumberOfTuples=\"2\" format=\"ascii\" RangeMin=\"%g\" RangeMax=\"%g\">\n",\
      Y0, Y0+L0);
  fprintf(fp, "\t\t\t\t%g %g",Y0, Y0+L0);
  fprintf(fp,"\n\t\t\t</DataArray>\n");	  
  fprintf(fp, "\t\t\t<DataArray type=\"Float64\" Name=\"ZCoordinates\" NumberOfTuples=\"2\" format=\"ascii\" RangeMin=\"%g\" RangeMax=\"%g\">\n",\
      X0, X0+L0);		
  fprintf(fp, "\t\t\t\t%g %g",X0, X0+L0);	  
#endif
  fprintf(fp,"\n\t\t\t</DataArray>\n");	  
  fprintf(fp,"\t\t</Grid>\n");
  fprintf(fp,"\t\t<Trees>\n");
  
  // Tree Begin
  unsigned int byte_offset = 0;
#if VTK_FILE_VERSION == 10
  fprintf(fp,"\t\t\t<Tree Index=\"0\" NumberOfLevels=\"%d\" NumberOfVertices=\"%u\">\n",\
      grid->maxdepth + 1, vertices_local);   	
  
  // Array Descriptor Bits    
  fprintf(fp,"\t\t\t\t<DataArray type=\"Bit\" Name=\"Descriptor\" NumberOfTuples=\"%u\" format=\"appended\" RangeMin=\"0\" RangeMax=\"1\" offset=\"%u\"/>\n",\
    descBits_local, byte_offset);        
  byte_offset += (descBits_local/8+1) * sizeof(u_int8_t) + sizeof(u_int32_t);
  
  // Array NbVerticesByLevel
  fprintf(fp,"\t\t\t\t<DataArray type=\"Int64\" Name=\"NbVerticesByLevel\" NumberOfTuples=\"%d\" format=\"ascii\" RangeMin=\"1\" RangeMax=\"%u\" >\n\t\t\t\t\t",\
      grid->maxdepth + 1, vertices_local_pL[grid->maxdepth]);
      
  for (int lvl = 0; lvl <= grid->maxdepth; lvl++) {
    fprintf(fp,"%u ", vertices_local_pL[lvl]);
  }		
  fprintf(fp,"\n\t\t\t\t</DataArray>\n");
#endif     

#if VTK_FILE_VERSION == 20
  // Array Descriptor Bits    
  fprintf(fp,"\t\t\t\t<DataArray type=\"Bit\" Name=\"Descriptors\" NumberOfTuples=\"%u\" format=\"appended\" RangeMin=\"0\" RangeMax=\"1\" offset=\"%u\"/>\n",\
    descBits_local, byte_offset);        
  byte_offset += (descBits_local/8+1) * sizeof(u_int8_t) + sizeof(u_int32_t);
  
  fprintf(fp,"\t\t\t\t<DataArray type=\"Int64\" Name=\"NumberOfVerticesPerDepth\" NumberOfTuples=\"%d\" format=\"ascii\" RangeMin=\"1\" RangeMax=\"%u\" >\n\t\t\t\t\t",\
      grid->maxdepth + 1, vertices_local_pL[grid->maxdepth]);
        
  for (int lvl = 0; lvl <= grid->maxdepth; lvl++) {
    fprintf(fp,"%u ", vertices_local_pL[lvl]);
  }		
  fprintf(fp,"\n\t\t\t\t</DataArray>\n");
  
  fprintf(fp,"\t\t\t\t<DataArray type=\"Int64\" Name=\"TreeIds\" NumberOfTuples=\"1\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\" >\n");
  fprintf(fp,"\t\t\t\t\t%i\n",0);
  fprintf(fp,"\t\t\t\t</DataArray>\n");
  
  fprintf(fp,"\t\t\t\t<DataArray type=\"UInt32\" Name=\"DepthPerTree\" NumberOfTuples=\"1\" format=\"ascii\" RangeMin=\"%i\" RangeMax=\"%i\" >\n",\
  grid->maxdepth+1, grid->maxdepth+1);
  fprintf(fp,"\t\t\t\t\t%i\n",grid->maxdepth+1);
  fprintf(fp,"\t\t\t\t</DataArray>\n");
          
  fprintf(fp,"\t\t</Trees>\n");
#endif  
  // Cell Data Begin
  fprintf(fp,"\t\t\t\t<CellData>\n");
  
#if WRITE_HTG_LEVEL    
  fprintf(fp,"\t\t\t\t\t<DataArray type=\"UInt8\" Name=\"Level\" NumberOfTuples=\"%u\" format=\"appended\" RangeMin=\"0\" RangeMax=\"%i\" offset=\"%u\"/>\n",\
      vertices_local, grid->maxdepth, byte_offset);		
  byte_offset += vertices_local * sizeof(u_int8_t) + sizeof(u_int32_t);
#endif    
  {
    int i = 0;
    for (scalar s in list)
    {			
      fprintf(fp,"\t\t\t\t\t<DataArray type=\"Float32\" Name=\"%s\" NumberOfTuples=\"%u\" format=\"appended\" RangeMin=\"%g\" RangeMax=\"%g\" offset=\"%u\"/>\n",\
          s.name, vertices_local, min_val[i], max_val[i], byte_offset);			    
      byte_offset += vertices_local * sizeof(float_t) + sizeof(u_int32_t);
      i++;
    }
  }
  {
    int i = 0;
    for (vector v in vlist)
    {
      char *vname = strtok(v.x.name, ".");
      fprintf(fp, "\t\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"%i\" Name=\"%s\" NumberOfTuples=\"%u\" format=\"appended\"  RangeMin=\"%g\" RangeMax=\"%g\" offset=\"%u\"/>\n",\
          3,vname, vertices_local, min_val_v[i],max_val_v[i],byte_offset);
      byte_offset += vertices_local * 3 * sizeof(float_t) + sizeof(u_int32_t);
      i++;
    } 
  }  
  fprintf(fp,"\t\t\t\t</CellData>\n");
#if VTK_FILE_VERSION == 10
  fprintf(fp,"\t\t\t</Tree>\n\t\t</Trees>\n");
#endif

  fprintf(fp,"\t</HyperTreeGrid>\n\t<AppendedData encoding=\"raw\">\n_");

#if HTG_SPEED_STATS    
  struct timespec t_header_stop;
  if( clock_gettime( CLOCK_REALTIME, &t_header_stop) == -1 ) {
    perror( "clock gettime" );
  } 
  
  double t_header = ( t_header_stop.tv_sec - t_header_start.tv_sec )
                  + (double)( t_header_stop.tv_nsec - t_header_start.tv_nsec )
                    / (double)1000000000L;
#endif
  /** Write Bitmask */
  int cell_size;

	cell_size=sizeof(u_int8_t);    
	  
  int vertices_local_corr = ((descBits_local/8)+1)*8;

  u_int32_t prepend_size = vertices_local_corr;
  fwrite(&prepend_size, sizeof(u_int32_t), 1, fp); 		
  u_int8_t* write_cache = (u_int8_t*)calloc(vertices_local_corr,cell_size);
  long index = 1;
  for(int lvl=0; lvl<grid->maxdepth;++lvl){ // Bitmask ist "0" auf grid->maxdepth	und muss nicht geschrieben werden
		foreach_level(lvl,serial)
			if (is_local(cell)) {
				if(is_leaf(cell)){					
					write_cache[index++] = 0; // ascii: 0					
				} else {					
					write_cache[index++] = 1; // ascii: 1
				}
      }
	}

  // Byte to bits
  for (int i = 0; i < vertices_local_corr/8 ; ++i){
    for ( int j = 0; j < 8; ++j ) {
      write_cache[i] |= write_cache[(8*i+j) + 1] << (7-j);        
      if ((j+1) == 8) {
         write_cache[i+1] = 0; // clear next byte         
      }
    }
  }
  
  fwrite(&write_cache[0], cell_size, vertices_local_corr/8, fp);
  free(write_cache);
  write_cache = NULL;

  /** Write LEVEL */
#if WRITE_HTG_LEVEL
 
 	cell_size=sizeof(u_int8_t);
	
  prepend_size = vertices_local * cell_size;
  fwrite(&prepend_size, sizeof(u_int32_t), 1, fp); 		

	for(int lvl=0; lvl<=grid->maxdepth;++lvl){
		u_int8_t * write_cache = (u_int8_t*) malloc(vertices_local_pL[lvl]*sizeof(u_int8_t));
		long index = 0;
		foreach_level(lvl,serial)
			if (is_local(cell))
				write_cache[index++] = (u_int8_t)lvl;

		fwrite(&write_cache[0], cell_size, vertices_local_pL[lvl], fp);
		free(write_cache);
		write_cache = NULL;
	}
	
#endif // end WRITE_HTG_LEVEL
  
	for (scalar s in list) {
		cell_size=sizeof(float_t);
		  
    u_int32_t prepend_size = vertices_local * cell_size;
    fwrite(&prepend_size, sizeof(u_int32_t), 1, fp); 		  
	  	
		for(int lvl=0; lvl < grid->maxdepth + 1; ++lvl){
      
			// offset per level
			float_t * write_cache = (float_t*) malloc(vertices_local_pL[lvl]*cell_size);
			long index = 0;
      
			foreach_level(lvl,serial)
				if (is_local(cell)) 			
					write_cache[index++] = val(s);
			
			fwrite(&write_cache[0], cell_size, vertices_local_pL[lvl], fp);      
			free(write_cache);
			write_cache = NULL;	        
    }
  }
  
	for (vector v in vlist) {
		cell_size = 3*sizeof(float_t);
		
    u_int32_t prepend_size = vertices_local * cell_size;
    fwrite(&prepend_size, sizeof(u_int32_t), 1, fp); 		  
			
		for(int lvl=0; lvl<=grid->maxdepth;++lvl){
			
			float_t * write_cache = (float_t*) malloc(vertices_local_pL[lvl]*cell_size);
			long index = 0;
			foreach_level(lvl,serial)
				if (is_local(cell)) {
					#if dimension == 2
			
					write_cache[index] = val(v.y);
					write_cache[index+1] = val(v.x);
					write_cache[index+2] = 0.;
					index += 3;
					#elif dimension == 3
			
					write_cache[index] = val(v.z);
					write_cache[index+1] = val(v.y);
					write_cache[index+2] = val(v.x);
					index += 3;
					#endif	
        }
			fwrite(&write_cache[0], cell_size, vertices_local_pL[lvl], fp);
			free(write_cache);
			write_cache = NULL;	
		}    	
	} 

	fprintf(fp, "ENDBINARY\n\t</AppendedData>\n</VTKFile>\n");    
  
#if HTG_SPEED_STATS  
 struct timespec stop;
 if( clock_gettime( CLOCK_REALTIME, &stop) == -1 ) {
    perror( "clock gettime" );    
  }
  
  long offset = ftell(fp);
  
  double t_file = ( stop.tv_sec - start.tv_sec )
                  + (double)( stop.tv_nsec - start.tv_nsec )
                    / (double)1000000000L;
	fprintf (stdout,"Write Time Taken: %lf ms, Throughput: %.2f MB/s, Stats %.2f%%, Header: %.2f%%\n",\
    t_file*1000., (double)offset/t_file/(double)sq(1024),\
    t_stats/t_file*100.,t_header/t_file*100.);
  fflush(stdout);
#endif
  fflush(fp);
  #if defined(_OPENMP)
    omp_set_num_threads(num_omp);
  #endif
  
}
#endif

