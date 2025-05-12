/**
#Read $\mu$-HH-formatted binary files into the Basilisk octree grid
This rather specefic function reads the binary files as outputted by the [$\mu$-HH](www.microhh.org) code into the Basilisk octree grid. Special care is taken to ensure that this function is parallized at read time when using MPI domain decomposition. This page discusses the steps taken to achieve this, and may serve others who wish to set-up a data input function themselves.  I have not considered the behaviour for the OpenMP-style parallelization. Since it would require special attention, the usage with OpenMP will not work. 

## Manual
This function assumes that there is an obvious 1:1 mapping of the Octree-grid cells to the micro-HH data points. Meaning that the octree grid in Basilisk has to be full at level $l$, and the Micro-HH output has a cubic, $2^l\times 2^l \times 2^l$ grid-structure. This is not such bad idea anyway, considering [FFTW's preferencies](http://www.fftw.org/speed/). There are two ways of dealing with the non-conforming [N-order MPI decomposition](The_Tree_Grid_Structure_in_Basilsisk) 'isue'. One approach exploits the special case where the number of MPI-threads is equal to a power of 8 (e.g. 1, 8, 64, 512, 4096). This way a cubic-block-structure decomposition is achieved. Another method is just to let the file pointer jump to a suitable location in the file.

## Helper function
First I need to remember myselft that this function should be upgraded in the future to deal with $2^l \times 2^l \times X$ Micro-HH grid data aswell. For the case where $2^l>X$, we need to fill in the remaining cells in the cubic octree. This function does not do much for now. 
*/
int fill_empty_cells(){
  return 1;
}

/**
## 1: the `readmuhh()` function
The shining star of this page is the `readmuhh()` function. You should provide it with a scalar field that you want to fill with the data that is stored in a file named `fname` char. Furthermore, it seems a good idea to specify the aforementioned level $l$ that dictates the grid resolution. As mentioned, when using MPI, this only works when the number of threads is a power of eight.   
*/

int readmuhh(scalar s,char * fname, int dlevel){
  /**
  First we check if the number of threads is a power of 8.     
  */
  double po8= log((double)npe())/log(8.);
  int ipo8 = (int)(po8+0.5); //Rounding error correction. 
  if (fabs(po8-(double)(ipo8))>1e-10){
    fprintf(stderr,"Error:\nThe reading of mu-hh files does not work with %d threads.\nUse a power of 8 in stead.\n",npe());
    return 2; //Error code 2;
  }
  /**
  The function assumes data to be stored in doubles. I think Micro-HH can also dump floats, that is for another day. 
  */
  int siz=8;
  int starti=1e8,endi=-(1e8),startj=1e8,endj=-(1e8),startk=1e8,endk=-(1e8);
  /**
  An over desinged piece of code finds the $\{x,y,z\}$ range within the domain that is alloted to each thread. 
  */
  foreach(){ 
    if (starti>point.i) starti=point.i;
    if (endi<point.i)   endi=point.i;
    if (startj>point.j) startj=point.j;
    if (endj<point.j)   endj=point.j;
    if (startk>point.k) startk=point.k;
    if (endk<point.k)   endk=point.k;
  }
#if !_MPI
  /**
  The file specified by `fname` is opened and the function exits by returning a non-zero error code and acompagnying message to `stderr` if this step has not gone as expected. Furthermore, a hint at the possible cause is provided by repeating the input. 
  */
  FILE * fpm;
  if (!(fpm = fopen(fname,"rb"))){
    fprintf(stderr,"Cound not open file with name '%s'\n",fname);
    return 1;
  }
    /**
When applying the MPI-domain decomposition we thinks a bit differently: We calculate the number of threads in each direction and declare and initialize *individual* file pointers.  
  */
#elif _MPI
  int nd = pow(2,ipo8);
  MPI_File fpm;
  MPI_File_open(MPI_COMM_WORLD,fname,MPI_MODE_RDONLY,MPI_INFO_NULL,&fpm);
  /**
  Each processor jumps its pointer to the location in the file that corresponds to the starting coordinates specified by $\{$`starti,startj,startk`$\}$. 
  */
  int o = -1-BGHOSTS;
  int set = (startk+o) + (N*(starti+o)) + (sq(N)*(startj+o)); 
  MPI_File_seek(fpm,set*siz,MPI_SEEK_CUR);
#endif
  /**
  Instead of using the trusty `foreach()` iterator, each thread loops over the cells within its range in a typical cartesian fashion, corresponding to the Micro-HH data format. Notice that we rotate the axis here.   
  */
  for (int j=startj;j<=endj;j++){      // z-dir in mu-HH -> y-dir in B
    for (int i=starti;i<=endi;i++){    // y-dir in mu-HH -> x-dir in B
      for (int k=startk;k<=endk;k++){  // x-dir in mu-HH -> z-dir in B
	/**
In grid units, $\{i,j,k\}$ represents the cartesian coordinates. We define a Point structure that points to the corresponding cell in the octree grid. 
	 */
	Point point;
	point.i=i; point.j=j; point.k=k;
	point.level=dlevel;
	/**
The algorithm is set-up such that we can now conviniently read a field value into the scalar field. We distinguish between the usage of a regular file pointer and the MPI-enabled individual version.
	 */
#if !_MPI
	fread(&s[],siz,1,fpm);
#elif _MPI
	MPI_Status status;
	MPI_File_read(fpm,&s[],1,MPI_DOUBLE,&status);
#endif
	/**
	   The parallel strategy also requires the file pointer of the various threads to make some jumps whilst reading the file.  
	*/
      }
#if _MPI
      MPI_File_seek(fpm,(nd-1)*N/nd*siz,MPI_SEEK_CUR);
#endif
    }
#if _MPI
    MPI_File_seek(fpm,N*(N*(nd-1))/nd*siz,MPI_SEEK_CUR);
#endif
  }
  /**
    As discussed above, we should fill the remaining cells in the octree with relevant values. This functionality is not implemented at the moment. So only produce and cubic data.    
  */
  fill_empty_cells();
  /**
We celebrate arriving at this point by returning an `int 0` error code.
   */
  return 0;
}
/**
## 2: the `readmuhh_with_any_number_of_threads()` function
This function does the same as above but also works with any number of threads. The price to pay is that for a grid with $N^3$ gridpoints, the file pointer makes $O(N^3)$ jumps, irrespective of the number of processors. Therefore, the file pionter makes more jumps compared to using `readmuhh`, but it may display different scaling properties when using more threads (better). Also it iterates the grid in a sequence the Basilisk data structure prefers.      
*/
int readmuhh_with_any_nr_of_threads(scalar s,char * fname, int dlevel){
#if !_MPI
  FILE * fpm;
  if (!(fpm = fopen(fname,"rb"))){
    fprintf(stderr,"Cound not open file with name '%s'\n",fname);
    return 1;
  }
#elif _MPI
  MPI_File fpm;
  MPI_File_open(MPI_COMM_WORLD,fname,MPI_MODE_RDONLY,MPI_INFO_NULL,&fpm);
#endif
  int siz=8;
  int o = -1-BGHOSTS;
  int g=0;
  /**
  Here we use the basilisk cell iterator sequence. Foreach data point, the file pointer is set to point at the desired location in the file.   
  */
  int CR = (1<<dlevel);
  foreach(){
    g = (point.k+o) + CR*(point.i+o) + sq(CR)*(point.j+o);
#if !_MPI
    fseek(fpm,g*siz,SEEK_SET);
    fread(&s[],siz,1,fpm);
#elif _MPI
    MPI_File_seek(fpm,g*siz,MPI_SEEK_SET);
    MPI_Status status;
    MPI_File_read(fpm,&s[],1,MPI_DOUBLE,&status);
#endif
  }
  fill_empty_cells();
  return 0;
}
/**
## Write a $\mu$-HH file. 
We also define a mirrored version that instead of reading a $\mu$-HH data file into Basilisk, writes Basilisk data into a $\mu$-HH$ start file. 
*/
int writemuhh_with_any_nr_of_threads(scalar s,char * fname, int dlevel){
#if !_MPI
  FILE * fpm;
  fpm = fopen(fname,"wb");
#elif _MPI
  MPI_File fpm;
  MPI_File_open(MPI_COMM_WORLD,fname,MPI_MODE_WRONLY,MPI_INFO_NULL,&fpm);
#endif
  int siz=8;
  int o = -1-BGHOSTS;
  int g=0;
  int CR = (1<<dlevel);
  foreach(){
    g = (point.k+o) + CR*(point.i+o) + sq(CR)*(point.j+o);
#if !_MPI
    fseek(fpm,g*siz,SEEK_SET);
    fwrite(&s[],siz,1,fpm);
#elif _MPI
    MPI_File_seek(fpm,g*siz,MPI_SEEK_SET);
    MPI_Status status;
    MPI_File_write(fpm,&s[],1,MPI_DOUBLE,&status);
#endif
  }
  if (pid()==0)
    fprintf(stdout,"Created a file with name %s\n",fname);
  return 0;
}

/**
## Read and dump
We also define a more user-friendly interfacing function that applies `readmuhh` for the usual suspect fields; $u_x,u_y,u_z$ and $\theta$. The idea is that you can read the data using the required number of threads and as a next step you may restore using any number of processes. You should provide it with a name for the file dumpfile, a file identifier *t*, the level of refinement, a vector field for the cell-centred volocity components, and a scalar field for the thermodynamic variable *th*. 

The velocity components in micro-HH are defined at faces, we will translate it into cell-centered values using a `low-storage' scheme with only a single scratch field. Notice that we also apply the rotation of coordinates for the velocity vector.  
*/

int dump_from_muhh(char * dname, double t, int level, vector u,scalar th){
  char name[100];
  sprintf(name,"v.%07g",t);
  readmuhh(u.x, name, level);
  sprintf(name,"w.%07g",t);
  readmuhh(u.y, name, level);
  sprintf(name,"u.%07g",t);
  readmuhh(u.z, name, level);
  sprintf(name,"th.%07g",t);
  readmuhh(th, name, level);
  dump(dname);
  return 1;
}



/**
## Performance.
Having tested the functions for importing a $64^3$ dataset: On a single thread, without MPI, the `readmuhh()` function is about 35 times faster than the `readmuhh_with_any_nr_of_threads()` alternative. With 8 threads, on a hyperthreaded quadcore CPU, the former outperforms the latter, because then the file pointer(s) in the `readmuhh` function are also making jumps now.  

###Performance data:
reading a $64^3$ file with doubles. A single field 5 times over.

#### 1 thread, no MPI:
`Readmuhh()`  : 21 msec  (!)  
`Any_nr_of()` : 721 msec  

#### 1 thread, with MPI:
`Readmuhh()`  : 603 msec ($\leftarrow$ Big MPI performance hit on a single thread)   
`any_nr_of()` : 904 msec

#### 8 threads
`Readmuhh()`  :117,117,117,117,117,114,117,118 msec ($\leftarrow$ Speed up!)   
`any_nr_of()` : 165,165,165,165,165,165,165,165 msec ($\leftarrow$ OK scaling)  

#### 64 threads (Something may be wrong with the way the time is diagnosed?)
`Readmuhh()`  : $\approx 20 $ msec   
`any_nr_of()` : $\approx 30$ msec 

Winner: `Readmuhh()`. But your performance may vary. 
*/