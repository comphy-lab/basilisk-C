/**
# Read X-Y-Z formatted binary data.

This function reades thee-dimensional (cubic $N^3$) data in file `fname` at level `dlev` ($N = 2^{\mathrm{dlev}}$) into the scalar field `s`. 

It works well and fast for large data sets and is compatible with MPI or OpenMP on Multigrid or full octrees  Unfortuantely, it is limited to run on a single node. (use `dump()`->`restore()`)

The first function reads in data stored in so-called single precision
*/
int read_xyz_float(char * fname, scalar s, int dlev){
  unsigned long long int size = (1 << (dimension*dlev));
  /**
The MPI parallel strategy requires special attention. We use the marvalous *shared-memory* functionality that is facilitated by modern MPI implementations.   
  */
#if _MPI 
  MPI_Win win;
  float * a;
  /**
  The root allocated an array and a MPI window is created for the other ranks.
  */
  if (pid() == 0){
    MPI_Win_allocate_shared (size*sizeof(float), sizeof(float), MPI_INFO_NULL,
			    MPI_COMM_WORLD, &a, &win);
  }
  else{ // Slaves obtain the location of the pid()=0 allocated array    
    int disp_unit;
    MPI_Aint  ssize; 
    MPI_Win_allocate_shared (0, sizeof(float), MPI_INFO_NULL,
			    MPI_COMM_WORLD, &a, &win);
    MPI_Win_shared_query (win, 0, &ssize, &disp_unit, &a);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  /**
  The root is also tasked with reading *all* the data. Notice that this is quite fast because it reads contiguous data and true parallel IO is, well... 
  */
  if (pid() == 0){ 
    MPI_Win_lock (MPI_LOCK_EXCLUSIVE, 0, MPI_MODE_NOCHECK, win);
    FILE * fp = fopen (fname, "rb");
    fread (a, sizeof(float), size,fp);
    MPI_Win_unlock (0, win);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  /**
  In serial, life is a bit easier.
  */
#else 
  float * a = (float*) malloc (sizeof(float)*size);
  FILE * fp = fopen (fname, "rb");
  fread (a, sizeof(float), size, fp);
#endif
/** We may need to take case of some specifics of MG parralelism and work with an offset. */
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
  unsigned long long int CR = (1 << dlev);
  int o = -BGHOSTS - 1;
  unsigned long long int index;
/** Loading the data itself is now straightforward*/
  foreach(){
    index = ((nxoffset + point.i + o) + (CR*(nyoffset + point.j + o)) + (sq(CR)*(nzoffset + point.k + o)));
    s[] = (double)a[index];
  }
  return 0;
}
/**
This simular function reads double-precision data:
*/
int read_xyz_double (char * fname,scalar s,int dlev){
  unsigned long long int size= (1<<(dimension*dlev));
#if _MPI
  MPI_Win win;
  double * a;
  if (pid() == 0){
    MPI_Win_allocate_shared(size*sizeof(double), sizeof(double), MPI_INFO_NULL,
			    MPI_COMM_WORLD, &a, &win);
  }
  else{ // Slave    
    int disp_unit;
    MPI_Aint  ssize; 
    MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL,
			    MPI_COMM_WORLD, &a, &win);
    MPI_Win_shared_query(win, 0, &ssize, &disp_unit, &a);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if (pid() == 0){
    MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, MPI_MODE_NOCHECK, win);
    FILE * fp = fopen(fname,"rb");
    fread(a,sizeof(double),size,fp);
    MPI_Win_unlock(0,win);
  }
  MPI_Barrier(MPI_COMM_WORLD);
#else //not _MPI, life is easier
  double * a= (double*) malloc(sizeof(double)*size);
  FILE * fp = fopen(fname,"rb");
  fread(a,sizeof(double),size,fp);
#endif
  
#if MULTIGRID_MPI
  int nz = pid() % (mpi_dims[2]);
  int ny = (pid()/mpi_dims[2]) % (mpi_dims[1]);
  int nx = pid()/(mpi_dims[2]*mpi_dims[1]);
  unsigned long long int nxoffset = ((1<<depth())*nx);
  unsigned long long int nyoffset = ((1<<depth())*ny);
  unsigned long long int nzoffset = ((1<<depth())*nz);
#else // Non MG-MPI, no offsets 
  int nxoffset = 0;
  int nyoffset = 0;
  int nzoffset = 0;
#endif
  unsigned long long int CR = (1<<dlev);
  int o = -BGHOSTS-1;
  unsigned long long int ind;
  
  foreach(){
    ind =((nxoffset+point.i+o)+(CR*(nyoffset+point.j+o)) + (sq(CR)*(nzoffset+point.k+o)));
    s[]=a[ind];
  }
  return 0;
}

