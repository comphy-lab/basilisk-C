/**
# Input functions

## *input_matrix()*: reads a matrix from a single precision binary file 

This function reads a matrix from a gnuplot-compatible binary file and stores it
in a scalar field. The binary file format is expected to be:

```
 <N+1>  <y0>   <y1>   <y2>  ...  <yN>
 <x0> <z0,0> <z0,1> <z0,2> ... <z0,N>
 <x1> <z1,0> <z1,1> <z1,2> ... <z1,N>
```

The matrix values are then assigned to the scalar field `s` based on the
coordinates and width provided.

The arguments and their default values are:

*s*
: scalar field where the matrix values will be assigned.

*fp*
: file pointer to the binary file.

*n*
: size of the matrix. Default is *N*.

*ox*
: x-coordinate offset. Default is *X0*.

*oz*
: z-coordinate offset. Default is *Z0*.

*width*
: width for mapping the matrix to the scalar field. Default is *L0*.

### Example Usage

```c
FILE *fp = fopen("matrix_single.bin", "rb");
if (fp == NULL) {
    perror("Failed to open file");
    exit(EXIT_FAILURE);
}
input_matrix(s, fp);
fclose(fp);
```
*/
void input_matrix(scalar s, FILE * fp, int n = N, double ox=X0, double oz=Z0, double width=L0){

  float nfile;
  int read_size;
  NOT_UNUSED(read_size);

  // Read the matrix size from file
  read_size = fread(&nfile, sizeof(float), 1, fp);
  n = (int)nfile;                                 

  // Allocate memory for coordinates and matrix
  float yp[n], xp[n];                                       
  float **v = matrix_new(n, n, sizeof(float));              

  // Read y coordinates from file
  read_size = fread(&yp, sizeof(float), n, fp);             
  for (int i = 0; i < n; i++){
    // Read x coordinate for each row
    read_size = fread(&xp[i], sizeof(float), 1, fp);        
    for (int j = 0; j < n; j++){
      // Read matrix values
      read_size = fread(&v[i][j], sizeof(float), 1, fp);    
    }
  }

  // Loop over the domain and assign matrix values to the scalar field
  foreach () {
    int i = (x - ox) * n / width;
#if dimension == 3
    int j = (z - oz) * n / width;
#elif dimension == 2
    int j = (y - oz) * n / width;
#endif
    if (i >= 0 && i < n && j >= 0 && j < n){
      s[] = v[i][j];
    }
    else{
      s[] = 0.0;
    }
  }
  matrix_free(v);
}

/** 
## *input_matrix_double()*: reads a matrix from a double precision binary file 

This function reads a matrix from a gnuplot-compatible binary file and stores it
in a scalar field. The binary file format is expected to be:

```
 <N+1>  <y0>   <y1>   <y2>  ...  <yN>
 <x0> <z0,0> <z0,1> <z0,2> ... <z0,N>
 <x1> <z1,0> <z1,1> <z1,2> ... <z1,N>
```

The matrix values are then assigned to the scalar field `s` based on the
coordinates and width provided.


The arguments and their default values are:

*s*
: scalar field where the matrix values will be assigned.

*fp*
: file pointer to the binary file.

*n*
: size of the matrix. Default is *N*.

*ox*
: x-coordinate offset. Default is *X0*.

*oz*
: z-coordinate offset. Default is *Z0*.

*width*
: width for mapping the matrix to the scalar field. Default is *L0*.

### Example Usage

```c
FILE *fp = fopen("matrix_double.bin", "rb");
if (fp == NULL) {
    perror("Failed to open file");
    exit(EXIT_FAILURE);
}
input_matrix_double(s, fp);
fclose(fp);
```

*/
void input_matrix_double(scalar s, FILE * fp, int n = N, double ox=X0, double oz=Z0, double width=L0){

  double nfile;
  int read_size;
  NOT_UNUSED(read_size);
  
  double *yp = NULL, *xp = NULL;
  double **v = NULL;

  if (pid() == 0) {
    // Process 0 reads the matrix size from the file
    read_size = fread(&nfile, sizeof(double), 1, fp);           
    n = (int)nfile;                                             

    // Allocate memory for y, x coordinates and matrix
    yp = (double *)malloc(n * sizeof(double));
    xp = (double *)malloc(n * sizeof(double));
    v = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++) {
      v[i] = (double *)malloc(n * sizeof(double));
    }

    // Read y coordinates from file
    read_size = fread(yp, sizeof(double), n, fp);              
    for (int i = 0; i < n; i++){
      // Read x coordinate for each row
      read_size = fread(&xp[i], sizeof(double), 1, fp);         
      for (int j = 0; j < n; j++){
        // Read matrix values
        read_size = fread(&v[i][j], sizeof(double), 1, fp);     
      }
    }
  }

#if _MPI
  /** Broadcast matrix size and allocate memory on all processes */ 
  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (pid() != 0) {
    yp = (double *)malloc(n * sizeof(double));
    xp = (double *)malloc(n * sizeof(double));
    v = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++) {
      v[i] = (double *)malloc(n * sizeof(double));
    }
  }

  /** Broadcast y and x coordinates */ 
  MPI_Bcast(yp, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(xp, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  /** Broadcast matrix values */ 
  for (int i = 0; i < n; i++) {
    MPI_Bcast(v[i], n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
#endif

  /** Loop over the domain and assign matrix values to the scalar field */ 
  foreach () {
    int i = (x - ox) * n / width;
#if dimension == 3
    int j = (z - oz) * n / width;
#elif dimension == 2
    int j = (y - oz) * n / width;
#endif
    if (i >= 0 && i < n && j >= 0 && j < n){
      s[] = v[i][j];
    }
    else{
      s[] = 0.0;
    } 
  }

  if (yp) free(yp);
  if (xp) free(xp);
  for (int i = 0; i < n; i++) {
    free(v[i]);
  }
  free(v);  
}

/** 
## *read_matrix()*: Reads and stores a matrix from a double precision binary file

This function constructs a filename using the provided prefix and suffix, opens the corresponding binary file, and reads the matrix data into a scalar field using the `input_matrix_double` function. The binary file is expected to contain double precision matrix data.

The arguments and their default values are:

*prefix*
: prefix for the filename.

*suffix*
: suffix for the filename.

*s*
: scalar field where the matrix values will be assigned.

### Example Usage

This function reads a file `test_input_matrix.bin` and stores the content into 
a scalar field `s`.
```c
read_matrix("test_", "input_matrix", s);
```
see, also [example use](test_input_matrix.c)

*/
void read_matrix(const char *prefix, const char *suffix, scalar s) {
  char filename[1030];
  snprintf(filename, sizeof(filename), "%s%s.bin", prefix, suffix);
  FILE *fp = fopen(filename, "r");
  if (!fp){
    printf("Binary file %s not found\n", filename);
    exit(EXIT_FAILURE);
  }
  input_matrix_double(s, fp);
  fclose(fp);
}
