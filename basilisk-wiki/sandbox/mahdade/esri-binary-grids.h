/**
## *output_flt()*: ESRI Binary Float Grid (.FLT,.HDR) format

The Esri Binary Float Grid format consists of 2 files: a binary .FLT image file and ASCII .HDR header file with the same file name but different file extension. For example, Test.FLT and Test.HDR [(Ref)](http://surferhelp.goldensoftware.com/subsys/subsys_ESRI_binary_float_grid_file_descr.htm).


The arguments are:

*f*
: a scalar field (compulsory).

*file*
: an array of char, path to the output file (with .FLT extension)

$\Delta$
: size of a grid element. Default is 1/N.

*linear*
: whether to use bilinear or first-order interpolation. Default is 
first-order.

*box*
: the lower-left and upper-right coordinates of the domain to consider.
 Default is the entire domain.

*mask*
: if set, this field will be used to mask out, the regions 
of the domain for which *mask* is negative. */

struct OutputFLT {
  scalar f;
  char * file;
  double Delta;
  bool linear;
  double box[2][2];
  scalar mask;
};

trace
int output_flt (struct OutputFLT p)
{
  FILE * flt , * hdr ;
  scalar input = p.f;
  char * hdrfile, * pch;
  char buffer[100];
  double * data_double;
  // default values
  if (p.box[0][0] == 0. && p.box[0][1] == 0. && 
      p.box[1][0] == 0. && p.box[1][1] == 0.) {
    p.box[0][0] = X0;      p.box[0][1] = Y0;
    p.box[1][0] = X0 + L0; p.box[1][1] = Y0 + L0;
    if (p.Delta == 0) p.Delta = L0/N;
  }

  double Delta = p.Delta;
  int nx = (p.box[1][0] - p.box[0][0])/Delta;
  int ny = (p.box[1][1] - p.box[0][1])/Delta;

  data_double = (double *)malloc(nx*ny*sizeof(double));  
  for(int i=0;i<(nx*ny);i++)
    data_double[i] = nodata;

  // data
  for (int j = ny-1; j >= 0; j--) {
    double yp = Delta*j + p.box[0][1] + Delta/2.;	// center of pixel
    for (int i = 0; i < nx; i++) {
      double xp = Delta*i + p.box[0][0] + Delta/2., v;
      float vf;
      if (p.mask.i) { // masking
	if (p.linear) {
	  double m = interpolate (p.mask, xp, yp);
	  if (m < 0.)
	    v = nodata;
	  else
    	    v = interpolate (p.f, xp, yp);
	}
	else {
	  Point point = locate (xp, yp);
	  if (point.level < 0 || val(p.mask) < 0.)
	    v = nodata;
	  else
	    v = val(p.f);
	}
      }
      else if (p.linear)
	v = interpolate (p.f, xp, yp);
      else {
	Point point = locate (xp, yp);
	v = point.level >= 0 ? val(p.f) : nodata;
      }
      if (v == nodata)
	data_double[i+(ny-1-j)*nx] = -9999.;
      else
	data_double[i+(ny-1-j)*nx] = v;

    }
  }

    if (pid() == 0) { // master writes .FLT and .HDR files
@if _MPI
    MPI_Reduce (MPI_IN_PLACE, data_double, nx*ny, MPI_DOUBLE, MPI_MAX, 0,
		MPI_COMM_WORLD);
@endif

// Set path to header file
  hdrfile = (char *)malloc(sizeof(char)*strlen(p.file));
  strcpy(hdrfile,p.file);
  
  pch = hdrfile + strlen(hdrfile)-3;
  strncpy (pch,"hdr",3);
 
  if(!(hdr = fopen (hdrfile, "w")))
  {
    printf("Failed to open header file %s\n",hdrfile);
    return -1;
  }

  // header
  fprintf (hdr, "ncols          %d\n", nx);
  fprintf (hdr, "nrows          %d\n", ny);
  fprintf (hdr, "xllcorner      %.8g\n", p.box[0][0]);
  fprintf (hdr, "yllcorner      %.8g\n", p.box[0][1]);
  fprintf (hdr, "cellsize       %.8g\n", Delta);
  fprintf (hdr, "nodata_value   -9999\n");
  fprintf (hdr, "byteorder   LSBFIRST\n");

  fflush(hdr);
  fclose(hdr);

  // open flt file

  bool opened = false;
  if( !(flt = fopen (p.file, "wb")) ) {
      perror (p.file);
      exit(1);
  }
  else
      opened = true;

  for(int i=0;i<(nx*ny);i++)
  {
    float vf = (float)data_double[i];
    fwrite(&vf,1,sizeof(float),flt);
  }
  fclose(flt);
}
@if _MPI
  else // slave does not write anything
    MPI_Reduce (data_double, NULL, nx*ny, MPI_DOUBLE, MPI_MAX, 0,
		MPI_COMM_WORLD);
@endif

  free(data_double);
  return(0);

}

int SYSTEM_IS_LITTLE_ENDIAN(void)
{
	// First, we declare a 16-bit integer (short int), which has the value 0x0001
        short int word = 0x0001;
	// then we get a pointer and dereference it
        char *b = (char *)&word;
	// If the LSB is stored at lower address (e.g. the value that pointer points to), then system is LSBFIRST i.e. Litte Endian
        return (b[0] ? 1 : 0);
}

float ReverseFloat( float inFloat )
{
   float retVal;
   char *floatToConvert = ( char* ) & inFloat;
   char *returnFloat = ( char* ) & retVal;

   // swap the bytes into a temporary buffer
   returnFloat[0] = floatToConvert[3];
   returnFloat[1] = floatToConvert[2];
   returnFloat[2] = floatToConvert[1];
   returnFloat[3] = floatToConvert[0];

   return retVal;
}

void lower_string(char s[]) {
   int c = 0;
   
   while (s[c] != '\0') {
      if (s[c] >= 'A' && s[c] <= 'Z') {
         s[c] = s[c] + 32;
      }
      c++;
   }
}

/**
# *input_flt()*: ESRI Binary Float Grid (.FLT,.HDR) format

The Esri Binary Float Grid format consists of 2 files: a binary .FLT image file and ASCII .HDR header file with the same file name but different file extension. For example, Test.FLT and Test.HDR [(Ref)](http://surferhelp.goldensoftware.com/subsys/subsys_ESRI_binary_float_grid_file_descr.htm).

This is the reciprocal function of [*output_flt()*](esri-binary-grids.h#output_flt).

The arguments and their default values are:

*s* 
: the scalar where the data will be stored. No default value. You
must specify this parameter

*file*
: the name of the file to read from.
The path to the header file is assumed to be the same as the binary file with the extension ".hdr"

*nodatavalue*
: the value of the NoDataValue. Default is the same as that defined in
the raster file. 

*linear*
: if true, the raster data is bilinearly interpolated. Default is false.
*/

struct InputFLT {
  scalar s;
  char * file;
  double nodatavalue;
  bool linear;
};

int input_flt (struct InputFLT p)
{
  FILE * flt , * hdr ;
  scalar input = p.s;
  char * hdrfile, * pch;
  char buffer[100];

  double DeltaGRD;
  int nx, ny,i,j;
  double XG0, YG0, ndv,center,dx,dy;
  bool lsbfirst,reverse_bytes;
  float f;

// Get path to header file
  hdrfile = (char *)malloc(sizeof(char)*strlen(p.file));
  strcpy(hdrfile,p.file);
  
  pch = hdrfile + strlen(hdrfile)-3;
  strncpy (pch,"hdr",3);
 
  if(!(hdr = fopen (hdrfile, "r")))
  {
    printf("Failed to open header file %s\n",hdrfile);
    return -1;
  }

// Parse header file 

  while(fscanf(hdr,"%s",buffer)!=EOF)
  {
    lower_string(buffer);

    if(!(strcmp(buffer,"ncols")))
      fscanf(hdr,"%d",& nx);

    if(!(strcmp(buffer,"nrows")))
      fscanf(hdr,"%d",& ny);

    if(!(strcmp(buffer,"cellsize")))
      fscanf(hdr,"%lf",& DeltaGRD);

    if(!(strcmp(buffer,"xllcenter")))
    {
      fscanf(hdr,"%lf",& XG0);
      center=1.0;
    }
 
    if(!(strcmp(buffer,"yllcenter")))
    {
      fscanf(hdr,"%lf",& YG0);
      center=1.0;
    }

    if(!(strcmp(buffer,"xllcorner")))
    {
      fscanf(hdr,"%lf",& XG0);
      center=0.0;
    }
 
    if(!(strcmp(buffer,"yllcorner")))
    {
      fscanf(hdr,"%lf",& YG0);
      center=0.0;
    }

    if(!(strcmp(buffer,"nodata_value")))
      fscanf(hdr,"%lf",& ndv);
   
    if(!(strcmp(buffer,"byteorder")))
    {
       fscanf(hdr,"%s",buffer);
       lower_string(buffer);
       lsbfirst = !strcmp(buffer,"lsbfirst");       		
    }
  }

  XG0 -= 0.5*DeltaGRD*center;
  YG0 -= 0.5*DeltaGRD*center;
   
  fclose(hdr);

  bool opened = false;
  if( !(flt = fopen (p.file, "rb")) ) {
      perror (p.file);
      exit(1);
  }
  else
      opened = true;
  
  //default value of NoData value
  if (!p.nodatavalue)
    p.nodatavalue = ndv;

  // Allocation of pointers to columns' start (size : nx)
  double ** value ;
  if( !(value = (double **) malloc(nx*sizeof(double *)) ) )
        return(-1);

  // Allocation of each column (each of size ny)
  for(i=0;i<nx;i++)
  {
    if ( !(value[i] = (double *) malloc(ny*sizeof(double)) ) )
      return(-1);
  }

  // read the data (float to double conversion and byte reverse if necessary)

  if ( ( lsbfirst && !SYSTEM_IS_LITTLE_ENDIAN() ) |
       ( !lsbfirst && SYSTEM_IS_LITTLE_ENDIAN() )      )
        reverse_bytes = true;
  else
        reverse_bytes = false;

  for(j = ny-1; j >= 0; j--){
  for(i = 0 ; i < nx; i++){
    fread(&f,sizeof(float),1,flt);
    if (reverse_bytes)
      f = ReverseFloat(f);
    value[i][j] = f ;
  }}
  
  printf(" %d values read\n",nx*ny);

  double LGx0 = nx*DeltaGRD;
  double LGy0 = ny*DeltaGRD;
  bool warning = false;
  bool internal,incl;
  double val;
  double onehalf=1./2.;
  double xllcenter = XG0+DeltaGRD/2., yllcenter = YG0+DeltaGRD/2.;

  foreach() {

    dx = (x - xllcenter)/DeltaGRD ; // relative offset in x w.r.t. center of LL pixel
    i = (int)floor(dx);               

    dy = (y - yllcenter)/DeltaGRD ; // relative offset in y w.r.t. center of LL pixel
    j = (int)floor(dy);             

    internal = (i >= 0 &&
	        i < nx-1 &&
		j >= 0 &&
		j < ny-1 );
    
    if (p.linear && internal ) { // bi-linear interpolation

      dx -= i;
      dy -= j;

      val = (1.-dx)*(1.-dy)*value[i][j]
	  + dx*(1.-dy)*value[i+1][j]
	  + (1.-dx)*dy*value[i][j+1]
          + dx*dy*value[i+1][j+1];
    }
    else {
    // Test if the point in the Basilisk area is included in the raster area 
      incl = (x >= XG0 &&
	      x <= XG0 + LGx0 &&
	      y >= YG0 &&
	      y <= YG0 + LGy0 );      
      if(incl) {
      // snap relative offsets dx and dy to nearest integers
        i = (int)(dx+onehalf);
	j = (int)(dy+onehalf);
        val = value[i][j];
      }
      else {
        val = p.nodatavalue;
        warning = true;
      }
      
    }

    input[] = val ;
  }
  
  // Deallocation of each column
  for(i=0;i<nx;i++)
   free(value[i]);

  // Deallocation of pointers to columns' start
  free(value);
  
  if (warning)
    fprintf (stderr,
	     "input_flt(): Warning: Raster data is not covering all"
	     " the simulation area\n");

  fclose (flt);
  return(0);
}
