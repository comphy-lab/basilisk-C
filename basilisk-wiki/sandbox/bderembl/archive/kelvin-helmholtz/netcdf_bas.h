/**
   Netcdf interface for basilisk
*/

#include <stdio.h>
#include <string.h>
#include <netcdf.h>

#define NDIMS 4
#define Y_NAME "y"
#define X_NAME "x"
#define REC_NAME "time"
#define LVL_NAME "level"

/* /\* For the units attributes. *\/ */
/* #define UNITS "units" */
/* #define PRES_UNITS "hPa" */
/* #define TEMP_UNITS "celsius" */
/* #define MAX_ATT_LEN 80 */

/* Handle errors by printing an error message and exiting with a
 * non-zero status. */
int nc_err;
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); return;}

/* IDs for the netCDF file, dimensions, and variables. */
int ncid;

// temporary
int nc_varid[1000];
char * nc_varname[1000];
int nvarout = 0;
int nc_rec = -1;
double loctime_nc = -1;
int rec_varid;

void create_nc(scalar * scalar_list, char* file_nc, int nl)
{
    if (pid() == 0) { // master

   /* LOCAL IDs for the netCDF file, dimensions, and variables. */
   int x_dimid, y_dimid, lvl_dimid, rec_dimid;
   int y_varid, x_varid;
   int dimids[NDIMS];
   
   /* Create the file. */
   if ((nc_err = nc_create(file_nc, NC_CLOBBER, &ncid)))
      ERR(nc_err);

   /* Define the dimensions. The record dimension is defined to have
    * unlimited length - it can grow as needed. In this example it is
    * the time dimension.*/
   if ((nc_err = nc_def_dim(ncid, LVL_NAME, nl, &lvl_dimid)))
      ERR(nc_err);
   if ((nc_err = nc_def_dim(ncid, Y_NAME, N, &y_dimid)))
      ERR(nc_err);
   if ((nc_err = nc_def_dim(ncid, X_NAME, N, &x_dimid)))
      ERR(nc_err);
   if ((nc_err = nc_def_dim(ncid, REC_NAME, NC_UNLIMITED, &rec_dimid)))
      ERR(nc_err);

   /* Define the coordinate variables. We will only define coordinate
      variables for lat and lon.  Ordinarily we would need to provide
      an array of dimension IDs for each variable's dimensions, but
      since coordinate variables only have one dimension, we can
      simply provide the address of that dimension ID (&y_dimid) and
      similarly for (&x_dimid). */
   if ((nc_err = nc_def_var(ncid, Y_NAME, NC_FLOAT, 1, &y_dimid, 
			    &y_varid)))
      ERR(nc_err);
   if ((nc_err = nc_def_var(ncid, X_NAME, NC_FLOAT, 1, &x_dimid, 
			    &x_varid)))
      ERR(nc_err);
  if ((nc_err = nc_def_var(ncid, REC_NAME, NC_FLOAT, 1, &rec_dimid, 
                           &rec_varid)))
    ERR(nc_err);
   /* The dimids array is used to pass the dimids of the dimensions of
      the netCDF variables. Both of the netCDF variables we are
      creating share the same four dimensions. In C, the
      unlimited dimension must come first on the list of dimids. */
   dimids[0] = rec_dimid;
   dimids[1] = lvl_dimid;
   dimids[2] = y_dimid;
   dimids[3] = x_dimid;

   /* Define the netCDF variables */
   char * str1 = "VariableName";

  for (scalar s in scalar_list){
          if (strcmp(str1,s.name) != 0) {
   
       if ((nc_err = nc_def_var(ncid, s.name, NC_FLOAT, NDIMS,
                                dimids, &nc_varid[nvarout])))
         ERR(nc_err);
       nc_varname[nvarout] = strdup(s.name);
       nvarout += 1;
       str1 = strdup(s.name);
     }
   }
   
   /* /\* Assign units attributes to the netCDF variables. *\/ */
   /* if ((nc_err = nc_put_att_text(ncid, pres_varid, UNITS,  */
   /*      			 strlen(PRES_UNITS), PRES_UNITS))) */
   /*    ERR(nc_err); */
   /* if ((nc_err = nc_put_att_text(ncid, temp_varid, UNITS,  */
   /*      			 strlen(TEMP_UNITS), TEMP_UNITS))) */
   /*    ERR(nc_err); */

   /* End define mode. */
   if ((nc_err = nc_enddef(ncid)))
      ERR(nc_err);

   /*  write coordinates*/
   float yc[N], xc[N];
   double Delta = L0*1.0/N;
   for (int i = 0; i < N; i++){
      yc[i] = Y0 + (i + 0.5)*Delta;
      xc[i] = X0 + (i + 0.5)*Delta;
   }

   if ((nc_err = nc_put_var_float(ncid, y_varid, &yc[0])))
      ERR(nc_err);
   if ((nc_err = nc_put_var_float(ncid, x_varid, &xc[0])))
      ERR(nc_err);

   /* Close the file. */
   if ((nc_err = nc_close(ncid)))
      ERR(nc_err);
   /* printf("*** SUCCESS creating example file %s!\n", file_nc); */
    }
}


struct OutputNetcdf {
  scalar * scalar_list;
  char * file_nc;
  int nl;
  double loctime;
  int n; 
  bool linear;
};

void write_nc(struct OutputNetcdf p) {
  if (p.n == 0) p.n = N;
  if (p.nl == 0) p.nl = 1;

  /**
     these line in case there are multiple calls to write_nc for the same record
   */
  if (loctime_nc != p.loctime){
    loctime_nc = p.loctime;
    nc_rec += 1;
  }
  
if (pid() == 0) {
  /* open file. */
  if ((nc_err = nc_open(p.file_nc, NC_WRITE, &ncid)))
    ERR(nc_err);

/*   /\* open file. *\/ */
/* #if _MPI */
/*   if ((nc_err = nc_open_par(file_nc, NC_WRITE, MPI_COMM_WORLD, MPI_INFO_NULL, &ncid))) */
/* #else */
/*   if ((nc_err = nc_open(file_nc, NC_WRITE, &ncid))) */
/* #endif */
/*     ERR(nc_err); */
  
/*   nc_rec_part += 1; */

/* #if _MPI */
/*   if (nc_err = nc_var_par_access(ncid, rec_varid_part, NC_COLLECTIVE)) */
/*     ERR(nc_err); */
/*   if (nc_err = nc_var_par_access(ncid, nc_varid_part, NC_COLLECTIVE)) */
/*     ERR(nc_err); */
/* #endif */

  size_t start[NDIMS], count[NDIMS];
  start[0] = nc_rec; //time
  count[0] = 1;
  
  float loctime = loctime_nc;
  if ((nc_err = nc_put_vara_float(ncid, rec_varid, start, count,
                                  &loctime)))
    ERR(nc_err);
 }


/*     start[0] = nc_rec_part; //time */
/*     start[1] =  pid(); */
/*     start[2] = 0;      //vars */
      
/*     count[0] = 1; */
/*     count[1] = 1; */
/*     count[2] = 1; */

/*     if ((nc_err = nc_put_vara_float(ncid, nc_varid_part, start, count, */
/*                                      &loctime))) */
/*        ERR(nc_err); */

/* #if _MPI */
/*   if (nc_err = nc_var_par_access(ncid, rec_varid_part, NC_INDEPENDENT)) */
/*     ERR(nc_err); */
/*   if (nc_err = nc_var_par_access(ncid, nc_varid_part, NC_INDEPENDENT)) */
/*     ERR(nc_err); */
/* #endif */

  float fn = p.n, Delta = L0/fn;
  float ** field = matrix_new (p.n, p.n, sizeof(float));
  
  
  /* The start and count arrays will tell the netCDF library where to
     write our data. */
  size_t start[NDIMS], count[NDIMS];
  
  
  /* These settings tell netcdf to write one timestep of data. (The
     setting of start[0] inside the loop below tells netCDF which
     timestep to write.) */
  start[0] = nc_rec; //time
  start[1] = -1;     //level
  start[2] = 0;      //y
  start[3] = 0;      //x
  
  count[0] = 1;
  count[1] = 1;
  count[2] = N;
  count[3] = N;
  
  
  /* char * str1; */
  for (scalar s in p.scalar_list){
    
    for (int j = 0; j < p.n; j++) {
      float yp = Delta*j + Y0 + Delta/2.;
      for (int i = 0; i < p.n; i++) {
        float xp = Delta*i + X0 + Delta/2.;
        if (p.linear) {
          field[j][i] = interpolate (s, xp, yp);
        }
        else {
          Point point = locate (xp, yp);
          field[j][i] = point.level >= 0 ? val(s) : nodata;
        }
      }
    }
    
    if (pid() == 0) { // master
@if _MPI
  MPI_Reduce (MPI_IN_PLACE, field[0], p.n*p.n, MPI_FLOAT, MPI_MIN, 0,MPI_COMM_WORLD);
@endif
  
      int nv; 
      for (nv = 0; nv < nvarout; nv ++)
        if (strcmp(s.name, nc_varname[nv]) == 0) {
          start[1] += 1;
          break;
        }
      
     if ((nc_err = nc_put_vara_float(ncid, nc_varid[nv], start, count,
                                     &field[0][0])))
       ERR(nc_err);
     
     if (start[1] == p.nl - 1)
       start[1] = -1;
  }
@if _MPI
  else // slave
  MPI_Reduce (field[0], NULL, p.n*p.n, MPI_FLOAT, MPI_MIN, 0,MPI_COMM_WORLD);
@endif
  }
  matrix_free (field);

if (pid() == 0) {
   /* Close the file. */
   if ((nc_err = nc_close(ncid)))
      ERR(nc_err);  
}
   /* printf("*** SUCCESS writing example file %s!\n", p.file_nc); */
}
