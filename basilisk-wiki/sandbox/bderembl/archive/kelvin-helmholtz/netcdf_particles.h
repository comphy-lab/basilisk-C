/**
   Netcdf interface for basilisk
*/

#include <stdio.h>
#include <string.h>
#include <netcdf.h>
#if _MPI
#include <netcdf_par.h>
#endif

#define FILE_NAME_NC_PART "particles.nc"
#define NDIMS_PART 3
#define TAG_NAME "tag"
#define P_NAME "variables"
#define REC_NAME "time"

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
int nc_varid_part;
int nvar_part = 8;
int rec_varid_part;

int nc_rec_part = -1;

void create_nc_particles(char* file_nc)
{
  /* LOCAL IDs for the netCDF file, dimensions, and variables. */
  int tag_dimid, p_dimid, rec_dimid;
  int tag_varid, p_varid;
  int dimids[NDIMS_PART];
   
  /* Create the file. */
#if _MPI
  if ((nc_err = nc_create_par(file_nc, NC_NETCDF4, MPI_COMM_WORLD, MPI_INFO_NULL, &ncid)))
#else
  if ((nc_err = nc_create(file_nc, NC_CLOBBER, &ncid)))
#endif
    ERR(nc_err);

  /* Define the dimensions. The record dimension is defined to have
   * unlimited length - it can grow as needed. In this example it is
   * the time dimension.*/
  if ((nc_err = nc_def_dim(ncid, TAG_NAME, npt*npt, &tag_dimid)))
    ERR(nc_err);
  if ((nc_err = nc_def_dim(ncid, P_NAME, nvar_part, &p_dimid)))
    ERR(nc_err);
  if ((nc_err = nc_def_dim(ncid, REC_NAME, NC_UNLIMITED, &rec_dimid)))
//  if ((nc_err = nc_def_dim(ncid, REC_NAME, 100, &rec_dimid)))
    ERR(nc_err);

  /* Define the coordinate variables. We will only define coordinate
     variables for lat and lon.  Ordinarily we would need to provide
     an array of dimension IDs for each variable's dimensions, but
     since coordinate variables only have one dimension, we can
     simply provide the address of that dimension ID (&y_dimid) and
     similarly for (&x_dimid). */
  if ((nc_err = nc_def_var(ncid, P_NAME, NC_FLOAT, 1, &p_dimid, 
                           &p_varid)))
    ERR(nc_err);
  if ((nc_err = nc_def_var(ncid, TAG_NAME, NC_FLOAT, 1, &tag_dimid, 
                           &tag_varid)))
    ERR(nc_err);
  if ((nc_err = nc_def_var(ncid, REC_NAME, NC_FLOAT, 1, &rec_dimid, 
                           &rec_varid_part)))
    ERR(nc_err);

  /* The dimids array is used to pass the dimids of the dimensions of
     the netCDF variables. Both of the netCDF variables we are
     creating share the same four dimensions. In C, the
     unlimited dimension must come first on the list of dimids. */
  dimids[0] = rec_dimid;
  dimids[1] = tag_dimid;
  dimids[2] = p_dimid;

  /* Define the netCDF variables */
  if ((nc_err = nc_def_var(ncid, "particles", NC_FLOAT, NDIMS_PART,
                           dimids, &nc_varid_part)))
    ERR(nc_err);
   
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
  float yc[nvar_part], xc[npt*npt];
  for (int i = 0; i < npt*npt; i++){
    xc[i] = i; 
  }
  for (int i = 0; i < nvar_part; i++){
    yc[i] = i;
  }

  if ((nc_err = nc_put_var_float(ncid, p_varid, &yc[0])))
    ERR(nc_err);
  if ((nc_err = nc_put_var_float(ncid, tag_varid, &xc[0])))
    ERR(nc_err);

  /* Close the file. */
  if ((nc_err = nc_close(ncid)))
    ERR(nc_err);
//  printf("*** SUCCESS creating example file %s!\n", file_nc);
}


void write_nc_particles(Particles ptr, char* file_nc) {
  
  /* open file. */
#if _MPI
  if ((nc_err = nc_open_par(file_nc, NC_WRITE, MPI_COMM_WORLD, MPI_INFO_NULL, &ncid)))
#else
  if ((nc_err = nc_open(file_nc, NC_WRITE, &ncid)))
#endif
    ERR(nc_err);
  
  nc_rec_part += 1;

#if _MPI
  if (nc_err = nc_var_par_access(ncid, rec_varid_part, NC_COLLECTIVE))
    ERR(nc_err);
  if (nc_err = nc_var_par_access(ncid, nc_varid_part, NC_COLLECTIVE))
    ERR(nc_err);
#endif

  size_t start[NDIMS_PART], count[NDIMS_PART];
  start[0] = nc_rec_part; //time
  count[0] = 1;
  
  float loctime;
  loctime = 1.*nc_rec_part;
  
  if ((nc_err = nc_put_vara_float(ncid, rec_varid_part, start, count,
                                  &loctime)))
    ERR(nc_err);

    start[0] = nc_rec_part; //time
    start[1] =  pid();
    start[2] = 0;      //vars
      
    count[0] = 1;
    count[1] = 1;
    count[2] = 1;

    if ((nc_err = nc_put_vara_float(ncid, nc_varid_part, start, count,
                                     &loctime)))
       ERR(nc_err);

#if _MPI
  if (nc_err = nc_var_par_access(ncid, rec_varid_part, NC_INDEPENDENT))
    ERR(nc_err);
  if (nc_err = nc_var_par_access(ncid, nc_varid_part, NC_INDEPENDENT))
    ERR(nc_err);
#endif


   particle_boundary(ptr); 

   scalar epsilon[]; 
   scalar ke[]; 
   scalar del2b[]; 
   foreach() { 
     del2b[] = laplacian(b); 
     ke[] = 0.5*(sq(u.x[]) + sq(ar*u.y[])); 
     //      viscous dissipation // FIXME: aspect ratio 
     epsilon[]= (sq(u.x[1] - u.x[-1]) + 
                 sq(u.x[0,1] - u.x[0,-1])/sq(ar) + 
                 sq(u.y[1] - u.y[-1])*sq(ar) + 
                 sq(u.y[0,1] - u.y[0,-1]) 
                 )/sq(2.*Delta); 
   } 
   boundary ({del2b,ke,epsilon}); 


   foreach_particle() { 

     float local_varout[nvar_part]; 
     size_t start[NDIMS_PART], count[NDIMS_PART]; 
      
      
     start[0] = nc_rec_part; //time 
     start[1] =  p().tag; 
     start[2] = 0;      //vars 
      
     count[0] = 1; 
     count[1] = 1; 
     count[2] = nvar_part; 

     local_varout[0] = x; 
     local_varout[1] = y; 
     local_varout[2] = interpolate_5 (b, x, y); 
     local_varout[3] = interpolate_5 (del2b, x, y); 
     local_varout[4] = interpolate_5 (epsilon, x, y); 
     local_varout[5] = interpolate_5 (ke, x, y); 
     local_varout[6] = p().u.x; 
     local_varout[7] = p().u.y; 

     if ((nc_err = nc_put_vara_float(ncid, nc_varid_part, start, count, 
                                     &local_varout[0]))) 
       ERR(nc_err); 

   } 

  /* Close the file. */
  if ((nc_err = nc_close(ncid)))
    ERR(nc_err);
   
//  printf("*** SUCCESS writing example file %s!\n", file_nc);
}
