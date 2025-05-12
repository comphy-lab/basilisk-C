/**
# Netcdf interface for basilisk

These input/output routines are meant to provide a simple way to store and read
netcdf files. It relies on the [netcdf](https://github.com/Unidata/netcdf-c)
library.

This file provides 3 routines:

~~~literatec
create_nc({scalar1, vector1}, "filename.nc");
~~~

which creates the structure of the netcdf file. This should be called only
once for instance in an init event. Then one can use

~~~literatec
write_nc();
~~~

to write a snapshot of the selected variables.

Last, the routine

~~~literatec
read_nc({scalar1, vector1}, "filename.nc");
~~~

reads the selected scalar from file "filename.nc".

Notes: This routine only works for multigrids. It should work for 2d, layered
variables

- TODO: add 3d
- TODO: add units
- TODO: for 2d fields, no need to have a z dimension
*/

#include <stdio.h>
#include <string.h>
#include <netcdf.h>


#define NDIMS 4
#define Y_NAME "y"
#define X_NAME "x"
#define REC_NAME "time"
#define LVL_NAME "level"

#if LAYERS == 0
int nl = 1;
int _layer = 0;
#endif

/* /\* For the units attributes. *\/ */
/* #define UNITS "units" */
/* #define PRES_UNITS "hPa" */
/* #define TEMP_UNITS "celsius" */
/* #define MAX_ATT_LEN 80 */

/**
 Handle errors by printing an error message and exiting with a non-zero status.
*/
int nc_err;
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); return;}

/** 
User global variables: IDs for the netCDF file, dimensions, and variables. 
*/

int ncid;

scalar * nc_scalar_list;
char nc_file[80];
int nc_varid[1000];
int nc_rec = -1;
int rec_varid;

/**
## Netcdf creation

create_nc is used to create the netcdf file
*/



void create_nc(scalar * list_out, char* file_out)
{


  // make it global variable
  sprintf (nc_file,"%s", file_out);
  nc_scalar_list = list_copy(list_out);

    if (pid() == 0) { // master

   /* LOCAL IDs for the netCDF file, dimensions, and variables. */
   int x_dimid, y_dimid, lvl_dimid, rec_dimid;
   int lvl_varid, y_varid, x_varid;
   int dimids[NDIMS];
   
   /* Create the file. */
   if ((nc_err = nc_create(nc_file, NC_CLOBBER, &ncid)))
      ERR(nc_err);

   /* Define the dimensions. The record dimension is defined to have
    * unlimited length - it can grow as needed. In this example it is
    * the time dimension.*/
#if _MPI
   int npx = mpi_dims[0];
   int npy = mpi_dims[1];
#else
   int npx = 1;
   int npy = 1;
#endif
   int Nloc = (1 << depth());
   int Nx = Nloc*npx;
   int Ny = Nloc*npy;

   if ((nc_err = nc_def_dim(ncid, REC_NAME, NC_UNLIMITED, &rec_dimid)))
      ERR(nc_err);
   if ((nc_err = nc_def_dim(ncid, LVL_NAME, nl, &lvl_dimid)))
      ERR(nc_err);
   if ((nc_err = nc_def_dim(ncid, Y_NAME, Ny, &y_dimid)))
      ERR(nc_err);
   if ((nc_err = nc_def_dim(ncid, X_NAME, Nx, &x_dimid)))
      ERR(nc_err);

   /* Define the coordinate variables. We will only define coordinate
      variables for lat and lon.  Ordinarily we would need to provide
      an array of dimension IDs for each variable's dimensions, but
      since coordinate variables only have one dimension, we can
      simply provide the address of that dimension ID (&y_dimid) and
      similarly for (&x_dimid). */
   if ((nc_err = nc_def_var(ncid, REC_NAME, NC_FLOAT, 1, &rec_dimid,
        		    &rec_varid)))
      ERR(nc_err);
   if ((nc_err = nc_def_var(ncid, LVL_NAME, NC_FLOAT, 1, &lvl_dimid,
        		    &lvl_varid)))
      ERR(nc_err);
   if ((nc_err = nc_def_var(ncid, Y_NAME, NC_FLOAT, 1, &y_dimid,
        		    &y_varid)))
      ERR(nc_err);
   if ((nc_err = nc_def_var(ncid, X_NAME, NC_FLOAT, 1, &x_dimid,
        		    &x_varid)))
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
   int nvarout = 0;
   for (scalar s in nc_scalar_list){
       if ((nc_err = nc_def_var(ncid, s.name, NC_FLOAT, NDIMS,
                                dimids, &nc_varid[nvarout])))
         ERR(nc_err);
       nvarout += 1;
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
   float yc[Ny], xc[Nx];
   double Delta = L0*1.0/Nx;
   for (int i = 0; i < Nx; i++)
     xc[i] = X0 + (i + 0.5)*Delta;
   
   for (int i = 0; i < Ny; i++)
      yc[i] = Y0 + (i + 0.5)*Delta;

   float zc[nl];
   for (int i = 0; i < nl; i++)
     zc[i] = i;

   if ((nc_err = nc_put_var_float(ncid, lvl_varid, &zc[0])))
      ERR(nc_err);

   if ((nc_err = nc_put_var_float(ncid, y_varid, &yc[0])))
      ERR(nc_err);
   if ((nc_err = nc_put_var_float(ncid, x_varid, &xc[0])))
      ERR(nc_err);

   /* Close the file. */
   if ((nc_err = nc_close(ncid)))
      ERR(nc_err);
   printf("*** SUCCESS creating example file %s!\n", nc_file);
}

}


/**
## Write in netcdf 

We use write_nc to dump a snapshot in the netcdf file at time t.
*/

void write_nc() {

#if _MPI
   int npx = mpi_dims[0];
   int npy = mpi_dims[1];
#else
   int npx = 1;
   int npy = 1;
#endif
   int Nloc = (1 << depth());
   int Nx = Nloc*npx;
   int Ny = Nloc*npy;


  bool nc_linear_interp = false;

  if (pid() == 0) { // master
  /* open file. */
  if ((nc_err = nc_open(nc_file, NC_WRITE, &ncid)))
    ERR(nc_err);
  }

  // write time
  nc_rec += 1;
  float loctime = t;

  size_t startt[1], countt[1];
  startt[0] = nc_rec; //time
  countt[0] = 1;
  if (pid() == 0) { // master
  if ((nc_err = nc_put_vara_float(ncid, rec_varid, startt, countt,
                                  &loctime)))
    ERR(nc_err);
  }


  float Delta = L0/Nx;
  float * field = (float *)malloc(Ny*Nx*nl*sizeof(float));

  
  /* The start and count arrays will tell the netCDF library where to
     write our data. */
  size_t start[NDIMS], count[NDIMS];
  
  
  /* These settings tell netcdf to write one timestep of data. (The
     setting of start[0] inside the loop below tells netCDF which
     timestep to write.) */
  start[0] = nc_rec; //time
  start[1] = 0;     //level
  start[2] = 0;      //y
  start[3] = 0;      //x
  
  count[0] = 1;
  count[1] = nl;
  count[2] = Ny;
  count[3] = Nx;

  
  int nv = -1;
  /* char * str1; */
  for (scalar s in nc_scalar_list){
    nv += 1;

    for (_layer = 0; _layer < nl; _layer++){
      for (int j = 0; j < Ny; j++) {
        float yp = Delta*j + Y0 + Delta/2.;
        for (int i = 0; i < Nx; i++) {
          float xp = Delta*i + X0 + Delta/2.;
          if (nc_linear_interp) {
            field[Ny*Nx*_layer + Nx*j + i] = interpolate (s, xp, yp);
          }
          else {
            Point point = locate (xp, yp);
            field[Ny*Nx*_layer + Nx*j + i] = point.level >= 0 ? val(s) : nodata;
          }
        }
      }
    }
    _layer = 0;

    
    if (pid() == 0) { // master
#if _MPI
        MPI_Reduce (MPI_IN_PLACE, &field[0], Ny*Nx*nl, MPI_FLOAT, MPI_MIN, 0,MPI_COMM_WORLD);
#endif
  

     if ((nc_err = nc_put_vara_float(ncid, nc_varid[nv], start, count,
        			      &field[0])))
         ERR(nc_err);

  }
#if _MPI
  else // slave
  MPI_Reduce (&field[0], NULL, Ny*Nx*nl, MPI_FLOAT, MPI_MIN, 0,MPI_COMM_WORLD);
#endif
  }
//  matrix_free (field);
  free(field);

   /* Close the file. */
  if (pid() == 0) { // master
   if ((nc_err = nc_close(ncid)))
      ERR(nc_err);
   }
//   printf("*** SUCCESS writing example file %s -- %d!\n", nc_file, nc_rec);
}

/**
## Read Netcdf file */

void read_nc(scalar * list_in, char* file_in, bool read_time = false)
{
  int i;
  int ncfile, ndims, nvars, ngatts, unlimited;
  int var_ndims, var_natts;
  int t_id;
  nc_type type;
  char varname[NC_MAX_NAME+1];
  int *dimids = NULL;

#if _MPI
  int npx = mpi_dims[0];
  int npy = mpi_dims[1];
#else
  int npx = 1;
  int npy = 1;
#endif
  int Nloc = (1 << depth());
  int Nx = Nloc*npx;
  int Ny = Nloc*npy;

  float * field = (float *)malloc(Ny*Nx*nl*sizeof(float));

  if ((nc_err = nc_open(file_in, NC_NOWRITE, &ncfile)))
    ERR(nc_err);

  if ((nc_err = nc_inq(ncfile, &ndims, &nvars, &ngatts, &unlimited)))
    ERR(nc_err);

  if (read_time) {
    if ((nc_err = nc_inq_varid(ncfile, "time", &t_id)))
      ERR(nc_err);
    size_t startt[1], countt[1];
    startt[0] = 0; //time
    countt[0] = 1;
    float loctime;
    
    if ((nc_err = nc_get_vara_float(ncfile, t_id, startt, countt,
                                    &loctime)))
      ERR(nc_err);
    t = loctime;
  }

  for (scalar s in list_in) {
    for (i = 0; i < nvars; i++) {

      if ((nc_err = nc_inq_var(ncfile, i, varname, &type, &var_ndims, dimids,
                               &var_natts)))
        ERR(nc_err);

      if (strcmp(varname,s.name) == 0) {
        fprintf(stdout,"Reading variable  %s!\n", s.name);

        // There is an issue with this strategy, because, there may be variables
        // with one layer... TO BE FIXED

        /* int nl_loc = _attribute[s.i].block; */

        /* if (nl_loc == 1){ // no layer */

        /*   size_t start[3], count[3]; */
        /*   start[0] = 0; //time */
        /*   start[1] = 0; */
        /*   start[2] = 0; */

        /*   count[0] = 1; */
        /*   count[1] = Ny; */
        /*   count[2] = Nx; */
        /*   if ((nc_err = nc_get_vara_float(ncfile, i, start, count, */
        /*                                   &field[0]))) */
        /*     ERR(nc_err); */

        /*   foreach(noauto) */
        /*     s[] = field[Nx*_J + _I]; */

        /* } else { // layer */

        size_t start[4], count[4];
        start[0] = 0; //time
        start[1] = 0;
        start[2] = 0;
        start[3] = 0;
        
        count[0] = 1;
        count[1] = nl;
        count[2] = Ny;
        count[3] = Nx;
        if ((nc_err = nc_get_vara_float(ncfile, i, start, count,
                                        &field[0])))
          ERR(nc_err);
        
        foreach_layer()
          foreach (noauto) // why noauto???
            s[] = field[Ny*Nx*_layer + Nx*_J + _I]; 
        /* } // end if layer */

      } // end if strcmp
    } // en loop nvars
  }// end scalar loop

  free (field);

  if ((nc_err = nc_close(ncfile)))
    ERR(nc_err);
}

event cleanup (t = end)
{
  free (nc_scalar_list);
}
