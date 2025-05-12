/**
   Netcdf interface for basilisk
*/

#include <stdio.h>
#include <string.h>
#include <netcdf.h>
#pragma autolink -lnetcdf

#define NDIMS 3
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

// User global variables
scalar * scalar_list_nc;
char file_nc[80];

/* IDs for the netCDF file, dimensions, and variables. */
int ncid;
int t_varid;

// temporary
int nc_varid[1000];
char * nc_varname[1000];
int nvarout = 0;
int nc_rec = -1;

// tempo
int nl_tmp = 1;

void create_nc()
{

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
   /* if ((nc_err = nc_def_dim(ncid, LVL_NAME, nl, &lvl_dimid))) */
   /*    ERR(nc_err); */
   if ((nc_err = nc_def_dim(ncid, Y_NAME, N+1, &y_dimid)))
      ERR(nc_err);
   if ((nc_err = nc_def_dim(ncid, X_NAME, N+1, &x_dimid)))
      ERR(nc_err);
   if ((nc_err = nc_def_dim(ncid, REC_NAME, NC_UNLIMITED, &rec_dimid)))
      ERR(nc_err);

   /* Define the coordinate variables. We will only define coordinate
      variables for lat and lon.  Ordinarily we would need to provide
      an array of dimension IDs for each variable's dimensions, but
      since coordinate variables only have one dimension, we can
      simply provide the address of that dimension ID (&y_dimid) and
      similarly for (&x_dimid). */
   if ((nc_err = nc_def_var(ncid, REC_NAME, NC_FLOAT, 1, &rec_dimid,
        		    &t_varid)))
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
   /* dimids[1] = lvl_dimid; */
   /* dimids[2] = y_dimid; */
   /* dimids[3] = x_dimid; */
   dimids[1] = y_dimid;
   dimids[2] = x_dimid;

   /* Define the netCDF variables */
//   char * str1;
   for (scalar s in scalar_list_nc){
     /* if (strcmp(str1,s.name) != 0) { */
       if ((nc_err = nc_def_var(ncid, s.name, NC_FLOAT, NDIMS,
                                dimids, &nc_varid[nvarout])))
         ERR(nc_err);
     /*   nc_varname[nvarout] = strdup(s.name); */
       nvarout += 1;
     /*   str1 = strdup(s.name); */
     /* } */
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
   float yc[N+1], xc[N+1];
   double Delta = L0*1.0/N;
   for (int i = 0; i < N+1; i++){
      yc[i] = Y0 + i*Delta;
      xc[i] = X0 + i*Delta;
   }

   if ((nc_err = nc_put_var_float(ncid, y_varid, &yc[0])))
      ERR(nc_err);
   if ((nc_err = nc_put_var_float(ncid, x_varid, &xc[0])))
      ERR(nc_err);

   /* Close the file. */
   if ((nc_err = nc_close(ncid)))
      ERR(nc_err);
   printf("*** SUCCESS creating example file %s!\n", file_nc);

}


struct OutputNetcdf {
  int it;
  int n; 
  bool linear;
} OutputNetcdf;

void write_nc(struct OutputNetcdf p) {
  if (p.n == 0) p.n = N + 1;

  if (pid() == 0) { // master
    /* open file. */
    if ((nc_err = nc_open(file_nc, NC_WRITE, &ncid)))
      ERR(nc_err);
  }

  // write time
  nc_rec += 1;
  float loctime = t;

  size_t startt[1], countt[1];
  startt[0] = nc_rec; //time
  countt[0] = 1;
  if (pid() == 0) { // master
    if ((nc_err = nc_put_vara_float(ncid, t_varid, startt, countt,
                                    &loctime)))
      ERR(nc_err);
  }



  float fn = p.n, Delta = L0/fn;
  float ** field = matrix_new (p.n, p.n, sizeof(float));
  
  /* The start and count arrays will tell the netCDF library where to
     write our data. */
  size_t start[NDIMS], count[NDIMS];
  
  
  /* These settings tell netcdf to write one timestep of data. (The
     setting of start[0] inside the loop below tells netCDF which
     timestep to write.) */
  /* start[0] = nc_rec; //time */
  /* start[1] = 0;     //level */
  /* start[2] = 0;      //y */
  /* start[3] = 0;      //x */
  
  /* count[0] = 1; */
  /* count[1] = 1; */
  /* count[2] = p.n; */
  /* count[3] = p.n; */
  start[0] = nc_rec; //time
  start[1] = 0;      //y
  start[2] = 0;      //x
  
  count[0] = 1;
  count[1] = p.n;
  count[2] = p.n;

  
  int nv = -1;
  /* char * str1; */
//  foreach_layer() {
  for (scalar s in scalar_list_nc){
    nv += 1;

    for (int j = 0; j < p.n; j++) {
      for (int i = 0; i < p.n; i++) {
        field[j][i] = nodata;
      }
    }

    foreach_vertex(noauto){
//      printf ("%d\t%d\t %g\n", point.i-GHOSTS, point.j-GHOSTS, s[]);
      field[_J][_I] = s[];
    }



    /* for (int j = 0; j < p.n; j++) { */
    /*   float yp = Delta*j + Y0 + Delta/2.; */
    /*   for (int i = 0; i < p.n; i++) { */
    /*     float xp = Delta*i + X0 + Delta/2.; */
    /*     if (p.linear) { */
    /*       field[j][i] = interpolate (s, xp, yp); */
    /*     } */
    /*     else { */
    /*       Point point = locate (xp, yp); */
    /*       field[j][i] = point.level >= 0 ? val(s) : nodata; */
    /*       printf ("%d\t%d\t %g\n", i, j, field[j][i]); */

    /*     } */
    /*   } */
    /* } */
    
    if (pid() == 0) { // master
@if _MPI
        MPI_Reduce (MPI_IN_PLACE, field[0], p.n*p.n, MPI_FLOAT, MPI_MIN, 0,MPI_COMM_WORLD);
@endif
  
 /*       int nv; */
 /* for (nv = 0; nv < nvarout; nv ++) */
 /*   if (strcmp(s.name, nc_varname[nv]) == 0) { */
 /*       start[1] += 1; */
 /*       break; */
 /*     } */

     if ((nc_err = nc_put_vara_float(ncid, nc_varid[nv], start, count,
        			      &field[0][0])))
         ERR(nc_err);

     /* if (start[1] == p.nl - 1) */
     /*   start[1] = -1; */
  }
@if _MPI
  else // slave
  MPI_Reduce (field[0], NULL, p.n*p.n, MPI_FLOAT, MPI_MIN, 0,MPI_COMM_WORLD);
@endif
//  }
  }
  matrix_free (field);


   /* Close the file. */
  if (pid() == 0) { // master
    if ((nc_err = nc_close(ncid)))
      ERR(nc_err);
  }
//   printf("*** SUCCESS writing example file %s -- %d!\n", file_nc, nc_rec);
}


/**
   Read NC
 */

void read_nc(scalar * list_in, char* file_in){

  int i, ret;
  int ncfile, ndims, nvars, ngatts, unlimited;
  int var_ndims, var_natts;
  nc_type type;
  char varname[NC_MAX_NAME+1];
  int *dimids=NULL;

  int Nloc = N+1;
  float ** field = matrix_new (Nloc, Nloc, sizeof(float));

  if ((nc_err = nc_open(file_in, NC_NOWRITE, &ncfile)))
    ERR(nc_err);

  if ((nc_err = nc_inq(ncfile, &ndims, &nvars, &ngatts, &unlimited)))
    ERR(nc_err);

  size_t start[NDIMS], count[NDIMS];
  start[0] = 0; //time
  start[1] = 0;
  start[2] = 0;

  count[0] = 1;
  count[1] = Nloc;
  count[2] = Nloc;


  for (scalar s in list_in){
    for(i=0; i<nvars; i++) {


  if ((nc_err = nc_inq_var(ncfile, i, varname, &type, &var_ndims, dimids,
                          &var_natts)))
    ERR(nc_err);


      if (strcmp(varname,s.name) == 0) {
        fprintf(stdout,"Reading variable  %s!\n", s.name);

          if ((nc_err = nc_get_vara_float(ncfile, i, start, count,
                                                 &field[0][0])))
            ERR(nc_err);

          foreach_vertex(noauto){
            s[] = field[_J][_I];
          }
          
        }


    }
  }

  matrix_free (field);

  if ((nc_err = nc_close(ncfile)))
    ERR(nc_err);

  boundary(list_in);

}

event cleanup (t = end)
{
  free(scalar_list_nc);
}
