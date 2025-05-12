/**


See [large.c]() for a more detailed description. */
#include "embed.h"

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "navier-stokes/conserving.h"
#include "reduced.h"

#include "tracer.h"
#include "navier-stokes/perfs.h"


const double h0 = 1;

#define Y_TOP 10

scalar dye[];
scalar * tracers = {dye};

face vector airforce[];



#define MAX_Y_BINS 1000
#define MAX_COLS 10

void output_along_profile (scalar F, FILE *fp, double time, double xmin, double xmax, int nx, double ymin, double ymax, int ny);

int main(int argc, char * argv[])
{
  //vertical steps

  
  //command simname Y_TOP ybins <cols>
  if (argc != 2){
    printf("Arguments: sim_folder  sim_dir\n");
    return -1;
  }
 

  char simname[1000];

  strcpy(simname, argv[1]);

  char fname[1000];
  sprintf(fname,"%s/shoal_wind_along_profiles.dat",simname);
  FILE * output = fopen(fname,"w");
  fprintf(output,"%%%% t x y ux\n");
  printf("Outputting to: %s\n",fname);

  int dumpfile = 0;
  sprintf(fname,"%s/shoal_dump/d%04d",simname,dumpfile);
  
  printf("First file to read: %s\n",fname);
  while(restore(file=fname)){
    printf("\rd%04d: t=%f        ",dumpfile,t);
    fflush(stdout);

    output_along_profile(u.x,output,t,1,59,59,1,8,8);

    dumpfile++;
    sprintf(fname,"%s/shoal_dump/d%04d",simname,dumpfile);
    //    printf("file = %s\n",fname);
  }
  printf("Done.                                          \n");
  fclose(output);
}


void output_along_profile (scalar F, FILE *fp, double time, double xmin, double xmax, int nx, double ymin, double ymax, int ny)
{


  scalar f = F;
  boundary ({f});  

  int nnx = nx;
  int nny = ny;  
  
  double dx = (double) (xmax-xmin)/(nx-1);
  double dy = (double) (ymax-ymin)/(ny-1);
  double xp, yp;  // values at the point
  double value;

  for (int j = 0; j < ny; j++) {
    yp = (double) (dy*j + ymin );
    for (int i = 0; i < nx; i++) {
       xp = (double) (dx*i + xmin);
       value = interpolate (f, xp, yp);
       fprintf(fp,"%f %f %f %f\n",time,xp, yp,value);
	 }
  }
  fflush (fp);
}

