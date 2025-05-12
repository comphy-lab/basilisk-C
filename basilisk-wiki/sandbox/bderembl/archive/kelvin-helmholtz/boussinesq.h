
/**
# boussinesq 2d
*/

double Ri = 1.;
double Re = 1.;
double Pr = 1.;
double N2b = 0.; // background stratification
double ar = 1; // aspect ratio

int N0 = 64;
int NMOUT = 5000;
double tend = 1.;
double dtout = 1.;
char dpath[80]  = "./";
double tol_T = 0;
double tol_u = 0;
timer tel;
double tel_m = 1e10;
int npt = 0;


//#include "centered-new.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"

scalar b[];
scalar * tracers = {b};
mgstats mgb;


// for mkdir
#include <sys/stat.h>
#include <sys/types.h>

#define laplacian(b) (b[1] + b[-1] + b[0,1] + b[0,-1] - 4*b[])/(sq(Delta))

event init (t=0) {

#if dimension == 2
  const face vector muc[] = {1/Re,1/(Re*sq(ar))};
#endif
#if dimension == 3
	const face vector muc[] = {1/Re,1/Re,1/Re};
#endif

  mu = muc;
  a = new face vector;

}

event tracer_diffusion (i++) {
#if dimension == 2
	const face vector D[] = {1./(Pr*Re), 1./(Pr*Re*sq(ar))};
#endif
#if dimension == 3
	const face vector D[] = {1./(Pr*Re), 1./(Pr*Re), 1./(Pr*Re)};
#endif
	mgb = diffusion (b, dt, D);
  boundary ({b});
}

event acceleration (i++) {
  face vector av = a;
  foreach_face(y)
    av.y[] += Ri*0.5*(b[] + b[0,-1])/ar;
}

#if TREE
event adapt (i+=2) {
  if (tol_b || tol_u) {

    double uemax = tol_u*normf(u.x).avg;
    int lmax = log2(N0);
    astats s = adapt_wavelet ((scalar *){b,u}, (double []){tol_b,uemax,uemax},
                              maxlevel=lmax, minlevel=6);
    fprintf (ferr, "# refined %d cells, coarsened %d cells, tol_u = %f\n", s.nf, s.nc, uemax);
  }
}
#endif


/**********************************************************************
*                       End of dynamical core                         *
***********************************************************************/

/**
   Read input parameters
 */

void read_params()
{
/**
   Read input parameters
 */

  FILE * fp;
  if ((fp = fopen("params.in", "rt"))) {
    char tempbuff[100];
    char tmps1[80];
    char tmps2[80];
    char tmps3[80];

    while(fgets(tempbuff,100,fp)) {
      sscanf(tempbuff, "%15s = %15s # %15s", tmps1, tmps2, tmps3);
      if      (strcmp(tmps1,"N")    ==0) { N0    = atoi(tmps2); }
      else if (strcmp(tmps1,"L0")   ==0) { L0    = atof(tmps2); }
      else if (strcmp(tmps1,"Re")   ==0) { Re    = atof(tmps2); }
      else if (strcmp(tmps1,"Pr")   ==0) { Pr    = atof(tmps2); }
      else if (strcmp(tmps1,"Ri")   ==0) { Ri    = atof(tmps2); }
      else if (strcmp(tmps1,"N2b")  ==0) { N2b   = atof(tmps2); }
      else if (strcmp(tmps1,"ar")   ==0) { ar    = atof(tmps2); }
      else if (strcmp(tmps1,"npt")  ==0) { npt   = atoi(tmps2); }
      else if (strcmp(tmps1,"DT")   ==0) { DT    = atof(tmps2); }
      else if (strcmp(tmps1,"CFL")  ==0) { CFL   = atof(tmps2); }
      else if (strcmp(tmps1,"tend") ==0) { tend  = atof(tmps2); }
      else if (strcmp(tmps1,"dtout")==0) { dtout = atof(tmps2); }
      else if (strcmp(tmps1,"NMOUT")==0) { NMOUT = atoi(tmps2); }
      else if (strcmp(tmps1,"tel_m")==0) { tel_m = atof(tmps2); }
      else if (strcmp(tmps1,"tol_T")==0) { tol_T = atof(tmps2); }
      else if (strcmp(tmps1,"tol_u")==0) { tol_u = atof(tmps2); }
    }
    fclose(fp);
  } else {
    fprintf(stdout, "file params.in not found\n");
    exit(0);
  }

  X0 = -L0/2;
  Y0 = -L0/2;
  N = N0;

  /**
     Viscosity CFL = 0.5
   */
  if (Re  != 0) DT = 0.5*min(DT,sq(L0/N)*Re/4.);

}

/**
   Create output directory and copy input parameter file for backup
*/
void create_outdir()
{
  if (pid() == 0) {
    for (int i=1; i<10000; i++) {
      sprintf(dpath, "outdir_%04d/", i);
      if (mkdir(dpath, 0777) == 0) {
        fprintf(stdout,"Writing output in %s\n",dpath);
        break;
      }
    }
  }
@if _MPI
  MPI_Bcast(&dpath, 80, MPI_CHAR, 0, MPI_COMM_WORLD);
@endif
}
void backup_config()
{
  char name[100];
  char ch;
  sprintf (name,"%sparams.in", dpath);
  FILE * source = fopen("params.in", "r");
  FILE * target = fopen(name, "w");
  while ((ch = fgetc(source)) != EOF)
    fputc(ch, target);
  fclose(source);
  fclose(target);

}
