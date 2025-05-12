
/* struct OutputMatrix {
  scalar f;
  FILE * fp;
  int n;
  bool linear;
  }; */



/**
 Iterate through dump files and extract the KE, PE, and dissipation rates in air & water
 */

#include "embed.h"

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "navier-stokes/conserving.h"
#include "reduced.h"

#include "tracer.h"
#include "navier-stokes/perfs.h"

#define FLOOR_BOUNDARY_EMBED y+max(0.37109375,1-max(0,0.0693*(x-30)))

#define SOLITON_C 1.2649110640673518
#define X0 15

#define RHOratio 0.000980392156862745
#define MUratio 0.018099999999999998

#define RE 40000

const double h0 = 1;

scalar dye[];
scalar * tracers = {dye};
scalar zeta[];

scalar mask[];

//face vector airforce[];
//int dissipation_rate (vector u, double* rates, scalar drate);

double rho1, rho2;
double mu1, mu2;

double *value;


void output_matrix_ff (scalar F, char *filename, double xmin, double xmax, int nx, double ymin, double ymax, int ny);

int main(int argc, char * argv[])
{

  char  dumpname[1000];  // dumpfile
  char  ctmp[1000];
  char  binaryname[1000];  // name of output binary file
  int dump_num_start, dump_num_end, dump_num_dt;
  double Lx, Ly;
  int NX, NY;
  char scalar_flag;
  double xsoliton;
  double xmin, xmax, ymin, ymax;  
  if (argc < 11){
    printf("Arguments: output_matrix  dump_num_start dump_num_end  dump_num_dt xmin xmax NX   ymin ymax NY  scalar\n");
    return 1;
  }
  dump_num_start = atoi(argv[1]);
  dump_num_end = atoi(argv[2]);
  dump_num_dt = atoi(argv[3]);    
  xmin = atof(argv[4]);
  xmax = atof(argv[5]);  
  NX = atoi(argv[6]);
  ymin = atof(argv[7]);
  ymax = atof(argv[8]);
  NY = atoi(argv[9]);
  scalar_flag = argv[10][0];

  
  strcpy(ctmp,"alldump");
  if (argc==12) {
    strcpy(ctmp,argv[11]);
    }
  value = malloc( sizeof(double)*NY);
  
  switch( scalar_flag) {
  case 'p':
  case 'u':
  case 'v':
  case 'f':
  case 'z':  // vorticity
  case 'a':  // vorticity    
    break;
  default: fprintf(stderr,"** Error: unknown scalar field\n");
    exit(-1);
  }	 
  printf("OUPUT MATRIX: dnum_start=%d  dnum_end=%d, dnum_dt=%d, xmin=%.3lf, xmax=%.3lf, NX=%d, ymin=%.3lf, ymax=%.3lf, NY=%d scalar=%c\n",
	 dump_num_start,dump_num_end,dump_num_dt, xmin,xmax, NX, ymin, ymax, NY, scalar_flag);

  int dumpnum = dump_num_start;
  sprintf(dumpname,"shoal_dump/d%04d",dumpnum);
  
  printf("First file to read: %s\n",dumpname);
  while(restore(file=dumpname))  {
    printf("\r %s: t=%f        \n",dumpname,t);
    //    xsoliton = X0 + (t-1001.0)*SOLITON_C;
    //xmin = xsoliton - 0.5*Lx;
    //xmax = xsoliton + 0.5*Lx;
    
    foreach() {
      if (FLOOR_BOUNDARY_EMBED > 0)
	mask[] = 1;
      else
	mask[] = 0;
    }


    sprintf(binaryname,"shoal_binary/%s_%c_d%04d.bdat",ctmp,scalar_flag,dumpnum);
    switch( scalar_flag) {
    case 'p':
      output_matrix_ff (p, binaryname,  xmin, xmax, NX, ymin, ymax, NY);
      break;
    case 'f':
      output_matrix_ff (f, binaryname,  xmin, xmax, NX, ymin, ymax, NY);
      break;
    case 'z':  // vorticity      
      vorticity(u,zeta);
      output_matrix_ff (zeta, binaryname,  xmin, xmax, NX, ymin, ymax, NY);
      break;
    case 'u':
      output_matrix_ff (u.x, binaryname,  xmin, xmax, NX, ymin, ymax, NY);    
      break;
    case 'v':
      output_matrix_ff (u.y, binaryname,  xmin, xmax, NX, ymin, ymax, NY);
    case 'a':
      sprintf(binaryname,"shoal_binary/alldump_f_d%04d.bdat",dumpnum);
      output_matrix_ff (f, binaryname,  xmin, xmax, NX, ymin, ymax, NY);      
      sprintf(binaryname,"shoal_binary/alldump_p_d%04d.bdat",dumpnum);
      output_matrix_ff (p, binaryname,  xmin, xmax, NX, ymin, ymax, NY);      
      // sprintf(binaryname,"shoal_binary/alldump_u_d%04d.bdat",dumpnum);
      //output_matrix_ff (u.x, binaryname,  xmin, xmax, NX, ymin, ymax, NY);            
      //sprintf(binaryname,"shoal_binary/alldump_v_d%04d.bdat",dumpnum);
      //output_matrix_ff (u.y, binaryname,  xmin, xmax, NX, ymin, ymax, NY);                  
      sprintf(binaryname,"shoal_binary/alldump_z_d%04d.bdat",dumpnum);
      vorticity(u,zeta);
      output_matrix_ff (zeta, binaryname,  xmin, xmax, NX, ymin, ymax, NY);                        
      break;
  default: fprintf(stderr,"** Error: unknown scalar field\n");
    exit(-1);
  }	 


    dumpnum += dump_num_dt;
    if (dumpnum>dump_num_end) {
      break;
    }
    sprintf(dumpname,"shoal_dump/d%04d",dumpnum);

  } 
  printf("Done.                                          \n");

}

void output_matrix_ff (scalar F, char *filename, double xmin, double xmax, int nx, double ymin, double ymax, int ny)
{

  FILE *fp = fopen(filename,"w");
  if (fp == NULL) {
    fprintf(stderr,"*ERROR null file pointer in output_matrix_ff");
    exit(0);
  }
  fprintf(stderr,"OPENED file %s for writing\n",filename);
  
  scalar f = F;
  boundary ({f});  /* what is this? */

  int nnx = nx;
  int nny = ny;  
  
  double dx = (double) (xmax-xmin)/nx;
  double dy = (double) (ymax-ymin)/ny;
  double xp, yp;  // values at the point


  fwrite (&t, sizeof(double),1, fp);   // time 8 bytes 
  fwrite (&nny, sizeof(int), 1, fp);  // Write nny 4 bytes
  fwrite (&nnx, sizeof(int), 1, fp);  // Write nnx 4 bytes

  for (int j = 0; j < ny; j++) {
    yp = (double) (dy*j + ymin + dy/2.0);
    fwrite (&yp, sizeof(double), 1, fp);  // write out y vector
  }
  for (int i = 0; i < nx; i++) {
    xp = (double) (dx*i + xmin + dx/2.0);
    fwrite (&xp, sizeof(double), 1, fp); // write out x vector
  }
  
  for (int i = 0; i < nx; i++) {
    xp = (double) (dx*i + xmin + dx/2.0);
    for (int j = 0; j < ny; j++) {
      yp = (double)(dy*j + ymin + dy/2.0);
      value[j] = interpolate (f, xp, yp);
    }
    fwrite (value, sizeof(double), ny, fp);    
  }
  fflush (fp);
}
