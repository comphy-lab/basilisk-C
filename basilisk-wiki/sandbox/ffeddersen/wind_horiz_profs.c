/**
# Large-amplitude standing wave

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


scalar dye[];
scalar * tracers = {dye};

face vector airforce[];



#define MAX_Y_BINS 1000
#define MAX_COLS 10

int main(int argc, char * argv[])
{
  //vertical steps
  int y_bins;
  double col_xmin[MAX_COLS];
  double col_xmax[MAX_COLS];
  double col_xmid[MAX_COLS];
  int num_cols = 1;
  col_xmin[0] = 0;
  col_xmax[0] = 2;
  col_xmid[0] = 1;
  double Y_TOP;
  
  //command simname Y_TOP ybins <cols>
  if (argc < 5){
    printf("Arguments: sim_folder dump_group Y_TOP y_bins <cols...>\n");
    return 1;
  }
  Y_TOP = atof(argv[3]);
  y_bins = atof(argv[4]);

  
  if (argc > 4){
    num_cols = 0;
    for(int i = 5; i+1 < argc;){
            col_xmin[num_cols] = atof(argv[i++]);
            col_xmax[num_cols] = atof(argv[i++]);
	    col_xmid[num_cols] = (col_xmin[num_cols] + col_xmax[num_cols])/2.0;
	    num_cols++;
	   
    }
  }
  printf("%d columns to read:\n",num_cols);
  for(int i = 0; i < num_cols; i++){
    printf("  Column %d: x in [%.4f,%.4f]\n",i,col_xmin[i],col_xmax[i]);
  }

  double ux_bins[MAX_Y_BINS*MAX_COLS];
  double y_step = Y_TOP / ((double) y_bins);


  char simname[1000];
  char * root_dir = "/home/kghanson/basilisk/shoaling/runs/test_2022_11_18";
  if(argc > 1){
    strcpy(simname, argv[1]);
  } 
  char fname[1000];
  sprintf(fname,"%s/%s_horiz_profs.dat",simname,argv[2]);
  FILE * output1 = fopen(fname,"w");
  fprintf(output1,"%%%% t x y ux\n");
  printf("Outputting to: %s\n",fname);

  int dumpfile = 0;
  sprintf(fname,"%s/%s_dump/d%04d",simname,argv[2],dumpfile);
  
  printf("First file to read: %s\n",fname);
  while(restore(file=fname)){
    printf("\rd%04d: t=%f        ",dumpfile,t);
    fflush(stdout);
    //initialize ux accumulator to zero:
    for(int bin = 0; bin < y_bins*num_cols; bin++)
      ux_bins[bin] = 0;

    foreach(reduction(+:ux_bins)){
//    foreach(){
      //accumulate if we are within a column
      double hw = Delta * 0.5; //half width

      for(int col = 0; col < num_cols; col++){
        double x_L = col_xmin[col], x_R = col_xmax[col];
        if((y + y_step * 0.5) >= 0 && (x - hw) < x_R && (x + hw) > x_L){
	
          //amount inside column on horizontal length
          double horiz_len = min(x + hw, x_R) - max(x - hw, x_L);
          int bin = max((int) floor( (y - hw)/y_step ), 0);
          //amount inside bin on vertical length
	  double vert_len = min(y + hw, (bin + 1)*y_step) - max(y - hw, bin*y_step);
	  do{
            // integral ( u.x[] * indicator(bin) * dv() )
            ux_bins[col*y_bins + bin] += u.x[] * (horiz_len * vert_len);
	    bin++;
	    vert_len = min(y + hw, (bin + 1)*y_step) - max(y - hw, bin*y_step);
	  }while(vert_len > 0 && bin < y_bins);

        }
      }
    }
    //output x-averaged bins
    for(int col = 0; col < num_cols; col++){
      for(int bin = 0; bin < y_bins; bin++){
        fprintf(output1,"%f %f %f %f\n",t,col_xmid[col], (bin + 0.5)*y_step,
	    ux_bins[col*y_bins + bin] / ((col_xmax[col]-col_xmin[col])*y_step));
      }
    }


    dumpfile++;
    sprintf(fname,"%s/%s_dump/d%04d",simname,argv[2],dumpfile);
    //    printf("file = %s\n",fname);
  }
  printf("Done.                                          \n");
  fclose(output1);
}
