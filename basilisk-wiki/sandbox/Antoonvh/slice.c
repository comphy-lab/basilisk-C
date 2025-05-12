/**
# Analysis of data slice
This code was used by Van Hooft et al. to obtain the results for the analysis of the adaptation function for for a slice of a 3D turbulent field. It was slightly modified here to read a $256\times256$ data slice instead of the original $1024\times1024$ data sample as presented in the paper.  
*/

#include "grid/quadtree.h"
#include "utils.h"
#define nn (256)
double zeta=0.00;
char name[100],name2[100],name3[100],name4[100];
int maxlevel=8;
scalar c[],err[];
int main(){
  fprintf(ferr,"Matrix initialization\n");
  double ** field = matrix_new(nn,nn,sizeof(double));

  fprintf(ferr,"Open file\n");
  FILE *fp;
  fp = fopen("slice.ux.dat", "r");
  for (int i=0;i<nn;i++){
    for (int j=0;j<nn;j++){
      field[i][j]=0.;
      fscanf(fp,"%lf",&field[i][j]);
    }
  }
  fprintf(ferr,"Done reading data\n");
  init_grid(1<<maxlevel);
  L0=nn;
  X0=Y0=-0.5;
  scalar c[];
  int cells=0;
  fprintf(ferr,"Done doing some basilisk stuff now loading field values\n");
  foreach(){
    c[]=field[(int)x][(int)y];
    cells++;
  }
  fprintf(ferr,"Done reading it into the grid\n");

  fprintf(ferr,"cells:%d\n",cells);

  scalar d[];
  foreach()
    d[]=c[];
  boundary({d});
  adapt_wavelet((scalar *){d},(double[]){9999},maxlevel);
  boundary({d}); //why?
  sprintf(name,"512x512");
  FILE * fp512 = fopen(name,"w");
  double g= pow(2,maxlevel);
  for (double m=0;m<g;m++){
    for (double n=0;n<g;n++){
      fprintf(fp512,"%g\t",interpolate(d,(double)m,(double)n));
    }
    fprintf(fp512,"\n");
    
    }
    fclose(fp512);
  
  for(int j=0;j<21;j++){
    init_grid(nn); 
    foreach()
      c[]=field[(int)x][(int)y];
    boundary({c});
    restriction({c});
    while(adapt_wavelet((scalar *){c},(double[]){zeta},maxlevel).nc){
    boundary({c});
    }
    
    cells=0;
    int cm=0;
    foreach(){
      cells++;
      if (level==maxlevel)
	cm++;
    }
    //fprintf(ferr,"%g %d %d\n",zeta,cells,cm);
    zeta=zeta+0.02;
    wavelet(c,err);
    boundary({c,err});
    sprintf(name2,"errj=%d",j);
    FILE * fpe = fopen(name2,"w");
    sprintf(name,"fieldj=%d",j);
    FILE * fpf = fopen(name,"w");
    double g= pow(2,maxlevel);
    output_field({err},fpe,g,linear=true);
    output_field({c},fpf,g,linear=true);
    fclose(fpf);
    fclose(fpe);
    sprintf(name4,"errcellsj=%d",j);
    FILE * fperr = fopen(name4,"w");
    foreach()
      fprintf(fperr,"%g\t",err[]);
    fclose(fperr);
  }
  matrix_free(field);
}

/**
## Results
let see the non adapted error field: 

 ~~~gnuplot
  set xr [1:256]
 set yr [1:256]
 set zr [-0.1:0.1]
 unset xtics
 unset ytics
 
 set size square
  set xlabel 'x'
  set ylabel 'y'
  set pm3d map
  splot 'errj=0' u 2:1:3
  ~~~
  
And now the adapted error field with $\zeta = 0.1$: 

 ~~~gnuplot
  splot 'errj=5' u 2:1:3
  ~~~
*/