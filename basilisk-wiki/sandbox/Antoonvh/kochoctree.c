/**
# The 1.26D Koch snowflake in 3D space. 

Very similar to what we did using a 2D quadtree [here](koch.c), we can try to do the same analysis employing a 3D -so-called- octree. 
*/

#include "grid/octree.h"
#include "navier-stokes/centered.h"

int main(){
  int jmax=8;
  int n=4;
  int nm=((n-1)*pow(4,jmax-1))+1;
  double xn[nm],yn[nm];
  double xk[nm];
  double yk[nm];
  xk[0]=0;
  xk[1]=1;
  xk[2]=0.5;
  xk[3]=0;
  yk[0]=0;
  yk[1]=0;
  yk[2]=-sqrt(3.)/2;
  yk[3]=0;
  fprintf(ferr,"\nyk[n-1]=%g\n\n",yk[n-1]);
  for (int j=1;j<jmax;j++){
    for (int k=0;k<(n-1);k++){
      xn[(k*4)]=xk[(k)];
      xn[(k*4)+1]=(2*xk[k]+xk[k+1])/3;
      xn[(k*4)+2]=((xk[k]+xk[k+1])/2)+((yk[k]-yk[k+1])*sqrt(3.)/6);      
      xn[(k*4)+3]=(xk[k]+2*xk[k+1])/3;
      yn[(k*4)]=yk[k];
      yn[(k*4)+1]=(2*yk[k]+yk[k+1])/3;
      yn[(k*4)+2]=((yk[k]+yk[k+1])/2)-((xk[k]-xk[k+1])*sqrt(3.)/6);
      yn[(k*4)+3]=(yk[k]+2*yk[k+1])/3;
    }
    yn[4*(n-1)]=yk[n-1];
    xn[4*(n-1)]=xk[n-1];
    FILE * fpkoch = fopen("koch.dat","w");
    for (int l=0;l<(4*(n-1)+1);l++){
      yk[l]=yn[l];
      xk[l]=xn[l];
    }
    if (j==jmax-1){
      for (int l=0;l<(4*(n-1)+1);l++){
	fprintf(fpkoch,"%g %g\n",xk[l],yk[l]);
      }
    }
    n=(n-1)*4+1;
  }
  
  FILE * fpit = fopen("cells.dat","w");
  X0=-0.5;
  Y0=-1.2;
  Z0=-1.0432; 
  L0=2.0;
  init_grid(32);
  
  scalar c[];
  for (int m=0;m<6;m++){
    foreach()
      c[]=0;
    for(int i=0;i<n;i+=64/pow(2,m)){
      /**
      We place the koch snowflake on a x-y plane at $z = 0.034124$
      */
      foreach_point (xn[i],yn[i],0.0343124)
        c[]+=1.0;
    }
    int gcc=0;
    int gcn=0;
    foreach(reduction(+:gcc) reduction(+:gcn)){
      gcn++;
      if(c[]>0.5)
	gcc++;
    }
    boundary({c});
    while(adapt_wavelet({c},(double[]){0.1},(7+m),4,{c}).nf);
    fprintf(fpit,"%d %d %d\n",m,gcc,gcn);
    fflush(fpit);
    scalar lev[];
    foreach()
      lev[]=level;
    
    static FILE * fp3 = NULL;
    if (!fp3) fp3 = popen ("gfsview-batch3D kochoctree.hop.gfv | ppm2gif > frac3d.gif --delay 100","w");
    output_gfs(fp3);
    fprintf (fp3, "Save stdout { format = PPM width = 600 height = 600}\n");
  }
}

/**
## Results
Below we see a visualization of the iterative refinement on the fractal curve in 3D space

![Visualization of the refinement iterations with the c=1 isosorfuce, colored with the level of refinement. The results form the last two iterations are not displayed as they require a sub-pixel resolution](kochoctree/frac3d.gif)

Now let us check if we can find the fractal dimension of this curve by plotting the number of grid cells against the refinement iteration.

~~~gnuplot
set xr [ -0.5:5.5]
set yr [100:100000]
set logscale y
set xlabel 'refinement iteration' 
set ylabel 'Grid cells on the curve'
plot "cells.dat" title "Number of cells on the curve" ,\
      (2**(1.26*x))*200 title "2^{1.26x}"
~~~

Apparently the fractal dimension is still $1.26$.  

Finally, we check if the total employed number of gridcells scales with the fractal dimension of this "$\text{problem}$".  

~~~gnuplot
set xr [ -0.5:5.5]
set yr [10000:3000000]
set logscale y
set xlabel 'refinement iteration'
set ylabel 'Total number of grid cells'
plot "cells.dat" using 1:3 title "Total number of cells" ,\
      (2**(1.26*x))*5000 title "2^{1.26x}"
~~~

After some initial hick-up, it does appear to be the case. Very beneficial compared to 3D scaling. Again, Well done adapt_wavelet() function! 

[But can we do it without the adapt_wavelet function as well?](kochrefine.c)
*/
