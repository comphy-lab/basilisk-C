/**
# The Fractal Dimension of the Koch Snowflake
This page presents an investigation on the Koch snowflake. You may want to  [watch this youtube video](https://www.youtube.com/watch?v=gB9n2gHsHN4) for a detailed, 20 minute introduction to the topic of fractal curves. Since we have the tools, let us use an adaptive quad-tree grid and the ppm2gif converter. 
*/
#include "grid/quadtree.h"
#include "utils.h"

int main(){

  /**
  ##The Koch Curve
  First, we need to obtain the snowflake curve. We determine it upto 8 fractal levels of refinement, resulting in $3 \cdot 4^{8} \approx 200.000$ points on the curve. This determines the maximum level of refinement that we can use in our grid based analysis.  
  */
  int jmax=9;
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
  for (int j=1;j<jmax;j++){//For every fractal level
    for (int k=0;k<(n-1);k++){//For every existing point add three new points. 
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
    for (int l=0;l<(4*(n-1)+1);l++){// Copy yn and xn into yk and xk. 
      yk[l]=yn[l];
      xk[l]=xn[l];
    }
    if (j==jmax-1){//Output the points of the curve after the last interation. 
      for (int l=0;l<(4*(n-1)+1);l++){
        fprintf(fpkoch,"%g %g\n",xk[l],yk[l]);
      }
    }
    n=(n-1)*4+1;
  }
  
  /** 
  Let us check-out the obtained snowflake. 
  
  ~~~gnuplot
  set xr [-0.2:1.2]
  set yr [-1.0:0.4]
  set xlabel 'x'
  set ylabel 'y'
  set size square
  plot 'koch.dat' title 'Koch snowflake' with lines
  ~~~
  
  Looks O.K. Note that the curve is defined at a much finer resolution than is displayed in the plot. 
  
  ## The Algorithm 
  To find the fractal dimension of this curve we iteratively refine the grid and log how many cells are located on the curve. 
  */
  FILE * fpit = fopen("cells.dat","w");
  X0=-0.5;
  Y0=-1.2;
  L0=2.0;
  init_grid(32);
  
  scalar c[];
  for (int m=0;m<7;m++){
    foreach()
      c[]=0;
    for(int i=0;i<n;i+=64/pow(2,m))
      foreach_point(xn[i], yn[i])
        c[]+=1.0;
    int gcc=0;
    int gcn=0;
    foreach(reduction(+:gcc) reduction(+:gcn)){
      gcn++;
      if(c[]>0.5)
	gcc++;
    }
    /**
    For the refinment we (once again) have faith in the wavelet based adaption algorithm to refine (and coarsen) the grid. 
    */
    boundary({c});
    while(adapt_wavelet({c},(double[]){0.1},(7+m),4,{c}).nf);
    
    fprintf(fpit,"%d %d %d\n",m,gcc,gcn);
    fflush(fpit);
    scalar lev[];
    foreach(){
      lev[]=level;
    }
    output_ppm (lev, file = "f.mp4", opt = "-r 1", 
                n = 512, min =4, max=15);
  }
}

/**
## Results
The adaptive grid that is used to refine the cells in the neighborhood of the curve is visualized below. More red colors represent a higher level of refinement. 

![The various stages of refinement](koch/f.mp4)

Again the displayed resolution ($512 \times 512$) is much less than the final resulution of our analysis (that corresponds to $4096 \times 4096$). 

Now we check if we can find the fractal dimension of this curve by plotting the number of grid cells against the refinement iteration.

~~~gnuplot
set xr [ -0.5:6.5]
set yr [100:100000]
set logscale y
set xlabel 'refinement iteration'
set ylabel 'Grid cells on the curve'
set key bottom right
plot "cells.dat" pt 3 title "Number of cells on the curve"  ,\
      (2**(1.26*x))*200 lw 2 title "2^{1.26x} "
~~~

Apparently the fractal dimension is 1.26 and not 1. The latter value would be typical for a 'regular' curve. As was explained in the aforementioned youtube video and also noted in [another source](https://en.wikipedia.org/wiki/Koch_snowflake), the fractal dimension can be expressed analytically as: $\frac[ln(4)][ln(3)] \approx 1.26$. That number is rather close to the one found with the present method.  

Finally, we check if the total employed number of gridcells scales with the fractal dimension of this "$\text{problem}$".  

~~~gnuplot
set xr [ -0.5:6.5]
set yr [500:1000000]
set logscale y
set xlabel 'refinement iteration'
set ylabel 'Total number of grid cells'
plot "cells.dat" using 1:3 pt 4 title "Total number of cells" ,\
      (2**(1.26*x))*1000 lw 3 title "2^{1.26x}"
~~~

It does appear to be the case. Well done adapt_wavelet() function! 

## Follow-up
There are two other, related cases to adress the following questions:

* [Can we do the same analysis without relying on the adapt_wavelet() function?](kochrefine.c)      
* [Does the algorith still work properly if we do our analysis in tree dimensional space?](kochoctree.c)
                                                                                                                                                    
*/
