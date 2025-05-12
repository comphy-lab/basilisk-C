/**
# The accuracy of various refinement/prolongation techniques in 2D
Here we test the accuracy of various techniques used for Prolongation and/or refinement. 

## The algorithm
We determine the accuracy of a a few interpolation techniques by applying them to a 2D Gaussian distribution. We therefore use the quadtree toolbox, the *statsf* function and two additional attributes.    
*/

#include "grid/quadtree.h"
#include "utils.h"
#include "additional_attributes.h"

scalar c[],err[];
FILE * fp1;
int cells = 0;
stats s;
int maxlevel = 11;

int main(){
  fp1=fopen("data.dat","w");
  L0=10.;
  X0=Y0=-L0/2;
  /**
  We start with initializing the Gaussian distribution on a fine-resolution grid corresponding to *maxlevel*.
  */
  init_grid(1<<maxlevel);
  foreach()
    c[]=exp(-(sq(x)+sq(y)));
   /**
  In this loop the analysis is done on succesively coarsened grid. 
  */
  
  for (int n=1;n<8;n++){
    cells=0;
    foreach()
      cells++;
    fprintf(ferr,"%d\t %d\n",n,cells);
    fprintf(fp1,"%d\t%g\t",n,L0/(pow(2.,(double)(maxlevel-n+1.))));
    /**
    For each resolution, five different interpolation techniques are tested.
    */
    for (int j=1;j<7;j++){
      if (j==1) c.prolongation=refine_injection;
      if (j==2) c.prolongation=refine_bilinear;
      if (j==3) c.prolongation=refine_linear;
      if (j==4) c.prolongation=refine_quadratic; //Third order?
      if (j==5) c.prolongation=refine_3x3_linear; // different from linear?
      if (j==6) c.prolongation=refine_biquadratic; //Third order?
      /**
      The error is evaluated using the *wavelet* methodology. 
      */
      wavelet(c,err);
      foreach()
	err[]=fabs(err[]);
      s=statsf(err);
      fprintf(fp1,"%g\t",s.sum);
    }
    fprintf(fp1,"\n");
    /**
    After the results are written to a file, the grid is coarsened such that the analysis can be repeated with a different $\Delta$. 
    */
    unrefine(level>=(maxlevel-n)); 
  }
}

/**
##Results 
~~~gnuplot Results from various interpolation/upsample techniques
  set xr [0.004:0.5]
  set yr [0.000001:1]
  set logscale y
  set logscale x
  set xlabel '{/Symbol D}'
  set ylabel 'Sum of Error'
  set key right bottom box 1
  set size square
  plot   (2*x**1) lw 3 lc rgb 'purple' title '{/Symbol \265}{/Symbol D}^{1}' ,\
   (2*x**2) lw 3 lc rgb 'blue' title '{/Symbol \265}{/Symbol D}^{2}' ,\
    (2*x**3) lw 3 lc rgb 'red' title '{/Symbol \265}{/Symbol D}^{3}' ,\
  'data.dat' using 2:3 title 'Injection',\
  'data.dat' using 2:4 title 'Bilinear' ,\
  'data.dat' using 2:5 title 'Linear' ,\
  'data.dat' using 2:6 title 'quadratic' ,\
  'data.dat' using 2:7 title '3x3-Linear' ,\
  'data.dat' using 2:8 title 'biquadratic'
  
~~~

Appearantly the *Injection* method is first-order accurate, *Linear* and *Bilinear* are both second-order accurate and the *3x3-linear* formulation is third-order accurate. The *quadratic* interpolation technique, that was designed to fall in 3-rd order catogory, turns out to be only second-order accurate. So something went wrong there.  

Also note that the *Injection, Linear* and *3x3-linear* techniques are conservative for the first-order moments. This could be favourable, as it may cause the numerical formulation of refinement to enherit some of the properties from the physical system of your interest (e.g. the refinement of a volicity field).  

## Todo

* Find out why the *quadratic* method is not third-order accurate
* Implement a *3x3(x3)-linear* analogue for 3D fields and test this if it is also third-order accurate. (DONE see elswere) 
*/
