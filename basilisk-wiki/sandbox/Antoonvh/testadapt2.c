/**
# Test adapt2.h
We test the *adapt2* function. This function can employ a different maximum level of refinment for each field. This can be interesting for case where multiple physical processes are present in a single domain. 

## Overview
After the siminal work of Stephane Popinet, introducing the *adapt_wavelet* function, some Basilisk users wanted a more flexible function. Most Notably, Ceasar Pairietti developed a variant of the *adapt_wavelet* function where the maximum allowed resolution could vary within the domain, [Look here](http://basilisk.fr/sandbox/pairetti/bag_mode/adapt_wavelet_limited.h). A similar approach is taken here, except that the maximum resolution is still global but now dependents on the field the algorithm is basing its adaptation upon. Allowing to distinguish between different physical processes and  to resolve them with different maximum resolutions. The algorithm should employ the 'highest resolution asked for'.      

## Test case
As a test, we initialize two droplets that will fall due to a body force. We assume (, without providing any motivation,) that the interfacial waves and interface advection requires a higher resolution than the wake behind the droplets. 
*/
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "adapt2.h" //<- Click here to see the code. 

int maxlevel = 11;
int gridm[12];
face vector grav[];
int m;
FILE * fp;

int main(){
  init_grid(256);
  L0=256;
  X0=Y0=-L0/2;
  rho1=1.;
  rho2=0.01;
  mu1=0.1;
  mu2=0.001;
  f.sigma = 0.05;
  a=grav;
  for (m=0;m<4;m++)
    run(); 
}

event init (t=0){
  char name[100];
  sprintf(name,"grid%d.dat",m);
  fp = fopen(name,"w");
  scalar f1[],f2[];
  refine(level<maxlevel && sq(x)+sq(y)<36.&&sq(x)+sq(y)>16.);
  refine(level<maxlevel && sq(x-20)+sq(y-20)<16&&sq(x-20)+sq(y-20)>4);
  fraction (f1,25.-sq(x)-sq(y));
  fraction (f2,9.-sq(x-20.)-sq(y-20.));
  foreach(){
    f[]=f1[]+f2[];
  }
  boundary(all);
}

event acceleration(i++){
  foreach_face(y)
    grav.y[]-=1.;
}
/**
## Methods
We use four different runs,

1. The default *adapt_wavelet* function
2. The new function *adapt_wavelet2* using the same maximum level for all fields. 
3. The new function *adapt_wavelet2* that uses a coarser maximum level of refinement for the velocity components. 
4. Same as 3, but now the list of fields (and corresponding lists of refinement criteria and maximum levels) are in a different order. 
*/
event adapt(i++){
  if (m==0)
    adapt_wavelet((scalar *){f,u.x,u.y},(double []){0.001,0.01,0.01},maxlevel,3);
  if (m==1)
    adapt_wavelet2((scalar *){f,u.x,u.y},(double []){0.001,0.01,0.01},(int []){maxlevel,maxlevel,maxlevel},3);
  if (m==2)
    adapt_wavelet2((scalar *){f,u.x,u.y},(double []){0.001,0.01,0.01},(int []){maxlevel,maxlevel-2,maxlevel-2},3);
  if (m==3)
    adapt_wavelet2((scalar *){u.x,u.y,f},(double []){0.01,0.01,0.001},(int []){maxlevel-2,maxlevel-2,maxlevel},3);
}

event gridmonitor(t+=0.1;t<=5){
  for (int gg=0;gg<=maxlevel;gg++)
    gridm[gg]=0;
  foreach(){
    gridm[level]++;
  }
  fprintf(fp,"%g\t",t);
  for (int gg=0;gg<=maxlevel;gg++)
    fprintf(fp,"%d\t",gridm[gg]);
  fprintf(fp,"\n");
  
  static FILE * fpmp4 =
    popen ("gfsview-batch2D testadapt2.gridbg.gfv | "
           "ppm2mp4 testing.mp4", "w");
  output_gfs (fpmp4);
  fprintf (fpmp4, "Save stdout { format = PPM width = 512 height = 512}\n");
}
/**
## Results
We check the corresponding grids using the default and the *adapt2* concept. 

<video width="512" height="512" controls>
<source src="testadapt2/testing.mp4" type="video/mp4">
</video> 
<br>
<video width="512" height="512" controls>
<source src="default.mp4" type="video/mp4">
</video> 
<br>
<video width="512" height="512" controls>
<source src="adapt2.mp4" type="video/mp4">
</video> 

As a first check, we look at the number of grid cells for the default and proposed approach, where the new approach is set to mimic the default one:

~~~gnuplot Results are close but not exactly the same
set yr [100:50000]
set logscale y
set xlabel 'time'
set ylabel 'Grid cells' 
set size square
plot 'grid0.dat' using 1:13 with line lw 4 title "Level = 11 Default" ,\
     'grid0.dat' using 1:11 with line lw 4 title "Level = 9 Default" ,\
     'grid1.dat' using 1:13 with line lw 2 title "Level = 11 Adapt 2 mimicking the Default" ,\
     'grid1.dat' using 1:11 with line lw 2 title "Level = 9 Adapt 2 micmicing the Default" 
~~~

To test an other aspect of the algorithm, we check the grid dependence on the order of scalar-field appearance in the list. It is not so obvious that the algorithm naively gives the same result. 

~~~gnuplot Results are close but not exactly the same
plot 'grid2.dat' using 1:13 with line lw 4 title "Level = 11, Adapt2" ,\
     'grid2.dat' using 1:11 with line lw 4 title "Level = 9, Adapt2" ,\
     'grid3.dat' using 1:13 with line lw 2 title "Level = 11, Adapt2 reversed list" ,\
     'grid3.dat' using 1:11 with line lw 2 title "Level = 9, Adapt2 reversed list" 
~~~

Note the difference between both approaches in the numbers of cells at $level=11$ (i.e. 10500 vs 700 at $t=5$). In conclusion the algorithm seems to do what is was supposed to do. The next question is if one can find good motivation to use it. 

*/