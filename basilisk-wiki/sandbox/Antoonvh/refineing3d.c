/**
# A Test of prolongation attributes in 3 dimensions

We test the accuracy of some prolongation attributes for a function $c(x,y,z)$ whos terms are linearly dependent on the $x,y$ and $z$ coordinate,  

$$c(x,y,z) = 1+x+y+z+xy+yz+xz+xyz. $$

*/

#include "grid/octree.h"
#include "utils.h"
#include "additional_attributes.h"

#define func (1.+x+y+z+x*y+x*z+y*z+x*y*z)

scalar c[],err[];
c[top]   = dirichlet(func); // On coarse grid we limit the effect of the boundaries 
c[bottom]= dirichlet(func);
c[left]  = dirichlet(func);
c[right] = dirichlet(func);
c[front] = dirichlet(func);
c[back]  = dirichlet(func);
FILE * fp1;
int cells = 0;
stats s;
int maxlevel = 8;
int main(){
  fp1=fopen("data.dat","w");
  L0=10.;
  X0=Y0=-L0/2;
  init_grid(1<<maxlevel);
  foreach()
    c[]=func;

  for (int n=1;n<5;n++){
    cells=0;
    foreach(reduction(+:cells))
      cells++;
    fprintf(ferr,"%d\t %d\n",n,cells);
    fprintf(fp1,"%d\t%g\t",n,L0/(pow(2.,(double)(maxlevel-n+1.))));
    for (int j=1;j<6;j++){
      if (j==1) c.prolongation=refine_injection;
      if (j==2) c.prolongation=refine_bilinear;
      if (j==3) c.prolongation=refine_linear;
      if (j==4) c.prolongation=refine_quadratic; //Third order?
      if (j==5) c.prolongation=refine_3x3x3_linear;
      wavelet(c,err);
      foreach(){// We are not interested in cells at the boundary. 
	if (abs(x)<3. && abs(y)<3. &&abs(z)<3.) 
	  err[]=fabs(err[]);
	else
	  err[]=0.;
      }
      s=statsf(err);
      fprintf(fp1,"%g\t",s.sum+0.00000001);
    }
    fprintf(fp1,"\n");
    unrefine(level>=(maxlevel-n)); 
  }
}
/**
~~~gnuplot Results from various interpolation/upsample techniques
  set xr [0.02:0.5]
  set yr [0.000000005:200]
  set logscale y
  set logscale x
  set xlabel '{/Symbol D}'
  set ylabel 'Sum of Error + 10^{-8}'
  set key center right box 1
  set size square
  plot   (200*x**1) lw 3 lc rgb 'purple' title '{/Symbol \265}{/Symbol D}^{1}' ,\
   (100*x**2) lw 3 lc rgb 'blue' title '{/Symbol \265}{/Symbol D}^{2}' ,\
   (50*x**3) lw 3 lc rgb 'blue' title '{/Symbol \265}{/Symbol D}^{3}' ,\
   0.00000001 lw 3 lc rgb 'red' title '{Exact!}' ,\
  'data.dat' using 2:3 title 'Injection',\
  'data.dat' using 2:4 title 'Bilinear' ,\
  'data.dat' using 2:5 title 'Linear' ,\
  'data.dat' using 2:6 title 'quadratic' ,\
  'data.dat' using 2:7 title '3x3x3-Linear'
  
~~~   

As was found in 2D, the *injection* method is only first order accurate. As expected the second-order accurate *bilinear* technique is exact for first-order functions. However, the *linear* technique does not seem to display this property. This is because it neglecs variations that are not grid alligned (i.e. introduced by the last four terms of the equation above). The *quadratic* furmulation seems to perform better, but is not exact. The *3x3x3 linear* that was designed to be exact for these functions is indeed exact. So apperently the weights are typed over correctly!    
*/