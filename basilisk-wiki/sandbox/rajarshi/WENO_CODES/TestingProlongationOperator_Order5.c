/**
#Testing the Order 5 prolongation operator convergence slope
*/

#include "grid/tree.h"
#include "Adaptivity_New.h"
#include "utils.h"

#define BGHOSTS 2
#define dimension 2


#if dimension == 1
double SolutionCompact1D (double x, double Delta){
  double xq1, xq2;
  xq1 = x - Delta/(2.*sqrt(3));
  xq2 = x + Delta/(2.*sqrt(3));
  if(sq(x) <=1)
    return ( 0.5 * ( pow((1-sq(xq1)),5) + pow((1-sq(xq2)),5) )  );
  else
    return (0);
}
#elif dimension == 2
double SolutionCompact2D (double x, double y, double Delta){
  double xq1, xq2, yq1, yq2;
  xq1 = x - Delta/(2.*sqrt(3));
  xq2 = x + Delta/(2.*sqrt(3));
  yq1 = y - Delta/(2.*sqrt(3));
  yq2 = y + Delta/(2.*sqrt(3));
  
  if (sq(x)<=1 && sq(y)<=1)
     return ( 0.25*( pow((1-sq(xq1)),5)*pow((1-sq(yq1)),5) + pow((1-sq(xq1)),5)*pow((1-sq(yq2)),5) + pow((1-sq(xq2)),5)*pow((1-sq(yq1)),5) + pow((1-sq(xq2)),5)*pow((1-sq(yq2)),5) ) );
  else
     return (0);
}
#endif


void convergence (int depth){

  FILE *fp;
  L0 = 4;
  origin (-2,-2);
  init_grid(1<<depth);
  scalar s[];
 
  foreach()
     s[] = SolutionCompact2D (x,y,Delta);

  if(depth==4){
     fp = fopen("BeforeRefine.dat","w");
     foreach()
        fprintf(fp,"%g %g %g\n",x,y,s[]);
     fclose(fp);
   }

  s.refine = refine_order5;
  s.prolongation = refine_order5;
  foreach_dimension(){
     s[left] = neumann(0);
     s[right] = neumann(0);
   }
  boundary({s});
  
  refine((level < depth+1 ) && (sq(x) + sq(y) <= 1) );
  refine((level < depth+2 ) && (sq(x) + sq(y) <= 0.25) );
 
  if(depth==4){
     fp = fopen("AfterRefine.dat","w");
     foreach()
       fprintf(fp,"%g %g %g \n",x,y,s[]);
     fclose(fp);
   }

  scalar S_An[], Error[];
  foreach(){
     S_An[] = SolutionCompact2D (x,y,Delta);
     Error[] = s[] - S_An[];
   }

  fp = fopen("ErrorvsGrid.dat","a");
  fprintf(fp,"%g %g \n",log10(sqrt(grid->tn)),log10(normf(Error).max));
  fclose(fp);
} 

int main(){
  system ("rm -f ErrorvsGrid.dat");
  Prolongation_Weight_Initialization();
  for (int depth = 5; depth<=8; depth++)
     convergence(depth); 
}

/**
## Results

~~~gnuplot Cell values - Before refine
set output 'Before_Refine.png'
set xlabel 'X'
set ylabel 'Y'
set zlabel 'Sol(X,Y)'
set grid
splot 'BeforeRefine.dat' u 1:2:3 w p t 'Before-Refine'
~~~

~~~gnuplot Cell values - After refine
set output 'After_Refine.png'
set xlabel 'X'
set ylabel 'Y'
set zlabel 'Sol(X,Y)'
set grid
splot 'AfterRefine.dat' u 1:2:3 w p t 'After-Refine'
~~~

~~~gnuplot Error Convergence of the 2D Prolongation operator
set output 'ErrorConvergence.png'
set xlabel 'Log(GridPoints)'
set ylabel 'Log(Max-Error)'
set grid
f(x) = a*x + b
fit f(x) 'ErrorvsGrid.dat' u 1:2 via a, b
title_f(a,b) = sprintf('f(x) = %.2fx + %.2f',a,b)
plot 'ErrorvsGrid.dat' u 1:2 w p t 'ErrorConvergence', f(x) t title_f(a,b)
~~~

*/