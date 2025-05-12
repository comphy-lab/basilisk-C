/**
# Poisson Solver on a Non Uniform Grid with Neumann Boundary Conditions and a Compact Polynomial Function

We solve the Poisson equation using and Order 2 Scheme.
*/

#include "grid/quadtree.h"
#include "poisson.h"
#include "Laplacian.h"
#include "utils.h"
#define BGHOSTS 1
#define Value 0

scalar a[], b[], res[], dp[], Error[];

void solve (int depth)
{ 
  origin (-2,-2);
  int nrelax = 4;
  L0 = 4;
  init_grid(1);
  refine(((level<depth)&&(x*x+y*y<=0.25)) || level<depth-1);

  foreach() {
    
    if(x*x<=1 && y*y <=1)

      b[] = (((-10.*pow(x+Delta/2.,9)+40.*pow(x+Delta/2.,7)-60.*pow(x+Delta/2.,5)+40.*pow(x+Delta/2.,3)-10.*(x+Delta/2.))-(-10.*pow(x-Delta/2.,9)+40.*pow(x-Delta/2.,7)-60.*pow(x-Delta/2.,5)+40.*pow(x-Delta/2.,3)-10.*(x-Delta/2.))) * ((-pow(y+Delta/2,11)/11. + 5.*pow(y+Delta/2.,9)/9. -10.*pow(y+Delta/2.,7)/7. + 10.*pow(y+Delta/2.,5)/5. - 5.*pow(y+Delta/2.,3)/3. + (y+Delta/2.)) - (-pow(y-Delta/2,11)/11. + 5.*pow(y-Delta/2.,9)/9. -10.*pow(y-Delta/2.,7)/7. + 10.*pow(y-Delta/2.,5)/5. - 5.*pow(y-Delta/2.,3)/3. + (y-Delta/2.))) + ((-10.*pow(y+Delta/2.,9)+40.*pow(y+Delta/2.,7)-60.*pow(y+Delta/2.,5)+40.*pow(y+Delta/2.,3)-10.*(y+Delta/2.))-(-10.*pow(y-Delta/2.,9)+40.*pow(y-Delta/2.,7)-60.*pow(y-Delta/2.,5)+40.*pow(y-Delta/2.,3)-10.*(y-Delta/2.))) * ((-pow(x+Delta/2,11)/11. + 5.*pow(x+Delta/2.,9)/9. -10.*pow(x+Delta/2.,7)/7. + 10.*pow(x+Delta/2.,5)/5. - 5.*pow(x+Delta/2.,3)/3. + (x+Delta/2.)) - (-pow(x-Delta/2,11)/11. + 5.*pow(x-Delta/2.,9)/9. -10.*pow(x-Delta/2.,7)/7. + 10.*pow(x-Delta/2.,5)/5. - 5.*pow(x-Delta/2.,3)/3. + (x-Delta/2.))))/sq(Delta);
     
    else
      b[] = 0.;

    a[] = 0;
  }
  
  
  stats s = statsf(b);
  foreach()
    b[] -= s.sum/s.volume;  
 
  b[left] = neumann(Value);
  b[right] = neumann(Value);
  b[top] = neumann(Value);
  b[bottom] = neumann(Value);
  b.prolongation = refine_biquartic;
  boundary ({b});

  #define NITER 30
  clock_t start = clock(), iter[NITER];
  double maxres[NITER];
  const scalar lambda[] = 0.;
  struct Poisson p;
  p.a = a; p.b = b; p.alpha = unityf; p.lambda = lambda;
  scalar * lres = {res};
  
  for (int i = 0; i < NITER; i++) {
    residual ({a}, {b}, lres, &p);
    mg_cycle ({a}, lres, {dp}, relax, &p, nrelax, 1, depth());
    double max = 0.;
    foreach(reduction(max:max))
      if (fabs(res[]) > max)
	max = fabs(res[]);
    iter[i] = clock();
    maxres[i] = max;
  }
         
  double max = 0;
  static FILE *fp, *fp1, *fp2;
  fp = fopen("ErrorDomain.dat","w");
  fp1 = fopen("Residual.dat","w");
  fp2 = fopen("ErrorvsGrid.dat","a");       
  foreach(reduction(max:max)) {
    if(x*x <=1 && y*y<=1)
     Error[] = a[] - (-pow(y,10)+5.*pow(y,8)-10.*pow(y,6)+10.*pow(y,4)-5.*pow(y,2)+1.)*(-pow(x,10)+5.*pow(x,8)-10.*pow(x,6)+10.*pow(x,4)-5.*pow(x,2)+1.);
    else
     Error[] = a[];
   }

  stats s2 = statsf(Error);
  foreach(){
    Error[] -= s2.sum/s2.volume;
    fprintf(fp,"%g %g %g \n",x,y,Error[]);
    fprintf(fp1,"%g %g %g \n",x,y,res[]);
    if (fabs(Error[]) > max) max = fabs(Error[]);
  }
  fprintf (fp2, "%g %g\n",log10(pow(2,depth)),log10(normf(Error).max));
  fclose(fp);
  fclose(fp1);
  fclose(fp2);

}

int main (int argc, char ** argv)
{
  system("rm -f ErrorvsGrid.dat");
  for (int depth = 6; depth <= 10; depth++)
    solve (depth);
}

/**
## Results

~~~gnuplot Error Distribution for Residual operator
set output 'ErrorDistribution.png'
splot 'ErrorDomain.dat'
~~~

~~~gnuplot Error Distribution for Residual operator
set output 'ResidualDistribution.png'
splot 'Residual.dat'
~~~

~~~gnuplot Error Scaling of the Poisson 4th order solver on a heterogenous grid

set output 'ErrorScaling.png'
set xlabel 'Log(GridPoints)'
set ylabel 'Log(Max-Error)'
set grid
f(x) = a*x + b
fit f(x) 'ErrorvsGrid.dat' u 1:2 via a, b
title_f(a,b) = sprintf('f(x) = %.2fx + %.2f',a,b)
plot 'ErrorvsGrid.dat' u 1:2 w p t 'ErrorScaling', f(x) t title_f(a,b)
~~~

*/

