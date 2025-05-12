/**
# Poisson Solver on a Uniform Grid with Neumann Boundary Conditions and a Compact Polynomial Function - ORDER 4

We solve the Poisson equation using an Order 4 Scheme.
*/

#include "grid/multigrid.h"
#include "Poisson-helmholtz4th.h"
#include "utils.h"
#define BGHOSTS 2
#define Value 0

void Solve(int depth){
  
   L0=4;
   origin(-2,-2);
   init_grid(1<<depth);

   scalar A[], B[], Res[], dA[];
   foreach(){
     if(x*x<=1 && y*y <=1)

      B[] = (((-10.*pow(x+Delta/2.,9)+40.*pow(x+Delta/2.,7)-60.*pow(x+Delta/2.,5)+40.*pow(x+Delta/2.,3)-10.*(x+Delta/2.))-(-10.*pow(x-Delta/2.,9)+40.*pow(x-Delta/2.,7)-60.*pow(x-Delta/2.,5)+40.*pow(x-Delta/2.,3)-10.*(x-Delta/2.))) * ((-pow(y+Delta/2,11)/11. + 5.*pow(y+Delta/2.,9)/9. -10.*pow(y+Delta/2.,7)/7. + 10.*pow(y+Delta/2.,5)/5. - 5.*pow(y+Delta/2.,3)/3. + (y+Delta/2.)) - (-pow(y-Delta/2,11)/11. + 5.*pow(y-Delta/2.,9)/9. -10.*pow(y-Delta/2.,7)/7. + 10.*pow(y-Delta/2.,5)/5. - 5.*pow(y-Delta/2.,3)/3. + (y-Delta/2.))) + ((-10.*pow(y+Delta/2.,9)+40.*pow(y+Delta/2.,7)-60.*pow(y+Delta/2.,5)+40.*pow(y+Delta/2.,3)-10.*(y+Delta/2.))-(-10.*pow(y-Delta/2.,9)+40.*pow(y-Delta/2.,7)-60.*pow(y-Delta/2.,5)+40.*pow(y-Delta/2.,3)-10.*(y-Delta/2.))) * ((-pow(x+Delta/2,11)/11. + 5.*pow(x+Delta/2.,9)/9. -10.*pow(x+Delta/2.,7)/7. + 10.*pow(x+Delta/2.,5)/5. - 5.*pow(x+Delta/2.,3)/3. + (x+Delta/2.)) - (-pow(x-Delta/2,11)/11. + 5.*pow(x-Delta/2.,9)/9. -10.*pow(x-Delta/2.,7)/7. + 10.*pow(x-Delta/2.,5)/5. - 5.*pow(x-Delta/2.,3)/3. + (x-Delta/2.))))/sq(Delta);
     
    else
        B[] = 0.;
       
     A[] = 0.;
    }
   
   stats s1 = statsf(B);
   foreach()
      B[] -= s1.sum/s1.volume;


      A[left] = neumann(Value);
      A[right] = neumann(Value);
      B[left] = neumann(Value);
      B[right] = neumann(Value);
      B[top] = neumann(Value);
      B[bottom] = neumann(Value);

   boundary({A,B});

   double maxres = 1E-08;
   mgstats Sol = poisson(A,B, tolerance = maxres);
   scalar Error[];
   foreach(){
     if(x*x<=1 && y*y<=1)  
       Error[] = A[] - (-pow(y,10)+5.*pow(y,8)-10.*pow(y,6)+10.*pow(y,4)-5.*pow(y,2)+1.)*(-pow(x,10)+5.*pow(x,8)-10.*pow(x,6)+10.*pow(x,4)-5.*pow(x,2)+1.);
     else
       Error[] = A[];   
   }

   stats s = statsf(Error);
   foreach()
     Error[] -= s.sum/s.volume;

   static FILE *fp,*fp2;
   fp2 = fopen("ErrorDomain.dat","w");
   foreach()
    fprintf(fp2,"%g %g %g\n",x,y,Error[]);
   fclose(fp2);

   fp = fopen("ErrorvsGrid.dat","a");
   fprintf(fp,"%d %g %g %g %g %g \n",depth,pow(2,depth),normf(Error).max,log10(pow(2,depth)),log10(normf(Error).max),s.sum/s.volume);
   fclose(fp);
}

int main(){
  
 system("rm -f ErrorvsGrid.dat");
 for(int depth=6;depth<=10;depth++)
    Solve(depth);  
}

/**
## Results

~~~gnuplot Error Distribution for Residual operator
set output 'ErrorDistribution.png'
splot 'ErrorDomain.dat'
~~~

~~~gnuplot Error Scaling of the Poisson 4th order solver
set output 'ErrorScaling.png'
set xlabel 'Log(GridPoints)'
set ylabel 'Log(Max-Error)'
set grid
f(x) = a*x + b
fit f(x) 'ErrorvsGrid.dat' u 4:5 via a, b
title_f(a,b) = sprintf('f(x) = %.2fx + %.2f',a,b)
plot 'ErrorvsGrid.dat' u 4:5 w p t 'ErrorScaling', f(x) t title_f(a,b)
~~~

*/