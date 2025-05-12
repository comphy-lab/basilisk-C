/**
## 1D Test case of 2nd Order Laplacian operator on a Non Uniform Grid

In this code the error scaling of the Laplacian operator has been demonstrated on a Non Uniform Grid. A compact function which includes a polynomial function of order 8 is used.
*/

#include "grid/bitree.h"
#include "utils.h"
#define BGHOST 1
#define Value 0


scalar A[],B[],Error[];

static inline double biquartic(Point point, scalar A){

 #if dimension == 1
   return ((105.*coarse(A,-2*child.x,0,0)-756*coarse(A,-child.x,0,0)+5670*coarse(A,0,0,0)+
            1260*coarse(A,child.x,0,0)-135*coarse(A,2*child.x,0,0))/(24*256));

 #elif dimension == 2
  return((105.*(105.*coarse(A,-2*child.x,-2*child.y,0)-756*coarse(A,-child.x,-2*child.y,0)
                +5670*coarse(A,0,-2*child.y,0)+1260*coarse(A,child.x,-2*child.y,0)-135*coarse(A,2*child.x,-2*child.y,0)) 
                -756.*(105.*coarse(A,-2*child.x,-child.y,0)-756*coarse(A,-child.x,-child.y,0)+5670*coarse(A,0,-child.y,0)
                +1260*coarse(A,child.x,-child.y,0)-135*coarse(A,2*child.x,-child.y,0)) + 5670.*(105.*coarse(A,-2*child.x,0,0)
                -756*coarse(A,-child.x,0,0)+5670*coarse(A,0,0,0)+1260*coarse(A,child.x,0,0)-135*coarse(A,2*child.x,0,0)) + 
                 1260.*(105.*coarse(A,-2*child.x,child.y,0)-756*coarse(A,-child.x,child.y,0)+5670*coarse(A,0,child.y,0)+
                 1260*coarse(A,child.x,child.y,0)-135*coarse(A,2*child.x,child.y,0)) - 135.*(105.*coarse(A,-2*child.x,2*child.y,0)
                 -756*coarse(A,-child.x,2*child.y,0)+5670*coarse(A,0,2*child.y,0)+1260*coarse(A,child.x,2*child.y,0)-
                  135*coarse(A,2*child.x,2*child.y,0)))/(256.*256.*24.*24.));
 #endif
}

static inline void refine_biquartic (Point point, scalar s)
{
  foreach_child()
    s[] = biquartic(point, s);
}

void solve(int depth){

   origin(-2);
   L0 = 4;
   init_grid(1);
   refine(((level<depth)&&(x*x<=0.25)) || (level<depth-1));

   foreach(){
      if(x*x<=1)
        A[] = pow(x,8) - 4*pow(x,6) + 6*pow(x,4) - 4*pow(x,2) + 1.;
      else
        A[] = 0.;
     
      B[] = 0.;
      Error[] = 0.;
     }
   A.prolongation = refine_biquartic;
   foreach_dimension(){
      A[left] = neumann(Value);
      B[left] = neumann(Value);
     }

   boundary({A,B});
   
#if TREE
   face vector g[];
   foreach_face()
      g.x[] = (A[] - A[-1])/(Delta);
   boundary_flux({g});
   foreach()
     B[] = (g.x[1]-g.x[])/(Delta);

#else

   foreach()
     B[] = (A[1]-2*A[0]+A[-1])/(sq(Delta));

#endif

   foreach(){ 
     if(x*x<=1)
       Error[] = B[] - 56*pow(x,6) + 120*pow(x,4) - 72*pow(x,2) + 8;
     else
       Error[] = B[];
    }

   FILE *fp1 = fopen("Error.dat","w");
   FILE *fp2 = fopen("Laplacian.dat","w");
   FILE *fp3 = fopen("ErrorvsGrid.dat","a");
   FILE *fp4 = fopen("G.dat","a");
   foreach(){
     fprintf(fp1,"%g %g %g %g\n",x,A[],B[],Error[]);
     fprintf(fp2,"%g %g %g \n",x,A[],B[]);
     fprintf(fp4,"%g %g \n",x,g.x[]);
   }
   fprintf(fp3,"%g %g \n",log10(pow(2,depth)),log10(normf(Error).max));
   fclose(fp1);
   fclose(fp2);
   fclose(fp3);
   fclose(fp4);

   printf("\n\n The Max Error is %g ",normf(Error).avg);
}

int main(){
  system("rm -f ErrorvsGrid.dat");
  for(int depth=6; depth<=8;depth++)
     solve(depth);
}

/**
## Results
~~~gnuplot Error Distribution of the 2nd order Laplacian scheme on a non uniform grid

set output 'ErrorDistribution.png'
set xlabel 'x'
set ylabel 'Error'
set grid
plot 'Error.dat' u 1:4 w p t 'Error Lap(g(x))'
~~~

~~~gnuplot  Computed Laplacian

set output 'Laplacian.png'
set xlabel 'x'
set ylabel 'Error'
set grid
plot 'Laplacian.dat' u 1:3 w p t 'Laplacian'
~~~

~~~gnuplot Error Scaling of the second order Laplacian Scheme
set output 'ErrorScaling.png'
set xlabel 'Log(GridPoints)'
set ylabel 'Log(Max-Error)'
set grid
f(x) = a*x + b
fit f(x) 'ErrorvsGrid.dat' u 1:2 via a, b
title_f(a,b) = sprintf('E(x) = %.2fx + %.2f',a,b)
plot 'ErrorvsGrid.dat' u 1:2 w p t 'ErrorScaling', f(x) t title_f(a,b)
~~~
*/