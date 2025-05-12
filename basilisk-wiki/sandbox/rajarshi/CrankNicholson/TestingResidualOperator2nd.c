/**
## Testing the Residual operator -  2nd order solver (Heterogenous Refinement / Quadtree)

This test case is performed to study the accuracy of the residual operator. The case is simplified, in that - the multigrid cycles have been left out. The results demonstrate that for the second order cases converged solution is obtained, hence we may conclude the residual & relaxation operators are fine. Refer now to the 4th order testing for the Residual and Relaxation functions. This is given in the file TestingResidualOperation4th.c

*/

#include "grid/tree.h"
#include "utils.h"
#define BGHOSTS 1
#define value 0

scalar A[],B[],Res[],dA[],Error[];
         
void relax2nd(scalar da, scalar residue, int NRelax){
   
   scalar c[];
   double n,d;
   for(int i=0;i<NRelax;i++){
     foreach(){
       n = -sq(Delta)*residue[];
       d = 0;
       foreach_dimension(){
          n += da[1]+da[-1];
          d += 2;
         }
       c[] = n/d;
      }
     foreach()
       da[] = (da[] + 2.*c[])/3.; 
     boundary({da}); 
   } 
}

double residual2nd (scalar a, scalar b, scalar residue)
{
  double maxres = 0.;

#if TREE
  face vector g[];
  foreach_face()
    g.x[] = (a[] - a[-1])/Delta;
  boundary_flux ({g});
  foreach () {
    residue[] = b[];
    foreach_dimension()
      residue[] += (g.x[] - g.x[1])/Delta;
    if (fabs (residue[]) > maxres)
      maxres = fabs (residue[]);
  }
#else
  foreach () {
    residue[] = b[];
    foreach_dimension()
      residue[] += (-a[1] + 2*a[] - a[-1])/sq(Delta);
    if (fabs (residue[]) > maxres)
      maxres = fabs (residue[]);
  }
#endif
  boundary ({residue});
  return maxres;
}


void solve(int depth){

   stats s,s1;
   origin(-0.5,-0.5);
   L0=1;
   double maxresidue=1;
   init_grid(1);
   refine (level < depth - 2 || level <= depth*(1. - sqrt(x*x + y*y)));

   foreach(){
      B[] = -8.*pi*pi*cos(2.*pi*x)*cos(2.*pi*y);
      A[] = 0;
      dA[] = 0;
      Error[] = 0;
    }    

/**
   As of this moment, the most suitable boundary conditions for both the cases is the Neumann boundary condition (with a slope zero) , as the dirichlet implementation for the the 4th order case has not been implemented yet and the periodic boundary conditions have not been implemented for a heterogenous refinement.
   
*/
  
   A[left] = neumann(value);
   A[right] = neumann(value);   
   A[top] = neumann(value);
   A[bottom] = neumann(value);
   B[left] = neumann(value);
   B[right] = neumann(value);
   B[top] = neumann(value);
   B[bottom] = neumann(value);  

/** 
This step is done to ensure that the volume average of the RHS of the equation is equal to zero to machine precision. This is an essential requirement to ensure compatibility of the solution. Without this step no solution can exist, given the boundary conditions
*/
  
   s = statsf(B);
   foreach()
     B[] -= s.sum/s.volume;  

   boundary({B});
   int NRelaxation = 4;
   int ctr=0;
   
   while(maxresidue > 1E-06 && ctr < 10000){
      ctr++;
      boundary({A,dA});
      maxresidue = residual2nd(A,B,Res);
      relax2nd(dA,Res,NRelaxation);
      foreach(){
        A[] += dA[];
        dA[] = 0;
       }
      printf("\n\n Iteration Number = %d & Maxresidue = %g ",ctr,maxresidue);
   }

/**
This step is done to remove the constant from the solution
*/  
  
   s1 = statsf(A);
   foreach()
    A[] = A[] - s1.sum/s1.volume;

   FILE *fp = fopen("ErrorDistribution.dat","w");
   foreach(){
     Error[] = A[]-cos(2.*pi*x)*cos(2.*pi*y);
     fprintf(fp,"%g %g %g\n",x,y,Error[]);
    }
   fclose(fp);
   
   FILE *fp1 = fopen("ErrorvsGrid.dat","a");
   fprintf(fp1,"%d %g %g %g %g\n",depth,pow(2,depth),normf(Error).max,log10(pow(2,depth)),log10(normf(Error).max));
   fclose(fp1);
   
   FILE *fp2 = fopen("Residual.dat","w");
   foreach()
      fprintf(fp2,"%g %g %g \n",x,y,Res[]);
   fclose(fp2);
}

int main(){
  system("rm -f ErrorvsGrid.dat");
  for(int depth=5; depth<=8;depth++)
     solve(depth);
}

/**
## Results

~~~gnuplot Error Distribution for Residual operator
set output 'ErrorDistribution.png
splot 'ErrorDistribution.dat'
~~~

~~~gnuplot Error Distribution for Residual operator
set output 'ResidualDistribution.png
splot 'Residual.dat'
~~~

~~~gnuplot Error Scaling of the Residual operator
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