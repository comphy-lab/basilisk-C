/**
#Shallow Water Test Case with Flux limiter order 2 reconstruction (Error-Scalability only in the hydraulic jump region) - Le Metayer Test Case
*/

#include "grid/cartesian1D.h"
#include "saint-venant.h"

void Compute_Error(double time){
 
  scalar u_analytical[], h_analytical[],u_error[],h_error[];
  double h_star,h_plus,h_minus,sigma,u_star,u_errorL2RF,h_errorL2RF;
  int ctr;  

  h_plus = 1.;
  h_minus = 1.8;
  h_star = 1.36898;  
  u_star = 2.*(sqrt(G*h_minus) - sqrt(G*h_star));
  sigma = (h_star*u_star)/(h_star - h_plus);

  foreach(){

    if(x < -1.*time*sqrt(G*h_minus)){
        u_analytical[] = 0.;
        h_analytical[] = h_minus;
      }
    else if( x > -1.*time*sqrt(G*h_minus) && x < time*(u_star - sqrt(G*h_star)) ){
        u_analytical[] = (2.*sqrt(G*h_minus) + 2.*x/time)/3.;
        h_analytical[] = sq(2.*sqrt(G*h_minus) - x/time)/(9.*G);
      }
    else if(x > time*(u_star - sqrt(G*h_star)) && x < time*sigma ){
        u_analytical[] = u_star;
        h_analytical[] = h_star;  
      }
    else if(x > time*sigma){
        u_analytical[] = 0.;
        h_analytical[] = h_plus;
      }
    u_error[] = u.x[] - u_analytical[];
    h_error[] = h[] - h_analytical[];
   }

   u_errorL2RF = 0.;
   h_errorL2RF = 0.;
   ctr = 0;
   foreach(){
     if(x>=-100 && x<=100){
        ctr++;
        u_errorL2RF += sq(u_error[]);
        h_errorL2RF += sq(h_error[]);
       }
     }
   u_errorL2RF = sqrt(u_errorL2RF/ctr);
   h_errorL2RF = sqrt(h_errorL2RF/ctr);

   FILE *fp1 = fopen("ErrorvsGrid.dat","a");
   fprintf(fp1,"%g %g %g \n",log10(N),log10(u_errorL2RF),log10(h_errorL2RF));  
   fclose(fp1);
   
}

int main()
{
  X0 = -300.;
  L0 = 600.;
  G = 10.;
  system("rm -f ErrorvsGrid.dat");
  for(N=64;N<=2048;N*=2)
    run();
}

event init (i = 0)
{
  foreach()
    h[] = x < 0. ? 1.8 : 1.;
}

event output (t = 48) {
  FILE *fp2 = fopen("Numerical_Solution.dat","w");
  foreach()
    fprintf (fp2, "%g %g %g\n", x, h[], u.x[]);
  fclose(fp2);
  Compute_Error(t);
}

/**

## Results

~~~gnuplot Error Scaling of the Weno 3rd order reconstruction (on velocity function)
set output 'ErrorScaling_velocity.png'
set xlabel 'Log(GridPoints)'
set ylabel 'Log(L2Norm-Error-velocity)'
set grid
f(x) = a*x + b
fit f(x) 'ErrorvsGrid.dat' u 1:2 via a, b
title_f(a,b) = sprintf('f(x) = %.2fx + %.2f',a,b)
plot 'ErrorvsGrid.dat' u 1:2 w p t 'Error-Scaling-of-velocity', f(x) t title_f(a,b)
~~~

~~~gnuplot Error Scaling of the Weno 3rd order reconstruction (on height function)
set output 'ErrorScaling_height.png'
set xlabel 'Log(GridPoints)'
set ylabel 'Log(L2Norm-Error-height)'
set grid
f(x) = a*x + b
fit f(x) 'ErrorvsGrid.dat' u 1:3 via a, b
title_f(a,b) = sprintf('f(x) = %.2fx + %.2f',a,b)
plot 'ErrorvsGrid.dat' u 1:3 w p t 'Error-Scaling-of-Height', f(x) t title_f(a,b)
~~~

*/

