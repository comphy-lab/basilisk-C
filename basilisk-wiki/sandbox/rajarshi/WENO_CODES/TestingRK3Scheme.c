/**
#Testing Order of RK3 scheme
*/

#define BGHOSTS 2
#include "grid/multigrid1D.h"
#include "utils.h"

void compute(double dt){

 double t; 
 t=0; 

 scalar h[], h1[], h2[], h_analytical[],Error[];
 foreach(){
    if(x<-1)
       h[]=0.;
    else if(x>=-1 && x<=1)
       h[]=pow((1-sq(x)),5);
    else if(x>=1)
       h[]=0.;
    }
 h[left] = neumann(0);
 h[right] = neumann(0);
 boundary({h});

 
 while(t<1){
   foreach()
     h1[] = h[] - h[]*dt/3.;
   h1[left] = neumann(0);
   h1[right] = neumann(0);
   boundary({h1});
   
   foreach()  
     h2[] = h[] - dt*h1[]/2.;
   h2[left] = neumann(0);
   h2[right] = neumann(0);
   boundary({h2});
   
   foreach()
     h[] = h[] - h2[]*dt;
   boundary({h});
  
   t +=dt;
  }
  
 foreach(){

    if(x<-1)
       h_analytical[]=0.;
    else if(x>=-1 && x<=1)
       h_analytical[]=pow((1-sq(x)),5)*exp(-1.*t);
    else if(x>1)
       h_analytical[]=0.;

    Error[] = h[] - h_analytical[];
   }
  
  FILE *fp = fopen("Error.dat","w");
  foreach()
    fprintf(fp,"%g %g %g %g \n",x,h[],h_analytical[],Error[]);
  fclose(fp);

  fp = fopen("ErrorvsGrid.dat","a");
  fprintf(fp,"%g %g \n",log10(dt),log10(normf(Error).rms));
  fclose(fp);
}

int main(){

 L0=4;
 origin(-2);
 init_grid(1<<12);
 system("rm -f ErrorvsGrid.dat");
 for(double dt=0.1; dt>=0.0002; dt=dt/2.)
    compute(dt);
}

/**
~~~gnuplot Analytical vs Numeric Solution
set output 'NumericVSAnalytic.png'
set xlabel 'X'
set ylabel 'H(X)'
set grid
plot 'Error.dat' u 1:2 w p t 'h-RK3', 'Error.dat' u 1:3 w p t 'h-analytical'
~~~

~~~gnuplot Error in Numerics
set output 'Error.png'
set xlabel 'X'
set ylabel 'Error(X)'
set grid
plot 'Error.dat' u 1:4 w p t 'h-RK3'
~~~

~~~gnuplot Error Scaling of the 3rd Order Runge Kutta Scheme
set output 'ErrorScaling_RK3.png'
set xlabel 'Log(dt)'
set ylabel 'Log(L2Norm-Error)'
set grid
f(x) = a*x + b
fit f(x) 'ErrorvsGrid.dat' u 1:2 via a, b
title_f(a,b) = sprintf('f(x) = %.2fx + %.2f',a,b)
plot 'ErrorvsGrid.dat' u 1:2 w p t 'Error-Scaling-RK3', f(x) t title_f(a,b)
~~~
*/
 