/**
# Resolution  
lubrication slump equation in 1D
it is exactly the same than [../Exemples/viscous_collapse_noSV.c]() but written with standard C

## Code
mandatory declarations: */
#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include <string.h>
/** definition of the field h, the flux Q, time step */
double*x=NULL,*h=NULL,*Q=NULL;
double dt,L0,Delta; 
double t;
int i,N;
/** Main with definition of parameters, `Delta=L0/N`  */
int main() {
  L0 = 5.;
  N = 128;    
  t=0;
  Delta = L0/N;
  dt =.0025;
/**  dynamic allocation  */
  x= (double*)calloc(N+1,sizeof(double));
  h= (double*)calloc(N+1,sizeof(double));
  Q= (double*)calloc(N+1,sizeof(double));
/** 
first cell between  `0-Delta/2` and `0+Delta/2`, centred in 0. The
ith cell is beween `(i-1/2) Delta` (left) and `(i +1/2) Delta`(right) centered in `(i)Delta`. Initial elevation: a "bump" */
  for(i=0;i<=N;i++)
    {  x[i]=0+(i)*Delta;  
       h[i] = (1)*(x[i]<1);}
/**  begin the time loop */  
   while(t<=100){ 
     t = t + dt;	
/**  print data  */
 for(i=0;i<=N;i++)
     fprintf (stdout, "%g %g %g  \n", x[i], h[i], t);
   fprintf (stdout, "\n");
/**      flux   between face `ì+1/2`, depending from either side `h[i]` and `h[i-1]`  */
 for(i=1;i<=N;i++)
    Q[i] =  - 1./3*pow(((h[i]+h[i-1])/2),3)*(h[i]-h[i-1])/Delta;
/**  explicit step  update and BC$$h_i^{n+1}=h_i^{n} -{\Delta t} \dfrac{Q_{i+1}-Q_{i}}{\Delta x}$$ */
 for(i=1;i<N-1;i++)
    h[i] +=  - dt* ( Q[i+1] - Q[i] )/Delta;
  h[0]=h[1];
  h[N]=h[N-1];
  }
/** clean */  
  free(h); 
  free(Q);
  free(x);
}

/**

## Run
Then compile and run:

~~~bash
 cc  -g -O2 -DTRASH=1 -Wall  vsicouscollapseC.c -o vsicouscollapseC ;./vsicouscollapseC > out
~~~

## Results
in gnuplot type

~~~gnuplot
set xlabel "x"
set ylabel "h"
 p[0:4]'out' u 1:2 ev :1000 w l
~~~
 
~~~gnuplot
 set xlabel "x"
 set ylabel "h"
 p[0:4][0:1.25](x<1?1:0) t 't=  0','out' u 1:($3==1? $2:NaN) t't=  1' w l,\
   ''u 1:($3==10? $2:NaN) t't= 10' w l,''u 1:($3==100? $2:NaN)  t't=100' w lp
~~~
 
 
 
~~~gnuplot
set xlabel "x/t**(1/5)"
set ylabel "h*t**(1/5)"
 p[0:2]'out' u ($1/$3**.2):($2*($3)**.2) ev :1000 w l
~~~
 
 
~~~gnuplot
 set key left
 set ytics 20
 set xlabel "x"
 set ylabel "t"
 set view 60,150
 set hidden3d
 sp[-4:4][][0:1]'out' u 1:3:2 ev :1000 t'h(x,t)' w l ,'' u (-$1):3:2 ev :1000  not  w l
~~~

## Exercice
 
Compare with the Basilisk [../Exemples/viscous_collapse_noSV.c]()

 * Lagrée  [M2EMN
Master 2 Ecoulements en Milieu Naturel](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/mainM2EMN.pdf)

*/
