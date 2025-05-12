/**
# Resolution of flood wave in 1D

This is a simple `C` code, not a `basilisk` one.

# Model equations
 
 
The flood wave or kinetic wave, see  Whitham p 82,  Fowler  p 76,  or Chanson, is a long wave approximation
 of shallow water equations were inertia and pressure gradient are neglected.
 There is only equilibrium between  slope and friction.
 The waves go only downstream.
We have to solve, with flux   $\bar F=\bar h^{3/2}$:
$$ \frac{\partial}{\partial \bar t} \bar h  + \frac{\partial}{\partial \bar x} \bar F = 0
 \text{ which is as well } \frac{\partial}{\partial \bar t} \bar h  + \bar c \frac{\partial}{\partial \bar x} \bar h = 0 \text{ with }  \bar c = \partial \bar F/\partial \bar h=\frac{3}{2}\sqrt{h}$$
 



~~~gnuplot
set samples 9 
set label "U i-1" at 1.5,3.1
set label "U i" at 2.5,3.15
set label "U i+1" at 3.5,2.5
set xtics ("i-2" 0.5, "i-1" 1.5, "i" 2.5,"i+1" 3.5,"i+2" 4.5,"i+3" 5.5)
set arrow from 2,1 to 2.5,1
set arrow from 3,1 to 3.5,1
set label "F i-1/2" at 2.1,1.25
set label "F i+1/2" at 3.1,1.25

set label "x i-1/2" at 1.5,0.25
set label "x i" at 2.4,0.25
set label "x i+1/2" at 3.,0.25

set label "x"  at 0.5,2+sin(0) 
set label "x"  at 1.5,2+sin(1)
set label "x"  at 2.5,2+sin(2) 
set label "x"  at 3.5,2+sin(3) 
set label "x"  at 4.5,2+sin(4) 
set label "x"  at 5.5,2+sin(5) 
p[-1:7][0:4] 2+sin(x) w steps not,2+sin(x) w impulse not linec 1
~~~

The numerical flux across face i+1/2 is denoted $F_{i+1/2}$ (or $F^n_{i+1/2}$ at time $n$), it is function (say $f$) of values 
before and after the face (i+1) which are  $U_i$ and $U_{i+1}$   
$$
 F_{i+1/2}=f(U_i,U_{i+1}).
$$
The position of the center of the cell is $x_{i}$. 


# Code
mandatory declarations:
*/
#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include <string.h>
/** definition of the field U, the flux F, time step 
*/
double*x=NULL,*h=NULL,*F=NULL;
double dt;  
double cDelta;

double L0,X0,Delta; 
double t;
int i,it=0;
int N;

/**

Here we code the flux $h^{3/3}$ as a function of the values of the left ($h_g$) and right ($h_d$) values of the face, with an approximation for the velocity $\frac{3}{2}h^{1/2}$, $c= \frac{3}{2} \frac{h_g^{1/2} + h_d^{1/2}}2$
so that
$$
F=\frac{h_g^{3/2} + h_d^{3/2}}2 - c  \frac{h_d -h_g}2
$$
  
*/  

double FR1(double hg, double hd){
    double c=1.5*(sqrt(hg)+sqrt(hd))/2;
    return (hg*sqrt(hg)+hd*sqrt(hd))*0.5-c*(hd-hg)*0.5;}


/**
Main with definition of parameters
*/
int main() {
  L0 = 12.;
  X0 = -L0/4;
  N = 256;
  t=0;
  dt = (L0/N)/2;
  Delta = L0/N;
/**
 dynamic allocation 
*/
  x= (double*)calloc(N+2,sizeof(double));
  h= (double*)calloc(N+2,sizeof(double));
  F= (double*)calloc(N+2,sizeof(double));
/** 

loop for initial `h(x,0)`: initial elevation: a "bump" 
 
 The celle `i=0` is a ghost cell/ 
The "first" (`i=1`) cell is  between  `X0` and `X0+Delta`, centred in  `Delta/2`
ith cell beween `(i-1) Delta` (left) and `i Delta`(right) centered in `(i-1/2)Delta`  

`Delta=L0/N`

*/
   for(i=0;i<N+2;i++)
    {  x[i]=X0+(i-1./2)*Delta;  
       h[i] = 0.5+exp(-x[i]*x[i]);
  }
/**
 begin the time loop
*/  
   while(t<=4){
   	
/** 
print data at t=1 t=2 t=3 
*/     
if((it == 1)||(it == (int)(1/dt))
   || (it == (int)(2/dt))
   || (it == (int)(3/dt))
   || (it == (int)(4/dt)) ){
 for(i=0;i<=N;i++)
    fprintf (stdout, "%g %g %g  \n", x[i], h[i], t);
  fprintf (stdout, "\n\n");
}
     t = t + dt;
     it++;
/** 
flux 
*/    
 for(i=1;i<=N+1;i++)
   {
    F[i] = FR1(h[i-1],h[i]);
    }
/** 
explicit step
update 
$$
 h_i^{n+1}=h_i^{n} -{\Delta t} \dfrac{F_{i+1/2}-F_{i-1/2}}{\Delta x}
$$*/ 
 for(i=1;i<=N;i++)
    h[i] +=  - dt* ( F[i+1] - F[i] )/Delta;
       
/**
     Boundary condition
*/
   h[0] = 0.5;
   h[N+1]=h[N];
 }
    
  free(h);
  free(F);
  free(x);
}
/**
## Run
Then compile and run:

~~~bash
 cc  floodwaveC.c -lm ; ./a.out> out
~~~
 

## Results
 
in gnuplot type

~~~bash
 reset
  set xlabel "x"
  set ylabel "h"
  h(x)=.5+exp(-x*x)
  c(x)=3./2*sqrt(h(x))
  set parametric
  set dummy x

  p[:][:]'out' w l , x,h(x) not w l lc -1,\
  x+1*c(x),h(x) not w l lc -1,\
  x+2*c(x),h(x) not w l lc -1,\
  x+3*c(x),h(x) not w l lc -1,\
  x+4*c(x),h(x) not w l lc -1,\
  x ,0 not w l lc -1
 
~~~

~~~gnuplot
 reset
  set xlabel "x"
  set ylabel "h"
  h(x)=.5+exp(-x*x)
  c(x)=3./2*sqrt(h(x))
  set parametric
  set dummy x

  p[:][:]'out' w l , x,h(x) not w l lc -1,\
  x+1*c(x),h(x) not w l lc -1,\
  x+2*c(x),h(x) not w l lc -1,\
  x+3*c(x),h(x) not w l lc -1,\
  x+4*c(x),h(x) not w l lc -1,\
  x ,0 not w l lc -1
 
~~~
 
 which gives $h(x,t)$ numerical and exact plotted here for t=0 1 2 3  4 and $-3<x<9$
 
 
## Note

the output is the standard one (terminal), we can save in a file and plot the generated file with `gnuplot`, see 
[http://basilisk.fr/sandbox/M1EMN/BASIC/gnuplot_examples.c]()


    FILE * fp = fopen ("out.txt", "w")
    for(i=0;i<=N;i++)
    fprintf (fp, "%g %g %g  \n", x[i], h[i], t);
    fprintf (fp, "\n\n");
    fclose(fp);
   

 
 
 
## Links
 
 * [http://basilisk.fr/sandbox/M1EMN/BASIC/advecte1c.c]()
 
 * [http://basilisk.fr/sandbox/M1EMN/BASIC/advecte1.c]()
 
 * [http://basilisk.fr/sandbox/M1EMN/Exemples/floodwave.c]()
 
 * [http://basilisk.fr/sandbox/M1EMN/BASIC/gnuplot_examples.c]()
 
## Bibliography
 
 * [Lagrée P-Y](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/MFEnv.pdf)
 "Equations de Saint Venant et application, Ecoulements en milieux naturels" Cours MSF12, M1 UPMC
 
  * [Lagrée P-Y](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/code_C_saintvenant.pdf)
  "Résolution numérique des équations de Saint-Venant,
mise en oeuvre en volumes finis par un solveur de Riemann bien balancé" Cours MSF12, M1 UPMC
 
 * G. B. Whitham "Linear and Nonlinear Waves" Wiley-Interscience, p 82
 
 * Fowler ["Mathematics of the environment"](https://people.maths.ox.ac.uk/fowler/courses/mathenvo/envonotes.pdf)
 
 
ready for new site  01/22
*/
