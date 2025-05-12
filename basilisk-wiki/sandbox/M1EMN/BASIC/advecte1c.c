/**
# Resolution of Advection equation in 1D

it is exactly the same than [advecte1.c]() but written with standard C
 
 
 we solve 
 $$\frac{\partial h}{\partial t}+\frac{\partial (a h)}{\partial x} = 0$$
 with $a$ constant. 
 we write it 
 $$\frac{\partial U}{\partial t}+\frac{\partial F(U)}{\partial x} = 0$$
 
# Resolution of Advection equation in 1D
 
All the problem consists to solve
$$\frac{\partial U}{\partial t}+\frac{\partial F(U)}{\partial x} = 0$$
the advection equation. We consider here the simple case $F=a U$
 

Space is decomposed in small segments of  a priori different length, but here we suppose that the length is a constant $\Delta x$.
The same for time increment: $\Delta t$ is constant.  
Those small segments are kind of  "volumes", as we are in  dimension one.


Consider the system 
$$\frac{\partial U}{\partial t}+\frac{\partial F(U)}{\partial x} = 0$$
and integrate in $x$ on a small interval between $x_i-\Delta x /2$ and $x_{i}+\Delta x /2$  and consider two times  $t^n$ and $t^{n+1}$.
 

A mean value of $U$ around $x_i$ between $x_i-\Delta x /2$ and 
 $x_i+\Delta x/2$ may be defined as 
  $U_i^n$ the integral between $x_i-\Delta x /2$ and
 $x_i+\Delta x/2$ is so by definition
$$
 U_i^n =\dfrac{1}{\Delta x} \int_{x_{i}-\Delta x /2}^{x_{i}+\Delta x /2} U(x,t_n)dx
$$
index $i$ is for the segment
$C_i=(x_{i}-\Delta x /2,x_{i}+\Delta x /2)$,  centered in $x_{i}$
 
index $n$ corresponds to time  $t_n$ with $t_{n+1}-t_{n}=\Delta t$. So that with the definition:
$$
 U_i^{n+1} =\dfrac{1}{\Delta x} \int_{x_{i}-\Delta x /2}^{x_{i}+\Delta x /2} U(x,t_{n+1})dx
$$


Hence by the integration on the "segment/ Volume" of the system,
we have if we integrate :

 $\int_{x_i-\Delta x /2}^{x_{i}+\Delta x/2} \frac{\partial U}{\partial t}  dx  =\frac{d}{dt}  \int_{x_{i}-\Delta x /2}^{x_{i}+\Delta x /2} U dx$
 and for the flux 
 $\int_{x_{i}-\Delta x /2}^{x_{i}+\Delta x /2} \partial_x F(U) dx = F^n_{i+1}-F^n_{i}.$
 This is  $\frac{d}{dt}  \int_{x_{i}-\Delta x/2}^{x_{i}+\Delta x/2} U dx = F^n_{i+1}-F^n_{i}.$
 We write it as the up to now exact expression:
 $$\frac{d}{dt}(\frac{1}{\Delta x} \int_{x_{i}-\Delta x/2}^{x_{i}+\Delta x/2} U dx )
 + \frac{F^n_{i+1}-F^n_{i}}{\Delta x}=0.$$
 We insisit, this is an "exact" integration from the flux on the  "volume" for the mean value. That is the finite volume method. 


 Taylor expansion:
$$ U(t + \Delta t) = U(t )+\Delta t \partial_t U  + (1/2)(\Delta t)^2 \partial^2_t U  +O(\Delta t)^3, $$
allows us to write the discrete finite volume problem (at first order in time) 
$$
 \dfrac{U_i^{n+1}-U_i^{n}}{\Delta t}+\dfrac{F^n_{i+1}-F^n_{i}}{\Delta x}=0,
$$

Numerical flux $F_{i+1}$ is an approximation of  $F(U)$ at right interface of segment  $C_i$ centered in $i$, the more pedagogical expression is
$F_{i+1/2}$, as it is a function of the value of  $U_i$ in the considered segment, which begins in  $i-1/2$ 
 and of the   value $U_{i+1}$ from the next one which begins in $i+1/2$:

~~~gnuplot
set samples 9 
set label "U i-1" at 1.5,3.1
set label "U i" at 2.5,3.15
set label "U i+1" at 3.5,2.5
set xtics ("i-2" 0.5, "i-1" 1.5, "i" 2.5,"i+1" 3.5,"i+2" 4.5,"i+3" 5.5)
set arrow from 2,1 to 2.5,1
set arrow from 3,1 to 3.5,1
set label "F i" at 2.1,1.25
set label "F i+1" at 3.1,1.25

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

The numerical flux across face i+1/2 is denoted $F_{i+1}$ (or $F^n_{i+1}$ at time $n$), it is function (say $f$) of values 
before and after the face (i+1) which are  $U_i$ and $U_{i+1}$   
$$
 F_{i+1}=f(U_i,U_{i+1}).
$$
The position of the center of the cell is $x_{i}$. 

So, if the length of domain is $L$, and if the domain starts in $x_0$, and if we take $N$ points, hence $\Delta x=L/N$
faces are $x_0 + (i-1) \Delta x$.

A the center of the cell, $x_{i}=x_0 + (i-1/2) (\Delta x)$, we have the mean value $U^n_i$. 


The finite volume method is 
$$
 \dfrac{U_i^{n+1}-U_i^{n}}{\Delta t}+\dfrac{F^n_{i+1}-F^n_{i}}{\Delta x}=0,
$$
It is explicit:  compute new values  $U_i^{n+1}$ as a function of the old ones $U_i^{n}$ .
 
Nota 1:

Do not confuse the $f(,)$ function, numerical flux across face $F_i$ and actual flux function $F$ comming from the physics.

Nota 2:

Attention,this is a simplified point of view for the faces, in  Basilisk, for  faces we should use 
"foreach_face()" and "face vector"

Nota 3:

Presentation may defer by changing from $x_i$ to $x_{i+1/2}$ 

### Flux  
 
The most simple expression for   $f$would be to take the mean value 
$$
 F_{i}=f(U_{i-1},U_{i})=\dfrac{F(U_{i-1})+F(U_i)}{2}
$$
so that 
$$
 U_i^{n+1}=U_i^{n} -{\Delta t} \dfrac{F(U_{i+1})-F(U_{i-1})}{2 \Delta x}
$$
but this is not a good idea...
It is unstable.  

 
### Flux  upwind
Here we propose to write the flux with the correction with $c_\Delta=dF/dU$
$$
F_i=f(U_{i-1},U_{i})=\dfrac{F(U_{i-1})+F(U_i)}{2} 
  -  c_\Delta \frac{({(U_{i})-(U_{i-1})})}{2}
$$  
 
In this case, as $F=a U$, $c_\Delta =a$. 
So that the flux on face centered in $i+1/2$ minus the flux on the face $i-1/2$ correspondig to the evolution in the volume $i$ is :
$$F_{i+1}-F_{i}= a\dfrac{(U_{i})+(U_{i+1})}{2}-a\dfrac{(U_{i-1})+(U_i)}{2} -  (a)( \frac{U_{i+1}-U_{i}}{2}-\frac{U_{i}-U_{i-1}}{2})=a(U_{i}-U_{i-1})$$
so that the temporal evolution in the volume $i$ du to the flux at his faces is
$$
 U_i^{n+1}=U_i^{n} -{\Delta t} \dfrac{U^{n}_{i}-U^{n}_{i-1}}{ \Delta x}
$$
This scheme is "consistant" (when ($\Delta t,\Delta x$, etc.) we reobtain the initial PDE). It is stable, and convergent.
  
 
 
 
 
## Code
mandatory declarations:
*/
#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include <string.h>
/** definition of the field U, the flux F, time step 
*/
double*x=NULL,*U=NULL,*F=NULL;
double dt;  
double cDelta;

double L0,X0,Delta; 
double t,a;
int i,it=0;
int N;

/**
Main with definition of parameters, note that CFL=1
*/
int main() {
  L0 = 12.;
  X0 = -L0/4;
  N = 128;    
  t=0;
  a=1.2;
  Delta = L0/N;
  dt = (L0/N)/a;
  
/**
 dynamic allocation of N cells + 2 ghosts 
*/
  x= (double*)calloc(N+2,sizeof(double));
  U= (double*)calloc(N+2,sizeof(double));
  F= (double*)calloc(N+2,sizeof(double));

  
    for(i=0;i<=N+1;i++)
/** 

ghost cell left i=0: between  `X0-Delta` and `X0`, centred in  `-Delta/2`
first cell i=1  between  `X0` and `X0+Delta`, centred in  `Delta/2`
ith cell beween `(i-1) Delta` (left) and `i Delta`(right) centered in `(i-1/2)Delta`  

`Delta=L0/N`

*/
    {  x[i]=X0+(i-1./2)*Delta;  
/**
initial elevation: a "bump", position 
*/
       U[i] = exp(-x[i]*x[i]);
  }
/**
 begin the time loop
*/  
   while(t<=3){  
/** 
print data
*/
// printdata at t = 1 2 and 3  
if( (it == (int)(1/dt)) || (it == (int)(2/dt)) || (it == (int)(3/dt))  ){
 for(i=0;i<N;i++)
    fprintf (stdout, "%g %g %g  \n", x[i], U[i], t);
  fprintf (stdout, "\n\n");
}    
     t = t + dt;
     it++;	
/** 
flux either $c=1$ upwind
$c=0$ unstable
$c=Delta/dt$ Lax Wendrof
*/    
 for(i=1;i<=N+1;i++)
   {
    //cDelta = Delta/dt;   // Lax Wendrof
   // cDelta = 0;          // Unstable
    cDelta = a;          // Upwind
    F[i] = a*(U[i]+U[i-1])/2.  - cDelta *(U[i]-U[i-1])/2;
    //F[i] = a*U[i-1]; // same!!!
    }
/** 
explicit step
update 
$$
 U_i^{n+1}=U_i^{n} -{\Delta t} \dfrac{F(U_{i+1})-F(U_{i})}{\Delta x}
$$*/ 
 for(i=0;i<=N;i++)
    U[i] +=  - dt* ( F[i+1] - F[i] )/Delta;
     
/**
 Boundary condition at the "ghost cells"
*/
  U[0] = U[1];
  U[N+1]=U[N];      
 }
   
 // clear memory (C instruction)

  free(U); 
  free(F);
  free(x);
}
/**
## Run
Then compile and run in a simple terminal:

~~~bash
 cc  -g -O2 -DTRASH=1 -Wall  advecte1c.c -o advecte1c ;./advecte1c > out
~~~

or better in basilisk point of view

~~~bash
 make advecte1c.tst;make advecte1c/plots    
 make advecte1c.c.html ; open advecte1c.c.html 
~~~


## Results
One analytical solution is 
$$U(x,t) =  exp(-(x-t)^2)$$
in gnuplot terminal type

~~~bash
 U(x,t)= exp(-(x-t)*(x-t))
 p'out' u ($1):($2)t'num'w l,'' u 1:(U($1,$3)) t'exact' w l
~~~
which gives $U(x,t)$ plotted here for t=0 1 2 3  and $-3<x<9$ 

note that if CFL=1, i.e. here $\Delta t=\Delta x$, the result is more precise...
~~~gnuplot
reset
 set xlabel "x"
 a=1.2
 U(x,t)= exp(-(x-a* t)*(x-a*t))
 p'out' u ($1):($2)t'num'w p,'' u 1:(U($1,$3)) t'exact' w l
~~~

## Exercice
Compare with the Basilisk [advecte1.c]()

Compare with [advecte.c]()



ready for new site 09/05/19
*/