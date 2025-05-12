/**
# Resolution of Advection equation in 1D
 
All the problem consists to solve
$$\frac{\partial U}{\partial t}+\frac{\partial F(U)}{\partial x} = 0$$
the advection equation. We consider here the simple case $F=U$
 

Space is decomposed in small sgments of  a priori different length, but here we suppose that the length is a constant $\Delta x$.
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


Hence by the integration on the "segment/ Volume" of 
$$\frac{\partial U}{\partial t}+\frac{\partial F(U)}{\partial x} = 0$$
we have if we integrate :

 $\int_{x_i-\Delta x /2}^{x_{i}+\Delta x/2} \frac{\partial U}{\partial t}  dx  =\frac{d}{dt}  \int_{x_{i}-\Delta x /2}^{x_{i}+\Delta x /2} U dx$
 and for the flux 
 $\int_{x_{i}-\Delta x /2}^{x_{i}+\Delta x /2} \partial_x F(U) dx = F^n_{i+1}-F^n_{i}.$
 This is  $\frac{d}{dt}  \int_{x_{i}-\Delta x/2}^{x_{i}+\Delta x/2} U dx = F^n_{i+1}-F^n_{i}.$
 We write it as:
 $$\frac{d}{dt}(\frac{1}{\Delta x} \int_{x_{i}-\Delta x/2}^{x_{i}+\Delta x/2} U dx )
 + \frac{F^n_{i+1}-F^n_{i}}{\Delta x}=0.$$
 
 This is an "exact" integration from the flux on the  "volume" for the mean value. That is the finite volume method. 


 Taylor expansion:
$$ U(t + \Delta t) = U(t )+\Delta t \partial_t U  + (1/2)(\Delta t)^2 \partial^2_t U  +O(\Delta t)^3, $$
allows us to write the approximation (at first order in time) 
$$
 \dfrac{U_i^{n+1}-U_i^{n}}{\Delta t}+\dfrac{F^n_{i+1}-F^n_{i}}{\Delta x}=0,
$$

Numerical flux $F_{i+1}$ is an approximation of  $F(U)$ at right interface of segment  $C_i$ centered in $i$.
 It is a function of the value of  $U_i$ in the considered segment, which begins in  $i-1/2$ 
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

Attention, cette présentation est simplifiée, dans la version complète de Basilisk, les faces sont traîtées avec 
"foreach_face()"  et les "face vector"

Nota 3:

Presentation may defer by changing from $x_i$ to $x_{i+1/2}$ 

### Flux avec la moyenne
 
Il faut maintenant trouver la bonne approximation  de $f$. Le plus simple serait de prendre la moyenne des valeurs (LeVeque)
$$
 F_{i}=f(U_{i-1},U_{i})=\dfrac{F(U_{i-1})+F(U_i)}{2}
$$
d'où
$$
 U_i^{n+1}=U_i^{n} -{\Delta t} \dfrac{F(U_{i+1})-F(U_{i-1})}{2 \Delta x}
$$
mais c'est une mauvaise idée...
car la méthode est instable (toute perturbation se trouve amplifiée)  .  

### Flux "Lax-Friedrich"
 
Une technique classique pour stabiliser est "Lax-Friedrich". Mais ensuite nous passerons à une méthode plus efficace (flux dépendants des valeurs propres du sytème, Rusanov et HLL). 
 
 En effet, on peut stabiliser le schéma instable
 $$
 U_i^{n+1}=U_i^{n} -{\Delta t} \dfrac{F(U_{i+1})-F(U_{i-1})}{2 \Delta x}
$$
si on remplace $U_i^{n}$ par une  moyenne $(U_{i+1}^{n} + U_{i-1}^{n})/2$, on obtient le schéma de Lax-Friedrich
$$
 U_i^{n+1}=\frac{1}{2}(U_{i+1}^{n} + U_{i-1}^{n}) - \dfrac{\Delta t}{2 \Delta x}({F(U_{i+1})-F(U_{i-1})})
$$
Cette expression a le mérite d'introduire une diffusion numérique qui tue l'instabilité, 
en effet en repassant au développement de Taylor:
$$\frac{1}{2}(U_{i+1}^{n} + U_{i-1}^{n}) -U_i^n = \dfrac{\Delta x^2}{2} \partial_x^2 U^n+...$$
si on repasse à  la description continue,  
$$
 U_i^{n+1}=U_i^{n} -{\Delta t} \dfrac{F(U_{i+1})-F(U_{i-1})}{2 \Delta x} + \dfrac{\Delta x^2}{2} \partial_x^2 U^n+...
$$
on obtient un terme de diffusion pour l'équation continue associée:
$$
 \partial_t U+\partial_x F(U)=\frac{\Delta x^2}{2 \Delta t}  \partial_x^2 U
$$

Ayant compris cet intéret stabilisateur, l'expression de Lax-Friedrich peut s'obtenir en choisissant le flux numerique de la forme suivante:
$$
F_i=f(U_{i-1},U_{i})=\dfrac{F(U_{i-1})+F(U_i)}{2} 
% - c_\Delta ({F(U_{i+1})-F(U_{i-1})})
  - c_\Delta \frac{({(U_{i})-(U_{i-1})})}{2}
$$  avec  $c_\Delta=\frac{\Delta x}{\Delta t}$.
En effet, substituant cette expression du flux numérique dans le schéma explicite des volumes finis:
$$  U_i^{n+1}=U_i^{n} -{\Delta t} \dfrac{F^n_{i+1}-F^n_{i}}{\Delta x}$$
  on retrouve bien l'expression avec la moyenne de Lax-Friedrich.
On va voir que cette forme est presque la bonne si on prend un  $c_\Delta$ mieux adapté.
Plusieurs flux, peuvent être employés, nous utiliserons les plus simples.

### Flux  upwind
Prenons la forme précédente avec $c_\Delta=dF/dU$
$$
F_i=f(U_{i-1},U_{i})=\dfrac{F(U_{i-1})+F(U_i)}{2} 
  -  c_\Delta \frac{({(U_{i})-(U_{i-1})})}{2}
$$  
nous allons voir ici que c'est le meilleur choix.

Le théorème de Lax dit que pour résoudre un problème d'EDP (comme ici)
 pour lequel on a posé un schéma numérique consistant (qui  retrouve bien l'équations aux dérivées partielles, 
 quand les pas de discrétisation ($\Delta t,\Delta x$, etc.) tendent tous vers 0), 
 la stabilité du schéma est une condition nécessaire et suffisante pour assurer sa convergence.
 
## Code
mandatory declarations:
*/
#include "grid/cartesian1D.h"
#include "run.h"
/** definition of the field h, the flux, its derivative, time step and 
*/
scalar U[];
scalar F[];
double dt;  
double cDelta;
/**

Boundary conditions

*/
U[left] = neumann(0);
U[right] = neumann(0);
/**
Main with definition of parameters
the cFL is 1
*/
int main() {
  L0 = 12.;
  X0 = -L0/4;
  N = 64;    
  DT = (L0/N);
  run();
}
/** 
initial elevation: an exponential "bump"
*/
event init (t = 0) {
  foreach()
    U[] = exp(-x*x);
  }
/** 
print data in console (standard output)

first point is in `X0+1/2*(L0/N))`, 
ith point is in `X0+(i-1/2)*(L0/N))`, 
last point is `X0+(N-1/2)*(L0/N))` which is `X0+L0-1/2*(L0/N))`
*/
event printdata (t += 1; t <=3) {
  foreach()
    fprintf (stdout, "%g %g %g  \n", x, U[], t);
  fprintf (stdout, "\n\n");
}
/** 
integration 
*/
event integration (i++) {
  double dt = DT;
/**
finding the good next time step
*/
  dt = dtnext (dt);
/**

 
Pour fixer les idées, dans le cas de l'équation d'advection simple
$$\partial_t U+\partial_x (U)=0$$
on a $F(U)=U$, la valeur propre est $dF/dU=1$ tout simplement, si on utilise le flux   
$$ 
 F_{i}=  \dfrac{U_{i-1}+U_i}{2} -c_\Delta (\dfrac{U_{i}-U_{i-1}}{2})
$$
avec  $c_\Delta=0$ valeurs moyenne pour le flux, ou $c_\Delta=\dfrac{\Delta x}{ \Delta t}$ cas Lax Wendrof,  ou  si on prend $c_{\Delta}=1$, c'est le flux  upwind:
et le nouvel $U$
$$
 U_i^{n+1}=U_i^{n} -{\Delta t} \dfrac{(F_{i+1}-F_{i})}{ \Delta x}
$$ 



  Expression of the flux  :
$$
F_i=f(U_{i-1},U_{i})=\dfrac{F(U_{i-1})+F(U_i)}{2} 
% - c_\Delta ({F(U_{i+1})-F(U_{i-1})})
  - c_\Delta \frac{({(U_{i})-(U_{i-1})})}{2}
$$ 
avec   
dans le cas Lax
 $c_\Delta=\Delta/\Delta t$
 dans le cas centré 
 $c=0$ 
 dans le cas upwind  $c=1$ 
*/ 
  foreach() {
   //cDelta = Delta/dt;
   // cDelta = 0;
    cDelta = 1;
    F[] = (U[0,0]+U[-1,0])/2.  - cDelta *(U[0,0]-U[-1,0])/2;}
  boundary ({F});
/** 
explicit step
update 
$$
 U_i^{n+1}=U_i^{n} -{\Delta t} \dfrac{F(U_{i+1})-F(U_{i})}{\Delta x}
$$*/
  foreach()
    U[] +=  - dt* ( F[1,0] - F[0,0] )/Delta;
  boundary ({U});
}
/**
## Run
Then compile and run:

~~~bash
 qcc  -g -O2 -DTRASH=1 -Wall  advecte1.c -o advecte1 ;./advecte1 > out
~~~

or better 

~~~bash
 ln -s ../../Makefile Makefile
 make advecte1.tst;make advecte1/plots    
 make advecte1.c.html ; open advecte1.c.html 
~~~

 



## Results
The analytical solution is 
$$U(x,t) =  exp(-(x-t)^2)$$
in gnuplot type

~~~bash
 U(x,t)= exp(-(x-t)*(x-t))
 p'out' u ($1):($2)t'num'w l,'' u 1:(U($1,$3)) t'exact' w l
~~~
which gives $U(x,t)$ plotted here for t=0 1 2 3   and $-3<x<9$  
$\Delta=L0/N=12/64=0.1875$
first point is $-3 + \Delta/2= -2.90625$ second point is in $-3 + \Delta/2+\Delta=-2.71875$ next is in -2.53125
up to the previous last one $-3 +1/2 \Delta +(N-2)\Delta = 8.71875$  and finally the 64th point ($N$) is  $-3 +1/2 \Delta + (N-1) \Delta =8.90625.$

(`x1=X0+(1-1/2) L0/N, x2=X0 +(2-1/2) L0/N, ...,xi= X0 +(i-1/2) L0/N, .. xN=X0 +(N-1/2) L0/N= X0 +L0-1/2 L0/N`)  

  




~~~gnuplot
 set output 'dessin.svg'; 
 reset
 set xlabel "x"
 U(x,t)= exp(-(x-t)*(x-t))
 p'out' u ($1):($2)t'num.'w p,'' u 1:(U($1,$3)) t'exact' w l
 
~~~


# Exercise

Change $c_\Delta$, 

if $c_\Delta = 0$ check that the amplitude increases (unstable), 

if $c_\Delta =\Delta/dt$ check that the amplitude decreases (stable), but decreases to much
    
if $c_\Delta =1$ check that the amplitude decreases (stable), but decreases less, increase number of points
  
Change `DT` to test stability

# Links
  
* [advecte1.c]() explains the notions of advection, testing the flux, coded with Basilisk 

* [advecte1c.c]() the same than the previous one but with standard C  

* [advecte.c]() advection with Basilisk

# Bibliography 
 
* [LeVeque](https://www.cambridge.org/core/books/finite-volume-methods-for-hyperbolic-problems/finite-volume-methods/CB7B0A27A6D37AE3B906D4AE7C60A87E) 
 chapitre 4, Finite Volume methods
* [Toro](https://www.springer.com/gp/book/9783540252023) "Riemann Solvers and Numerical Methods"  Chap 5 Springer
* [PYL](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/code_C_saintvenant.pdf)   cours M1, 
"Résolution numérique des  équations de Saint-Venant,
mise en oeuvre en volumes finis par un solveur de Riemann bien balancé"


Version 1 Feb 2015, 30 June 2018, ready for new site 09/05/19
*/


