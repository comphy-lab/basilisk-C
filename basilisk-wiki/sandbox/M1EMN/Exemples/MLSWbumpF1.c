/** 
#  Interactive flow  


Supercritical and subcritical flow over a bump in a thin layer flow.



## Code

We use the Multilayer solver
*/
#include "grid/cartesian1D.h"
#include "saint-venant.h" 

double tmax,Fr,hout;
double alpha;
scalar tau[],delta1[],delta2[];
FILE * fq;
 
/**
position of domain `X0=0`, length `L0`, no dimension `G=1` 
run with 64 points. We define `hout` from the analytical solution, 
it will be usefull for the boundary condition.
Note that the velocity is $Fr$ the Froude number. 
The Reynolds is $Re=500$ ($\nu=0.001$ and $Fr=0.5$)
*/

int main() {
  
  double web=1;
  X0 = 0;
  L0 = 1.-X0;
  G = 1;
  N = 128*2/web;
  tmax=5; 
  nu=0.0025;
  nl = 200;

#if 1 
  Fr=0.5;
  alpha=0.09/web;   
  hout = 1+1.73*sqrt(L0*nu/Fr)*(Fr*Fr)/(Fr*Fr-1);
  run();
  system("cp   xproff.txt xproffSUB.txt");
#endif
  
#if 0  
  Fr=1.5;
  alpha=0.04/web;
  hout = 1+1.73*sqrt(L0*nu/Fr)*(Fr*Fr)/(Fr*Fr-1);
  run();
  system("cp   xproff.txt xproffSUP.txt");
#endif

}
/** 
Initialisation,  we initialize *h* and *hc*.  
We define a scalar field *hc* to check for convergence on $h$. */

scalar hc[];

event init (i = 0)
{
  foreach(){
    zb[] = 0  ;
    delta1[]=1.7*sqrt(x*nu/Fr);
    h[] = hc[] = 1+(delta1[]+zb[])*(Fr*Fr)/(Fr*Fr-1); 
   }
   boundary({zb,h,hc});

 for (vector u in ul) {
  foreach()
  u.x[] =Fr*(1+(delta1[]+zb[])/(1-Fr*Fr));

 }

/**
Boundary condtion, neumann for $h$ and impose  flux  and a Rieman invariant 
$u - 2 * \sqrt{gh} $ is constant  around  $h=1$ 
neumann for $u$ and impose  the height from the analytical linearized solution for $Fr<1$.
But for $Fr>1$, given values at the entrance and neumann at the exit.
*/
  
  h[left]   = (Fr<1? neumann(0) : dirichlet(1));
  eta[left] = (Fr<1? neumann(0) : dirichlet(1));
 
  for (vector u in ul) {
    u.n[left] = dirichlet(h[left] ? Fr /h[left] : 0.)  - (2*sqrt(G*h[]) -  2*sqrt(G*1.))*(Fr<1);  
    u.n[right] = neumann(0.);
  }
 eta[right] = (Fr<1?  dirichlet(hout) : neumann(0)); 
}
/**
 Compute friction:
*/ 
event friction (i++) {
/**
smooth update of the bump
*/

   foreach()
     zb[] = 0 + alpha*exp(-800*sq(x-L0/2))*(1-exp(-t/10.));
/**
 Poiseuille friction: implicit scheme (time splitting, here commented)
 $$\frac{\partial u}{\partial t} = - C_f \frac{u}{h^2}$$
 hence, if implicited for $u$, explicited for $h$:
 $$ u^{n+1} = u^n -C_f \Delta t \frac{u^{n+1}}{(h^{n})^2}$$
*/ 
  
 if(nl==1){ 
  foreach()
  {    double ff = h[] < dry ? HUGE : (1. +   3*nu*dt/h[]/h[]);
    u.x[] /= ff;
    }

   foreach()
    tau[]=3*u.x[]/(h[]); 
  }
/**
Multilayer case, solve full friction across tthe layers
 $$\frac{\partial u}{\partial t}+...  = ...+ \frac{\partial^2 u}{\partial z^2}$$

Tne value at the wall of the velocity is the friction $\frac{\partial u}{\partial z}$,  
this derivative is twice the value in the first layer, 
*/   
else{
  vector u0 = ul[0];
  foreach()
    tau[]=2*u0.x[]/(h[]/nl);
   }
}

/**
We check for convergence. */
event logfile (t += 0.1; i <= 75000) {
	
/*
	FILE *f; 
	f = fopen("F1.IN", "r"); 
fscanf(f,"alpha=%lf     \n",&alpha); 
fclose(f); */

	//alpha = .11*(1-exp(-t/10));
	//foreach()     zb[] = 0 + alpha*exp(-800*sq(x-L0/2));

  double dh = change (h, hc);
  fprintf (stderr,"coucou t= %g  -- var=  %g, Re = %lf alpha = %lf\n", t, dh,  Fr/nu, alpha);
  if (i > 0 && dh < 10e-5)
    return 1;

/** }
And save fields

event sauv (t<=tmax;t+=.1) { */


 foreach() { 
      double z = zb[],psi=0;
      vector ue;
      ue=ul[nl-1];
      int l = 0;
      z -= (layer[0]/2)*h[];
      for (vector uz in ul) {
  z += (layer[l++])*h[];
  psi += (1-uz.x[]/ue.x[])*(layer[l])*h[];
      }
  delta1[]=  psi; 
}

 vector ue;
 ue=ul[nl-1];
 
 FILE *fq =  fopen("xproff.txt","w");  
 foreach() { 
      fprintf(fq,"%g %g %g %g %g %g %g %g %g\n", x, h[], tau[], ue.x[],t,Fr,nu,delta1[],zb[]);}
  fprintf(fq," \n");

  fflush(fq);
  fclose(fq);
}


/**
 
## Run
To compile and run  :

~~~bash
   qcc -O2  -o MLSWbumpF1  MLSWbumpF1.c  -lm;  ./MLSWbumpF1
~~~
 
using `Makefile`

~~~bash
  make MLSWbumpF1.tst
  make MLSWbumpF1/plots
  make MLSWbumpF1.c.html 
~~~  

## Plots


### Supercritical flow


Note upstream influence (perturbation before the bump) and boundary layer separation before the bump.

Shape of the bump, friction divided by gauge of Blasius friction $\tau /Fr/(.32 \sqrt{Re})$,
reduced Blasius friction $1/\sqrt{x}$

~~~gnuplot
 set xlabel "x"
 p[0:1][:4]'xproffSUP.txt'u 1:($9/.04) t'bump' w lp ,\
 '' u 1 ($3/sqrt(600)/.32/1.5) t'tau comp' w lp,1/sqrt(x) t 'tau Blasius'
~~~

Shape of the bump, displacement thickness divided by gauge of Blasius  displacement 
 $\delta_1/(1.7 \sqrt{Re})$,
reduced Blasius displacement thickness $\sqrt{x}$.

~~~gnuplot
 set xlabel "x"
 p[:][:2]'xproffSUP.txt'u 1:($9/.04) t'bump' w lp ,'' u 1:($8*sqrt(600))/1.7 t'delta_1',sqrt(x) t'delta_1 Blasius',''u 1:4 t'ue' w l,'' u 1:(1.5*(1-$8/(1.5*1.5-1))) t 'lin' w l
~~~


### Subcritical flow


Note no upstream influence (perturbation appear at the bump foot), boundary layer separation is after the bump.

Shape of the bump, friction divided by gauge of Blasius friction $\tau /Fr/(.32 \sqrt{Re})$,
reduced Blasius friction $1/\sqrt{x}$

~~~gnuplot
 set xlabel "x"
 p[:][:5]'xproffSUB.txt'u 1:($9/.1) t'bump' w lp ,\
 '' u 1:($3/sqrt(200)/.32/.5) t'tau comp' w lp,1/sqrt(x) t 'tau Blasius'
 #p[:][:4]'xproffSUB.txt'u 1:($9/.1)  t'bump' w lp ,'' u 1:($8*sqrt(100))/1.7 t'delta_1',sqrt(x) t'delta_1 Blasius',''u 1:4 t'ue' w l,'' u 1:(.5*(1-$8/(.5*.5-1))) t 'lin' w l
~~~

Shape of the bump, displacement thickness divided by gauge of Blasius  displacement 
 $\delta_1/(1.7 \sqrt{Re})$,
reduced Blasius displacement thickness $\sqrt{x}$.

~~~gnuplot
 set xlabel "x"
 p[:][:2]'xproffSUB.txt'u 1:($9/.1) t'bump' w lp ,'' u 1:($8*sqrt(200))/1.7 t'delta_1',\
 sqrt(x) t'delta_1 Blasius',''u 1:4 t'ue' w l,'' u 1:(.5*(1-$8/(.5*.5-1))) t 'lin' w l
~~~







## Bibliography


*  Lagree IBL CISM
 
* [Lagree P-Y](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/MFEnv.pdf)
"Boundary Layer" Cours MSF12, M2 SU

* [James Lagree Le Legrand ]()



 
*/