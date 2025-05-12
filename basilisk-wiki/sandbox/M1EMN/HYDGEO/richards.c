/**
# Resolution of Richard infiltration problem  (with variable diffusion coefficient!!!)

## Problem:

Finding the position of the front of water being imbibited in a porous media, such as coffee in sugar (*"canard"*, but with gravity up) or water in the porous soil during rain ( gravity down). Here we consider an horizontal infiltration.

## Equations

 $\theta$ volume fraction saturation, volumetric water content (or moisture content :the quantity of water contained in the material) (*teneur en eau volumétrique*)

$$ \frac{\partial \theta}{\partial t } =
\frac{\partial }{\partial x } (K \frac{\partial (h+z) } {\partial x } )$$

 $h$ (charge hyrdaulique)
 $K(h)$ hydraulic conductivity (*perméabilité*)
 $h=p/(\rho g )$ capillary potential

 
 $H=p/(\rho g ) +z$ hydraulic head (*charge hydraulique*),
 
 $$ C(h) \frac{\partial h}{\partial t } =
 \frac{\partial }{\partial x } (K \frac{\partial h } {\partial x } -0)$$
 with $C(h)$ specific moisture capacity
 and $K(h)$ hydraulic conductivity,
 
 
 Saturation-based
 $$  \frac{\partial \theta }{\partial t } =
 \frac{\partial }{\partial x } (K \frac{\partial H } {\partial x } )\text{ is written }
 \frac{\partial \theta }{\partial t } =
 \frac{\partial }{\partial x } (K\frac{\partial H } {\partial \theta } \frac{\partial \theta} {\partial x } )$$
 
 
 $$  \frac{\partial \theta }{\partial t } =
 \frac{\partial }{\partial x } (D \frac{\partial \theta } {\partial x } )$$
 
 
 $D(\eta)$ is the hydraulic diffusivity, we take the following $D= {(0.05/(1.05 - \theta))}^2$ proposed by Caputo.
 
 
The surface is in $x=0$, we impose $\theta=\theta_S$ (we put 1, but it should be .4), we impose a Neumann at depth $x=L_0$.
This is the imbibition problem.


it is   written with standard C
 
## Code
mandatory declarations in C:
*/
#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include <string.h>


/** definition of the field U, the flux F, time step 
*/
double*x=NULL, *theta=NULL;
double*Fx=NULL, *D=NULL;
double dt;  

double L0,X0,Delta; 
double t;
int i,it=0;
int N;

/**
Main with definition of parameters
*/
int main() {
  L0 = 15.;
  X0 = 0;
  N = 128*2;
  t=0;
  Delta = L0/N;
  dt = 0.0001;
  
/**
 dynamic allocation 
*/
  x = (double*)malloc((N+2)*sizeof(double));
  theta = (double*)malloc((N+2)*sizeof(double));
  D = (double*)malloc((N+2)*sizeof(double));
  Fx= (double*)malloc((N+2)*sizeof(double));
  
    for(i=0;i<=N+1;i++)
/** 
~~~bash
  |   |   |   |   |
0 | 1 | 2 | 3 | N | N+1
  --|---|---|---|---
 0  x1  x2      L
  <--------------->
      <--> 
N cells, +2 ghost cells, (N+2) data

   |           |
   |->Fi       |->Fi+1
   |  thetai   |
xi+Delta/2
       xi   xi+Delta/2
~~~

Fi is the flux across interface between hi-1 and hi  (sometimes Fi-1/2)

Fi+1 is the flux across interface between hi and hi+1  (sometimes Fi+1/2)




first cell, index 1,  between  `X0` and `X0+Delta`, centred in  `Delta/2`
ith cell beween `(i-1) Delta` (left) and `i Delta`(right) centered in `(i-1/2)Delta`  

`Delta=L0/N`

*/
    {  x[i]=X0+(i+1-1./2)*Delta;  
/**
initial   
*/
       theta[i] = 0;
  }
/**
 begin the time loop
*/  
   while(t<=9){
     t = t + dt;
     it++;  
/** 
print data
*/
//event printdata (t += 1; t <=3)  
if(   (it == (int)(1/dt))
   || (it == (int)(3/dt))
   || (it == (int)(5/dt))
   || (it == (int)(7/dt))
   || (it == (int)(9/dt)) ){
 for(i=1;i<=N;i++)
    fprintf (stdout, "%g %g %g\n", x[i], theta[i], t);
  fprintf (stdout, "\n\n");
}
/** 
Richards equation
$$C \frac{\partial h}{\partial t } =
\frac{\partial }{\partial x } (K \frac{\partial  } {\partial x }h -0)$$
is a conservative equation 
$$  \frac{\partial \theta}{\partial t } = -
\frac{\partial }{\partial x } F$$
with the flux
$$F = - (D \frac{\partial  } {\partial x }\theta -0)$$
*/
       for(i=1;i<=N+1;i++)
       { double thetam;
           thetam =(theta[i]+theta[i-1])/2;
           D[i] =  sq(0.05/(1.05 - thetam));
       }
       /**
*/    
 for(i=1;i<=N+1;i++)
   {
    Fx[i] = - D[i]*((theta[i]-theta[i-1])/Delta -0.00);
    }
/** 
explicit step
update 
$$
 theta_i^{n+1}=theta_i^{n} -{\Delta t} \dfrac{F(U_{i+1})-F(U_{i})}{\Delta x}
$$*/ 
 for(i=1;i<N;i++)
    theta[i] +=  - dt* ( Fx[i+1] - Fx[i] )/Delta;
/**
 Boundary condition, imposed value at the surface
*/
    
   theta[0] = - theta[1] + 2*(1);
   theta[N] = theta[N-1];
 }
  
/**
 final purge
*/
  free(theta);
  free(Fx);
  free(D);
  free(x);
}
/**
## Run
Then compile and run:

~~~bash
 cc  -g -O2 -DTRASH=1 -Wall  richards.c -o richards ;./richards > out
~~~

or better 

~~~bash
 make richards.tst;make richards/plots    
 make richards.c.html ;
 
 
 source ../Exemples/c2html.sh richards
 open richards.c.html
~~~


## Results
One analytical solution is 
$$U(x,t) =  erfc(\frac{x}{2\sqrt{t}})$$
in gnuplot type

~~~bash
 h(x,t)= erfc(x/\sqrt(4 t))
 p'out' u ($1):($2)t'num'w l
~~~
which gives $h(x,t)$ plotted here for t=0 1 2 3  and $0<x<5$

We see the   the selfsimilar solution $erfc$ corresponding to constant coefficient in green

~~~gnuplot
 set xlabel "x"
 h(x,t)= erfc(x/2./sqrt(t))
 p[:]'out' u ($1):($2)t'num'w lp,'' u ($1):(h($1,$3)) t'erf' w l
~~~

 
 Self similar solution
 $\eta=x/\sqrt{t}$ and
$$-\frac{\eta}{2} \frac{d \theta}{d \eta} = \frac{d  }{d \eta}\left( D\frac{d \theta}{d \eta}\right)$$
 
 
 Solution is indeed self similar, note that the front is defined.
~~~gnuplot
 set xlabel "x/sqrt(t)"
 p[:1]'out' u ($1/sqrt($3)):($2)t'num'w l
~~~

## Exercice
Compare with heat equation


## bibliography

* Caputo et al. 
["Experimental and numerical study of a transient, twodimensional unsaturatedsaturated water table recharge problem"](http://lmi2.insa-rouen.fr/~jgcaputo/cs07.pdf)
 
* Hogg et al. [Numerical Solution of Richards Equation  Review of Advances and Challenges](https://www.researchgate.net/publication/320515457_Numerical_Solution_of_Richards%27_Equation_A_Review_of_Advances_and_Challenges)


* [PYL erf solution](http://www.lmm.jussieu.fr/~lagree/COURS/ENSTA/PC1.ENSTA.pdf)

 * [Lagrée P-Y](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/MFEnv.pdf)
 "Equations de Saint Venant et application, Ecoulements en milieux naturels" Cours MSF12, M1 UPMC


update Minorque 03/20
*/
