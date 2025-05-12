/**
##Problem

In this section, we solve the collapse of a heap over a slope using a Saint-Venant solver. Previously, we addressed this problem using the flood-wave approximation and a simplified code that neglected inertia effects (see [Chezy](http://basilisk.fr/sandbox/mjafari/TurbulentC/Chezy.c)). This approach assumed the flow was dominated by gravity and friction, omitting inertial and pressure variation terms for simplicity. 


###Governing Equation
In the current script, we solve the classical (1D) Saint-Venant equations while incorporating a turbulent basal friction term. This friction is modeled using empirical formulations, such as those proposed by Chezy or Manning, to account for the interaction between the flow and the bed roughness. The governing equations are:

$$
 \left\{\begin{array}{l}
         \partial_t h+\partial_x (hu)=0\\
  \partial_t hu+ \partial_x hu^2  
= -g h \partial_x \eta -C_f u^2 ,\;\;\;\eta=h + z_b  
        \end{array}\right. 
$$

Here, $h$ represents the flow depth, $u$ is the flow velocity, $z_b$ is the bed elevation, and $C_f$ denotes the friction coefficient.

The basal friction term $\tau_0$ (last term in the momuntum equation), proportional to the square of the flow velocity ($u^2$), is incorporated into the acceleration calculation manually. This is done outside the Saint-Venant solver using a simple implicit scheme to capture the effect of turbulent friction on the flow dynamics.

###Codes
*/

#include "grid/multigrid1D.h"
#include "saint-venant.h"

/** Definition of the Parameters
 */

double DTprt = 5., tmax = 50.;
double S, Cf;
char s[80];

u.n[left]  = neumann(0);
u.n[right] = neumann(0.); 

/**
 Main with definition of parameters, note that $S$ is the normalized slope and $C$ the friction coefficient.
 */
int main() {
  X0 = -1.;
  L0 = 20;
  G = 1;
  N = 2040;
  //nl = 2;
  S = 1;
  Cf = 0.5;
  run();
}

/**
 Initial elevation: a constant rectangle that we release over the slope.
 */
event init (i = 0)
{

  foreach() {

    zb[] = -(x-X0)*S;
    h[] = (x>=0 && x<=.5)*.5 + dry;
    u.x[] = 0;
  }

}

/**
Adding acceleration due to the frition using the splitting methode:

$$ \frac{u^{n+1}-u^n}{dt} = -C_f \lvert u \lvert \frac{u^{n+1}}{h} $$
*/

event friction(i++){
  
  foreach() {
    double a = h[] < dry ? HUGE : 1. + Cf*dt*norm(u)/h[];
    foreach_dimension()
      u.x[] /= a;
  }
}    
 
/**
Saving data;
*/
event printdata (t += DTprt; t <= tmax){
    sprintf (s, "shape.txt");
    FILE * fp = fopen (s, "a");
    foreach(){
      fprintf (fp, "%g %g %g %g %g %g %g\n", x, h[], t, u.x[], eta[],
               eta[]-zb[], zb[]);
    }
    fprintf (fp, "\n");
    fclose(fp);
    }

/**
####Self-Similarity
By studying the governing equations, we find stationnary solutions in the self-similar domain which are expressed as:

$$  \hat{h} \sim \hat{t}^{-2/3} \mathcal{H}(\eta), \quad \mathrm{where} \quad \eta = V_0 \hat{x}\hat{t}^{-2/3}, \quad \mathcal{H}(\eta) = \left( \frac{2 \eta}{3 \beta V_0} \right)^2.$$
*/

/**
#plot
     

~~~gnuplot  

 set xlabel "x " 
 set ylabel "h " 
 
 p  'shape.txt'  w l
~~~


~~~gnuplot  

 set xlabel "x " 
 set ylabel "h, \hat{h}" 
 
 S=1.
 C=0.5
 V0=0.25
 n = 3./2.
 beta = (S/C)**(1./2.)   
 H(x)=  (x/(beta*n*V0))**(1/(n-1))
 
 p[0:][:0.15]  'shape.txt'  w l, '' u 1:(H(V0*$1/($3**(1/n))))*($3**(-1/n)) w l
~~~


~~~gnuplot  

 set xlabel "xt^{-2/3}" 
 set ylabel "ht^{2/3}" 
 tt=20
 S=1.
 C=0.5
 V0=0.25
 n = 3./2.
 beta = (S/C)**(1./2.)   
 etaf = V0*((beta*n)**(1/n))*((V0*n/(n-1))**((n-1)/n))
 H(x)=  (x/(beta*n*V0))**(1/(n-1))
 Hmax=H(etaf)
 set arrow from etaf, graph 0 to etaf, graph 1 nohead lc rgb "red"
 p[-0.1:.5][0:.6] 'shape.txt' u (V0*$1/($3**(1/n))):($3>tt?($2*($3**(1/n))):NaN) t'h' w l,\
  H(x) t ' ss' w l,Hmax
~~~

*/
