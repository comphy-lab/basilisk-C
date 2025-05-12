/**
# 1D Shallow Water Saint Venant Von K&aacute;rm&aacute;n

 see
["svvk"](http://basilisk.fr/sandbox/M1EMN/Exemples/svvk.c) for the explanation of the equations

### SVVK system
The final SVVK (Saint-Venant Von K&aacute;rm&aacute;n) equations are
$$
 \left\{\begin{array}{l}
         \partial_t h+\partial_x hu_e= \partial_x q_1\\
	\partial_t (hu_e)+ \partial_x \left[hu_e^2+ \dfrac{h^2}{2}\right]
= - gh \partial_x z_b + u_e \partial_x q_1 - \tau_\pi\\
\partial_t q_1 + \partial_x (u_e q_1 ( 1+\frac{1}H )) = u_e \partial_x q_1  -\tau_\pi  + \tau_0 \\
        \end{array}\right. 
$$
at time  t=0, $h=1 $, $u_e=u_0$, $q_1=q_{10}$ and a given bottom $z_b$  

note that $H$, $\tau_\pi$ and $\tau_0$ need a closure

note that the source terms are topography $\partial_x z_b$ and equivalent topography $\partial_x (\delta_1 u_e)$.

note that there is an exchange of momentum $u_e \partial_x q_1$ and friction $-\tau_\pi$ between the two layers.

# Code
## Declarations

Before including the conservation solver, we need to overload the
default *update* function of the predictor-corrector scheme in order
to add our source term. 
*/
#include "grid/cartesian1D.h"
#include "predictor-corrector.h"
static double momentum_source (scalar * current, scalar * updates, double dtmax);
event defaults (i = 0)
  update = momentum_source;
#include "conservation.h"
/**
## Declarations of Variables
We define the conserved scalar fields $h$ $Q$ and $q_1$ which are passed to the
generic solver through the *scalars* list. We do not have any conserved
vector field, we should do this for 2D. */
scalar h[], Q[], q1[], zb[], dxq1[],tauPi[],tau0[],d1[],coef[];
scalar * scalars = {h,Q,q1};
vector * vectors = NULL;
/**
The other parameters are specific to the example. */
double sec =1e-4,tmax;
double theta,alpha,h0,Q0;
double H=2.59,f2=1./2.59;//.23
/**
## Functions

We define the *flux* function required by the [generic
solver](/src/conservation.h). */
void flux (const double * s, double * f, double e[3])
{  
  double h = s[0], Q = s[1], q1 = s[2],  u ;
/**
The solution is done by spliting, first :  
  $$   \dfrac{h^*-h}{\Delta t} +\partial_x hu_e = 0$$
  $$ \dfrac{(h^*u_e^*) -(hu_e)}{\Delta t}+ \partial_x \left[hu_e^2+ \dfrac{h^2}{2}\right]
= 0 $$
$$ \dfrac{ q_1^* -q_1}{\Delta t} + \partial_x F =0$$
$$F = (u_e q_1 ( 1+\frac{1}H )), $$ 
are solved by `"conservation.h"`  
We have to furnish the expression of the flux
*/   
  u = (h > sec ? Q/h: 0 );
  f[0] = Q;
  f[1] = u*u*h + h*h/2.;
  f[2] = (1+1./H)*u*q1; 
/**
and the explicit eigenvalues
*/  
  // min/max eigenvalues
  double c = sqrt(h);
  e[0] = fmin(fmin(u - c,u + c),(1+1./H)*u); // min
  e[1] = fmax(fmax(u - c,u + c),(1+1./H)*u); // max
}

/**
We need to add the source term of the momentum equation. We define a
function which, given the current states, fills the *updates* with the
source terms for each conserved quantity. */

static double momentum_source (scalar * current, scalar * updates, double dtmax)
{
  /**
  We first compute the updates from the system of conservation laws. */
  double dt = update_conservation (current, updates, dtmax);

  /**
  We recover the current fields and their variations from the lists... */

  scalar h = current[0], dh = updates[0];
  scalar Q = current[1], dQ = updates[1], q1 = current[2], dq1 = updates[2];

  /**
 
  
  Next the viscous sources for Saint-Venant and Von-K&aacute;rm&aacute;n equations allow to update the estimated $h^*,u_e^*,q_1^*$ to the finals $h,u_e,q_1$:
  $$ \dfrac{h-h^*}{\Delta t} = \partial_x q_1$$
  $$\dfrac{(hu_e)-(h^*u_e^*)}{\Delta t}  = 
  - gh^* \partial_x z_b + 
  u_e^* \partial_x q_1/h^*  - \tau_\pi /h^*. $$
 $$\dfrac{ q_1 -q_1^*}{\Delta t} =  u_e \partial_x q_1  -\tau_\pi  $$	
the update of $q_1$ is then without basal friction (to avoid the division by a 0 ),  
   We add the source term for all ths variables. 
 */
  foreach(){
    dh[] += dxq1[];	
    dQ[] += -h[] * (zb[]-zb[-1])/Delta;
    dQ[] += (h[]> sec ? dxq1[]*Q[]/h[] : 0) -  tauPi[];
    dq1[]+= (h[]> sec ? dxq1[]*Q[]/h[] : 0) -  tauPi[];
  }
  return dt;
}
/**
## Subroutines
The viscous correction is computed here
*/
event boundary_layer(i++){
/** 
Upstream differentiation to evaluate $\partial_x q_1$ */
  foreach()
       dxq1[]= (Q[]>0 ? (q1[1] - q1[0])/Delta : (q1[0] - q1[-1])/Delta ) ; 
  boundary ({dxq1});
/** 
Definition of the friction functions
 $\tau_\pi$ and $\tau_0$, as part of the closure, so they depend on the other variables
*/  
 foreach(){
     coef[] = (h[]>sec && Q[] >0 ? pow(sin(pi/2*min((q1[]/(Q[]/h[])/(h[]/3)),1)),16) : 0);    
     tauPi[] = (h[]> sec ?  3*Q[]/h[]/h[] * f2*H * coef[]:0) ;
     tau0[] =  (h[]> sec && fabs(q1[]) >0 ? (f2*H)*sq(Q[]/h[])/q1[] :0) ; 
 }
  boundary ({tauPi,tau0});  
/**
and the final split $\partial_t q_1 =  \tau_0$ with $\tau_0=  (f_2H) u_e^2/q_1$ 
 final evaluation of $q_1$ with  friction $\tau_0$, so that 
$$\partial_t (q_1^2)/2 =   (f_2 H)u_e^2$$
hence we integrate explicitely to avoid the division by a 0 in the friction in $1/q_1$, hence 
$$q_1^2(t+\Delta t)  =  q_1^2(t )  +  2 (f_2H)u_e^2 \Delta t$$
*/   
  foreach(){
   q1[] = (h[]>sec ?   sqrt(q1[]*q1[] + 2* dt * (f2*H)*sq(Q[]/h[])): 0 );
   q1[] = (Q[]>0 ? q1[] : -q1[]);
   } 
  boundary ({q1});
}
/**
## Boundary conditions
We impose   $Q(t)$ at the left of the domain.  
A given flat profile with thin boundary layer at the left
*/
Q[left] = dirichlet(Q0);
h[left] = dirichlet(h0);//neumann(0);//
q1[left] = 0;
/**
right, output
*/
Q[right] = neumann(0);;
h[right] = neumann(0);
q1[right] = neumann(0);
/**
## The Main
The domains starts from `X0`, length `L0`, without dimension   the velocity is measured by $\sqrt{gh_0}$),
were  $\alpha$ is the bump height, $-\theta$ the slope angle. 
*/
int main() {
  X0=0;
  L0=1./2;        //1.
  N = 1024*2;  //1024*4
  tmax = 5.;   // 5 
  theta = 0.5 ;
  Q0 = 0.5;
  h0 = 1;
  alpha = 0.2;
  run();
}
/**
## Initial conditions 
The initial conditions are $h=h0$ and $Q=Q0$. 
and  bumps on a mild slope
*/
event init (i = 0) {
  double xi,xii,xiii,x0=0.1;	
  foreach(){
    xi   = (x - 0.05)/.001;
    xii  = (x - 0.4*5)/.1;
    xiii = (x - 0.35)/.01;	
    zb[] =  - theta *(x - x0)*(x>x0) + alpha * exp(- xi*xi)    + 2*alpha * exp(- xii*xii)   +  alpha * exp(- xiii*xiii)  ;}
	
  foreach(){
   h[] = h0;
   q1[] = 1*(x<.3)*(x>.1)*0;
   Q[] = Q0;
  // zb[] = (0.2*exp(-10000*sq(x-.5)) + 0.*exp(-10*sq(x-.75)) - theta*x) ;
  }
} 
/**
and this is done for the computations.

## Outputs
*/
event field_short (t<.1;t+=.001) {	 	
  foreach()
    fprintf (stderr, "%g %g %.6f %.6f %.6f %.6f %.6f \n", t, x, h[], Q[], q1[], zb[],coef[]);
  fprintf (stderr, "\n\n");
}

//log
event  field_long (t += 0.5; t < tmax) {
  foreach()
    fprintf (stderr, "%g %g %.6f %.6f %.6f %.6f %.6f \n", t, x, h[], Q[], q1[], zb[],coef[]);
  fprintf (stderr, "\n\n");
}
/**
gnuplot output for SVVK solution at long time
 */
 /**
End
*/
event end (t = tmax) {
    double fun; 
    static int i=0;	
 
   if(i==0){
   	i++;
    foreach(){
    	fun = 3*f2*H*(Q[]/h[])/h[]/h[]/theta;
    	fun =   3*q1[]/(Q[]/h[])/h[];
    	fun = coef[];
      fprintf (stdout," %g %g %g %g %g %g %g %g %g \n", x, h[], (Q[]/h[]),  q1[], zb[], tauPi[],tau0[],(f2*H)*sq(Q[]/h[])/q1[] , fun );}
  }   
}
/**
fin des procédures.

# Run and Results
## Run
To compile and run (save in `out` and `log`) :

~~~bash
qcc -g -O2 -Wall -o svvk_3 svvk_3.c -lm ; ./svvk_3 > out 2> log
~~~
  
 


## Plot of main fields
Plot of shear $\tau_0$ and $\tau_{\pi}$, free surface (blue) and bottom (black). In red the $\tau_\pi$, and in green the $\tau_0$, they differ when the flow is dominated by ideal fluid (at the entrance, and accelerated on the bumps). They merge at the end when the flow is a Nusselt: when the displacement thickness $\delta_1= q_1/u_e$  (in cyan) is equal to $h/3$.


~~~gnuplot  result, free surface (blue), displacement thickness (cyan) and bottom (black)
set xlabel 'x'
p[:][:2]'out'u 1:6 w l t'tau_{Pi}' ,'' u 1:7 t'tau_0'w l,''u 1:(($5+$1*.5))t 'z_b +\theta*x'w l linec -1,'' u 1:2 t'h' w l linec 3,'' u 1:($4/$3) t'delta_1' w l
~~~
 

## Small time/ space
At small time, one recovers the Blasius solution and Rayleigh

at small $x$, we have the Blasius problem $\partial_x (u_e^2 \dfrac{\delta_1}{H})  = f_2 H u_e / \delta_1$, hence
$q_1=\delta_1 u_e = \sqrt{2 f_2 H^2 u_e x}$
(comparison is not so good here, we need more points?).

~~~gnuplot  result small $x$
p[:.1]'out' u 1:4 t'q_1'w l,''u 1:(sqrt(2* 2.59*$3*$1)) t'q_1 est'w lp
~~~ 



at small $t$, we have the Rayleigh problem 
$\partial_t (u_e \delta_1 )  = f_2 H u_e / \delta_1$, hence
$q_1= \delta_1 u_e = \sqrt{2 f_2 H u_e^2 t}$
Il y a un problème de $sqrt(2)$???

~~~gnuplot q_1(x,t)
set xlabel "x" 
set xlabel "t"
set zlabel "q_1"  
sp[:.01][0:.05]'log'u 2:1:5 w l,sqrt(2*2.59*x*.5)*(x<y/2.55? 1 :NaN) ,sqrt(2*y*.5*.5)*(x>y/2.55? 1 :NaN) 
~~~

  


## Large time/ space
far from the bump, the solution should be invariant in $x$ so that $\partial_x=0$, hence 

$u_e(h- \delta_1)$ is constant , 
$0 = -   (-\theta h)  - \tau_\pi$ and $0= \tau_{\pi}- \tau_0$

with $\tau_{\pi}= 3 u_e /h$ and $\tau_0 f_2 H  u_e/\delta_1$, hence $\delta_1=h/3$ and 
$u_e=  (h^2/3)\theta$


we have a region were the Nusselt (half Poiseuille) solution is valid.




~~~gnuplot  result far from the entrance: the Nusselt film is recovered $q_1=h^3\theta/9$
set xlabel "x"   
   p'out' u 1:($2*$2*$2*.5/3/3) t 'theta h^3/9',''u 1:4 t'q_1'
~~~ 

~~~gnuplot  result far from the entrance: the Nusselt film is recovered: $\tau_{Pi}= \theta h = \tau_0$...
set xlabel "x"   
  p[.1:]'out' u 1:($6) t'tauPi' w l,''u 1:7 t 'tau0' w l,''u 1:($2*.5) t'theta h'w l,'' u 1:(3*$3/$2) w l,''u 1:($3*$3/$4)
~~~ 

 

## Conclusion
seems to work, need extra work to check and improve the closures.
Check signe for reverse flow...



V.1 Noeux les Mines, 04 Juillet 15

#Bibliography


* see
["svvk"](http://basilisk.fr/sandbox/M1EMN/Exemples/svvk.c) 

* Schlichting, 1968 (with Gersten 00) ["Boundary Layer theory"](https://books.google.fr/books?id=8YugVtom1y4C&printsec=frontcover&dq=schlichting+boundary+layer+theory&hl=fr&sa=X&ei=gwOZVe_JPMmtU83jgIAE&ved=0CCEQ6AEwAA#v=onepage&q=schlichting%20boundary%20layer%20theory&f=false)

* P.-Y. Lagrée:
["A rem(a)inder on Ideal Fluid - Boundary Layer decomposition"](http://www.lmm.jussieu.fr/~lagree/COURS/CISM/blasius_CISM.pdf)
UPMC Master curse


* P.-Y. Lagrée & S. Lorthois (2005):
["The RNS/Prandtl equations and their link with other asymptotic descriptions. Application to the computation of the maximum value of the Wall Shear Stress in a pipe"](http://www.lmm.jussieu.fr/~lagree/TEXTES/PDF/lagreelorthois05.pdf)
Int. J. Eng Sci., Vol 43/3-4 pp 352-378.
 
* P.-Y. Lagrée, E. Berger, M. Deverge, C. Vilain & A. Hirschberg (2005):
["Characterization of the pressure drop in a 2D symmetrical pipe: some asymptotical, numerical and experimental comparisons"](http://www.lmm.jussieu.fr/~lagree/TEXTES/GLOTTE/1190.pdf),
ZAMM: Z. Angew. Math. Mech. 85, No. 2, pp. 141-146. 

* F. James P.-Y. Lagrée, M. Legrand, Draft

* P.-Y. Lagrée (2000): 
[An inverse technique to deduce the elasticity of a large artery](http://www.lmm.jussieu.fr/~lagree/TEXTES/ARTERE/ap9016.pdf)
European Physical Journal, Applied Physics 9, pp. 153-163 

*/

 
   
