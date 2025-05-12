/**
# Saint-Venant Von K&aacute;rm&aacute;n: Interactive Boundary Layer Shallow Water equations (SVVK/IBLSW)
## Scope
In this file, we present an extension of the Shallow Water equations (Saint-Venant) with an extra layer of viscous flow (Von K&aacute;rm&aacute;n). Hence, the basic profile is no more an half Poiseuille, but may be flat in the core flow with a boundary layer near the wall. Those two layers interact through friction and exchange of momentum. This extension allows to better estimate the friction at the bottom, especially when starting from a flat profile (boundary layer growth), and on short bumps (acceleration in fluvial flow).



## Long Wave approximation of Navier Stokes Eq. (RNSP)
Starting from the Long Wave approximation of Navier Stokes (RNSP see Lagrée &amp; Lorthois 05), $O<y<h(x,t)$
$$
\partial_t  u+
u 
\partial_x  u + v \partial_y u = - \partial_x p + \partial_y^2u,
$$
$$ 0 = - \partial_y p -1
$$
$$
\partial_x  u + \partial_y v=0
$$$$
u(0,y<h(0))=u_0,\;\;u(x,z_b(x))=0,\;\;v(x,z_b(x))=0,\;$$
$$\partial_y u(x,y=h(x))=0,\;\;v(x,h(x))=\partial_t h + u(x,h(x))\partial_xh \;$$ 
these equations are without dimension.

## Usual Shallow Water Eq. (Laminar flow)


The classical way consists to integrate equations across the depth $h$, giving the usual Shallow Water Eq.
("Saint-Venant" in french), comming back with dimensions:
$$
 \left\{\begin{array}{l}
         \partial_t h+\partial_x h Q= 0\\
	\partial_t Q+ \partial_x \left[\alpha \dfrac{Q^2}{h}+ g\dfrac{h^2}{2}\right]
= - gh \partial_x z_b  + \tau_p \\
        \end{array}\right. 
$$
with $\tau_p=3  \nu Q/h^2$, $\alpha=6/5$, if one supposes that the velocity profile is always close to an half Poiseuille
(though in practice $\alpha=1$).
There are drawbacks of this equation, the most important is that shear stress is badly evaluated in accelerating flow...

Note that the scale is $h Re$ longitudinaly and $h$ transversally, with $Re=U_0 h_0/\nu$ (were $u_0=\sqrt{gh_0}$)
the `invsqrtRe` variable is reminiscent of this if the longitudinal scale is a given $L$ which is not $h (U_0 h_0/\nu)$, and $Re=U_0L/\nu$, and `invsqrtRe` is $1/\sqrt{U_0L/\nu}$

Note that $g$ is the `G` in `Basilisk`.

## SVVK Eq. (Laminar flow)
### Von K&aacute;rm&aacute;n and Saint-Venant
To better evaluate the shear, we propose boundary layer equations at the bottom of the flow,
we define $u_e$ the velocity at the surface, substract it from the RNSP equations, 
the equations are again integrated across the depth $h$.
Let define $q_1= \int_{z_b}^h (u_e-u)dy$  this is $u_e \delta_1$, were $\delta_1$ is the displacement thickness, 
$q_2= \int_{z_b}^h (u/u_e) (u_e-u)dy$ (or $u_e \delta_2$ with $\delta_2$ momentum thickness)
and $\tau_0=\partial_y u|_0$ the basal shear stress, let define as usual
$H = q_1/q_2$ the shape factor so that by integration  
the Von K&aacute;rm&aacute;n equation is recovered (Schlichting): 
$$
\partial_t q_1 + \partial_x (u_e q_2)+ q_1 \partial_x u_e = 
h [-\partial_x p + \partial_t u_e + u_e \partial_x u_e]  + \tau_0 
$$
note that the term in braces is not null as it is usual in Von Karman equation, but here, it is equal to
a quantity we call $-\tau_\pi/h$. Hence the Saint-Venant equation is recovered:
$$[-\partial_x p + \partial_t u_e + u_e \partial_x u_e]  = 
-\tau_\pi/h.$$
We have now an equation for the boundary layer (equation for $q_1=u_e \delta_1$) and an equation for the ideal fluid layer (in fact the free surface velocity $u_e$).

### Closure
This $-\tau_\pi/h$ needs a closure, just like $H$ and $\tau_0$, it is 0 for $\delta_1 \ll h$,
when there is a very thin boundary layer, it is equal to $3 u_e/h$ for 
 $\delta_1 = h/3$ when the boundary layer has merged with the layer.
 

So, we have to close the equations with $H$, $\tau_0$ and $\tau_{\pi}$ given and depending on the chosen velocity profile familly.
Here, $\tau_0=f u_e/\delta_1$, and  $\tau_{\pi} = 3 \dfrac{u_e}{h}\Pi(\dfrac{3 \delta_1}{h})$, with $\Pi(0)=0$, $\Pi(1)=1$.

with $f = f_2H$ (with $f_2$ and  $H$ constants, here, but variables in fact) for small $\delta_1$ and $f =1$ for  $\delta_1=h/3$.

Variable $f_2$, $H$ are given in Lagrée et al 05, Lagrée Lorthois 05  in the case $\tau_\pi=0$ 

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
On utilise le solveur de Saint-Venant et on utilise une grille simple cartésienne: 
*/
#include "grid/cartesian1D.h"
#include "saint-venant.h"
/**   declarations ... 
*/
scalar q1[],dxq1[],F[],tauPi[],tau0[],d1[];
double h0,alpha,tmax,q10,ue0,theta; 
double H=2.59,f2=1./2.59;//.23
double invsqrtRe=1 ; 
/**
given flat profile with thin boundary layer at the left
*/
u.n[left] = ue0; 
h[left] = 1;
q1[left] = q10;
F[left] = 0;
/**
right, output
*/
u.n[right] = neumann(0);
h[right] =neumann(0);
q1[right] = neumann(0);
F[right] = neumann(0);
/**
## The Main
The domains starts from `X0`, length `L0`, without dimension `G=1`, (it means that the velocity is measured by $\sqrt{gh_0}$),
were  $\alpha$ is the bump height, $-\theta$ the slope angle. 
*/
int main()
{  
  X0 = 0.0;
  L0 = 10;
  G = 1;
  N = 1024*4; 
  ue0 = 0.5;
  q10=0.01;
  // DT = (L0/N)/50;
  DT = HUGE[0];
  tmax=100;
  alpha=.1;
  theta=.5;

  /* L0=1;
    N=1024*16;
    invsqrtRe=1;
    tmax=1;  
    alpha=0;
    theta=0;*/
  run();
}
/**
the subroutines, 
the bottom: we put three bumps on a flat plate followed by a constant slope 
*/ 
event init (i = 0) {
  double xi,xii,xiii;	
  foreach(){
    xi   = (x - 0.05)/.01;
    xii  = (x - 0.4)/.1;
    xiii = (x - 0.35)/.01;	
    zb[] =  - theta *(x -0.1)*(x>0.1) + alpha * exp(- xi*xi)  + 2*alpha * exp(- xii*xii)   +  alpha * exp(- xiii*xiii);
    u.x[] = ue0  ;
    h[]= 1 ;
    q1[] =   q10;
    }
}
/**
## Subroutines
The viscous correction is computed here
*/
event boundary_layer(i++){
 double cDelta,coef;	
/** 
Upstream differentiation to evaluate $\partial_x q_1$ */
  foreach()
      dxq1[]= (u.x[]>=0 ? (q1[1] - q1[0])/Delta : (q1[0] - q1[-1])/Delta ) ;
  
/** 
Definition of the friction functions
 $\tau_\pi$ and $\tau_0$, as part of the closure, so they depend on the other variables
*/  
 foreach(){
     coef = pow(sin(pi/2*min((q1[]/u.x[])/(h[]/3),1)),16);    
     tauPi[] = 3*u.x[]/h[] * f2*H * coef  ;
     //tau0[] =  (f2*H)*u.x[]*u.x[]/q1[];  
 }
  
/**
The solution is done by spliting, first :  
  $$   \dfrac{h^*-h}{\Delta t} +\partial_x hu_e = 0$$
  $$ \dfrac{(h^*u_e^*) -(hu_e)}{\Delta t}+ \partial_x \left[hu_e^2+ \dfrac{h^2}{2}\right]
= - gh \partial_x z_b $$
are solved by `"saint-venant.h"` (well balanced). 
Next the viscous sources for Saint-Venant allow to update the estimated $h^*,u_e^*$ to the finals $h,u_e$:
  $$ \dfrac{h-h^*}{\Delta t} = \partial_x q_1$$
  
  $$\dfrac{u_e-u_e^*}{\Delta t}  = \frac{u_e^*}{h^*}\partial_x q_1  - \frac{\tau_{\pi }}{h^*}.$$
  
*/
  foreach(){
    h[] += dt * ( dxq1[] ) * invsqrtRe  ;	
    u.x[] +=  dt * (u.x[]+u.x[1])/2 * dxq1[] / h[] * invsqrtRe - dt * tauPi[]/h[] ; 
    }
 	
/**
Finaly the Von-K&aacute;rm&aacute;n equation, is written as
$$\partial_t q_1 + \partial_x F = u_e \partial_x q_1  -\tau_\pi  + \tau_0$$
were the flux is defined as
$$F = (u_e q_1 ( 1+\frac{1}H )), $$ it is computed here:
*/
  foreach() {	
   //cDelta = Delta/dt;
    cDelta = 0;
   // cDelta = 1; // upwind
    F[] = (1.+1./H) * ((q1[0]*u.x[0]+q1[-1]*u.x[-1])/2.  - cDelta*u.x[] *(q1[0]-q1[-1])/2);}

/**
we split in $$ \dfrac{ q_1^* -q_1}{\Delta t} + \partial_x F = u_e \partial_x q_1  -\tau_\pi$$
$$\dfrac{ q_1 -q_1^*}{\Delta t} =  \tau_0$$
 the update of $q_1$ is then without basal friction
*/
  foreach()
    q1[] += - dt* ( F[1,0] - F[0,0] )/Delta + dt * (u.x[]+u.x[1])/2 * dxq1[] - dt * tauPi[] + 0.00*dt*tau0[];
 
/**
and the final split $\partial_t q_1 =  \tau_0$ with $\tau_0=  (f_2H) u_e^2/q_1$ 
 final evaluation of $q_1$ with  friction $\tau_0$, so that 
$$\partial_t (q_1^2)/2 =   (f_2 H)u_e^2$$
hence we integrate explicitely to avoid the division by a 0 $q_1$:
$$q_1^2(t+/\Delta t)  =  q_1^2(t )  +  2 (f_2H)u_e^2 \Delta t$$
*/   
  foreach(){
   q1[] = (h[]>dry ?   sqrt(q1[]*q1[] + 2* dt * (f2*H)*u.x[]*u.x[]): 0 );
   q1[] = (u.x[]>0 ? q1[] : -q1[]);
   }	 	
}
/**
and this is done for the computations.

## Outputs
gnuplot output for Blasius/ Rayleigh solution at short time
*/
event field_short (t<.1;t+=.001) {	 	
    foreach()
    printf (" %g %g %g %g %g %g %g \n", x, h[], u.x[], zb[], q1[], dxq1[], t);
    printf ("\n\n");
  double fun;  
  FILE *fp = fopen("shorttime.txt", "w"); 
    foreach(){
    	fun = 3*f2*H*u.x[]/h[]/h[]/theta;
    	fun = 3*q1[]/u.x[]/h[];
    fprintf (fp," %g %g %g %g %g %g %g %g \n", x, h[], u.x[], zb[], q1[], tauPi[],(f2*H)*u.x[]*u.x[]/q1[] , fun );}
  fclose(fp);
}
/**
gnuplot output for SVVK solution at long time
*/
event field_long (t<tmax;t+=1) {	 	
    foreach()
    printf (" %g %g %g %g %g %g %g \n", x, h[], u.x[], zb[], q1[], dxq1[], t);
    printf ("\n\n");
      
  double fun;  
  FILE *fp = fopen("longtime.txt", "w"); 
    foreach(){
    	fun = 3*f2*H*u.x[]/h[]/h[]/theta;
    	fun =   3*q1[]/u.x[]/h[];
    fprintf (fp," %g %g %g %g %g %g %g %g \n", x, h[], u.x[], zb[], q1[], tauPi[],(f2*H)*u.x[]*u.x[]/q1[] , fun );}
    fclose(fp);
}
/**
End
*/
event end (t = tmax) {
    double fun; 	
 
    foreach(){
    	fun = 3*f2*H*u.x[]/h[]/h[]/theta;
    	fun =   3*q1[]/u.x[]/h[];
  //  fprintf (out," %g %g %g %g %g %g %g %g \n", x, h[], u.x[], zb[], q1[], tauPi[],(f2*H)*u.x[]*u.x[]/q1[] , fun );
  }   
}
/**
fin des procédures.

# Run and Results
## Run
To compile and run:


~~~bash  
qcc -g -O2 -Wall -o svvk svvk.c -lm ; ./svvk > out
~~~


## Plot of main fields
Plot of shear $\tau_0$ and $\tau_{\pi}$, free surface (blue) and bottom (black). In red the $\tau_\pi$, and in green the $\tau_0$, they differ when the flow is dominated by ideal fluid (at the entrance, and accelerated on the bumps). They merge at the end when the flow is a Nusselt: when the displacement thickness $\delta_1= q_1/u_e$  (in cyan) is equal to $h/3$.
~~~gnuplot  result, free surface (blue), displacement thickness (cyan) and bottom (black)
set xlabel "x"   
p[:1][0:2]'longtime.txt' u 1:6 t'tau_{Pi}' w l,'' u 1:7 t'tau_0'w l,''u 1:(5*($4+$1*.5))t 'z_b +\theta*x'w l linec -1,'' u 1:2 t'h' w l linec 3,'' u 1:($5/$3) t'delta_1' w l
~~~


## Small time/ space
At small time, one recovers the Blasius solution and Rayleigh

at small $x$, we have the Blasius problem $\partial_x (u_e^2 \dfrac{\delta_1}{H})  = f_2 H u_e / \delta_1$, hence
$q_1=\delta_1 u_e = \sqrt{2 f_2 H^2 u_e x}$
(comparison is not so good here, we need more points?).

at small $t$, we have the Rayleigh problem 
$\partial_t (u_e \delta_1 )  = f_2 H u_e / \delta_1$, hence
$q_1= \delta_1 u_e = \sqrt{2 f_2 H u_e^2 t}$

 
~~~gnuplot q_1(x,t)
set xlabel "x" 
set xlabel "t"
set zlabel "q_1"  
sp[:.01][0:.05]'out'u 1:7:5 w l,sqrt(2*2.59*x*.5)*(x<y/2.55? 1 :NaN) ,sqrt(2*y*.5*.5)*(x>y/2.55? 1 :NaN) 
~~~



## Large time/ space
far from the bump, the solution should be invariant in $x$ so that $\partial_x=0$, hence 

$u_e(h- \delta_1)$ is constant , 
$0 = -   (-\theta)  - \tau_\pi$ and $0= \tau_{\pi}- \tau_0$

with $\tau_{\pi}= 3 u_e /h = \tau_0 f_2 H  u_e/\delta_1$, hence $\delta_1=h/3$ and 
$u_e=  (h/3)\theta$


we have a region were the Nusselt (half Poiseuille) solution is valid.


~~~gnuplot  result far from the entrance: the Nusselt film is recovered
set xlabel "x"   
 p[1:2][0:1]'longtime.txt' u 1:6 t'tau_{pi}' w l,'' u 1:7 t'tau_0'w l,'' u 1:($5/$3) t'd1' w l,'' u 1:($2/3) t' h/3'w l linec 3,'' u ($3*($2-$5/$3)) t'flux' w l
~~~ 

Note there is a problem at the output (due to boundary conditions?) a kind of increasing exponential.

~~~gnuplot  result far from the entrance: the Nusselt film, normalized value
set xlabel "x"     
 p[:][0:2]'longtime.txt' u 1:($6/$7) t'tau_{pi}/tau_0' w l, '' u 1:(3*$5/$3/$2) t'3 d1/h' w l,1
~~~


## Conclusion
seems to work, need extra work to check and improve the closures.
Check signe for reverse flow...



V.1 Noeux les Mines, 04 Juillet 15

#Bibliography
  
* see
["svvk"](http://basilisk.fr/sandbox/M1EMN/Exemples/svvk_3.c) 

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

* Toro 
[finite Volumes](http://marian.fsik.cvut.cz/~bodnar/PragueSum_2012/Toro_1-FORCE.pdf)
*/

 
   

 
   
