/**
#   Dam Break  on inclined plane, with Chézy turbulent friction (toward kinetik wave) 

A reservoir held back by a dam suddenly flows down a mountain, 
at early stages we see the discharge of the dam, the turbulent friction slows down the flow. Finaly a self similar wave (solution of the associated kinetik wave) appears.


Early time: 
 ![animation  ](damb_slope_turb/animatezoom.gif)

at longer times we see the flow on the slope (the slope has been diminished, so not to scale), a self similar solution develops.
 ![animation  ](damb_slope_turb/animate.gif)

We solve the Saint-Venant equation with turbulent friction in 1D on a slope :
$$\left\{\begin{array}{l}
         \partial_t h+\partial_x (hu)=0\\
         \partial_t (hu)+ \partial_x (hu^2)= - g h\partial_x h     - g h\partial_x z_b -\frac{c_f}{2}u^2\\
  \end{array}\right.$$
With at  $t=0$, a given triangular heap as  $h(0<x<\sqrt{2},t=0)=x$ of unit surface.

We compare the numerical solution with the analytical self similar solution of the kinematic wave

## Kinematik wave 

For long times inertia can be neglected, so that only slope and friction remain, in Chézy'case:
$$ Q  = h \sqrt{\frac{2  g \alpha h}{c_f}}$$
Mass conservation gives the kinetik wave (flood wave)
$$
\frac{\partial h}{\partial t} + \frac{\partial Q}{\partial x} = 0, \text{ with } Q = \beta h^n 
$$
with 
$\beta = \sqrt{\frac{2 g  \alpha }{c_f} }$ and $n=3/2$, 
with initial  distribution of height: 
    $\int_{x_0}^{x_f} h \, dx = V_0$.


 

The solution is self similar with $\eta=xt^{-2/3}$ and $h=t^{-2/3}\mathcal{H}(\eta)$, etc.

we find $\mathcal{H}(\eta) = \left( \frac{\eta}{ n} \right)^{\frac{1}{n-1}}$ so that:
$$h =     \left( \frac{ 2 x}{  3 \beta   t} \right)^2.$$
for $x<x_f$ position of the front 
$$x_f=\frac{3}{2^{2/3}}   V_0^{\frac{1}{3}} (\beta  t)^{\frac{2}{3}}= 
  \left(3 V_0 \right)^{\frac{1}{3}} \left(  \frac{3 \beta t}{2}\right)^{\frac{2}{3}}.
$$
at front for $x_f(t)$ or $t(x_f)$:
$$h_{max}(t(x_f))=\frac{2^{2/3} V_0^{2/3}}{(\beta t(x_f))^{2/3}}$$
note that we have 
$$ h_{max}(x_f)=\frac{3}{V_0},\text{ and }
 h/{h}_{\mathrm{max}} = (x/x_f)^{2}$$


## Code

We use here the [Multilayer Euler Lagrange solver](layered/hydro.h).
Of course this can be done with [shallow water solver](saint-venant.h) in 1D as well.
*/
#include "grid/multigrid1D.h"
#include "layered/hydro.h"
//#include "saint-venant.h"
// declare parameters

double tmax,Lbump,hbump,LLake,Cfs2,e;
/**
 position of domain `X0`, length `L0`, no dimension `G=1`
run with 1024 points  
*/

int main()
{
    X0 = -2;
    L0 = 150; 
    N =  1024*2*2; 
    nl = 1;
    G = 1.0;
    Cfs2 = 1;
    DT = .01/4;
    LLake = sqrt(2);
    tmax=600;  
    run();
}
/** start by a horizontal lake on constant unit negative slope  
 */
event init (i = 0)
{
    foreach(){
        zb[] = - x;
        h[]=  (x)*(x>0)*(x<=LLake);
    }
}

/** 
We have  Turbulent friction (implicit scheme) treated as a split
$$
  \partial_t (u) =  
   -\frac{c_f (u^2)}{2}
$$
*/  

event friction (i++) {
  foreach()
  {       double ff = h[] < dry ? HUGE : (1. + Cfs2*fabs(u.x[])*dt/h[]);
    u.x[] /= ff;}
}

/**
Error between the self similar solution and the computed one

*/
event errorS (t+=2.0) {
  double H,eta,pt23;
  e = 0;
  pt23 = pow(t+.0000001,2./3);
foreach()
  {  eta = (x)/pt23;  
     H = 4./9*(eta*eta)*(eta<3./pow(2,2./3))*(eta>0);
     e = e+fabs(h[]*pow(t,2./3) - H)*Delta/pt23;
    } 
    e=e/L0;
    fprintf (stdout,"%f %f \n",t,e);
}

/**
## Outputs

Output for plot at at $t<2$ and larger

*/
event outputfile1 (t+=.2;t <= 2) {
  foreach()
  fprintf (stderr,"%g %g %g  \n", x, h[],t); 
} 
 
event outputfile (t+=10;t <= tmax) {
  foreach()
  fprintf (stderr,"%g %g %g \n", x, h[],t); 
}
/** animation
*/
event animatedplot (t+=1.0) {
    static FILE * fp = popen ("gnuplot -persist 2> /dev/null", "w");
    if(t==0) fprintf (fp,"set term gif animate;set output 'animate.gif';set size ratio .333333\n");
    fprintf (fp,"\nset grid\n");
    fprintf (fp,"set title 'Turbulent flow on a slope 1D --- t= %.2lf '\n"
      "t= %.2lf  ; "
      "p[%g:%g][-1:.5]  '-' u 1:($2+$4) t'free surface' w l lt 3,'-' u 1:4 t'topo/L0' w l\\\n",
           t,t,X0,X0+L0);
    foreach()
    fprintf (fp,"%g %g %g %g %g\n", x, h[], u.x[], zb[]/L0, t);
    fprintf (fp,"e\n\n");
    fflush (fp);
}

event animatedplotzoom (t+=.01,t<2) {
    static FILE * fp = popen ("gnuplot -persist 2> /dev/null", "w");
    if(t==0) fprintf (fp,"set term gif animate;set output 'animatezoom.gif';set size ratio 1\n");
    fprintf (fp,"\nset grid\n");
    fprintf (fp,"set title 'Dam-Break 1D, early stages t= %.2lf <2'\n"
      "t= %.2lf  ; "
      "p[-1.25:2.75][-2.75:1.25] '-' u 1:($2+$3) t'free surface' w l lt 3,'-' u 1:3 not w l -1 \\\n",
           t,t);
    foreach()
    fprintf (fp,"%g %g %g  \n", x, h[], zb[]);
    fprintf (fp,"e\n\n");
    fflush (fp);
}
 
/**
end of subroutines

 
##  Run

To compile and run by hand or  using `make`
 
~~~bash

qcc damb_slope_turb.c  -o damb_slope_turb; ./damb_slope_turb  2> log  1>out


make damb_slope_turb.tst
make damb_slope_turb/plots
make damb_slope_turb.c.html
source ../c2html.sh damb_slope_turb
open damb_slope_turb.c.html

~~~
 
 
##   Results
 
 
An example of collapse along a unit slope with Chézy friction.
At early times:
 
~~~gnuplot collapse
set xlabel 'x '
set ylabel 'h(x,t)'
p [-1:4]'log' u 1:($3<5.?($2-$1):NaN) not w l
~~~

flow a  larger time, the wave is longer and longer, but of decreasing height.
  It looks a self similar solution...
  
  
~~~gnuplot collapse
set xlabel 'x '
set ylabel 'h(x,t)'
p [-1:][0:.25]'log' u 1:($3==10||$3==20||$3==40||$3==75 ?$2:NaN) t'num t=1,20,40,80,150,300,600' w l lc 1,\
'log' u 1:($3==150||$3==300||$3==600?$2:NaN) not w l lc 1
~~~




Computation of the error (approaching the self similar solution)
$$e(t)= \int (\mathcal{H} - h t^{2/3})^2 d \eta$$
so that for $t >300$ the self similar solution is obtained, this is for $x \sim 90$ as 
(3^(1./3)   (3./2 300)^(2./3)=85)
  
 
~~~gnuplot time to selfsimilar
set xlabel 't '
set ylabel 'error'
p [:]'out' 
~~~




The solution is indeed selfsimilar

~~~gnuplot collapse self similar 
reset
set xlabel 'eta'
H(x)=4./9*(x)**2*(x>0)*(x<3/2**(2./3))
p [-1:3][0:2]'log' u ($1/$3**(2./3)):($3>100?$2*$3**(2./3):NaN) w l t 'cal', H(x) t'self',2**(2./3) not
~~~

## Conclusion


We obtain the self similar solution of the kinematik wave equation as a solution at long time of the Shallow-water equation.
Height is:
$$h =     \left( \frac{ 2 x }{  3 \beta   t} \right)^2.$$
for $x<x_f$ position of the front 
$$x_f=\left(3 V_0 \right)^{\frac{1}{3}} \left(  \frac{3 \beta t}{2}\right)^{\frac{2}{3}}.$$


Exactly the same analysis can be done with Manning's friction where 
$\tau_0 = \rho g n_{GM}^2 Q^2/h^{7/3}$, this time $Q=\beta h^n$ with $n= 5/3$, so $1/(n-1)= 3/2$ and  $\beta = \frac{\sqrt{\alpha}}{n_{GM}}$ so that:
  $$Q=\frac{\sqrt{\alpha}}{n_{GM}} h^{5/3}. $$  
  etc. 

 
 
## Links

* Ideal fluid solution [Ancey-Dressler](/sandbox/M1EMN/Exemples/damb_slope.c)

* see classical Dam-Break on horizontal bottom (no slope)
  [Ritter solution](/sandbox/M1EMN/Exemples/damb.c) 
  
* see non viscous dam break with [standard C](/sandbox/M1EMN/Exemples/svdb.c)
and with [Basilisk](/sandbox/M1EMN/Exemples/damb.c)  
                    
* see Dam-Break on horizontal bottom (no slope) Chanson's solution
  [with turbulent friction](/sandbox/M1EMN/Exemples/damb_dressler.c)
  this is another Dressler solution
 
* see all the viscous collapse examples ([with laminar friction](/sandbox/M1EMN/Exemples/viscous_collapse.c))

## Bibliography

* [Lagrée P-Y](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/MFEnv.pdf)
"Equations de Saint Venant et application, Ecoulements en milieux naturels" Cours MSF12, M1 Sorbonne 



 */
