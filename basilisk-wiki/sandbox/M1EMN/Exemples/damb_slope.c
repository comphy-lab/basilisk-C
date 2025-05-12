/**
#   Dam Break  on inclined plane, Ancey's solution

A reservoir held back by a dam suddenly flows down a mountain, at early stages we see the discharge of the dam

 ![animation  ](damb_slope/animatezoom.gif)

at longer times we see the flow on the slope (the slope has been diminished, so not to scale)

 ![animation  ](damb_slope/animate.gif)

We solve the Saint-Venant equation with no friction in 1D on a slope :
$$\left\{\begin{array}{l}
         \partial_t h+\partial_x (hu)=0\\
         \partial_t (hu)+ \partial_x (hu^2)= - g h\partial_x h     - g h\partial_x z_b\\
  \end{array}\right.$$
With at  $t=0$, a given triangular heap as  $h(-1<x<0,t=0)=1+x$ .

We compare with the analytical solution of Ancey et al. 2008  and Dressler 1958. 

## Code

We use here the [Multilayer Euler Lagrange solver](layered/hydro.h).
Of course this can be done with [shallow water solver](saint-venant.h) in 1D as well.
*/
#include "grid/multigrid1D.h"
#include "layered/hydro.h"
//#include "saint-venant.h"
// declare parameters

double tmax,LLake;
/**
 position of domain `X0`, length `L0`, no dimension `G=1`
run with 1024 points  
*/

int main()
{
    X0 = -2;
    L0 = 60; 
    N =  1024; 
    nl = 1;
    G = 1.0;
    DT = .001;
    LLake=1;
    tmax=8;  
    run();
}
/** start by a horizontal lake on constant unit negative slope  
 */
event init (i = 0)
{
    foreach(){
        zb[] = - x;
        h[]= (x+LLake)*(x<=0)*(x>=-LLake);
    }
}

/** 
If 1, this is the  second Huppert problem. 
We have  Poiseuille friction (implicit scheme) treated as a split
$$
  \partial_t (hu) =  
   + \frac{- 3 ( h u )}{
   h^2  }
$$
*/  
#if 0
event friction (i++) {
  foreach()
  {       double ff = h[] < dry ? HUGE : (1. + 3.*dt/h[]/h[]);
    u.x[] /= ff;}
}
#endif

/**
## Outputs

Output for plot at at $t<1$ and larger

*/
event outputfile1 (t+=.2;t <= 1) {
  foreach()
  fprintf (stderr,"%g %g %g  \n", x, h[],t); 
} 
 
event outputfile (t+=1;t <= tmax) {
  foreach()
  fprintf (stderr,"%g %g %g \n", x, h[],t); 
}
/** animation
*/
event animatedplot (t+=.1) {
    static FILE * fp = popen ("gnuplot -persist 2> /dev/null", "w");
    if(t==0) fprintf (fp,"set term gif animate;set output 'animate.gif';set size ratio .333333\n");
    fprintf (fp,"\nset grid\n");
    fprintf (fp,"set title 'Dam-Break 1D --- t= %.2lf '\n"
      "t= %.2lf  ; "
      "p[%g:%g][-1:1]  '-' u 1:($2+$4) t'free surface' w l lt 3,'-' u 1:4 t'topo/L0' w l\\\n",
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

To compile and run  using `make`
 
~~~bash
make damb_slope.tst
make damb_slope/plots
make damb_slope.c.html
source ../c2html.sh dambslope
open idambslope.c.html

~~~
 
 
##   Results
 
 
An example of collapse along a slope
 
~~~gnuplot collapse
set xlabel 'x '
set ylabel 'h(x,t)'
p [-1:]'log' not w l
~~~

Compare to Ancey analytical solution for free surface at small time 
(almost figure 7)
~~~gnuplot collapse
set xlabel 'x '
set ylabel 'h(x,t)'
p [-2:2]'log' u 1:($3<=1?$2:NaN) t'num t=0.2 0.4 0.6 0.8 1.0' w l,\
'../REFCASES/ex.2.txt' t'Ancey sol t=0.2'w l,\
'../REFCASES/ex.4.txt' t'Ancey sol t=0.4'w l,\
'../REFCASES/ex.6.txt' t'Ancey sol t=0.6'w l,\
'../REFCASES/ex.8.txt' t'Ancey sol t=0.8'w l,\
'../REFCASES/ex1.txt' t'Ancey sol t=1.0' w l

~~~



Compare to Ancey analytical solution for free surface figure 6 of the paper 

~~~gnuplot collapse
set xlabel 'x '
set ylabel 'h(x,t)'
p [-5:50][0:.3]'log'  u 1:($3==1||$3==2||$3==4||$3==8?$2:NaN) t'num t=1,2,4,8' w l lc 1,'../REFCASES/ex1.txt' t'Ancey sol t=1'w l,\
'../REFCASES/ex2.txt' t'Ancey sol t=2'w l,\
'../REFCASES/ex4.txt' t'Ancey sol t=4'w l,\
'../REFCASES/ex8.txt' t'Ancey sol t=8'w l 
~~~



Compare to Ancey analytical solution for free surface at intermediate time 

~~~gnuplot collapse
set xlabel 'x '
set ylabel 'h(x,t)'
p [-1:16][:.1]'log'  u 1:($3>=1?$2:NaN) t 'num t=1, 2, 3, 4,5..' w l,\
'../REFCASES/ex1.txt' t'Ancey sol t=1' w l,\
'../REFCASES/ex2.txt' t'Ancey sol t=2' w l,\
'../REFCASES/ex4.txt' t'Ancey sol t=4' w l 

~~~
  
Compare to Ancey analytical solution for free surface at larger time $t=8$

~~~gnuplot collapse
set xlabel 'x '
set ylabel 'h(x,t=8)'
p [-1:][:.05]'log'  u 1:($3==8?$2:NaN)t 'num t=8' w l,\
'../REFCASES/ex8.txt' t'Ancey sol t=8' w l 
~~~

 
## Links

* see classical Dam-Break on horizontal bottom (no slope)
  [Ritter solution](http://basilisk.fr/sandbox/M1EMN/Exemples/damb.c) 
  
* see non viscous dam break with [standard C](http://basilisk.fr/sandbox/M1EMN/Exemples/svdb.c)
and with [Basilisk](http://basilisk.fr/sandbox/M1EMN/Exemples/damb.c) (this file)
                    
* see Dam-Break on horizontal bottom (no slope)
  [with turbulent friction](http://basilisk.fr/sandbox/M1EMN/Exemples/damb_dressler.c)
  this is another Dressler solution
 
* see all the viscous collapse examples ([with laminar friction](http://basilisk.fr/sandbox/M1EMN/Exemples/viscous_collapse.c))

## Bibliography

 * [Ancey](https://lhe.epfl.ch/articles/2007wrr.pdf)  Rentschler, R.M. Iverson, R.P. Denlinger "An Exact Solution for Ideal Dam-Break Floods on Steep Slopes"
 Water Resources Research, 2008, 44(1): 1-10. 
 
* [Dresssler](https://royalsocietypublishing.org/doi/10.1098/rspa.1958.0177)
   1958 Unsteady non-linear waves in sloping channels, Proc. R. Soc. Lond. A247186–198

* [Lagrée P-Y](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/MFEnv.pdf)
"Equations de Saint Venant et application, Ecoulements en milieux naturels" Cours MSF12, M1 UPMC 

*  C. Ancey 2024 Private communication 

 */
