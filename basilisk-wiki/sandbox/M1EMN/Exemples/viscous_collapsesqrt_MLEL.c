/**
#    Collapse of a rectangular viscous column along a slope 

Or collapse of a viscous fluid along a slope

From the paper by Huppert 82  "Flow and instability of a viscous current along a slope"


The 1D theory is presented  [with shallow water](http://basilisk.fr/sandbox/M1EMN/Exemples/viscolsqrt.c), here we solve the full profiles

We solve shallow water in 1D using [Multilayer horizontaly Euleur vertically Lagrange](http://basilisk.fr/src/layered/hydro.h) of Popinet 2020 with one layer only. 
$$
 \left\{\begin{array}{l}
         \partial_t h+\partial_x (hu)=0\\
	\partial_t hu+ \partial_x hu^2  
= -g h \partial_x \eta  -  C_f Q h^{-2},\;\;\; \eta = h+Z
        \end{array}\right. 
$$
 where the friction is laminar $C_f=3$ and
with at  t=0, a given heap $h(|x|<.5,t=0)=1$ so that $\int_{-\infty}^{\infty} h(x,0)dx =1$.

as the flow is very viscous, the convective term is small, 
and there only is a balance between pslope and viscosity, this equation has a self similar solution, we compare to it.

*/


//#include "grid/cartesian1D.h"
//#include "saint-venant.h"
/**
This is NOT by solved with [http://basilisk.fr/src/saint-venant.h]() but by [http://basilisk.fr/src/layered/hydro.h]() */
#include "grid/multigrid1D.h"
#include "layered/hydro.h"

 
// declare parameters
double Cf,alpha;
double tmax;


/**
position of domain `X0`, length `L0`, no dimension `G=1`
run with 1024*32 points to fix the tail of the falling heap (less is not enough!!!).
In the dam break case it was OK was far more  few points  
[http://basilisk.fr/sandbox/M1EMN/Exemples/damb_MLEL.c]() 
*/

int main() {
  X0 = -5;
  L0 = 30;
  N =  1024*4;
  G = 1.;
  alpha=0.5;
  Cf = 3.;
  tmax=1000;
  
  run();
}
 
/**
Boundary conditions (not used)

*/

u.n[left]   = neumann (0); 
u.n[right]  = neumann (0); 
h[left] = neumann (0);
eta[left] = neumann (0);
h[right] = neumann (0);

/** start by a rectangular column of surface 1

*/
event init (i = 0)
{
  foreach(){
     h[] =   (fabs(x)<.5);
    zb[] = -(x-X0)*alpha;}
}

/** split for friction 
 $$\frac{\partial u}{\partial t} = -C_f \frac{u}{h^2} \text { so that a semi implicit discretization gives }
 \frac{  u^{n+1} - u^n }{\Delta t} = -C_f \frac{u^{u+1}}{h^{n2}}$$
 we impose $u^{n+1}=0$ when $h^n$ very small (parameter `dry`).
*/
event friction (i++) {
  // Poiseuille friction (implicit scheme)
  foreach()
  {  	 double ff = h[] < dry ? HUGE : (1. + Cf*dt/h[]/h[]);
    u.x[] /= ff;
    }
}

/** saving  
*/
event outputfile (t <= tmax;t+=50) {
  foreach()
  fprintf (stderr,"%g %g %g %g \n", x, h[], zb[], t );
  fprintf (stderr,"\n");
}

/**
end of subroutines
*/


/**
#    Run

To compile and run

~~~bash  
qcc -g -O2 -DTRASH=1 -Wall -o viscous_collapsesqrt_MLEL viscous_collapsesqrt_MLEL.c -lm
./viscous_collapsesqrt_MLEL > vviscous_collapsesqrt_MLEL.out 2> viscous_collapsesqrt_MLEL.out2 
~~~

using `make`
 
~~~bash
make viscous_collapsesqrt_MLEL.tst
make viscous_collapsesqrt_MLEL/plots
make viscous_collapsesqrt_MLEL.c.html
source ../c2html.sh viscous_collapsesqrt_MLEL
open viscous_collapsesqrt_MLEL.c.html
~~~
 
#    Results 
 
  
~~~gnuplot collapse
 set xlabel 'x'
 set ylabel 'h(x,t)'
 set arrow from 0,.47 to 0,0.08
 set label "t" at 0.2,0.1
 set arrow from 3,.05 to 7,0.05
 set label "t" at 8,0.06
 p [-5:][:1] 'log'   w l not
~~~
 
gives an example of collapse

 
 
 
~~~gnuplot self similar collapse
 reset
 set title 'Self Similar Solution for a viscous collapse, t=200, 400...1000'
 set key top left
 set xlabel "x/t^{1/3}"
 set ylabel "h(x,t) t^{1/3}"
 n=1./3
 p [-0.5:1.5]'log' u ($1/($4**n)):($4>200? $2*($4**n): NaN) t'comp.' w l,sqrt(x/(.5)) t'anal''
~~~

 
 

plus on augmente la résolution, moins la partie gauche va bouger: ici on voit que le bloc a glissé, ce qui est une erreur numérique.
Ce glissement est dû à l'étape visqueuse 



# Links
* same example   [with classical shallow water](http://basilisk.fr/sandbox/M1EMN/Exemples/viscolsqrt.c)
* same example with [Multilayer Saint-Venant with mass exchange](http://basilisk.fr/sandbox/M1EMN/Exemples/viscous_collapsesqrt_ML.c)

 
 with Bingham

 * [Bingham 1D collapse on a incline](bingham_collapse_noSV.c)


# Bibliography
* [Huppert H.](http://www.itg.cam.ac.uk/people/heh/Paper49.pdf)
 "Flow and instability of a viscous current along a slope"
 Nature volume 30 1982 p 427   
* Lagrée  [M1EMN
Master 1 Ecoulements en Milieu Naturel](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/MFEnv.pdf) 



 */

 
