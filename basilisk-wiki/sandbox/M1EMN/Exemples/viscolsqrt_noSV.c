/**
# "full Huppert problem" 1D
## Equations 

We look at the following problem the
 the viscous collapse with hydrostatic pressure gradient of a lubricated film on a inclined plate of contant angle. 
This is as well a diffusive wave problem


$$ \frac{\partial}{\partial t}  h  + \frac{\partial}{\partial x}  Q = 0,\text{ with } Q= \frac{h^3}3 (1 - \frac{\partial h}{\partial x})$$
 
 
 
## "Huppert's problems" 

This is a mix of (Huppert first problem on a horizontal plate and second problem on an inclined plate:
 
Huppert's first problem is 
$$\frac{\partial}{\partial   t}   h = \frac{\partial}{\partial   x}(\frac{h^3}{3}    \frac{\partial}{\partial   x}   h )$$
whereas Huppert's second one is 
$$\frac{\partial}{\partial t} h + h^2 \frac{\partial}{\partial   x} h =0 $$

so that we have both effects here. This is refered as "Spreading of a shallow mass on an incline" by Mei. 
See bibliography


## Numerics 

We do not use here the shallow water solver (as in the two links below).
We do a splitting to deal with the two parts of $Q$.

~~~gnuplot
set size ratio .3
set samples 9 
set label "h i-1" at 1.5,3.1
set label "h i" at 2.5,3.15
set label "h i+1" at 3.5,2.5
set xtics ("i-2" 0.5, "i-1" 1.5, "i" 2.5,"i+1" 3.5,"i+2" 4.5,"i+3" 5.5)
set arrow from 2,1 to 2.5,1
set arrow from 3,1 to 3.5,1
set label "Q i" at 2.1,1.25
set label "Q i+1" at 3.1,1.25

set label "x i-1/2" at 1.5,0.25
set label "x i" at 2.4,0.25
set label "x i+1/2" at 3.,0.25

set label "x"  at 0.5,2+sin(0) 
set label "x"  at 1.5,2+sin(1)
set label "x"  at 2.5,2+sin(2) 
set label "x"  at 3.5,2+sin(3) 
set label "x"  at 4.5,2+sin(4)  
f(x) = (2+sin(x))*(x<=4)
p[-1:7][0:4] f(x) w steps not,f(x) w impulse not linec 1
~~~

The full problem is 
solved with a splitting $Q=Q_c+ Q_d$,  each part has a different behavior.
$$ Q_c=\frac{h^3}3  \;\text{ and }\; Q_d= -\frac{h^3}3  \frac{\partial h}{\partial x}$$
where subscript "c" refers to the "advective" part with a celerity $c_0=h^2$:

$$\frac{\partial}{\partial t} h + c_0 \frac{\partial}{\partial x} h = 0 \;\;\text{ with }\;\; c_0 = \partial   Q_c/\partial h$$


and subscript "d" refers to the "diffusive" part of non linear viscosity $\frac{h^3}{3}$:
$$\frac{\partial}{\partial   t}   h  = -   \frac{\partial}{\partial   x}   Q_d $$
 
 The first leads to shocks, the second is diffusive



## Code

Mandatory declarations, libraries and tables. Note we use the face vectors.
 
*/
#include "grid/cartesian1D.h"
#include "run.h"

scalar h[];
scalar zb[];

face vector  hx[];
face vector Q[];
face vector c0[];
face vector Qc[];
face vector Qd[];

double alpha;
double X0,L0;
double DT;
double DD;

/**
Main with definition of parameters :


Note that we run the code twice, the first run is without diffusion `DD=0`
*/


int main() 
 {
  L0 = 15;
  X0 = -1;
  N = 200;
  DT = (L0/N)*(L0/N)/20;
  DD = 0;
  run();
  DD = 1 ;
  run();
}
/**
Boundary conditions

*/

h[left] =  neumann (0);
hx.n[left] = neumann (0);
Q.x[left] = dirichlet(0);

/**
Initial elevation: a " quare" of surface $A=2$  
*/
event init (t = 0) 
 {
  foreach()
    h[] = (fabs(x-1)<1);
    
  foreach_face()
   {   
    Q.x[]=0;
    c0.x[]=1;
  }
}
/**
print data in stdout

*/

event printdata ( t+=5; t <= 100)
 {
  foreach()
    fprintf (stdout, "%g %g %g %g %g \n", x, h[], Q.x[], t, DD);
  fprintf (stdout, "\n");
}

/**
Integrating the domaine over time, updating the elevation at each time step.

*/
event integration (i++)
 {
    double dt = DT;
   
    
/**
finding the good next time step

*/
    dt = dtnext (dt);
/**
  the flux $Q(h)$ : split in two fluxes. 
  
 1-  The avective part
 
  First $Q_c$  the advective part gives a estimation of the height say   $h^*$, 
 
$$ \frac{h^* - h^n}{\Delta t} = -\frac{1}{\Delta x} (Q^n_{ai} - Q^n_{ai-1})$$  with $Q_a$ obtained from $Q_c$ with a correction

$$ Q_{ai} = \frac{Q_{ci} + Q_{ci-1}}{2}  - c \frac{h_j - h_{j-1}}{2}\text{  with }  c_0 = \frac{\partial Q_{c}}{\partial h}$$

*/
    foreach_face()
      {
      double hm = ((h[0,0] + h[-1,0])/2);
      Q.x[] =  (pow(hm,3)/3) ;
      c0.x[] = pow(hm,2) *(1-DD);
      Qc.x[] = (Q.x[]+Q.x[-1])/2.  - (c0.x[] *(h[]-h[-1])/2);
      }
  
    foreach()
        h[] = h[] - dt*( Qc.x[1,0] - Qc.x[0,0] )/Delta;
       

  
/**

2- The diffusion part : calculate the next step elevation $h^{n+1}$,
`n`: the temporal index.

$$ \frac{h^{n+1} - h^*}{\Delta t} = -\frac{1}{\Delta x} (Q^*_{di} - Q^*_{di-1})$$

 
*/

    
    foreach_face()
     {
      hx.x[] = ((h[0,0] - h[-1,0] )/Delta);
      double hm = ((h[0,0] + h[-1,0])/2);
      Qd.x[] = -(pow(hm,3)/3) *  hx.x[];
    }
    
/**

final upadte, note the `DD` switch.  
the first run is without diffusion `DD=0`, the second is with: `DD=1`
*/
    foreach()
      h[] = h[] - dt*DD*( Qd.x[1,0] - Qd.x[0,0] )/Delta;
    
}




/**
## Plots 

Solution of the full problem (this is as well a diffusive wave problem). 

~~~gnuplot
reset
set xlabel "x"
set size ratio .25
p'out' u 1:($5==1?$2:NaN)w l
~~~


Solution of only the advected problem (second problem of Huppert, inclined plate, this is as well a  kinetic wave problem); swith `DD=0`.


~~~gnuplot
reset
set xlabel "x"
set size ratio .25
p'out' u 1:($5==0?$2:NaN)w l
~~~


Both superposed to see differences
 
~~~gnuplot
reset
set xlabel "x"
set size ratio .25
p'out' u 1:($5==0?$2:NaN) t 'adv' w l,'' u 1:($5==1?$2:NaN)t 'full' w l
~~~





## Plots of self similar solution

Wave solution of form : $\sqrt{x/t}=t^{-1/3} \sqrt{x/t^{1/3}}$
the position of the front is $x_f= (3  A /2)^{2/3} t^{1/3}$, where $A* is the inition surface, here $A=2$ hence $\eta_f=2.088$

the solution is self similar at long times,


~~~gnuplot
reset
set xlabel "x"
set size ratio .25
xf=(3.*2/2)**(2./3)
p'out' u 1:($5==0?$2:NaN)w l,'' u 1:($5==1?$2:NaN)w l,\
'' u 1:(sqrt($1/$4)*($1<xf*($4**(1./3)))) w l
 
~~~


 


the solution is self similar at long times, we plot in the selfsimilar space, 
we shift of a value of $x$ equal to 1 the solution of the full problem (diffusive wave) to see better the differences induced by the addition of pressure gradient in the kinetic wave 


~~~gnuplot
reset
set xlabel "x"
set size ratio .25
xf=(3.*2/2)**(2./3)
p[-1:4]'out' u ($1/$4**.33333333+$5):($2*$4**.3333333) not w l,\
sqrt(x)*(x<xf),sqrt(x-1)*(x<xf+1)
 
~~~


 

# Links
 * [http://basilisk.fr/sandbox/M1EMN/Exemples/viscolsqrt.c]()
 * [http://basilisk.fr/sandbox/M1EMN/Exemples/viscous_collapse.c]()
 
 
# Bibliography
 * [Huppert H.](http://www.itg.cam.ac.uk/people/heh/Paper47.pdf)
 ”The propagation of two-dimensional and axisymmetric viscous gravity currents over a rigid horizontal surface” J . Fluid Mech. (1982), vol. 121, p p . 43-58
 * [Huppert H.](http://www.itg.cam.ac.uk/people/heh/Paper49.pdf)
 "Flow and instability of a viscous current along a slope"
 Nature volume 30 1982 p 427  
 * Chiang C. Mei,  2007  [Spreading of a shallow mass on an incline](https://web.mit.edu/1.63/www/Lec-notes/chap2_slow/2-4spread-mud.pdf)
 * Lagrée  [M1EMN
Master 1 Ecoulements en Milieu Naturel](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/MFEnv.pdf)
 * Lagrée  [M2EMN
Master 2 Ecoulements en Milieu Naturel](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/mainM2EMN.pdf)



*/
