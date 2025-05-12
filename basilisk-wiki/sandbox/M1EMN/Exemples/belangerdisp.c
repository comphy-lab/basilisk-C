/** 
# Saint Venant Shallow Water with surface tension effects

## Scope
 
We test the Bélanger relation for velocity and wtaer height before and after a steady standing jump.
 Here we add Laplace pressue due to surface tension.
 
 ~~~gnuplot  configuration
 set arrow from 4,0 to 4,1.2 front
 set arrow from 4,1 to 4,0 front
 set arrow from 8,0 to 8,2. front
 set arrow from 0.1,.1 to 9.9,0.1 front
 set arrow from 9.1,.1 to 0.1,0.1 front
 set label  "L0" at 6,.15 front
 set label  "depth h(x,t)" at 2.9,1.4 front
 set label  "depth h2" at 8.2,2.2 front
 set label  "water" at 1.,.8 front
 set label  "air" at 1.4,2.05 front
 set xlabel "x"
 set ylabel "h"
 p [0:10][0:3]0 not,(x<5? 1+1*exp(1*(x-5))*cos(6*(x-5)) :2) w filledcurves x1 linec 3 t'free surface'
 reset
 ~~~
 
 
  
## Cases
### ideal fluid, no surface tension, no viscosity: the jump
Conservation de la quantité de mouvement
$$U_1^2h_1+g\frac{h_1^2}{2}= U_2^2h_2+g\frac{h_2^2}{2}
$$ 
 Conservation de la masse:
$$ U_1h_1=U_2h_2 $$
 
 
 $$\big(\frac{h_2}{h_1}\big)= \big(\frac{U_1}{U_2}\big)=
 \big(\frac{F_1}{F_2}\big)^{2/3}
 = \frac{-1+\sqrt{1+8 F_1^2}}{2}$$
 
### viscous  fluid, no surface tension: Burgers
 
 bla bla
 
### slightly viscous  fluid, with surface tension: dispersive jump
 That is the case we look at
 
# Code 
 */
#include "grid/cartesian1D.h"
#include "saint-venant.h"

double tmax,Fr,deltah,x_rs,W,sigma,nu;
scalar K[],hxc[],hxcm[];
int main(){
  X0 = 0;
  L0 = 256;
  X0 = -200;
  N = 1024*2;
  Fr = 1.5;    // Nombre de Froude initial
  G = 1;
  tmax = 150;
  
  W = 0;       //vitesse du ressaut, ici 0
  sigma = 0.5;  // surface tension
  nu = 0.0625;      // additional viscosity
  x_rs=0;    // position nitiale du saut dans [X0;X0+L0]
  DT = 1.e-3;   // 0.25 5.e-3;
  run();
}
/** 
sortie libre
*/

u.n[left]  = dirichlet(1.5);
h[left]    = dirichlet(1);

h[right]   = neumann(0);
u.n[right] = neumann(0);//radiation(Fr/((-1 + sqrt(1+8*Fr*Fr))/2));

K[right]   = dirichlet(0);//neumann(0);
hxc[right]  = dirichlet(0);//neumann(0);
hxcm[right] = dirichlet(0);//neumann(0);

K[left]    = neumann(0);
hxc[left]   = neumann(0);//dirichlet(0);
hxcm[right] = neumann(0);//dirichlet(0);
 /**
 initialisation: Belanger with $h_1=1$ and $h_2$ computed
*/
event init (i = 0)
{
  foreach(){
    zb[]=0.;
    //  Bélanger formula with h1=1;
    deltah=(((-1 + sqrt(1+8*Fr*Fr))/2)-1);
     h[]=1+(((-1 + sqrt(1+8*Fr*Fr))/2)-1)*(1+tanh((x-x_rs)/2))/2;
    // vitesse associée par conservation du flux Fr/h + translation
     u.x[]=Fr/h[]+W;
    }
  boundary({zb,h,u});
#ifdef gnuX
  printf("\nset grid\n");
#endif
}
/**
### Compute the approximate curvature:
  $$K = \frac{\partial^2 h}{\partial x^2}$$
wigh gives Laplace pressure:
 $$p = -\sigma \frac{\partial^2 h}{\partial x^2}$$
*/
event courv(i++)
{
/**
 Compute the a centered derivative:
*/
foreach()
    hxc[] = (eta[+1]-eta[-1])/2/Delta;
boundary ({hxc});

/**
 mean value gives an estimation of the derivative with an error $O(\Delta^2)$
 $$
 (2 \frac{h(1) - h(-1)}{2 \Delta} +
 \frac{h(0) - h(-2)}{2 \Delta}+
 \frac{h(2) - h(0)}{2 \Delta})/4 = 
 h_x +  \frac{5 \Delta^2}{12}h_{xxx}$$
*/
foreach()
    hxcm[] =  (hxc[1] + 2* hxc[]+ hxcm[-1])/4 ;
boundary ({hxcm});
/** derivative
 $$ h_{xx} + 1/2 \Delta    h_{xxx} + (7/12)  \Delta^2  h_{xxxx}+   (1/4)  \Delta^3  h_{xxxxx}+... $$
 the curvature, second order derivative is:
 */
foreach()
    K[] = (hxcm[0] -  hxcm[-1])/(Delta);
boundary ({K});
/**
 Curvature contribution to pressure
 $$\frac{\partial u}{\partial t}=
 + \sigma \frac{\partial K}{\partial x} \text{ or }
 \frac{\partial u}{\partial t}=
 + \sigma \frac{\partial^3h}{\partial x^3}
 $$
 */
foreach()
    u.x[] += dt*(sigma*(K[1] - K[0])/Delta) ;
 boundary ({u.x});
/**
 all the differences all together give:
 $$h_{xxx}+ (1/2)  \Delta^2  h_{xxxxx}+...$$
 precision is $O(\Delta^2)$.
*/
}
/**
 Some viscous dissipation; at the entrance and at the output....
 $$\frac{\partial u}{\partial t}=
 + \nu \frac{\partial^2 u}{\partial x^2} $$ */
event visq(i++)
{
foreach()
    K[] = ((x<X0+8)||(x>X0+L0-32))*nu*(u.x[0]-u.x[-1])/Delta;
boundary ({K});
    foreach()
    u.x[] += dt*((K[1] - K[0])/Delta) ;
 boundary ({u.x});
}
/**
 
 Plots
 */
#ifdef gnuX
event plot (t<tmax;t+=0.1) {
    printf("set title 'Ressaut dispersif 1D  --- t= %.4lf '\n"
	 "p[%g:%g][-.25:2]  '-' u 1:($2+$4) t'free surface' w l lt 3,"
	 "'' u 1:3 t'velocity' w l lt 4,\\\n"
	 "'' u 1:4 t'topo' w l lt -1,\\\n"
     "'' u 1:(sqrt($3*$3/($2*%g))) t'Froude' w l lt 2,\\\n"
     "'+' u (%lf):(1+%lf) t'jump theo' \n",
           t,X0,X0+L0,G,x_rs+W*t,deltah/2);
    foreach()
    printf ("%g %g %g %g %g\n", x, h[], u.x[], zb[], t);
    printf ("e\n\n");
}
#endif
event end(t=tmax ) {
    foreach()
    fprintf (stderr,"%g %g %g %g %g\n", x, h[], u.x[], zb[], t);
}


/**
# Run and Results
## Run
To compile and run:

~~~bash
 qcc  -DgnuX=1 belangerdisp.c -o belangerdisp -lm; ./belangerdisp  2>log | gnuplot
 
 make belangerdisp.tst
 make belangerdisp/plots;
 source  c2html.sh belangerdisp
 
~~~

## Plot of jump
 
 
~~~gnuplot  result, free surface (blue)
set xlabel "x"
 set size ratio .33333
 p [:56][0:3] 'log' t'free surface' w l lc 3
 
~~~

## Plot of jump
 
 Plot of jump comparision with analytical result for `Fr=1.5` and `h1=1` we have
 `h2=(((-1 + sqrt(1+8*Fr*Fr))/2))`
 ~~~gnuplot  result, free surface (blue) and bottom (black)
 set xlabel "x"
 Fr=1.5
 h2=(((-1 + sqrt(1+8*Fr*Fr))/2))
  set size ratio .4
 p [-150:56][-1:3] 'log' t'free surface'w l lc 3,'' u 1:3 w l t'speed',\
'' u 1:4 t'zb' w l lc -1,(x>0?h2:NaN) t'h2'  w l lc 1,(x<0?1:NaN)  t'h1' w l lc 1
 
 ~~~
 
# Liens
 
 * Cas dispersif [http://basilisk.fr/sandbox/M1EMN/BASIC/disperse.c]()
 
 * [http://basilisk.fr/sandbox/M1EMN/Exemples/airy_watertrainfront.c]()
 
 * [http://basilisk.fr/sandbox/crobert/1_Layered/tension.h]()
 
 * [http://basilisk.fr/sandbox/M1EMN/Exemples/ressaut_mascaret.c]()
 
# Bibliographie
 
* [Lagrée P-Y](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/MFEnv.pdf)
"Equations de Saint Venant et application, Ecoulements en milieux naturels" Cours MSF12, M1 UPMC
 
 
*/


