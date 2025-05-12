/**
# Viscous hydraulic Jump

We want to reproduce the hydraulic jump of Higuera (1994) using Audusse et al (2011) multilayer algorithm. */

#include "grid/cartesian1D.h"
#include "saint-venant.h"
double Cf;

int main() {
  X0 = 0.14;
  L0 = 1. - X0;
  G  = 1.;
  N  = 128*2;
/**
one layer Poiseuille
*/
  nl = 1;
  Cf = 3;  
  run();
/**
one layer Watson
*/
  Cf = 2.2799;
  run();
/**
multi layer Higuera
*/
  nl = 25;
  nu = 1.;
  run();
}
/**
We impose boundary condition for $h$ and $\eta$. */

h[left] = dirichlet (.2);
eta[left] = dirichlet (.2);

h[right] = dirichlet (0.0);
eta[right] = dirichlet (0.0);

/**
## Initialization

We define a scalar field *hc* to check for convergence on $h$. */

scalar hc[];

event init (i = 0) {

  /**
  We set a constant velocity at the inlet and a free outlet. */

  for (vector u in ul) {
    u.n[left] = 5.90;
    u.n[right] = neumann(0.);
  }
  
  /**
  We initialize *h* and *hc*. */
  
  foreach()
    hc[] = h[] = 0.2;
}


/**
 friction, only in the case of one layer
*/
event frotte (i++) {
  if(nl==1){
   foreach()
  {    
    double ff = h[] < dry ? HUGE : (1. + Cf*dt/h[]/h[]);
    u.x[] /= ff;
    }
  }
}
/**
We check for convergence. */

event logfile (t += 0.1; i <= 50000) {
  double dh = change (h, hc);
  printf ("%g %g\n", t, dh);
  if (i > 0 && dh < 1e-5)
    return 1;
}

/**
## Output

We print the elevation and the stress, $C_f u/h$ or $\partial u/\partial z|_0$. */

event output (t = end) {
  vector u0 = ul[0];
  double tau ;
 
  foreach(){
     if(nl==1){tau = Cf*(u.x[]/h[]);}else{tau = 2.000*u0.x[]/(h[]/nl);}
    fprintf (stderr, "%g %g %g\n", x, eta[], tau);
  }
  fprintf (stderr, "\n\n");  
}

/**

## Run

compile and run:

~~~bash
 
qcc higuera.c -o higuera -lm
./higuera 2> log

make higuera.tst  
make higuera/plots
make higuera.c.html 

~~~


## Results

~~~gnuplot raw plot.
 
plot[][:2]  'log' u 1:2 t 'h' w l,'' u 1:3 t 'tau' w l,0 not w l linec -1

~~~



~~~gnuplot Comparison with Figure 2 of Higuera (1994).
set output 'f.png'
X0 = 169
X1 = 604
Y0 = 222.24
Y1 = 528
unset tics
plot [0:][0:605] '../Img/higueragraph.png' binary filetype=png with rgbimage not, \
  'log' i 2 u (X0+$1*(X1-X0)):($2/2*(Y1-Y0)+Y0) t 'h' w l,\
  '' i 2 u (X0+$1*(X1-X0)):($3/15*(Y1-Y0)+Y0) t 'tau' w l
~~~
 
 
 
We compare several solutions, the jump appears sooner with Poiseuille than with Watson. 
 
~~~gnuplot Comparison  with one layer (two closures, Red Poiseuille, Green Watson) and multilayer (blue)
reset
set xlabel "x"
set ylabel "h"
set y2label "{/Symbol t}_b"
set y2range [-1:30]
set y2tics border out 
H= .0

  p [0:][-.1:3] 0 not w l linec -1,1.813*(x<.15? x:NaN) not linec -1,\
            (x>.35? (H+(12*(1-x))**.25):NaN)  not linec -1,\
              'log' i 0 u 1:2 t'P' w l ls 1,'' i 0 u 1:($3) not w l ls 1 axes x1y2,\
                   '' i 1 u 1:2 t'W' w l ls 2,'' i 1 u 1:($3) not w l ls 2 axes x1y2,\
                 '' i 2 u 1:2 t'ML' w l ls 3,'' i 2 u 1:($3) not w l ls 3 axes x1y2
  
~~~

The same in two graphs 
 
~~~gnuplot Comparison  with one layer (two closures, Red Poiseuille, Green Watson) and multilayer (blue)
reset
set xlabel "x"
set ylabel "h" 
H= .0
 set multiplot layout 2,1 title "Jump"
  p [0:][-.1:3.5] 0 not w l linec -1,1.813*(x<.15? x:NaN) not linec -1,\
              'log' i 0 u 1:2 t'P' w l ls 1,\
                   '' i 1 u 1:2 t'W' w l ls 2,\
                 '' i 2 u 1:2 t'ML' w l ls 3
 set ylabel "{/Symbol t}_b"                
   p [0:][-1:25]0 not w l linec -1,\
              'log'   i 0 u 1:($3) t'P' w l ls 1,\
                   '' i 1 u 1:($3) t'W' w l ls 2  ,\
                 ''   i 2 u 1:($3) t'ML' w l ls 3  
  
  unset multiplot
~~~


## Links

  * [http://basilisk.fr/src/test/higuera.c]() with Basilisk (final source case)
  
  * [http://basilisk.fr/sandbox/M1EMN/Exemples/svdbvismult_hydrojump.c]() std C, no Basilisk

* [http://basilisk.fr/sandbox/M1EMN/Exemples/higuera_vs_SWradial.c]() axi case


## Bibliography

* Higuera, F. 1994. [The hydraulic jump in a viscous laminar
flow](http://dx.doi.org/10.1017/S0022112094002041). J. Fluid
Mech. 274, 69-92.

* Audusse Sainte-Marie 2011

* Francesco De Vita, Pierre-Yves Lagrée, Sergio Chibbaro, Stéphane Popinet 
[Beyond Shallow Water: appraisal of a numerical approach to hydraulic jumps based upon the Boundary Layer Theory.](http://www.lmm.jussieu.fr/~lagree/TEXTES/PDF/devita20.pdf) 
Volume 79, January–February 2020, Pages 233-246
European Journal of Mechanics - B/Fluids
https://doi.org/10.1016/j.euromechflu.2019.09.010

*/
