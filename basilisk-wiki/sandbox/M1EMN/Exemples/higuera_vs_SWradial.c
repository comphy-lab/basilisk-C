/**
# Viscous hydraulic Jump

We want to reproduce the circular hydraulic jump of Higuera (1997) and compare with Saint Venant (one layer!). */

#include "grid/cartesian1D.h"
#include "radial.h"
#include "saint-venant.h"
double Cf;


int main() {
  X0 = 0.14;
  L0 = 1. - X0;
  G  = 1;
  N  = 128*2;
  nl = 1;
 // Cf = 0;
 // run();

  Cf = 3;
  run();

  Cf = 2.2799;
  run();


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


//h[right] = dirichlet (t>.25? 1.46:0);
//eta[right] = dirichlet (t>.25? 1.46:0);

/**
## Initialization

We define a scalar field *hc* to check for convergence on $h$. */

scalar hc[];

event init (i = 0) {

  /**
  We set a constant velocity at the inlet and a free outlet. */

  for (vector u in ul) {
    u.n[left] = 5.9;
    u.n[right] = neumann(0.);
  }
  
  /**
  We initialize *h* and *hc*. */
  
  foreach()
    hc[] = h[] = 0.2;
}

event frotte (i++) {

  /**
 friction
  */
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

We print the elevation and the stress. */

event output (t = end) {
  vector u0 = ul[0];
  double tau ;
  foreach(){
     if(nl==1){tau = Cf*u.x[]/h[];}else{tau = 2.*u0.x[]/(h[]/nl);}
    fprintf (stderr, "%g %g %g \n", x, eta[], tau);
  }
   fprintf (stderr, "\n\n"); 
}


event profiles (i += 50) {
  static FILE * fp = popen ("gnuplot", "w");
  if (i == 0)
    fprintf (fp,
       "set term x11 noraise\n"
       "set grid\n"
       "set xrange [%g:%g]\n", 0.0, L0);
  fprintf (fp,
     "set title 't = %.2f'\n"
     "p '-' u 1:3:2 w filledcu lc 3 t ''\n", t);
  foreach()
      fprintf (fp, "%g %g %g\n", x, eta[], zb[]);
  fprintf (fp, "e\n\n");
  fflush (fp);
}


/**

## Run

compile and run

~~~bash
 
qcc higuera_vs_SWradial.c -o higuera_vs_SWradial -lm
./higuera_vs_SWradial 2> log


make higuera_vs_SWradial.tst
make higuera_vs_SWradial/plots
make higuera_vs_SWradial.c.html

~~~


## Results

~~~gnuplot raw plot.
 
plot[][:2]  'log' u 1:2 t 'h' w l,'' u 1:3 t 'tau' w l,0 not w l linec -1

~~~


  Comparison with Figure 2 of Higuera (1997) to be done...
 
 
 
 
  
 
~~~gnuplot Comparison with  one layer (two different closures) and multilayer

reset
H=0.339092
set ylabel "h"
set y2label "{/Symbol t}_b"
set y2range [-1:30]
set y2tics border out 

p [0:][-.1:3] 0 not w l linec -1,1.813*(x<.15? x:NaN) not linec -1 axes x1y1,(H**4-12/2/pi*log(x))**.25  not linec -1 axes x1y1,\
            'log' i 0 u 1:2 t'P' w l ls 1 axes x1y1,'' i 0 u 1:($3) not w l ls 1 axes x1y2,\
                 '' i 1 u 1:2 t'W' w l ls 2 axes x1y1,'' i 1 u 1:($3) not w l ls 2 axes x1y2,\
               '' i 2 u 1:2 t'ML' w l ls 3 axes x1y1,'' i 2 u 1:($3) not w l ls 3 axes x1y2

~~~

 

## Links
  
  * [http://basilisk.fr/sandbox/M1EMN/Exemples/higuera.c]() same in 2D

  * [http://basilisk.fr/src/examples/swasi.c]() reverse jump (injection at the external radius)
   

## Bibliography

* Higuera, F. 1994. [The hydraulic jump in a viscous laminar
flow](http://dx.doi.org/10.1017/S0022112094002041). J. Fluid
Mech. 274, 69-92.

* J. Watson The radial spread of a liquid jet over a horizontal plane
 J.Fluid Mech. (1964), uol. 20, part 3,pp. 481-499

* F. J. Higuera The circular hydraulic jump, 
Phys. Fluids 9 (5), May 1997 1070-6631/97/9(5)/1476/3


*/
