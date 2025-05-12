/**
# Dune in a shear flow

 

The shear flow over a dune is solved using the [hydrostatic
multilayer](/src/layered/hydro.h) approximation. 

The dune moves due to [erosion/deposition](erosion.h).
 
*/

#include "grid/multigrid1D.h"
#include "layered/hydro.h"
#include "layered/implicit.h"
#include "layered/remap.h"

/**
We add erosion/deposition to the bathymetry. */

#include "erosion.h"

/**
These are the dune size parameters. */

double eps,a;
FILE * gp, * ggp,* ffp;

int main()
{
  periodic (right);
  L0 = 20.;
  X0 = - 5;
  G = 10;
  nu = .01;
  //DT = 0.0001;
  /**
  We drive the flow with $\partial_z u = 1$ on the free-surface.
  *G=10* to impose sub critical flow 
  */
  
  dut.x = unity;

  /**
  Numerical parameters: */
  
  N = 256*2;
  nl = 50;
 
  theta_H = 1.; // more free-surface damping
   
  gp = fopen ("tmv.txt", "w");
  fclose (gp);
  ggp = fopen ("mvc.txt", "w");
  fclose (ggp);
  ffp = fopen ("finalprof.txt", "w");
  fclose (ffp);
  /**
  These are the parameters for the dune material. */
  
  
    l_s = 1/2.5, tau_s = .9, E = 1;

  /**
  The non-erodible "bedrock" is at $z_b = 0$. */
  
  z_br = zeroc;

  /**
We run several cases to illustrate the different dynamics.
   
## Check Double Deck scaling

We show here an asymptotic approximation of friction $\tau=\frac{\partial u}{\partial z}|_0$ in the case of large $Re$ in a pure shear flow.

The perturbation induced by a small bump on a shear flow  $u=U'_0 z$, is solved with the "double deck" theory.
The logitudinal scale is the same, the transverse scale is $z=Re^{-1/3}\tilde z$, in this scale the bump is
$\tilde z_b = \alpha f(x)$ with say $f(x)=1/\sqrt(\pi) e^{-x^2}$.

Then we linearise for $\alpha \rightarrow 0$ to obtain
the double deck solution
$$u = Re^{-1/3}\tilde u   = U'_0 Re^{-1/3}\tilde z + \alpha Re^{-1/3}\tilde u_1$$

$$\tilde{\tau}_{DD} = \frac{\partial \tilde u}{\partial \tilde z}|_0= U'_0  + \alpha \tilde{\tau}_1, \text{ with  }\tilde \tau_1= (3 AiryAi(0))  U'_0 TF^{-1}((-i k  U'_0)^{1/3} TF (f))$$

one has to compare $\tilde \tau_1$ versus $(\frac{\partial u}{\partial z}|_0- U'_0)/\alpha$
   
   
the physical bump is $z_b = \alpha  Re^{-1/3} f(x)$ in the code (note the $Re^{-1/3}$ which is `pow(nu,1./3))`in the code,
   for numerical test we take $\alpha=0.2$ and $Re=1/0.01$, we plot $(\tau_{layered}-U'_0)/\alpha$ compared to $\tilde \tau_1$

(see Chouly Lagrée)

   we start by a bump
   */
    
#if 0
    E=0;
    eps = .2*pow(nu,1./3), a = 1;
    gp = fopen ("tmv.txt", "a");
    fprintf (gp, " \n" );
    fclose(gp);
    run();
    
// save manually "finalprof.txt" in "DD/finalprof1.txt"
    //system(" mkdir DD;cp finalprof.txt DD/finalprof1.txt" );
#endif
/**

~~~gnuplot
p [-4:8]'DD/finalprof1.txt' u 1:($3/.2) w l t'bump 1/\sqrt(\pi) e^{-x^2}',''u 1:(($4-1)/.2) w l t'layered','../REFCASES/refDD.txt' u 1:($2) w l t'tau DD'
~~~

*/

/**
 
 same shape but longer
 */
    
#if 0
    E=0;
    eps = .2*pow(nu,1./3), a = 1/2.;
    gp = fopen ("tmv.txt", "a");
    fprintf (gp, " \n" );
    fclose(gp);
    run();
// save manually "finalprof.txt" in "DD/finalprof2.txt
    system(" mkdir DD;cp finalprof.txt DD/finalprof2.txt" );
#endif
/**
 
 We check that when the scale of a bump is reduced by a factor, the shear in increased by this factor to the power 4/3
 
 ~~~gnuplot
 p [-4:8]'DD/finalprof2.txt' u ($1*.5):(($4-1)/(.5)**(4./3)/.2) w l t'layered long','DD/finalprof1.txt' u 1:(($4-1)/.2) w l t'layered short','../REFCASES/refDD.txt' u 1:($2) w l t'tau DD'
 ~~~


## Erosion, two initial heaps with same mass give the same final dune
    
Erosion is there, we change the shape

 */

#if 0
  E = 1;
  eps = 1*pow(nu,1./3), a = .5;
  gp = fopen ("tmv.txt", "a");
  fprintf (gp, " \n" );
  fclose(gp);
  run();

  /**
  Same mass, but smaller. */

  eps = 1*pow(nu,1./3), a = sqrt(2);
  gp = fopen ("tmv.txt", "a");
  fprintf (gp, " \n" );
  fclose(gp);
  run();
    
    // save manually "finalprof.txt" in "REF/finalprof.txt
  system("mkdir REF; cp finalprof.txt REF/finalprof.txt;cp out REF/out");
 
#endif
/**

we see the evolution of the bump, as the domain is periodic we add the size of the domain to see the evolution
~~~gnuplot
set term @SVG size 1024,200
set xlabel 'x'
set ylabel 'z'
#set yrange [0:0.16]
p[:32][:0.2]'REF/out' u 1:3 index 0:20:1 w l t 't=0,5,10,15...95,100','' u 1:3 i 30:30 w l t't=150',''u ($1+20):3 index 40:60:10 w l t 't 200, 250,300 ','' u ($1+20):3 i 70 w l t 't=350'
~~~


With the same mass but half the initial length (short),

~~~gnuplot Dune evolution
p[:32][:0.2]'REF/out' u 1:3 index 70:90:1 w l t 't=0,5,10,15...95,100','' u 1:3 i 100:100 w l t't=150',''u ($1+20):3 index 110:130:10 w l t 't 200, 250,300 ','' u ($1+20):3 i 140 w l t 't=350'
~~~

the dune shape
converges toward the same solution.


~~~gnuplot Dune evolution
set ylabel 'z_b'
p[0:][:0.2]'REF/out' u ($1):3 i 140 w l t 't=350 short',''u ($1):3 i 70 w l t 't=350 long'
~~~

the dune shape
converges toward the same solution here in the refernce frame of the moving dune.



~~~gnuplot plot in the frame of the bump
set xlabel 'x-xmax'
set ylabel 'z_b'
p[-5:5]'REF/finalprof.txt' u ($1-$8):5 w l
~~~

The two dunes and the shear stress associated

~~~gnuplot The corresponding final frictions.
reset
set xlabel 'x'
set ylabel 'dudz, z_b'
set key top left
plot [0:]'REF/out' u  1:3 index 70 t ' initial long bump ', '' u 1:6 index 70 w l t ' frcition long bump ', \
'REF/out' u ($1):3 index 140 w l t ' initial short bump ',''u 1:6 index 140 w l t ' friction short bump ',
#, \
#    'out' u 1:6 index 32 w l t 'eps =  , a = 0.25'
~~~


## explore the change of saturation length
 

 the sediment law $l_s \partial_x q +q = (\tau-\tau_s)_+$
 it has a relaxaton length $l_s$, we vary it
 
*/
#if 1
    E = 1;
    eps = 1*pow(nu,1./3), a = .5;
    l_s = .5;
    gp = fopen ("tmv.txt", "a");
    fprintf (gp, " \n" );
    fclose(gp);
    run();
    
    eps = 1*pow(nu,1./3), a = .5;
    l_s = .45;
    gp = fopen ("tmv.txt", "a");
    fprintf (gp, " \n" );
    fclose(gp);
    run();

    eps = 1*pow(nu,1./3), a = .5;
    l_s = .42;
    gp = fopen ("tmv.txt", "a");
    fprintf (gp, " \n" );
    fclose(gp);
    run();
    
    eps = 1*pow(nu,1./3), a = .5;
    l_s = 0.40;
    gp = fopen ("tmv.txt", "a");
    fprintf (gp, " \n" );
    fclose(gp);
    run();
    
    eps = 1*pow(nu,1./3), a = .5;
    l_s = .38;
    gp = fopen ("tmv.txt", "a");
    fprintf (gp, " \n" );
    fclose(gp);
    run();
   /** 
    eps = 1*pow(nu,1./3), a = .5;
    l_s = .36;
    gp = fopen ("tmv.txt", "a");
    fprintf (gp, " \n" );
    fclose(gp);
    run();
    */
   system("mkdir L_S; cp finalprof.txt L_S/finalprof.txt" );
   system("cp mvc.txt L_S/mvc.txt" );
   system("cp tmv.txt L_S/tmv.txt" );
#endif
    
/**

 several different $L_s$, the smaller this length, the thinner the bump 
 
~~~gnuplot ls=.38 0.4  0.42 .45 .5
set xlabel 'x'
set ylabel 'q'
set arrow from 1,.1 to 2,.15
p[-5:5]'REF/finalprof.txt' u ($1-$8):5 w l,'L_S/finalprof.txt' u ($1-$8):($11>0?$5:NaN)  w l
~~~
    


~~~gnuplot
set xlabel 'x'
set ylabel 'q'
p[-5:5]'REF/finalprof.txt' u ($1-$8):5 w l,'L_S/finalprof.txt' u ($1-$8):($11<=.4?$5:NaN)  w l
~~~


plot of height  as function of time

~~~gnuplot
reset
set xlabel 't'
set ylabel 'm'
p[][0:]'M/tmv.txt' u 1:4 w l 
~~~



## Erosion, variation of the mass

We consider different masses, the hihgher the slower 
*/      

 
#if 0
    E = 1;
    l_s = 1/2.5;
    
    eps = 2*pow(nu,1./3), a = sqrt(2);
    gp = fopen ("tmv.txt", "a");
    fprintf (gp, " \n" );
    fclose(gp);
    run();
    
    eps = 1.75*pow(nu,1./3), a = sqrt(2);
    gp = fopen ("tmv.txt", "a");
    fprintf (gp, " \n" );
    fclose(gp);
    run();
    
    eps = 1.5*pow(nu,1./3), a = sqrt(2);
    gp = fopen ("tmv.txt", "a");
    fprintf (gp, " \n" );
    fclose(gp);
    run();
    
    eps = 1.25*pow(nu,1./3), a = sqrt(2);
    gp = fopen ("tmv.txt", "a");
    fprintf (gp, " \n" );
    fclose(gp);
    run();
    
    eps = 1.1*pow(nu,1./3), a = sqrt(2);
    gp = fopen ("tmv.txt", "a");
    fprintf (gp, " \n" );
    fclose(gp);
    run();
    
    eps = 1.*pow(nu,1./3), a = sqrt(2);
    gp = fopen ("tmv.txt", "a");
    fprintf (gp, " \n" );
    fclose(gp);
    run();
    
    eps = .95*pow(nu,1./3), a = sqrt(2);
    gp = fopen ("tmv.txt", "a");
    fprintf (gp, " \n" );
    fclose(gp);
    run();
    
    eps = .9*pow(nu,1./3), a = sqrt(2);
    gp = fopen ("tmv.txt", "a");
    fprintf (gp, " \n" );
    fclose(gp);
    run();
    
    eps = .8*pow(nu,1./3), a = sqrt(2);
    gp = fopen ("tmv.txt", "a");
    fprintf (gp, " \n" );
    fclose(gp);
    run();
    
   // eps = .75*pow(nu,1./3), a = sqrt(2);
   // gp = fopen ("tmv.txt", "a");
   //fprintf (gp, " \n" );
   // fclose(gp);
   // run();
   system("mkdir M; cp finalprof.txt M/finalprof.txt" );
   system("cp mvc.txt M/mvc.txt" );
   system("cp tmv.txt M/tmv.txt" );
#endif
    
   /**

variable mass

~~~gnuplot
set xlabel 'x'
set ylabel 'z_b'
m=(0.01)**(1./3)
p[-6:3]'REF/finalprof.txt' u ($1-$8):3 w l,'M/finalprof.txt' u ($1-$8):($11>0?$3:NaN)  w l
~~~

not loosing mass (mass enough large):

~~~gnuplot
set xlabel 'x'
set ylabel 'q'
m=(0.01)**(1./3)
p[-5:5]'REF/finalprof.txt' u ($1-$8):5 t'ref' w l,'M/finalprof.txt' u ($1-$8):($7>=0.9*m?$5:NaN) t' m variable' w l
~~~

## Velocity of the dune

We measure the velocity with the position of the max as function of time.
    We measure it also with $q_{max}/z_{bmax}$
    
As $$\frac{\partial z_b}{\partial t }+ \frac{\partial q}{\partial x }=0$$

we define $c= \frac{\partial q}{\partial z_b }$ as the dune velocity
indeed $q/z_b/(q_{max}/z_{bmax})$ is almost constant
    
    
~~~gnuplot velocity
p[0:20][0:2]'M/finalprof.txt' u 1:(($5/$3)/$12) w l,1
~~~

..... velocity as function of $m$

~~~gnuplot
reset
set xlabel 'mass'
set ylabel 'velocity'
m=(0.01)**(1./3)
 p'M/mvc.txt' u ($2/m):7 w lp,'' u ($2/m):($5/.02) w lp,0
~~~ 
 
 speed decreases with $l_$ and $m$
 
~~~gnuplot
  set ylabel 'speed'
  set xlabel 'l_s or mass '
  p'M/mvc.txt' u ($2/m):7 w lp,\
   'L_S/mvc.txt' u 6:7  w lp
~~~

Self similar curve from KL 06 valid for small dunes?

~~~gnuplot
  set ylabel 'reduced c: cm^{1/4}'
  set xlabel 'mass^{3/4}/l_s'
  p'M/mvc.txt' u ($2>.9*m?((1/$6)*($2/m)**.75):NaN):($7)*($2**.25) w lp,\
   'L_S/mvc.txt' u ($6>0?((1/$6)*($2/m)**.75):NaN):($7)*($2**.25) w lp
~~~
 
 
*/
    
}


/**

## init of program and save
We initialise the topography and the initial thickness of each layer
*h*. */

event init (i = 0)
{
  fprintf (stderr, "# eps = %g, a = %g\n", eps, a);
  printf ("# eps = %g, a = %g\n", eps, a);

  foreach() {
    zb[] = eps*a*exp(- sq(a*x))/sqrt(3.141692);
    double z = 0.;
    foreach_layer() {
      h[] = (1. - zb[])/nl;
      z += h[]/2.;
      u.x[] = z;
      z += h[]/2.;
    }
  }
  boundary (all);
}

/**
Uncomment this part if you want on-the-fly animation. */

#if 0
event gnuplot (i += 10)
{
  static FILE * fp = popen ("gnuplot", "w");
  if (i == 0)
    fprintf (fp, "set term x11\n"
	     // "set size ratio -1\n"
	     );
  fprintf (fp,
           "set title 'nl=%d, t=%f'\n"
           "set xl 'x'\nset yl 'h'\n"
           "plot [-5:]'-' u 1:2 w l t 'eta', '-' u 1:3 w l t 'zb', '-' u 1:4 w lp t 'tau', "
	   "'-' u 1:($5 - $6) w l t 'q - q_e'\n",
	   nl, t); 
  foreach()
  fprintf (fp, "%g %g %g %g %g %g\n", x, eta[], zb[], frott[], q[], q_e[]);
  fprintf (fp, "e\n");
  fflush (fp);
  //  usleep (10000);
}
#endif


/**
## Outputs

We output the profiles at regular intervals. */

event profiles (t += 5; t <= 350)
{
  foreach (serial)
    printf ("%g %g %g %g %g %g %g\n",
	    x, eta[], zb[], q[], q_e[], dudz(u), t);
  printf ("\n\n");
  fflush (stdout);
    
    
    gp = fopen ("tmv.txt", "a");
    fprintf (gp, "%g %g %g %g %g %g %g\n",t,m,xmax,zbmax,cdune,l_s,qmax/zbmax);
    fclose(gp);
    
     cdune=(xmax-xmax0)/5;
     xmax0=xmax;
}

event end (t = 350)
{
    ggp = fopen ("mvc.txt", "a");
    fprintf (ggp, "%g %g %g %g %g %g %g \n",t,m,xmax,zbmax,cdune,l_s,qmax/zbmax);
    fclose(ggp);
    
    ffp = fopen ("finalprof.txt", "a");
    foreach()
    fprintf (ffp, "%g %g %g %g %g %g %g %g %g %g %g %g \n", x, eta[], zb[], frott[], q[], q_e[],m,xmax,zbmax,cdune,l_s,qmax/zbmax);
    fprintf (ffp, "\n");
}


/**



 

## compilation 


~~~bash
 
make dune13.tst
source ../../c2html.sh dune13
~~~
 
 
 
 
 







## References

~~~bib
@article{kouakou2006,
title={Evolution of a model dune in a shear flow},
author={Kouakou, Kouam{\'e} Kan Jacques and Lagr{\'e}e, Pierre-Yves},
journal={European Journal of Mechanics-B/Fluids},
volume={25},
number={3},
pages={348--359},
year={2006},
publisher={Elsevier},
pdf = {http://www.lmm.jussieu.fr/~lagree/TEXTES/PDF/kouakoulagree06.pdf}
}
~~~

 * Chouly Lagrée
 [http://www.lmm.jussieu.fr/~lagree/TEXTES/PDF/choulylagree12.pdf]()
 
 */
