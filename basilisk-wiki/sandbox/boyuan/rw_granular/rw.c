/**

**Inspired by [the sandbox example](http://basilisk.fr/sandbox/M1EMN/Exemples/front_poul_ed.c).**

# Periodic and Non-periodic granular roll waves
Reproduction of the case in Gray & Edwards (2014),  Forterre & Pouliquen (2003).
 
## Problem
Frontal roll wave and trailing waves for $Fr>Fr_c=2/3$.
 
## Equations
Simulation in 1D of the fundamental experiments by Pouliquen 99 revisited by the longitudinal viscosity by Edwards & Gray 14

We solve the Savage-Hutter-Depth-Averaged-Shallow-Water-Saint-Venant equations in 1D
 
$$\frac{\partial }{\partial t}  h \; +\; \frac{\partial }{\partial x} uh=0$$
$$ \frac{\partial }{\partial t} hu +  
 \frac{\partial }{\partial x}  \dfrac{(hu)^2}{h} +
 \frac{\partial }{\partial x}g\dfrac{h^2}{2}
= - gh \frac{\partial }{\partial x} Z-\mu(I) g h\frac{u}{|u|}
+ \frac{\partial }{\partial x}(\nu_e h^{3/2}\frac{\partial }{\partial x}u)
$$
 over a rough inclined plane $Z=-\tan \theta x$.
 
Basal  friction law  based on $\mu(I)$ rheology ($\mu$):
 $$\mu=\mu_0 + \frac{ \Delta \mu}{1+I_0/I}$$
More directly:
 $$\mu=\mu_1+\frac{\mu_2-\mu_1}{\frac{\beta h}{\mathcal{L} Fr}+1}$$
where $I=\frac{5 (d_g)|u|}{2\sqrt{gh}h}$ is the mean inertial/ Froude / Savage number, where $d_g$ is the grain diameter.
Note that the problem is without dimension, if space is measured by $H$, the speeds are measured by 
$\sqrt{g H}$  (see remark on t). So, there is no $g$ and $\rho$ in $I$.
Note that in the fundamental papers of Pouliquen, the function $\mu$ is different (and the value of the parameters of the friction law as well). 

$\mu_1=\tan(\zeta_1)$, $\mu_2=\tan(\zeta_2)$: two critical angles proposed by Pouliquen & Forterre 2002. $\beta$ is a dimensionless constant, $\mathcal{L}$ is a constant with dimension of length.

An extra viscous term has been added by Edwards and Gray corresponding to a longitudinal viscosity
depending on the height. This new viscous term $\nu_e h^{3/2}$ will be explained later. (thickness three halves power viscosity term)

## Code
*/ 
#include "grid/multigrid1D.h"
#include "saint-venant.h" 
/**The problem is solved in one dimension, but can be extended to two. 

Problem-specific parameters consistent with Table 1 and Table 2 of Gray & Edwards (2014).

Perturbation frequency follow Edwards & Gray(2015).

*Note that the parameters are in SI units.*
*/
double zeta1 = 0.36477;
double zeta2 = 0.5718;
double zeta = 0.506145;
double ell = 0.000825;
double beta = 0.136;
double mu = 0.0; //friction-coeff. depending on h and Fr
double u0, h0; //stead-uniform flow, to be calculated later
double pertFreq = 1.50;
double pertAmp = 0.10;
double piConst = 3.1415926536;
double nue;
double tantheta; //define this variable maybe more convient
double simTime;

/**
Main and parameters
*/
int main()
{
/**
  The domain is 10 long (follow Edwards & Gray(2015)'s granular wave sim), 
  
  tantheta is $\tan \theta$, 
  
  Initial shape and velocity */
  X0 = 0.;
  L0 = 120.0;
  /*steady-uniform flow depth used in  Forterre & Pouliquen (2003) experiment
  */
  h0 = 0.0031989;

  G = 9.81;

  N = 1024*4*2*4*2;
  simTime = 104.0;
  //use CFL number for $dt$ control
  CFL = 0.25;
  theta = 1.50; // use Sweby to slightly increase resolution

/**
   Bagnold velocity and depth for steady-uniform flow
*/   

  u0 = beta/ell*((tan(zeta)-tan(zeta1))/(tan(zeta2)-tan(zeta)))*(sqrt(G*cos(zeta)))*(pow(h0, 3.0/2.0));
  // u0 = 0.16898;
/**
Edwards and Gray 14 proposed a viscosity comming from the $\partial_x \tau_{xx}$ term, 
through integration, this gives 
$$ \frac{\partial }{\partial x}(\nu_e h^{3/2}\frac{\partial }{\partial x}u)$$
where 
$\nu_e=5 d sin \theta /(9 I_\theta \sqrt{cos \theta}) $
or
$$\nu_e=\frac{2}{9}\frac{\mathcal{L} \sqrt{g}}{\beta}\frac{\sin \zeta}{\sqrt{\cos \zeta}}(\frac{\tan \zeta_2-\tan \zeta}{\tan \zeta-\tan \zeta_1})$$
that we approximate as 
*/
  nue =  (2.0/9.0)*(ell*sqrt(G)/beta)*(sin(zeta)/sqrt(cos(zeta)))*((tan(zeta2)-tan(zeta))/(tan(zeta)-tan(zeta1)));
  fprintf (stderr,"u0=%lf  h0=%lf nu=%lf\n",u0,h0,nue);  

  tantheta = tan(zeta);
/** 
If viscosity $\nu_e$  is not zero, the maximal time step is defined due to this viscosity:

*$dt$ will be updated later??*
*/
 DT = (L0/N)*(L0/N)/2/nue/2.0;
  run();
}

/**
The initial conditions */
event init (i = 0){

/**
BCs: since the case is supercritical, we use a Dirichlet BC for inlet and Neumann-0 Bc for outlet.
*/
    u.n[left] = dirichlet(u0);
    // u.n[left] = neumann(0.0);
    h[left] =  dirichlet(h0*(1.0+pertAmp*(sin(2.0*piConst*pertFreq*t))));
    u.n[right] =  neumann(0.0);
    h[right] = neumann(0.0);//h[];
/**
the inclined plane is the topography
*/   
   foreach(){
    // zb[] =  -tantheta*x*cos(zeta);
    zb[] =  0.0;
    }
/**
Initial condition: steady-uniform flow
*/
   foreach(){
    h[]= h0;
    u.x[]= u0;
    }
    boundary ({u.x,h});
}


/** bed-slope and bed-friction: use a simple forward Euler scheme for now
*/
event coulomb_friction (i++) {  
  double Froude, ff;
/**
We use a simple implicit scheme to implement coulomb bottom friction i.e.
(note the simplification by $h$) 
$$\frac{d\mathbf{u}}{dt} = -\mu g \frac{\mathbf{u}}{|\mathbf{u}|}$$
  with $\mu$ fonction of $I$. 
  Note that the good implementation preserving equilibrium balance is in Bouchut's book
*/  
  foreach() {
    Froude = u.x[]/sqrt(G*cos(zeta)*h[]); 
    // ff = norm(u) > 0 ? (1. +  dt *(tan(zeta1)+(tan(zeta2)-tan(zeta1))/(beta*h[]/(ell*Froude)+1.0))*G*cos(zeta)) : HUGE ;

  foreach_dimension()
      // u.x[] /= ff;
      u.x[] = u.x[]+dt *(-1.0)*(tan(zeta1)+(tan(zeta2)-tan(zeta1))/(beta*h[]/(ell*Froude)+1.0))*G*cos(zeta)+dt*(tan(zeta)*G*cos(zeta));
  }
  boundary ({u.x,h});
}
/**
 The new viscous term from Gray & Edwards
$$ \frac{\partial }{\partial x}(\nu_e h^{3/2}\frac{\partial }{\partial x}u)$$
*/ 
event friction_long (i++) {
 foreach_dimension() {
    face vector g[];
    scalar a = u.x;
    foreach_face()
      g.x[] = nue*(a[] - a[-1,0])/Delta*pow((h[0,0] + h[-1,0])/2.,3./2);
    foreach ()
      u.x[] += dt/Delta*(g.x[1,0] - g.x[]);
      // u.x[] += dt/Delta*(g.x[1,0] - g.x[] + g.y[0,1] - g.y[]);
  }
  boundary ((scalar *){u});
}

/**
# Output Control 
*/

void plot_profile(double t, FILE *fp)
{
     fprintf(fp,
             "set term pngcairo enhanced size 800,600 font \",10\"\n"
             "set output 't%.0f.png'\n"
             "set title 't = %.2f'\n"
             "set xlabel 'x(m)'\n"
             "set ylabel 'eta(m)'\n"
             "plot [0:][0:]'-' u 1:2 w l lw 2\n",
             t, t);
     foreach ()
     {
          fprintf(fp, "%g %g\n", x, h[]);
     }
     fprintf(fp, "e\n\n");
     fflush(fp);
}

event gnuplot(t = 0; t <= simTime; t += 2.0)
{
     static FILE *fp = popen("gnuplot 2> /dev/null", "w");
     plot_profile(t, fp);
     // fprintf(fp,
     //         "set term pngcairo enhanced size 800,600 font \",10\"\n"
     //         "set output 't%.0f.png'\n"
     //         "set title 't = %.2f'\n"
     //         "set xrange [0:40]\n"
     //         "plot u 1:2 w l t\n",
     //         t, t);
     // fprintf(fp, "\n");
     // foreach ()
     //      fprintf(fp, "%g %g\n", x, h[]);
     // fprintf(fp, "e\n\n");
     // fflush(fp);
     // fprintf(stderr, "%.3f %.3f\n", t, statsf(h).max); // uncomment if needed
}

event output(t = 0; t <= simTime; t += 4.0)
{
     char name[80];
     sprintf(name, "out-%.1f", t);
     FILE *fp = fopen(name, "w");
     foreach ()
          fprintf(fp, "%g %g %g \n", x, h[], u.x[]);
     fprintf(fp, "\n");
     fclose(fp);
     fprintf(stderr, "%g\n", dt);
}

event outputMax(t = 0; t <= simTime; t += 0.1)
{
  double maxDepth = 0.0;
  double maxDepthLocX = 0.0;
  double maxDepthVel = 0.0;

  FILE *fp2 = fopen("maxDepth", "a+");


    foreach ()
  {
    if (h[]>maxDepth){
        maxDepth = h[];
        maxDepthLocX = x;
        maxDepthVel = u.x[];
    }
  }
  
  fprintf(fp2, "%g %g %g %g\n", t, maxDepthLocX, maxDepth, maxDepthVel);

  fclose(fp2);
}

/**
# Run

To run and plot with gnuplot through a pipe 

~~~bash
qcc -g -O2 -DTRASH=1 -Wall -DgnuX=1 -o front_poul_ed front_poul_ed.c -lm
./front_poul_ed | gnuplot 
~~~

To run and make a film (a gif animated)

~~~bash
 qcc -g -O2 -DTRASH=1 -Wall -DgnuX=1 -o front_poul_ed front_poul_ed.c -lm
  ./front_poul_ed > v.out
 echo "set term gif animate; set output 'front_poul_ed.gif'" > dmp
 cat v.out >> dmp
 mv dmp v.out 
 cat v.out | gnuplot
~~~

~~~bash 
  gifsicle --colors 256 --optimize --delay 1   --loopcount=0 front_poul_ed.gif > front_poul_edo.gif 
~~~ 

To run like on the web

~~~bash
qcc -g -O2 -DTRASH=1 -Wall   -o front_poul_ed front_poul_ed.c -lm
./front_poul_ed > out 
~~~

 
## Links
 
 * savagestaron.c
 
## Bibliography

 * Guillaume Saingier, Stéphanie Deboeuf, and Pierre-Yves Lagrée
 "On the front shape of an inertial granular flow down a rough incline"

 * O. Pouliquen
["On the shape of granular fronts down rough inclined planes"](http://iusti.polytech.univ-mrs.fr/~pouliquen/publiperso/PoFfront99.pdf)
 Phys of Fluids Vol 11, Nbr 7 July 1999

 * O. Pouliquen 
 ["Scaling laws in granular flows down rough inclined planes"](http://iusti.polytech.univ-mrs.fr/~pouliquen/publiperso/PoFscaling99.pdf) Phys.  Fluids, 11, 542-548 (1999) .

 * A. Edwards N. Gray
"A depth average $\mu(I)$-rheology for shallow granular free-surface flows"
JFM 2014

 * F. Bouchut
4.12 Coulomb Friction
p 97 


*/

 
