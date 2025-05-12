/**
#  Periodic free surface flow of a Plastic fluid
*/

#include "navier-stokes/centered.h"
//Borrowed from acastillo's sandbox
#include "output_vtu_foreach.h"

#define LEVEL 5
double mumax,dg,Co=0.125,angle=0.6,mus=0.38,dmu=0.26,I0=0.3,lambdaR=1e-3,oldErreur=HUGE;
long int itmax=1000000;
scalar mu_eq[],foo[];

double error_rms_vit(void);


/**
now, come back to the velocity shear which  is
  $$\frac{\partial u}{\partial y }= K_B \sqrt{ (h-y)}
 - t_c \frac{1}{\sqrt{ (h-y)}}$$
 with  $\frac{\partial u}{\partial y }=0$  for $y>Y$ with $Y= h -\frac{t_c}{K}$
 $K_B =(\frac{(\tan \alpha    - \mu_0)\sqrt{g \cos \alpha}}{a d })$ and $t_c=
 - \frac{\tau_c }{\rho g h}  ( \frac{\sqrt{g }}{d})\frac{1}{a \sqrt{\cos \alpha}}$
 */
double dUbc( double y){
    double dU;
    double a = dmu/I0;
    double h = 1.;
    double KB = ((tan(angle) - mus)* sqrt(cos(angle))/(a*dg));
    double tc = Co/(sqrt(cos(angle))*(a*dg));
    double Y = h -tc/KB;
    dU=(y<Y ? KB*sqrt(h - y)-tc/sqrt(h - y):0);
    return dU;
}
/**
This gives the Bagnold's solution if $\tau_c=0$
$$ u = \frac{2}{3} \frac{I_\alpha}{d_g} \sqrt{g h^3 \cos(\alpha)}   (1 - (1 - \frac{y}{h})^{3/2})$$
 (note the  numerical value  for alpha=0.43, dg=0.04 U0 = 2.06631)
*/
double Ub( double y){
    double U0=(sqrt(cos(angle))*(-I0*mus + I0*tan(angle))/(dg*(mus+dmu - tan(angle))));
    return ((1. - pow(1. - y, 1.5))*2./3.*U0) ;
}
/**
 This gives the Cohesive Bagnold's solution if $\tau_c \ne0$
 $$ u =  \frac{(2(h^{3/2})K_B-3\sqrt{h} t_c-h K_B \sqrt{(h-y)}+3 t_c\sqrt{h-y}+K_B \sqrt{h-y}y)))}{3}$$
 
 $$U =\frac{2}{3} ((h^{3/2} K_B - 3 t_c \sqrt(h)+ 2 t_c \sqrt{(t_c/K_B)}))$$
 */
double Ubc( double y){
    double a = dmu/I0;
    double h = 1;
    double KB = ((tan(angle) - mus)* sqrt(cos(angle))/(a*dg));
    double tc = Co/(sqrt(cos(angle))*(a*dg));
    double Y = h -tc/KB;
    double U ,Uy;
    Uy=(2*(pow(h,1.5)*KB-3*sqrt(h)*tc-h*KB*sqrt(h-y)+3*tc*sqrt(h-y)+KB*sqrt(h-y)*y))/3.;
    U = 2*((pow(h,1.5)*KB - 3*tc *sqrt(h)+ 2*tc*sqrt(tc/KB)))/3.;
    U = (y<Y ? Uy : U);
    return (U) ;
}

int main() {
  L0 = 1.;
  origin (0., 0);
  dg = 0.04;
  N=1<<5;

  mumax=1000;

  periodic (right);
  u.t[top] = neumann(0);
  u.n[top] = dirichlet(0);
  u.n[bottom] = dirichlet(0);
  u.t[bottom] = dirichlet(0);
  p[top] = dirichlet(0);

  DT = 0.001;
  //stokes = true;
  
  run(); 
}

face vector muv[];

event init (t = 0) {

  mu = muv;
/**
  pressure gradient, gravity acceleartion `mdpdx`
 $$-\frac{\partial p}{\partial x} = \sin(alpha) $$
 $$-\frac{\partial p}{\partial y} = -\cos(alpha) $$
*/
  const face vector mdpdx[] = {sin(angle),-cos(angle)};

  a = mdpdx;

  foreach() {
    u.x[] = Ub(y);
    u.y[] = 0;
    p[]=cos(angle)*(1-y);
  }
}

/**
We check the number of iterations of the Poisson and viscous
problems. */
event logfile (i++) 
  fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);



scalar un[];
event init_un (i = 0) {
    foreach()
    un[] = u.x[];
}

event conv (t += 1; i <= itmax) {
  double du = change (u.x, un);
  double err=error_rms_vit(),errDu=fabs(err-oldErreur)/oldErreur;
  oldErreur=err;

  fprintf(stdout,"%g %g %g %g %g %g %g\n",t,interpolate (u.x, L0/2, .999)/Ubc(1),interpolate (u.x, L0/2, .999),Ubc(1),du,err,errDu);
  fflush(stdout);
  if (i > 0 && du < 5e-5) {//was -6
    scalar solEx[];
    foreach() {
      solEx[]=Ubc(y);
    }
    
    char name[80];
    sprintf(name, "./sol_at_%g",t);
    output_vtu((scalar *) {p,solEx}, (vector *) {u}, name);
  
    return 1; /* stop */
  }
}

/**
Not parallel yet (with MPI) 
*/
double error_rms_vit(void) {
  double rms=0;
  int cpt=0;
  
  foreach(reduction(+:rms) reduction(+:cpt)) {
    double val,valEx;
    val=u.x[];
    valEx=Ubc(y);
    rms+=(val-valEx)*(val-valEx);
    //printf("%g %d %g %g\n",t,cpt,val,valEx);
    cpt++;
  }
  return sqrt(rms)/cpt;
}


/**
## Implementation of the Bagnold cohesive viscosity
*/
#define regul 1
event bagnold(i++) {
  foreach() {
    #if regul
      double valRegul;
    #endif
    mu_eq[] = mumax;
    if (p[] > 0.) {
      double D2 = 0.,In,muI;
      foreach_dimension() {
        double dxx = u.x[1,0] - u.x[-1,0];
        double dxy = (u.x[0,1] - u.x[0,-1] + u.y[1,0] - u.y[-1,0])/2.;
        D2 += sq(dxx) + sq(dxy);
      }
      #if regul
        D2 = sqrt(2.)*sqrt(D2)/(2.*Delta) ;
	valRegul=sqrt(D2*D2+lambdaR*lambdaR);
	In = dg*D2/sqrt(fabs(p[]));
	muI = mus + (dmu)*In/(I0 + In);
	muI = mus + (dmu)*In/(I0 );
	foo[] = D2;
	mu_eq[] =  (Co+muI*p[])/valRegul;
      #else
        if (D2 > 0.) {
           D2 = sqrt(2.)*sqrt(D2)/(2.*Delta) ;
           In = dg*D2/sqrt(fabs(p[]));
           muI = mus + dmu*In/(I0 + In);
           muI = mus + dmu*In/(I0 );
           foo[] = D2;
           mu_eq[] =  min( (Co + muI*fabs(p[]))/D2 , mumax );
        } else{
           mu_eq[] =  mumax;
        }
      #endif
    }
  }
  boundary ({mu_eq});
  foreach_face() {
    muv.x[] = (mu_eq[] + mu_eq[-1,0])/2.;
  }
  boundary ((scalar *){muv});
}
/**
  Save profiles
*/
event profiles (t += 1)
{
    FILE * fp = fopen("xprof", "w");
    scalar shear[];
    foreach()
    shear[] = (u.x[0,1] - u.x[0,-1])/(2.*Delta);
    boundary ({shear});
    for (double y = 0.; y < 1.0; y += 1./pow(2.,LEVEL))
        fprintf (fp, "%g %g %g %g %g %g %g %g %g \n", y, interpolate (u.x, L0/2, y), interpolate (shear, L0/2, y),
                  Ubc(y), cos(angle)*(1-y),interpolate (p, L0/2, y),
                 interpolate (mu_eq, L0/2, y), interpolate (foo, L0/2, y),dUbc(y) );
    fclose (fp);
}

/**
We adapt according to the error on the velocity field. 
*/
event adapt (i++) {
 // adapt_wavelet ({u}, (double[]){3e-3,3e-3}, 8, 6);
}
/**
## Results and plots
 
To run the program
 
~~~bash
 qcc -g -O3 -o bagnold_periodic_cohesif bagnold_periodic_cohesif.c -lm
 ./bagnold_periodic_cohesif
 
 lldb bagnold_periodic
 
 make bagnold_periodic_cohesif.tst;make bagnold_periodic_cohesif/plots;
 make bagnold_periodic_cohesif.c.html;open bagnold_periodic_cohesif.c.html;
 
~~~

 bug from MacOS
~~~bash
 cp bagnold_periodic_cohesif.c.html XXX.c.html
 
 sed -i - 's/\\)//g' XXX.c.html;sed -i - 's/\\(//g' XXX.c.html
 sed -i - 's/\\\[//g' XXX.c.html
 sed -i - 's/\\\]//g' XXX.c.html;
 mv XXX.c.html bagnold_periodic_cohesif.c.html
 open bagnold_periodic_cohesif.c.html
 
~~~

Plots of the velocity,  $\tau$ and $p$, the two last are linear as expected:

 
~~~gnuplot Velocity, pressure and $\tau= \mu_{eq} du/dy$ profiles computed
 set xlabel "y"
 set xlabel "U, tau, p"
 p'xprof' u 1:2 t'U computed' ,''u 1:($7*$3) t'tau comp.','' u 1:6 t'p'
~~~

We check the pressure $p$ is $cos(alpha)*(1-y)$
 
 
~~~gnuplot pressure profiles compared to the lithostatic for Cohesive Bagnold flow
 set xlabel "y"
 set xlabel "p"
 p[][0:]'xprof' u 1:6 t'Pression',''u 1:($5) t'cos(alpha)*(1-y)' w l
~~~
 
 we verify that the computed velocity is the Bagnold one.
 
 
 
~~~gnuplot Velocity and stress profiles computed and exact for Bagnold flow
 set xlabel "y"
 p[][0:2.5]'xprof' u 1:($2*1.) t'U computed',''u 1:4 t'U Cohesive Bagnold exact' w l linec 1,'' u 1:($6) t'pressure' w l
~~~

We check that $(\partial u/ \partial y)$ is $\sqrt{2}D_2$
 and that the computed gradient of velocity is the Bagnold one.
 
 
~~~gnuplot shear velocity profiles for Cohesive Bagnold flow
 set xlabel "y"
  p[][0:]'xprof' u 1:3 t'Computed dUbc/dy',''u 1:($8*sqrt(2))t'Sqrt(2)*D2' ,''u 1:($9) t' d Ubc/dy'w l linec 1
~~~
 
 
## Related examples, links
 

 * [The basic Bagnold flow](http://basilisk.fr/sandbox/M1EMN/Exemples/bagnold_periodic.c)
 
 * [http://basilisk.fr/sandbox/M1EMN/Exemples/bagnold_periodic_cohesif.c]()
 
 * [related example in Basilisk for pure Bagnold](http://basilisk.fr/sandbox/M1EMN/Exemples/bagnold_periodic.c)
 
 * [related example in Gerris](http://gerris.dalembert.upmc.fr/gerris/tests/tests/poiseuille.html#bagnold)
 
 
 
## Bibliography
 
 * Lagrée, Staron, Popinet
 ["The granular column collapse as a
 continuum: validity of a two–dimensional
 Navier–Stokes model with a μ(I)-rheology"](http://www.lmm.jussieu.fr/~lagree/TEXTES/PDF/JFMcollapsePYLLSSP11.pdf) J. Fluid Mech. 2011
 
 * Pierre Jop, Yoël Forterre & Olivier Pouliquen
 "A constitutive law for dense granular flows", Vol 441 8 June 2006 doi:10.1038/nature04801


 

 
 Paris  dec 2019
 
*/

