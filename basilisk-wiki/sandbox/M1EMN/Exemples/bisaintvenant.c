/**
# Two Saint-Venant all together : a "Bi Saint-Venant" model

## Scope
 
This is the problem of the wave generation by a granular heap collapsing in a layer of water at rest.
 
 

## Model equations
 
 Two fluids interact, a standard fluid of depth $h$, and of velocity $u=q/h$, and a granular fluid
 of depth $h_s$, and of velocity $u_s=q_s/h_s$,
 we use depth averaged Saint Venant Equations written in a conservative way:
 $$
 \frac{\partial }{\partial t}\left(\begin{array}{c}
h\\
h_s\\
q\\
q_s\\
\end{array}\right) +
 \frac{\partial }{\partial x}\left(\begin{array}{c}
 q  \\
 q_s\\
 \frac{q^2}{h} + \frac{gh^2}{2}  \\
\frac{q_s^2}{h} + \frac{gh_s^2}{2}
 \end{array}\right) =
 \left(\begin{array}{c}
 0\\
 0\\
 \alpha h_m ( \frac{q_s}{h_q} -  \frac{q}{h}) -  C_{\nu}\frac{q}{h^2})\\
 - \alpha h_m \frac{\rho_f}{\rho_s}( \frac{q_s}{h_q} -  \frac{q}{h}) - \mu g h
  \end{array}\right)
 $$
 with $g$ the acceleration of gravity .
 
 An inter exchange force proportional to the difference of velocity is in source term:
 $\pm \alpha ( u_s -  u)$, as we integrate over the layers, the force plays a role only for the smaller height of fluid $min(h,h_s)$,
 that is why we define $h_m=min(h,h_s)$.
 
 The friction is a laminar one for the fluid, and a Coulomb friction for the solid.
 
This is a first test, many other sophistications are possible...

 
We use the generic solver [http://basilisk.fr/src/conservation.h]() to solve the system.
 The system has
 4 eigen values
 $u \pm \sqrt{gh}$ and  $u_s \pm \sqrt{gh_s}$ (no coupling in the conservative part here, density are supposed the same).
 
## Code
Code is strongly inspired from example `artery1D` of `conservation.h`.
*/
#include "grid/cartesian1D.h"
#include "predictor-corrector.h"

static double momentum_source (scalar * current, scalar * updates, double dtmax);

event defaults (i = 0)
update = momentum_source;

#include "conservation.h"

/**
## Variables

We define the conserved scalar fields $h$ (fluid) $h_s$ (second fluid), $q$ and $q_s$ which are passed to the
 generic solver through the `scalars` list.
 */

scalar h[],q[],hs[],qs[],u[],us[];
scalar * scalars = {h,hs,q,qs};
vector * vectors = NULL;

double g,alpha,Cnu,tmax;

/**
## Flux Functions
 
 We define the `flux` function required by the [generic
 solver](/src/conservation.h). */

void flux (const double * s, double * f, double e[2])
{
    double h = s[0], hs = s[1], q = s[2], qs = s[3];
    double   u = q/h, us = qs/hs;
    f[0] = q;
    f[1] = qs ;
    f[2] = q*q/h    + g*h*h/2;
    f[3] = qs*qs/hs + g*hs*hs/2;
    // min/max eigenvalues
    double c = sqrt(g*h), cs = sqrt(g*hs);
    e[0] =  min(u-c,us-cs);
    e[1] =  max(u+c,us+cs);
}

/**
 We need to add the source term of the momentum equation. We define a
 function which, given the current states, fills the *updates* with the
 source terms for each conserved quantity.



note that the friction is splitted and computed after. */

static double momentum_source (scalar * current, scalar * updates, double dtmax)
{
    /**
     We first compute the updates from the system of conservation laws. */
    double dt = update_conservation (current, updates, dtmax);
    
    /**
     We recover the current fields and their variations from the lists... */
    
    scalar h = current[0], hs = current[1], q = current[2], qs = current[3];
    scalar dq  = updates[2] ;
    scalar dqs = updates[3] ;
    
    /**
     We add the source term for `q` and next for `qs`:
     it is an interforce proportional to the difference of velocities.
     Note $\frac{\rho_f}{\rho_s}=1/2.5$
     */
    double Finter;
    foreach(){
      u[] = q[]/h[];
      us[]=qs[]/hs[];
      Finter=alpha*(u[]-us[]);
      dq[] += - min(h[],hs[])*Finter ;
      dqs[]+=   min(h[],hs[])*Finter/2.5 ;
    }
    return dt;
}
/**
## Boundary conditions
No penetration at the left. */
q[left] = dirichlet(0);
qs[left] = dirichlet(0);
h[left] = neumann(0);
hs[left] = neumann(0);
u[left] = dirichlet(0);
us[left] = dirichlet(0);

/**
## Parameters
 Initial values, grid, parameters, etc. */

int main() {
    init_grid (256*2);
    L0 = 40;
    X0 = 0;
    g = 1 ;
    alpha = 1;
    Cnu = 0.0;
    tmax = 100;
    run();
}

/**
## Initial conditions
 
The initial conditions, a heap of granular media and water at rest. */
event init (i = 0) {
   // theta = 1.3; // tune limiting from the default minmod
    foreach(){
        h[] = 1. ;
        q[] = 0;
        hs[] = 2*(x<8) +.001;}
}
/**
## Frictions
 Viscous Friction for the fluid
*/
event viscous_frcition(i++){
    foreach()
    q[]=q[]/(1+Cnu*dt/h[]/h[]);
}
/**
Coulomb Friction for granular flow
*/ 
event coulomb_friction (i++) {
 double In,nbrg=32.,mu,Qs;
 foreach() {
   Qs=fabs(qs[]);
   In=(1./nbrg)*5./2.*Qs/pow(g*hs[],0.5);
   mu = (0.4 + 0.28*In/(0.4+In)) ;
   // mu = (0.4 + 0.26*exp(-0.136/In));
   if(Qs>0){
     qs[] = max(Qs -dt * mu * g * hs[],0)*qs[]/Qs;}
  }
 boundary ({qs});
 }

/**
## Outputs
We print to standard error the depths, we may plot with `gnuplot` as well. */
event printdata (t += 1; t <= tmax) {
//event printdata (i++,t<=2){
    foreach()
    fprintf (stderr, "%g %.6f %g  %g \n", x, hs[], h[],t);
    fprintf (stderr, "\n");
}

#ifdef gnuX
event plot (t<tmax;t+=0.05) {
    printf("set title '  --- t= %.4lf '\n"
           "p[%g:%g][-.25:2]  '-' u 1:($2) t'free surface' w l lt 3,"
           "'' u 1:3 t'hs' w l lt 4\n" ,
           t,X0,X0+L0);
    foreach()
    printf ("%g %g %g %g\n", x, h[], hs[], t);
    printf ("e\n\n");
}
#endif
/**
 
## compilation
 
~~~bash
 qcc -DgnuX=1 bisaintvenant.c -lm; ./a.out 2> log | gnuplot
 
 make bisaintvenant.tst;  make bisaintvenant/plots ;  source  c2html.sh bisaintvenant
 
~~~

## Results
 
Plots of the wave generated by the collapse of the granular heap,
plot at time 0, 1, 5 and 10
 
~~~gnuplot
 set multiplot layout 4,1
 set xlabel 'x'
 plot [][:3]'log' every :::0::0 w l lc 2 t 'granular','' u 1:3 every :::0::0 t'water'w l lc 3
 plot [][:3]'log' every :::1::1 w l lc 2 t 'granular','' u 1:3 every :::1::1 t'water'w l lc 3
 plot [][:3]'log' every :::5::5 w l lc 2 t 'granular','' u 1:3 every :::5::5 t'water'w l lc 3
 plot [][:3]'log' every :::10::10 w l lc 2 t 'granular','' u 1:3 every :::10::10 t'water'w l lc 3
 unset multiplot
 
~~~

 
## To do
 
 Put dispersion!
 
 improve dry / wet transition...
 
 improve with density different...
 
 
## Links
 
 * [http://basilisk.fr/src/test/artery1D.c]()
 
 * etc
 
## Bibliography
 
 * Manon Robbe-Saule, Cyprien Morize, Robin Henaff, Yann Bertho, Alban Sauret, and Philippe Gondret
 "Experimental investigation of tsunami waves generated by granular collapse into water"
 
 * S. Viroulet A. Sauret O. Kimmoun  C. Kharif
 "Granular collapse into water: toward tsunami
 landslides"
 J Vis (2013) 16:189–191
 DOI 10.1007/s12650-013-0171-4
 
 * etc
 
 */
