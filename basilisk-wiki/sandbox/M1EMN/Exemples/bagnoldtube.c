/** 

#  Bagnold flow in a pipe 

We condider the flow of granular media in a horizontal tube 

$$\tau = 2 \eta_{eq} D$$

Flow under gravity (in $-e_x$), with confining pressure  $P_0$.


Un écoulement "bouchon" à vitesse constante se produit dans la majeure partie du tube, 
une bande de cisaillement près des parois permettent à la  vitesse de passer à zéro à la paroi. 

 

Equations d'équilibre à accélération négligée (stationnaire):
 
Prenons $P=P_0$ constant
$$0=-\frac{\partial P}{\partial r}\;   
   et\;   0=   \frac{\partial }{r \partial r} ( r \tau_{rx})
-\rho g $$
donc, tant que $\tau_{rx}  < \mu_s P_0$ on est dans le bouchon solide, la contrainte n'est pas assez forte pour créer un mouvement relatif. 
Par symétrie au centre $\tau_{rx}$ est nul, et on peut intégrer
$\tau_{rx} = \rho g \frac{r}2,$ [remarque, en 2D plan, il n'y a pas de 1/2, et le reste de l'analyse est le même].

Pour une certaine position $r=r^*$ on est tout juste au seuil  $\tau_{rx} (r=r^*) = \mu_s P_0$.
Ensuite, pour $r>r^*$, $\tau_{rx}$ reste linéaire, mais on l'écrit:
$$\tau_{xy}   = \frac{1}2 \rho g (r-r^{*}) +  \mu_s P_0$$
On a calculé la contrainte, pour $r<r^*$ on est dans l'écoulement solide en bloc (le bouchon), 
pour $r>r^*$ on est dans l'écoulement fluide donc  on a $\tau_{rx}=\mu(I) P$ 
et donc   :
$$    \frac{\Delta \mu}{1+I_0/I} =  \frac{1}{2 P_0} \rho g (r-r^{*})   $$
ou  après inversion de $\mu(I)$ et intégration
$$ u(r) = \left(\frac{I_0 R}{ d } \sqrt{\frac{P_0}{\rho }} \right)
\frac{( R-r)   \frac{1}{2} \rho g    +    \Delta \mu P_0 Log[\frac{ \Delta \mu P_0 +  
 \frac{1}{2 } \rho g (r^* -R)}{ \Delta \mu P_0 +   \frac{1}{2} \rho g (r^*-r)}]}{  \frac{1}{2 } \rho g R}$$

vitesse au centre 
$$ u(r) = \left(\frac{I_0 R}{ d } \sqrt{\frac{P_0}{\rho }} \right)
\frac{( R-r^*)   \frac{1}{2} \rho g    +    \Delta \mu P_0 Log[\frac{ \Delta \mu P_0 +  
 \frac{1}{2 } \rho g (r^* -R)}{ \Delta \mu P_0 }]}{  \frac{1}{2 } \rho g R}$$
 


*/

#include "axi.h"
#include "navier-stokes/centered.h"
#define LEVEL 6
double mumax,dg,P0;
scalar mu_eq[];

/** 
Bagnold solution for comparison, 

*/
double Uba( double r){
    double dmu=.26,mu0=.38,I0=.3;
    //double rgs2=1;//
    double rgs2=1./2;
    double rs=mu0*P0/rgs2;
    double uu = 30*2*I0*((1 - r)*rgs2 + dmu*P0*log(rgs2*(rs-1) + dmu*P0 ) - 
       dmu*P0*log(dmu*P0 + rgs2*(rs-r) ));

    double uup = 30*2*I0*((1 - rs)*rgs2 + dmu*P0*log(rgs2*(rs-1) + dmu*P0 ) - 
       dmu*P0*log(dmu*P0 + rgs2*(rs-rs) ));
     return (r < rs? uup:uu ); }
/**
Main with parameters
*/
 int main() {
  L0 = 1.;  
  P0 = 1;
  DT = 0.05;

  // P0 =2;2D 
/**
  the regularisation value of viscosity 
*/
  mumax=5000;
/**
 Boundary conditions are periodic
*/
    periodic (right);
/**
  no slip at the top, confinment pressure $P_0$ 
*/
    u.t[top] = dirichlet(0);
    u.n[top] = dirichlet(0);
    p[top] = dirichlet(P0);
/**
  symmetry at the bottom
*/  
    u.n[bottom] = dirichlet(0);
    u.t[bottom] = neumann(0); 
  //  stokes = true; // because U=u(y)e_x. stokes true > no CFL condition
  run(); 
}


face vector muv[];

event init (t = 0) {
/** 
 prepare viscosity
*/
  mu = muv;
/**
  minus pressure gradient, or equivalent gravity acceleration `mdpdx`, this is equivalent to add the 
  gravity 
 $$-\frac{\partial p}{\partial x} = -1 $$
 $$-\frac{\partial p}{\partial y} = 0 $$
*/
    const face vector mdpdx[] = {-1,0};
/**
 note that in "accceleration" in "navier-stokes/centered.h" there is the `fm`metric term in front.

 `event acceleration (i++,last)`

  `uf.x[] = fm.x[]*(face_value (u.x, 0) + dt*a.x[]);`

  will be the same for `alphav.x[] = fm.x[]/rho(ff); ` next...
*/
  a = mdpdx;
/**
 Initialy at rest
*/
  foreach() {
    u.x[] = 0;
    u.y[] = 0;
    p[]=1;
  }
}

/**
We check the number of iterations of the Poisson and viscous
problems. */
//event logfile (i++)
// fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);
/**
 old value of the velocity is saved
*/
scalar un[];
event init_un (i = 0) {
    foreach()
    un[] = u.x[];
}
/**
 so that when it does not more change we are converged
*/
event conv (t += 1; t < 150) {
    double du = change (u.x, un);
    fprintf(stdout,"t= %g %g %g %g \n",t,interpolate (u.x, L0/2, 0),interpolate (p, L0/2, 0),du);
    if (i > 0 && du < 1.0e-6)
        return 1; /* stop */
}  

/**
## Implementation of the $\mu(I)$ viscosity
*/

event nonnewviscosity(i++) {
    scalar eta_eq[]; 

/** computation of the second invariant as defined by Darby $-II_2 = 2 D:D$  and $D_2=\sqrt{D:D}$
$$ 2 D:D = (2 [(\frac{\partial v}{\partial y})^2  + (\frac{ v}{ y})^2) +(\frac{\partial u}{\partial x})^2] +
   [\frac{\partial v}{\partial x} + \frac{\partial u}{\partial y}]^2) $$
   Note that $y$ is $r$

so viscosity is 
$$
\eta_{eq} = \mu(I)P/(\sqrt(2.)D2)
$$
with regularisation
*/
  
    scalar shear[];
    foreach()
    shear[] = fabs((u.x[0,1] - u.x[0,-1])/(2.*Delta));
    boundary ({shear});

 foreach() { 
      double mI2 = 0.,D2 = 0,In = 0, muI = 0;
      dg = 1./30;
      double duxx = (u.x[1,0] - u.x[-1,0])/(2 * Delta);
      double duxy = (u.x[0,1] - u.x[0,-1])/(2 * Delta);
      double duyx = (u.y[1,0] - u.y[-1,0])/(2 * Delta);
      double duyy = (u.y[0,1] - u.y[0,-1])/(2 * Delta); 
	    mI2 =  sq(duyx+duxy) + 2*(sq(duyy) + sq(duxx) + sq(u.y[]/ max(y, 1e-20)));
      D2 = sqrt(mI2/2.);
      In = sqrt(2.)*dg*D2/sqrt(fabs(p[]));
      //In =  dg*shear[]/sqrt(fabs(P0));
      muI = .38 + (.26)*In/(.3 + In);
      if(D2>0){
        eta_eq[] = min(muI*fabs(p[])/(sqrt(2.)*D2) , mumax );}
       else {
        eta_eq[]=mumax;	
       }
    }   
      boundary ({eta_eq});

 
    boundary ({mu_eq});
    foreach_face() {
        muv.x[] = fm.x[]*(eta_eq[] + eta_eq[-1,0])/2.;
    }
    boundary ((scalar *){muv});
}
/**
  Save profiles computed, shear and exact
*/
event profiles (t += 1)
{
    FILE * fp = fopen("xprof", "w");
    scalar shear[];
    foreach()
    shear[] = (u.x[0,1] - u.x[0,-1])/(2.*Delta);
    boundary ({shear});
    for (double y = 0.; y < 1.0; y += 1./pow(2.,LEVEL))
        fprintf (fp, "%g %g %g %g %g  %g \n", y, interpolate (u.x, L0/2, y), interpolate (u.y, L0/2, y),
        	   interpolate (shear, L0/2, y),interpolate (p, L0/2, y),
                  Uba(y));
    fclose (fp);
}
/**
We adapt according to the error on the velocity field. 
*/
event adapt (i++) {
	 fprintf (stderr," %g \n",t);
	  // adapt_wavelet ({u}, (double[]){3e-3,3e-3}, 8, 6);
}

event profile (t = end) {
  foreach()
    printf ("%g %g %g %g %g\n", x, y, u.x[], u.y[], p[]);
}

/**

## Compilation 

~~~bash
make bagnoldtube.tst;
make bagnoldtube/plots;
make bagnoldtube.c.html;
~~~



## Results and plots

Plot of the exact and computed velocities.

~~~gnuplot profiles
set ylabel "u(y)";set xlabel "y" 
p[][-1:2]'xprof' t'u comp',''u 1:5 t'p','' u 1:($6/1.) t'exact' w l      
~~~



## Bibliography

* R. Darby Viscoelastic fluids, Dekker ed. (1976) p 223-225, p 194 

* [same with Bingham](http://basilisk.fr/sandbox/M1EMN/Exemples/nonnewtube.c)


*/

