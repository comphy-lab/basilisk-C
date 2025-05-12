/**
 
#  Dicharge of a cylindrical silo 
 
## Problem

 We look at an axi flow of cohesive grains going out a silo: How a cohesive granualar media flows in a cylindrical silo?

The gravity is from right to left, the bottom of the image is the symmetry line. The hole is at the left at the bottom corner...

Yellow grains, blue air. 

![Snapshot of silo. t=0](siloaxiCo/f0.png)
![Snapshot of silo. t=1](siloaxiCo/f1.png)
![Snapshot of silo. t=2](siloaxiCo/f2.png)

## Equations

 
 We solve the granular flow in a silo of cylindrical shape  (using `mask`) of radius $R$ chosen as the unit.
 
 
 We change the radius $R_0$ of the hole to test Beverloo's law, the flow rate is
 compared to the scaling law
 $$Q \sim \sqrt{g} R_0^{5/2}$$
for dry granular flow. Then  we change the value of the cohesion yeld stress $\tau_0$
in a $\mu(I)$ rheology with  $\tau = \tau_0 + \mu(I)p$
to  see when the flow stops.
 We use cohesion without dimension
 $Co=\tau_0/(\rho g R)$, we will see that the pertinent parameter is $\tau_0/(\rho g R_0)$,
 
 
~~~gnuplot silo
 set arrow 1 from 1.6,.2 to .3,.2 front
 set arrow 2 from 0.05,0 to 0.05,.2 front
 set arrow 3 from 0.15,0 to .15,1 front
 set arrow 4 nohead lw 3 from 0.,0.2 to .0,1 front
 set arrow 5 nohead lw 3 from 0.,1 to 2,1 front
 set arrow 6 from 0.05,.2  to 0.05,0 front

 set label 1 "gravity" at 1,.25 front
 set label 2 "R0" at 0.1,.15 front
 set label 3 "\"top\"" at 1,1.1 front
 set label 4 "grains" at 1.2,.5 front
 set label 5 "air" at 1.8,.5 front
 set label 6 "R" at .15,1.05 front
 set xlabel "axis of cylindrical symmetry    x"
 set ylabel "r"
 
p [0:2][0:2]0 not,(x<1.8) w filledcurves x1 linec rgb "yellow" not,(x>=1.75)  w filledcurves x1  lc rgb "cyan" not
 unset arrow 1; unset arrow 2; unset arrow 3
 unset arrow 4; unset arrow 5; unset arrow 6
 unset label 1; unset label 2; unset label 3
 unset label 4; unset label 5; unset label 6
 
~~~

*/
#include "axi.h"
#include "navier-stokes/centered.h"
#include "vof.h"
#define LEVEL 7
#define RHOF 1e-4
#define mug  1e-5

double mumax,Co,dg,P0,D,Q,tmax, h0,R0;
scalar mu_eq[];
scalar f[];
scalar * interfaces = {f};
face vector alphav[];  
face vector muv[];
scalar rhov[];


/**
 Boundary conditions: we have to remember that $y$ is the radial coordinate, $x$ is along the axis of the silo.
 The gravity will be horizontal along $x$. Everything is turned.
 
 no slip at the top (latteral wall of the silo)
 */
u.t[top] = dirichlet(0);
u.n[top] = dirichlet(0);
/**
 pressure 0 at the right and at the hole
 */
p[right] =  dirichlet( 0);
p[left] = (y <= R0 ? dirichlet(0) : neumann(0) );
u.n[left] =  (y<= R0 ? neumann(0) :dirichlet(0));
u.t[left] =  (y <= R0 ? neumann(0) :dirichlet(0));
u.n[right] =  neumann(0);
u.t[right] = neumann(0);
/**
 symmetry at the bottom (which is the axis of the silo and which is alined with gravity)
 */
u.n[bottom] = dirichlet(0);
u.t[bottom] = neumann(0);

/**
 Main with parameters, height of the silo is twice of the radius; diameter of grain is $R/90$. Iinital relative filling height 1.75.
 */
int main() {
    
    L0 = 2.; 
    DT = 0.015;  // 8 0.025;
    dg = 1./90;
    h0 = 1.75;
    /**
     the regularisation value of viscosity
     */
    mumax=1000;
    /**
     equivalent gravity acceleration (horizontal and toward the left), viscosity density,
     */
    const face vector grav[] = {-1,0};
    a = grav;
    mu = muv;
    alpha = alphav;
    rho = rhov;
/**   
 Loop in $R_0$ and $Co$
*/ 
    FILE * fp3 = fopen("Qout", "w");
    fclose(fp3);
    FILE  *fwq = fopen ("outQCo", "w");
    fclose(fwq);
    for(R0 = 0.2;  R0 <= .45      ; R0+= .05  ){
    for(Co = 0.  ;  Co <= .5*R0+.025; Co+= .025 ){
        tmax = 4;
        run();
        fwq = fopen ("outQCo", "a");
        fprintf(fwq,"  %lf %lf %lf \n",  Q, Co ,R0);
        fclose (fwq);
    }
        fwq = fopen ("outQCo", "a");
        fprintf(fwq,"\n");
        fclose (fwq);
    }
    fclose(fwq);

}

/**
Initalisation, relative filling height $h_0$
*/
event init (t = 0) { 

    fraction (f, (h0-x));
/**
  domain is 2x1
*/
    mask(y>1 ? top : none);
/**
     Initialy at rest
*/
    foreach() {
        u.x[] = 0;
        u.y[] = 0;
        p[]=(y<D && fabs(x) <= .1) ?
        0 : max(h0 - x,0);
    }
}

/**
 We check the number of iterations of the Poisson and viscous
 problems. */
event logfile (i++)
 fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);

/**
## Implementation of the $\mu(I)$ viscosity
 */
#define rho(f) ((f) + RHOF*(1. - (f)))
/**
 prepare viscosity
 */
event nonnewviscosity(i++) {
    scalar eta_eq[];
    
 /** computation of the second invariant $-II_2$ as defined by Darby $-II_2 = 2 D:D$  and $D_2=\sqrt{D:D}$
$$2 D:D = (2 [(\frac{\partial v}{\partial y})^2  + (\frac{ v}{ y})^2) +(\frac{\partial u}{\partial x})^2] +
\frac{\partial v}{\partial x} + \frac{\partial u}{\partial y}]^2) $$
     Note that $y$ is $r$.
     So viscosity is
     
$$\eta_{eq} = min(\frac{(Co+\mu(I)p)}{(\sqrt{2}D_2)},\eta_{max})$$
     with regularisation
 */
    
    scalar shear[];
    foreach()
    shear[] = fabs((u.x[0,1] - u.x[0,-1])/(2.*Delta));
    boundary ({shear});
     
    foreach() {
        double mI2 = 0.,D2 = 0,In = 0, muI = 0;
        double duxx = (u.x[1,0] - u.x[-1,0])/(2 * Delta);
        double duxy = (u.x[0,1] - u.x[0,-1])/(2 * Delta);
        double duyx = (u.y[1,0] - u.y[-1,0])/(2 * Delta);
        double duyy = (u.y[0,1] - u.y[0,-1])/(2 * Delta);
        mI2 =  sq(duyx+duxy) + 2*(sq(duyy) + sq(duxx) + sq(u.y[]/ max(y, 1e-20)));
        D2 = sqrt(mI2/2.);
        In = sqrt(2.)*dg*D2/sqrt(fabs(p[])+1e-10); 
        muI = .4 + (.28)*In/(.4 + In);
        if(D2>0){
            eta_eq[] = min((Co + muI*fabs(p[]))/(sqrt(2.)*D2) , mumax );}
        else {
            eta_eq[]=mumax;
        }
    }
 
    boundary ({eta_eq});
    boundary ({mu_eq});
    
    scalar fa[];
    foreach()
    fa[] = (4.*f[] +
            2.*(f[-1,0] + f[1,0] + f[0,-1] + f[0,1]) +
            f[1,1] + f[-1,1] + f[1,-1] + f[-1,-1])/16.;
    boundary ({fa});
    
    foreach()
    rhov[] = rho(fa[]); 
    boundary ({rhov});

    /**
     note that in "accceleration" in "navier-stokes/centered.h" there is the `fm`metric term in front.
     
     `event acceleration (i++,last)`
     
     `uf.x[] = fm.x[]*(face_value (u.x, 0) + dt*a.x[]);`
     
     will be the same for `alphav.x[] = fm.x[]/rho(ff); ` next...
     */
    
    foreach_face() {
        double fM = (fa[] + fa[-1])/2.;
        muv.x[] = (fM*(eta_eq[] + eta_eq[-1])/2. + (1. - fM)*mug);
        alphav.x[] =   fm.x[]/rho(fM);
    }
    boundary ((scalar *){muv,alphav});
}

event interface (  t += 1  ) {
    char s[80];
    sprintf (s, "field-%g.txt", t);
    FILE * fp = fopen (s, "w");
    output_field ({p,u,uf,pf}, fp, linear = true);
    fclose (fp);
}

/**
 
 Here we compute the flow rate $\int_0^{R_0} r u(r)dr$
 */
event flowrate (t += 0.25; t < tmax ) {
   // static double Vold,V=1,Qinst=0;
   // Vold=V;
   // V=0;
   // foreach()
   // V = V + f[]*y*Delta*Delta;
   // Qinst = -(V-Vold)/.25;
 
    FILE * fp3 = fopen("Qout", "a");
    Q=0;

    //foreach_boundary (left)
    // Q -= (Delta)*y*u.x[] ;
    
    double dy=L0/pow(2.,LEVEL);
    for (double y = 0.; y < 1; y += dy)
        Q-=interpolate (u.x, 0, y)*y;
     Q=Q*dy;

    fprintf (fp3, "%6.4g %g \n", t, Q);
    fclose (fp3);
}

/**
## films
 
*/

#if 0
event movie (t += 0.05) {
    
    scalar l[];
    foreach()
    l[] = f[]*(1+sqrt(sq(u.x[]) + sq(u.y[])));
    boundary ({l});
    static FILE * fp2 = popen ("ppm2mpeg > velo.mpg", "w");
    output_ppm (l, fp2, min = 0, max = 2., linear = true,
                n = 512, box = {{0,-.5},{2,1}});
    foreach()
    l[] = f[]*p[];
    boundary ({l});
    static FILE * fp3 = popen ("ppm2mpeg > p.mpg", "w");
    output_ppm (l, fp3, min = 0, max = 2., linear = true,
                n = 512, box = {{0,0},{L0,L0}});
    
}
#else
event movie (t += 0.05) {
    
    scalar l[];
    foreach()
    l[] = f[]*(1+sqrt(sq(u.x[]) + sq(u.y[])));
    boundary ({l});
    output_ppm (l, file =  "velo.mp4", min = 0, max = 2., linear = true,
                n = 512, box = {{0,-.5},{2,1}});
    foreach()
    l[] = f[]*p[];
    boundary ({l});
    output_ppm (l,  file =  "p.mp4", min = 0, max = 2., linear = true,
                n = 512, box = {{0,0},{L0,L0}});
    
}

#endif

event pictures (t = {0, 1. , 2., 3., 4.} ) {
    if(Co<0.0000001){
    scalar l[];
    foreach()
      l[] = f[]*p[];
    boundary ({l});
    if(t==0.)
    output_ppm (f, file = "f0.png", min=-1, max = 2,  spread = 2, n = 256, linear = true,
                box = {{0,0},{4,2}});
    if(t==1.)
    output_ppm (f, file = "f1.png", min=-1, max = 2,  spread = 2, n = 256, linear = true,
                box = {{0,0},{4,2}});
    if(t==2.)
    output_ppm (f, file = "f2.png", min=-1, max = 2,  spread = 2, n = 256, linear = true,
                box = {{0,0},{4,2}});
    if(t==3.)
    output_ppm (f, file = "f3.png", min=0, max = 2,  spread = 2, n = 256, 
linear = true,
                box = {{0,0},{4,2}});}

}
/**
 
## Compilation
 
~~~bash
 make siloaxiCo.tst;
 make siloaxiCo/plots;
 make siloaxiCo.c.html;
~~~
 
~~~bash
 qcc -g -O2 -Wall -o siloaxiCo siloaxiCo.c -lm
 ./siloaxiCo
~~~
 
 
## Results and plots

 
 Plot of Hagen-Beverloo law : for $Co=0$ we reobtain $Q \sim \sqrt{g} R_0^{5/2}$.
 
 ~~~gnuplot
 reset
 set xlabel "R0"
 set ylabel "Q"
 p[0:][0:]'outQCo' u 3:($2==0?$1:NaN)t'Co=.00'  w lp,.8*x**2.5,'' u 3:($2==0.050000?$1:NaN) w lp t'Co=.05', ''u 3:($2==0.100000?$1:NaN) w lp t'Co=.10'
 ~~~
 
 
 
 
 The silo has a radius $R$, a simple model is to consider the equilibrium of a cylinder over the hole of radius $R_0$.
 This cylinder experiences cohesion and weight:
 $$\tau_0 (2 \pi R_0)= \rho g (\pi R_0^2), \text{  hence the no flow condition is }  \frac{\tau_0}{\rho g R_0}  = \frac{1}{2}.$$
 Remember that cohesion  $Co=\frac{\tau_0}{\rho g R}$
 is measured in the code with $R$, so that 
 $\frac{\tau_0}{\rho g R_0} = \frac{C_0}{(R_0/R)}$.
 
 
 We plot here $Q/R_0^{5/2}$ as a function of 
 $\frac{C_0}{(R_0/R)}$ (note that in the code $R=1$ so that the $R_0$ of the code is in fact $R_0/R$.
 
~~~gnuplot
 reset
 set ylabel "Q/(R0^2(gR0)^{1/2})";set xlabel "tau0/(rho g R_0)"
 p[:.8][0:1]'outQCo' u ($2/$3):($3>0.?$1/($3**2.5):NaN)w lp not
~~~
 
 
 
 ![Animation  ](./siloaxiCo/velo.mp4)

 ![Animation  ](./siloaxiCo/velo.mpg)
 
 
[![](./siloaxiCo/f0.png)](./siloaxiCo/velo.mp4)

  velocity  (click on image for animation)
 
## Links
 
 * [http://basilisk.fr/sandbox/M1EMN/Exemples/granular.h]()  $\mu(I)$ and silo and Co 
 
 * same in 2D, to come....
 
 * [http://basilisk.fr/sandbox/M1EMN/Exemples/granular_column_cohesif.c]() collapse granular 
 
## Bibliography
 
 * cohesive collapse granular
 
 */
