/** 

#  Bagnold flow in cylindrical silo


 


*/

#include "axi.h"
#include "navier-stokes/centered.h"
#define LEVEL 8
double mumax,dg,P0,R0,Q;
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
  /**
'Jansen' pressure, but in fact any pressure gives almost the same 'Q'
  */
  P0 = 1./2/.38;
  DT = 0.01/2;
  R0 = 0.1 ; 

  // P0 =2;2D 
/**
  the regularisation value of viscosity 
*/
  mumax=1000;
/**
 Boundary conditions are periodic
*/
     
/**
  no slip at the top, confinment pressure $P_0$ 
*/
    u.t[top] = dirichlet(0);
    u.n[top] = dirichlet(0);
    p[left] =  dirichlet( (t < 1 )? P0 : P0);
    p[right] = (y <= R0 ? dirichlet(0) : neumann(0) );
    u.n[right] =  (y<= R0 ? neumann(0) :dirichlet(0));
    u.t[right] =  (y <= R0 ? neumann(0) :dirichlet(0));
    u.n[left] =  neumann(0);
    u.t[left] = neumann(0); 
/**
  symmetry at the bottom
*/  
    u.n[bottom] = dirichlet(0);
    u.t[bottom] = neumann(0); 
 
  run(); 
}


face vector muv[];

event init (t = 0) {
/** 
 prepare viscosity
*/
  mu = muv;
/**
   equivalent gravity acceleration 
*/
    const face vector grav[] = {1,0};
/**
 note that in "accceleration" in "navier-stokes/centered.h" there is the `fm`metric term in front.

 `event acceleration (i++,last)`

  `uf.x[] = fm.x[]*(face_value (u.x, 0) + dt*a.x[]);`

  will be the same for `alphav.x[] = fm.x[]/rho(ff); ` next...
*/
  a = grav;
/**
 Initialy at rest
*/
  foreach() {
    u.x[] = 0;
    u.y[] = 0;
    p[]=P0;
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
      In = sqrt(2.)*dg*D2/sqrt(fabs(p[])+1e-10);
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

 event interface (  t += .25  ) {
  char s[80];
  //sprintf (s, "field-%g.txt", t);
  sprintf (s, "field.txt");
  FILE * fp = fopen (s, "w");
  output_field ({p,u,uf,pf}, fp, linear = true);
  fclose (fp);
}

/**
We adapt according to the error on the velocity field. 
*/
event adapt (i++) {
      Q=0;
    double dy=1./pow(2.,LEVEL);
    for (double y = 0.; y < 1.0; y += dy)
        Q+=interpolate (u.x, L0-dy, y);
    Q=Q*dy;
	 fprintf (stderr," %6.4g  %g \n",t,Q);
	  // adapt_wavelet ({u}, (double[]){3e-3,3e-3}, 8, 6);
}

event profile (t = end) {
  foreach()
    printf ("%g %g %g %g %g\n", x, y, u.x[], u.y[], p[]);
}

/**

## Compilation 

~~~bash
ln -s siloaxi.c.page siloaxi.c              
make siloaxi.tst;
make siloaxi/plots;
make siloaxi.c.html;
~~~



## Results and plots

Plot of the exact and computed velocities.

~~~gnuplot profiles
set ylabel "u(y)";set xlabel "y" 
p'field.txt' u ($4+$1*3):2 w l
~~~



~~~gnuplot 
reset
set pm3d; set palette rgbformulae 22,13,-31;unset surface;
set ticslevel 0;
unset border;
unset xtics;
unset ytics;
unset ztics;
unset colorbox; 
set view 0,0
sp'field.txt' u 1:2:3 not
~~~  

## Bibliography

*  mu(I) and silo and Co


*/

