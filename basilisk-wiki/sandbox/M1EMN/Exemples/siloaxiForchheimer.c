/**
# Darcy Forchheimer flow in axicylindrical granular silo

*/
#include "axi.h"
#include "navier-stokes/centered.h"
#include "vof.h"
#define LEVEL 6
#define RHOF 1e-4
#define mug  1e-5

double mumax,dg,P0,D,Q,tmax, h0;
scalar mu_eq[];
scalar f[];
scalar * interfaces = {f};
face vector alphav[];
face vector muv[];
scalar rhov[];


face vector  ug[];
scalar In[];
scalar phi[];
face vector av[];
double ugranular;
double betam = 100;
scalar pDarcy[],sourceDarcy[];
face vector beta[];
face vector gradpDarcy[];
mgstats mgpDarcy;

pDarcy[top]   = neumann(0);
pDarcy[bottom]  = neumann(0);
pDarcy[right]    = dirichlet(P0);//neumann((Qair / (LDOMAIN * (1 - maskr) ) - ugranular) / betam );
pDarcy[left] =  (y<= D ? dirichlet(0) :neumann(0));


ug.t[top] = dirichlet(0);
ug.n[top] = dirichlet(0);
ug.n[left] =  (y<= D ? neumann(0) :dirichlet(0));
ug.t[left] =  (y <= D ? neumann(0) :dirichlet(0));
ug.n[right] =  neumann(0);
ug.t[right] = neumann(0);
ug.n[bottom] = dirichlet(0);
ug.t[bottom] = neumann(0);

In[left] = neumann(0.);
In[top] = neumann(0.);
In[bottom] = neumann(0.);
In[right] = neumann(0.);

phi[left] = neumann(0.);
phi[top] = neumann(0.);
phi[bottom] = neumann(0.);
phi[right] = neumann(0.);


u.t[top] = dirichlet(0);
u.n[top] = dirichlet(0);
p[right] =  dirichlet( 0);
p[left] = (y <= D ? dirichlet(0) : neumann(0) );
u.n[left] =  (y<= D ? neumann(0) :dirichlet(0));
u.t[left] =  (y <= D ? neumann(0) :dirichlet(0));
u.n[right] =  neumann(0);
u.t[right] = neumann(0);
u.n[bottom] = dirichlet(0);
u.t[bottom] = neumann(0);




/**
 Main with parameters
 */
int main() {
    L0 = 4.;
    /**
     'Jansen' pressure, but in fact any pressure gives almost the same 'Q'
     */
    P0 = 1.;//1./2/.4;
    DT = 0.01/2;
    D = 0.25 ;
    tmax=15;
    dg = 1./90;
    h0 = 3.8; 
    
    a = av;
    /**
     the regularisation value of viscosity
     */
    mumax=1000;
    /**
     Boundary conditions are periodic
     */
    alpha = alphav;
    mu = muv;
    rho = rhov;
    run();
}




event init (t = 0) {
/** silo of lenght (or hight!) 4 and radius 1
*/
    mask (y > L0/4. ? top : none);
    
    fraction (f, (h0 - x));
    scalar phi[];
    foreach_vertex()
    phi[] = (h0 - x);
    fractions (phi, f);
    
    foreach() {
        u.x[] = 0;
        u.y[] = 0;
        p[]=(y<D && fabs(x) <= .05) ?0 : max(h0 - x,0);//P0;
        pDarcy[] = P0;
        sourceDarcy[] = 0.;
    }
    boundary ({pDarcy, sourceDarcy});
    
    
    foreach()
    In[]=0;
    foreach()
    phi[]=0;
    
    
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
event conv (t += 1; t < tmax) {
    double du = change (u.x, un);
    fprintf(stdout,"t= %g %g %g %g \n",t,interpolate (u.x, L0/2, 0),interpolate (p, L0/2, 0),du);
    if (i > 0 && du < 1.0e-6)
        return 1; /* stop */
}

/**
## Implementation of the $\mu(I)$ viscosity
 */

#define rho(f) ((f) + RHOF*(1. - (f)))

event nonnewviscosity(i++) {
    scalar eta_eq[];
    trash ({alphav});
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
        double duxx = (u.x[1,0] - u.x[-1,0])/(2 * Delta);
        double duxy = (u.x[0,1] - u.x[0,-1])/(2 * Delta);
        double duyx = (u.y[1,0] - u.y[-1,0])/(2 * Delta);
        double duyy = (u.y[0,1] - u.y[0,-1])/(2 * Delta);
        mI2 =  sq(duyx+duxy) + 2*(sq(duyy) + sq(duxx) + sq(u.y[]/ max(y, 1e-20)));
        D2 = sqrt(mI2/2.);
        In = sqrt(2.)*dg*D2/sqrt(fabs(p[])+1e-10);
        //In =  dg*shear[]/sqrt(fabs(P0));
        muI = .4 + (.28)*In/(.4 + In);
        if(D2>0){
            eta_eq[] = min(muI*fabs(p[])/(sqrt(2.)*D2) , mumax );}
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
    
    
    foreach_face() {
        double fM = (fa[] + fa[-1])/2.;
        muv.x[] = (fM*(eta_eq[] + eta_eq[-1])/2. + (1. - fM)*mug);
 //       double coeff=In[]>0.6 ? .33 *.6 : 0.33 * In[];
 //       alphav.x[] = 1./( fm * (1. - coeff) + RHOF*(1. - fm));
        alphav.x[] = fm.x[]/rho(fM);
    }
    foreach()
    rhov[] = rho(fa[]);
    boundary ({muv,alphav,rhov});
    //boundary ((scalar *){muv,alphav});
}

 

event pressureDEF (i++)
{
    
    /**
     'betam' is a penalisation parameter, so that pressure is constant ouside
     the grains and as
     $(\beta dp/dn)$ is constant across the interface
     and $\beta_m dp/dy = Qair$ at the top */
    double B_l=1.0;
    double B_i=1.0;
    double R_f=1.0;//4.0*1.2/2500;
    
    //double beta_isl = 0.1;
    scalar un[];
    vector uold[];
    face vector ur[];
    double err=1;
    
    foreach_face(){
        ur.x[]=ug.x[]-u.x[];
    }
    boundary ((scalar *){ur});
    // while (err > 0.01){
/**

  note that we have to solve a problem like 
  $$\nabla \cdot \beta \nabla p = s$$
  where the source term is as well a divergence of a vector: $s=\nabla \cdot U_s$
  
  in two D it is 
  $$\frac{\partial }{\partial x}[ \beta_x \frac{\partial }{\partial x}p]+
  \frac{\partial }{\partial y}[ \beta_y \frac{\partial }{\partial y}p] = 
  \frac{\partial }{\partial x}U_{sx}+ \frac{\partial }{\partial y}U_{sy}.$$

In axi $y$ is $r$ and :

$$\frac{\partial }{\partial x}[ y \beta_x \frac{\partial }{\partial x}p]+
  \frac{\partial }{\partial y}[ y \beta_y \frac{\partial }{\partial y}p] = 
  \frac{\partial }{\partial x}(yU_{sx})+ \frac{\partial }{\partial y}(yU_{sy}).$$

that the reason why we multiply by `y`  the vector `beta` and others.
*/ 
  
  
  
  
    foreach_face(){
        uold.x[] = ug.x[];
        //ur.x[]=ug.x[]-u.x[];
        un[] = norm(ur);
        beta.x[] =1.0/(R_f+dt*(B_l+B_i*un[]))*f[]  + (betam/(R_f+dt))*(1. - f[]);
        beta.x[] =  beta.x[] *y;
    }
    boundary ((scalar *){beta});
    //beta.x[] = 1.0/(1.0+beta_isl*un[])*f[]  + betam*(1. - f[]);}
    //beta.x[] = (f[] + betam*(1. - f[]));
    vector us[];
    foreach_face(){
        //us.x[]=R_f*ug.x[]/(R_f+dt*(B_l+B_i*un[]));
        us.x[]=(R_f*ug.x[]+dt*(B_l*u.x[]+B_i*u.x[]*un[]))/(R_f+dt*(B_l+B_i*un[]));
        us.x[]= us.x[]*y;
    }
    boundary ((scalar *){us});
    scalar div[];
    foreach(){
        div[]=0.0;
        foreach_dimension()
        div[] +=us.x[1]-us.x[];
        div[] /=dt*Delta;
    }
    
    sourceDarcy=div;
    
    
    
    
    mgpDarcy = poisson (pDarcy, sourceDarcy, beta);
    // betam-big?
    face vector gradp[];
    foreach() {
        //ur.x[]=ug.x[]-u.x[];
        un[] = norm(ur);
        foreach_dimension(){
            gradp.x[] =(pDarcy[0,0] - pDarcy[-1,0])/Delta;
            ug.x[]= (R_f*uold.x[] +(-gradp.x[]+B_l*u.x[]+B_i*u.x[]*un[])*dt)/(R_f +dt*( B_i * un[]+B_l))*f[]+
            (R_f/betam*uold.x[] -gradp.x[]*dt) / ((R_f +dt)/betam)*(1-f[]) ;
        }
    }
    boundary ((scalar *){ug});
    
    double max=0;
    foreach() {
        double du = fabs(un[] -norm(ug));
        if(du > max) max = du;
    }
    err=max;
    fprintf (stderr, "%g   du= %g\n", t, err);
    
}

event acceleration (i++)
{
    double eps = 1.0;
    
    /**
     Coupling gravity : Grad(p) added to gravity */
    
    foreach_face(y)
    av.y[] +=  - eps*(pDarcy[0,0] - pDarcy[0,-1])/Delta;
    foreach_face(x)
    av.x[] += - 1.- eps*(pDarcy[0,0] - pDarcy[-1,0])/Delta ;
    boundary ((scalar *){av});
}
 



/**
 Save profiles computed, shear and exact
 */

event interface (  t += 1  ) {
    char s[80];
    sprintf (s, "field-%g.txt", t);
    //sprintf (s, "field.txt");
    foreach_face()
    gradpDarcy.x[] = (pDarcy[0] - pDarcy[-1])/Delta ;
    boundary ((scalar *){gradpDarcy});
    
    FILE * fp = fopen (s, "w");
    output_field ({f,p,u,uf,pf,pDarcy,gradpDarcy,In,ug,alphav}, fp, linear = true);
    fclose (fp);
}

/**
 We adapt according to the error on the velocity field.
 */
event adapt (t += 0.25 ) {
    FILE * fp3 = fopen("Qout", "a");
    Q=0;
    double dy=1./pow(2.,LEVEL);
    for (double y = 0.; y < 1.; y += dy)
        Q+=interpolate (u.x, dy, y);
    Q=Q*dy;
    //fprintf (stderr," %6.4g  %g \n",t,Q);
    fprintf (fp3, "%6.4g %g \n", t, Q);
    fclose (fp3);
    // adapt_wavelet ({u}, (double[]){3e-3,3e-3}, 8, 6);
}

/**
 event profile (t = end) {
 foreach()
 printf ("%g %g %g %g %g\n", x, y, u.x[], u.y[], p[]);
 }
 */

#if 1
event movie (t += 0.05) {
    
    scalar l[];
    foreach()
    l[] = f[]*(1+sqrt(sq(u.x[]) + sq(u.y[])));
    boundary ({l});
    static FILE * fp2 = popen ("ppm2mpeg > velo.mpg", "w");
    output_ppm (l, fp2, min = 0, max = 2., linear = true,
                n = 512, box = {{0,-2},{L0,L0}});
    
    foreach()
    l[] = f[]*p[];
    boundary ({l});
    static FILE * fp3 = popen ("ppm2mpeg > p.mpg", "w");
    output_ppm (l, fp3, min = 0, max = 2., linear = true,
                n = 512, box = {{0,-2},{L0,L0}});
    
    foreach()
    l[] = f[]*pDarcy[];
    boundary ({l});
    static FILE * fp4 = popen ("ppm2mpeg > pDarcy.mpg", "w");
    output_ppm (l, fp4, min = 0, max = 2., linear = true,
                n = 512, box = {{0,-2},{L0,L0}});
    
}

#endif



/**
 
## Compilation
 
~~~bash
 ln -s siloaxiForchheimer.c.page siloaxiForchheimer.c
 make siloaxiForchheimer.tst;
 make siloaxiForchheimer/plots;
 make siloaxiForchheimer.c.html;
~~~
 
~~~bash
 qcc -g -O2 -Wall -o siloaxiForchheimer siloaxiForchheimer.c -lm
 ./siloaxiForchheimer
~~~
 
 
## Results and plots
 
 Plot of the Darcy Forchheimer pressure
 
~~~gnuplot profiles
 set ylabel "p(x)";set xlabel "x"
 set key left
 p'field-14.txt' u 1:($2==0? $10:NaN) t'p Darcy' w lp,''u 1:($2==0? $4:NaN) t'p grains' w lp
~~~
 
 
 Granular pressure 
 
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
 sp 'field-14.txt' u 1:($2<1?$2:NaN):($4*$3) not
~~~
 
 
Darcy pressure 
 
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
 sp 'field-14.txt' u 1:($2<1?$2:NaN):($10) not
~~~


## Bibliography
 
*  mu(I) and silo and Co
 
 
*/
