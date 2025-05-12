/**
# Flow in a Sandglass / Hourglass with a bottom orifice
 
 We propose an implementation of the Jop Pouliquen Forterre $\mu(I)$
 rheology for the flow in a hourglass.

Coupled trough Jackson method with a gaz with Darcy Forchheimer
  
 
# Code
 
 Includes and definitions
 */

#include "navier-stokes/centered.h"
#include "vof.h"
// Domain extent
#define LDOMAIN 4.
// heap definition
double  H0,R0,D,W,tmax,Q,Qair,Wmin,muwall,WDOMAIN,Hmask,maskr,P1,position,epai;
//film size
double  Lz;

/**
 passive fluid small density to preserve 0 pressure
 and small viscocity */

#define RHOF 1e-4
#define mug  1e-5
// Maximum refinement
#define LEVEL 7
// ratio of mask section
#define maskr 0.75
//#define MINLEVEL 7
char s[80];
FILE * fpf,*fwq;
scalar f[];
face vector  ug[];
scalar In[];
scalar phi[];

scalar * interfaces = {f};
face vector alphav[];
face vector muv[];
face vector av[];

//velocity of granular materials
double ugranular;

/**
 Darcy fluid */

double betam = 100;
scalar pDarcy[],sourceDarcy[];
face vector beta[];
face vector gradpDarcy[];
mgstats mgpDarcy;

/**
 Boundary conditions for air flow, impose pressure at the top $P_1$ 
 */

pDarcy[left]   = neumann(0);
pDarcy[right]  = neumann(0);
pDarcy[top]    = dirichlet(P1); 
pDarcy[bottom] = dirichlet(0);
pDarcy[right]= neumann(0);
pDarcy[left]=neumann(0);

ug.n[top]    = neumann(0);
ug.n[bottom] = neumann(0);
ug.n[right] = dirichlet(0);
ug.n[left] = dirichlet(0);



/**
 Boundary conditions for granular flow, pressure must be zero at the
 surface.  The pressure is zero in the hole $|x|<=W/2$, but the
 lithostatic gradient is given elsewhere on the bottom wall.  No slip
 boundary conditions on the other walls. */
#define test 1

#if test
p[top]      = dirichlet(0);
u.n[top]    = neumann(0);
u.t[top]    = dirichlet(0);
u.t[bottom] = dirichlet(0);
u.n[bottom] = dirichlet(0);
u.n[right]  = dirichlet(0);
u.n[left]   = dirichlet(0);

/*
p[bottom]   = neumann(0);
p[right]= neumann(0);
p[left]=neumann(0);
*/

 

In[left] = neumann(0.);
In[top] = neumann(0.);
In[bottom] = neumann(0.);
In[right] = neumann(0.);

phi[left] = neumann(0.);
phi[top] = neumann(0.);
phi[bottom] = neumann(0.);
phi[right] = neumann(0.);

#else
p[top]      = dirichlet(0);
u.n[top]    = neumann(0);
u.t[bottom] = neumann(0);
u.n[bottom] = neumann(0);
p[bottom]   = dirichlet(0);
p[right]=(fabs(y)>= 0. && fabs(y)<= (LDOMAIN*position-epai)) ? dirichlet(0):  neumann(0);
p[left]=(fabs(y)>= 0. && fabs(y)<= (LDOMAIN*position-epai)) ? dirichlet(0):  neumann(0);


u.n[right] = (fabs(y)>= 0. && fabs(y)<= (LDOMAIN*position-epai)) ? neumann(0):  dirichlet(0);
u.t[right] = (fabs(y)>= 0. && fabs(y)<= (LDOMAIN*position-epai)) ? neumann(0):  dirichlet(0);
u.n[left] = (fabs(y)>= 0. && fabs(y)<= (LDOMAIN*position-epai)) ? neumann(0):  dirichlet(0);
u.t[left] = (fabs(y)>= 0. && fabs(y)<= (LDOMAIN*position-epai)) ? neumann(0):  dirichlet(0);
bid plaque;
u.t[plaque] = dirichlet(0.);
u.n[plaque] = dirichlet(0.);

p[plaque] = neumann(0.);

In[left] = neumann(0.);
In[top] = neumann(0.);
In[bottom] = neumann(0.);
In[right] = neumann(0.);

phi[left] = neumann(0.);
phi[top] = neumann(0.);
phi[bottom] = neumann(0.);
phi[right] = neumann(0.);

#endif

int main(int argc, char ** argv) {
    L0 = LDOMAIN;
    // number of grid points
    N = 1 << LEVEL;
    // maximum timestep
    DT = 0.01;
    // coefficient of friction of wall
    muwall = 0.1;
    TOLERANCE = 1e-3;
    H0 = 2.5;
    R0 = 20.000;
    // Grain size
    D = 1. / 90.;
    fwq = fopen("outWQ", "w"); // ????
    fclose(fwq);               // ????
    system("rm uprof.txt");
    Lz = LDOMAIN;
    // size of the hole
    W = 0.25;
    //Qair = 1.0;
    P1=5.;
    WDOMAIN = 2.;
    a = av;
    alpha = alphav;
    mu = muv;
    epai=0.1;
    position=1./3.;
    Q = 0;
    tmax = 2.;
    fpf = fopen("interface.txt", "w");
    run();
    fclose(fpf);
    fprintf(stdout, "\n");
    fwq = fopen("outWQ", "a");
    fprintf(fwq, " %lf %lf \n", W, Q);
    fclose(fwq);
}

#if 0
int adapt() {
#if TREE
    astats s = adapt_wavelet ({f}, (double[]){5e-3}, LEVEL, MINLEVEL);
    return s.nf;
#else // Cartesian
    return 0;
#endif
}
#endif

/**
 initial heap, a rectangle */

// note the v
event init (t = 0) {
      
    scalar phiinit[];
    foreach_vertex()
     {
        phiinit[] = min(H0 - y, R0 - x);
    }
    fractions (phiinit, f);
    
    /**
     lithostatic pressure, with a zero pressure near the hole
     to help */
    
    foreach()
     p[] = max(H0- y,0) ;
    // the boundary conditions for the pressure need to be handled by
    // the Navier--Stokes solver
    foreach()
    In[]=0;
    
    foreach()
    phi[]=0;
    /**
     initial for pressure darcy */
    
    foreach(){
        pDarcy[] = P1;
        sourceDarcy[] = 0.;
        ug.x[]=0.;
        ug.y[]=0.;
    }
    boundary ({pDarcy, sourceDarcy});
}

/**
 total density */

//#define rho(f) ((f) + RHOF*(1. - (f)))

/**
 Viscosity computing $D_2=D_{ij}D_{ji}$;
 
 In the pure shear flow
 $D_{11}=D_{22}=0$ et $D_{12}=D_{21}=(1/2)(\partial u/\partial y)$,
 so that
 $D_2=\sqrt{D_{ij}D_{ij}} =\sqrt{ 2 D_{12}^2} = \frac{\partial u}{ \sqrt{2}\partial y}$.
 In a pure shear flow, $\partial u/\partial y= \sqrt{2} D_2$.
 The inertial number $I$ is $D \sqrt{2} D_2/\sqrt(p)$
 and $\mu = \mu_s+ \frac{\Delta \mu}{1+I/I_0}$
 the viscosity is $\eta = \mu(I)p/(\sqrt{2}D_2}$:
 
 note that if $\eta$ is too small an artificial small viscosity $\rho D \sqrt{gD}$
 is taken see Lagrée et al. 11 § 2.3.1. */

event properties (i++) {
    trash ({alphav});
    scalar eta[];
    foreach() {
        eta[] = mug;
        if (p[] > 0.) {
            double D2 = 0.;
            foreach_dimension() {
                double dxx = u.x[1,0] - u.x[-1,0];
                double dxy = (u.x[0,1] - u.x[0,-1] + u.y[1,0] - u.y[-1,0])/2.;
                D2 += sq(dxx) + sq(dxy);
            }
            if (D2 > 0.) {
                D2 = sqrt(2.*D2)/(2.*Delta);
                In[] = D2*D/sqrt(p[]);  // this D2 is sqrt(2) D2
                double muI = .4 + .28*In[]/(.4 + In[]);
                double etamin = sqrt(D*D*D);
                eta[] = max((muI*p[])/D2, etamin);
                eta[] = min(eta[],100);
            }
        }
    }
    boundary ({eta});
    scalar fa[];
    foreach()
    fa[] = (4.*f[] +
            2.*(f[-1,0] + f[1,0] + f[0,-1] + f[0,1]) +
            f[1,1] + f[-1,1] + f[1,-1] + f[-1,-1])/16.;
    boundary ({fa});
    foreach_face() {
        double fm = (fa[] + fa[-1])/2.;
        muv.x[] = (fm*(eta[] + eta[-1])/2. + (1. - fm)*mug);
        alphav.x[] = 1./( fm * (1. - 0. * In[]) + RHOF*(1. - fm));
    }
    boundary ((scalar *){muv,alphav});
}

/**
 convergence outputs */

void mg_print (mgstats mg)
{
    if (mg.i > 0 && mg.resa > 0.)
        fprintf (stderr, "#   %d %g %g %g\n", mg.i, mg.resb, mg.resa,
                 exp (log (mg.resb/mg.resa)/mg.i));
}

/**
 convergence stats */

event logfile (i += 1) {
    stats s = statsf (f);
    fprintf (stderr, "%g %d %g %g %g %g\n", t, i, dt, s.sum, s.min, s.max - 1.);
    mg_print (mgp);
    mg_print (mgpf);
    mg_print (mgu);
    fflush (stderr);
}
/**
Computation of Darcy Forchheimer pressure at each time step 
*/
event pressureDEF (i++)
{
    double B_l=1.0;
    double B_i=1.0;
    double R_f=1.0;//4.0*1.2/2500;
    scalar un[];
    vector uold[];
    double err=1; 
/** 
 Darcy Forchheimer with penalization to have a constant pressure in the gas ($1-f$). 
 grains correspond to   $f=1$.
 
With 
$R_f= {\rho_f} (1 + C \frac{\phi}{1-\phi} )* (f + (1-f)\frac{1}{\beta_m}) $, 
$B_l=[ \eta \beta_l (1-\phi)] f + [\frac{1}{\beta_m}](1-f)$
and
$B_i= [\eta \beta_l (1-\phi)] f$, 

and $\beta_m \gg 1$ the penalization parameter.
We have  to solve by projection
$$ R_f\frac{\partial \overrightarrow{u}}{\partial t}= -  \overrightarrow{\nabla} p - B_l  \overrightarrow{u} -B_i |\overrightarrow{u}|\overrightarrow{u}
\mbox{  avec } \overrightarrow{\nabla} \cdot \overrightarrow{u} = 0$$ 
 
 semi implicit formulation 

 $$ R_f \frac{\overrightarrow{u}^{n+1} -  \overrightarrow{u}^{n}}
 {\Delta t}= -  \overrightarrow{\nabla} p^{n+1} - B_l \overrightarrow{u}^{n+1} -
B_i |\overrightarrow{u}^{n}|\overrightarrow{u}^{n+1}
$$
the new velocity is 
 $$\overrightarrow{u}^{n+1}  =  \frac{ R_f \overrightarrow{u}^{n}- \Delta t \overrightarrow{\nabla} p^{n+1}}{ 
R_f +  \Delta t ( B_l  +  
B_i |\overrightarrow{u}^{n}| 
)}$$
with
$$\overrightarrow{\nabla} \cdot \overrightarrow{u}^{n+1} = 0$$ 
the problem reduces to the computation of 
 
  $$\overrightarrow{\nabla} \cdot (\beta p^{n+1} ) = S_{darcy}$$ 
  with 
$$
\beta = \frac{  1}{ R_f+  \Delta t ( B_l  +  
B_i |\overrightarrow{u}^{n}| ) }f
+   (\frac{\beta_m}{(R_f+ \Delta t)})(1 - f);
$$ 



*/
    foreach_face(){
        uold.x[] = ug.x[];
        un[] = norm(ug);
        beta.x[] =1.0/(R_f+dt*(B_l+B_i*un[]))*f[]  + (betam/(R_f+dt))*(1. - f[]);
    }
    boundary ((scalar *){beta});
/**
computation of source term
$$
\overrightarrow{u}_s =  (
\frac{R_f \overrightarrow{u}^{n} }{ 
R_f +  \Delta t ( B_l  +  
B_i |\overrightarrow{u}^{n}| )})f+
(\frac{R_f \overrightarrow{u}^{n} }{ 
R_f +  \Delta t    })(1-f)$$

*/
    vector us[];
    foreach_face(){
        us.x[]=R_f*ug.x[]/(R_f+dt*(B_l+B_i*un[]))*f[] +   R_f*ug.x[]/(R_f+dt)*(1-f[]);
    }
    boundary ((scalar *){us});
 /**
 $$  S_{darcy}= \overrightarrow{\nabla} \cdot \overrightarrow{u}_s/\Delta t$$
 */
    scalar div[];
    foreach(){
        div[]=0.0;
        foreach_dimension()
           div[] +=us.x[1]-us.x[];
        div[] /=dt*Delta;
        }     
    sourceDarcy=div;            
 /**
solve $\nabla \cdot (\beta \nabla p_{darcy} )= S_{darcy}$ 
with [Poisson solver](http://basilisk.fr/src/poisson.h)
 */   
    mgpDarcy = poisson (pDarcy, sourceDarcy, beta);
  
/**

update
 $$\overrightarrow{u}^{n+1}  =  \frac{ R_f \overrightarrow{u}^{n}- \Delta t \overrightarrow{\nabla} p^{n+1}}{ 
R_f +  \Delta t ( B_l  +  
B_i |\overrightarrow{u}^{n}| 
)}$$

*/
    foreach() {
        un[] = norm(ug);
        foreach_dimension(){
           gradpDarcy.x[] =(pDarcy[0,0] - pDarcy[-1,0])/Delta;
           ug.x[]= (R_f*uold.x[] - gradpDarcy.x[]*dt)/(R_f +dt*( B_i * un[]+B_l))*f[]+
                    (R_f/betam*uold.x[] -gradpDarcy.x[]*dt) / ((R_f +dt)/betam)*(1-f[]) ;
        }
    }
    boundary ((scalar *){ug,gradpDarcy});
 
 /** monitoring   unsteadiness */   
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
    av.y[] += - 1. - eps*(pDarcy[0,0] - pDarcy[0,-1])/Delta;
    foreach_face(x)
    av.x[] += - eps*(pDarcy[0,0] - pDarcy[-1,0])/Delta ;
    boundary ((scalar *){av});
}

/**
 save velocity for analytical solution */

event sauf (t += 0.1) {

  FILE * fp = fopen("uprof.txt","a");
  fprintf (fp,"%g %g   \n", t,interpolate(ug.y, 3.5, 2));
}      


event interface (t = 0 ; t += 1. ; t <= tmax) {
#if dimension == 2
    output_facets (f, fpf);
#endif
    char s[80];
    sprintf (s, "field-%g.txt", t);
    
  

    
    FILE * fp = fopen (s, "w");
    output_field ({f,p,u,uf,pf,pDarcy,gradpDarcy,In,ug}, fp, linear = true);
    fclose (fp);
}

/**
 Rate of flowing materials across the hole */

event debit (t = 0 ; t += 0.1 ; t <= tmax) {
    static double V = 1;
    V = 0;
    foreach()
    V = V + f[]* Delta * Delta;
    if (t >= 0.) fprintf (stdout,"%lf %lf %lf\n",t,V,ugranular);
    fflush (stdout);
}


/**
 Calculate velocity of flowing materials in the bulk and add it to Qair at the top
 to simulate the counter flow */

event velocitybulk (i++) {
    static double VVold, VV = 3.9,Qinst = 0;
    VVold = VV;
    VV = 0;
    foreach()
    VV = VV + f[]* Delta * Delta;
    Qinst = -(VV-VVold)/dt;
    ugranular = Qinst/(LDOMAIN*(1-maskr));
}

/**
 film output */

#if 1
event movie (t += 0.05) {
    static FILE * fp1 = popen ("ppm2mpeg > level.mpg", "w");
    scalar l[];
    foreach()
    l[] = level;
    boundary ({l});
    output_ppm (l, fp1, min = 0, max = LEVEL,
                n = 512, box = {{0,0},{Lz,Lz}});
    
    foreach()
    l[] = f[]*(1 + sqrt(sq(u.x[]) + sq(u.y[])));
    boundary ({l});
    
    static FILE * fp2 = popen ("ppm2mpeg > velo.mpg", "w");
    output_ppm (l, fp2, min = 0, max = 2., linear = true,
                n = 512, box = {{0,0},{Lz,Lz}});
    
    foreach()
      l[] =  sqrt(sq(ug.x[]) + sq(ug.y[]));
    boundary ({l});
    static FILE * fp4 = popen ("ppm2mpeg > velo_f.mpg", "w");
    output_ppm (l, fp4, min = 0, max = 2., linear = true,
                n = 512, box = {{0,0},{Lz,Lz}});
    
    static FILE * fp3 = popen ("ppm2mpeg > pDarcy.mpg", "w");
    foreach()
    l[] = pDarcy[];
    boundary ({l});
    output_ppm (l, fp3, min = 0, linear = true,
                n = 512, box = {{0,0},{Lz,Lz}});
}

event pictures (t==3) {
    output_ppm (f, file = "f.png",
                min = 0, max = 2,  spread = 2, n = 512, linear = true,
                box = {{0,0},{2,2}});
}
#endif


#if 0
event gfsview (i++)
{
    static FILE * fp = popen ("gfsview2D -s", "w");
    output_gfs (fp, t = t);
}
#endif

/**
# Run
 
 to run
 
~~~bash
 qcc -g -O2 -Wall -o forchheimer forchheimer.c -lm
 ./forchheimer> out
~~~


# Results

~~~gnuplot   pressure
set pm3d map
set palette rgbformulae 22,13,-31;
unset colorbox
set xlabel "x  iso p"
splot [][:] 'field-1.txt' u 1:2:10   not
reset
~~~
 

~~~gnuplot   pressure
  set ylabel "p_darcy"
set xlabel "y"
 p 'field-1.txt' u 2:10 
~~~
 

Analytical solution of the problem, $\frac{\partial p}{\partial x} = 2 $

$$ \frac{\partial u}{\partial t}=  -2 - u + u^2$$
  $$u(t)= \frac{2 (exp(3 t )-1)}{(1+2 exp(3 t))}$$

~~~gnuplot 
set xlabel "t"
set ylabel "u(t)"
 p'uprof.txt' t'num.' ,-2*(exp(3*x)-1)/(1+2*exp(3*x)) t'analytic'

~~~

 
# Bibliography 
 

*/
