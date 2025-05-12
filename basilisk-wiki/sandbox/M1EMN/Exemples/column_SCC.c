/**

# Collapse of a Bingham flow


## Physical problem
First, we consider simple viscous collapse of a Nextonian fluid
 and find the $t^{1/5}$ selfsimilar solution. To be compared with all the 1D models...
 
 
 Then wew consider   collapse of a Non-Nextonian fluid : the Bingham fluid.
 
 
Application to mud flows,debris flows etc
application to wet concrete flow: collapse of columns (Abrahams cone test).



## equations
 We propose an implementation of the Bingham rheology.
 For those flows, when stress is larger than a yield value, the media flows.
 For non Newtonian fluids, the "viscosity" is at least a function of  the second principal invariant of the shear strain rate tensor that we define here as (other definitions proportional to this one are possible):
 $$D_2=\sqrt{\sum_{i,j}D_{ij}D_{ij}}$$
 In littérature two norms are used, the Euclidian:
 $$||D||=\sqrt{\frac{1}{2}\sum_{i,j}D_{ij}D_{ij}}$$
 and the Frobenius one:
 $$|D|=\sqrt{\sum_{i,j}D_{ij}D_{ij}}$$
 Obviously  $$||D||=\sqrt{\frac{1}{2}} |D|= \frac{D_2}{\sqrt{2}}$$
 The strain rate tensor is $D_{ij}=(\partial_iu_j+\partial_ju_i)/2$, it has unit of $\text{s}^{-1}$.
 the components in 2D:
 
 $D_{11}=\frac{\partial u}{\partial x}$,
 $D_{12} =\frac{1}{2}( \frac{\partial u}{\partial y}+ \frac{\partial v}{\partial x})$,
 $D_{21} =D_{12} =\frac{1}{2}( \frac{\partial u}{\partial y}+ \frac{\partial v}{\partial x})$,
 $D_{22}=\frac{\partial v}{\partial y}$
 
 And where the second invariant is $D_2=\sqrt{D_{ij}D_{ij}}$
 $$D_2^2= D_{ij}D_{ij}= D_{11}D_{11} + D_{12}D_{21} + D_{21}D_{12} + D_{22}D_{22}$$
 
 
 hence, as $D_{12}D_{21} = D_{21}D_{12}$:
 $$D_2^2= D_{ij}D_{ij}= ( \frac{\partial u}{\partial x})^2 + (\frac{\partial v}{\partial y})^2 +  \frac{1}{2}( \frac{\partial u}{\partial y}+ \frac{\partial v}{\partial x})^2$$
 
 We have by definition of the Bingham rheology
 $$\tau_{ij} = 2 \mu_1 D_{ij} + \tau_y \frac{D_{ij}}{||D||}$$
 (hence if $||\tau_{ij}|| > \tau_y$ there is a flow).
 In order to identify this as an effective viscosity, $\tau_{ij} = 2 \eta_{eft} D_{ij}$,
 
 $$\tau_{ij} = 2 \mu D_{ij} + \tau_y \frac{D}{||D||} = 2 ( \mu +  \frac{\tau_y}{2||D||})D_{ij}$$
 then the equivalent, or effective, of apparent viscosity is with $D_2= \sqrt{2}||D||$
 $$\mu_{eq}= \mu_1 + \frac{\tau_y}{\sqrt{2} D_2 }$$
 
 
 We can use a general Herschel-Bulkley formulation, with here $N=1$,
  $$ \tau = \tau_y + \mu_N \dot{\gamma}^{N}$$
 which is coded with an effective  viscosity of the form:
 
 
 $$\mu(||D||)= \mu_N(2||D||)^{N-1} + {\frac{\tau_y}{2||D||}}.$$
 
 $$ (|| {\dot  \gamma} ||^2)^{(N-1)/2}  
 =  (\sqrt{2} D_2)^{(N-1)}  =
 =  (({2} ||D||)^{(N-1)} )  $$
 
 Where $\tau_y$ is the  yield stress and $\mu_N$ the generalized viscosity.
 
 

check with [bingham simple](http://basilisk.fr/sandbox/M1EMN/Exemples/bingham_simple.c), it should be the same.
 
 
 
## Adimensionalisation
 
 The Navier Stokes equations with Herschel-Bulkley/ Bingham are:
 
 $$\rho _0 \frac{d  u }{d t } =  -    \nabla p -  \rho_0   g  +   \nabla \cdot    \tau $$
 
 we use the effective viscosity
 $\tau_{ij} =  2\mu_{eff} D_{ij}$
 with
 $$ 2\mu_{eff}  D_{ij}   =  (\frac{\tau_y}{ ||D||}+ 2 \mu_N||D||^{N-1})  D_{ij}$$
 so
 $$  \tau_{ij} =  2  (\frac{\tau_y}{2||D||}+\mu_N||D||^{N-1}) D_{ij}  $$
 Using $x=L \bar x$, $y=L \bar y$, and $(u,v)=U_0(\bar u,\bar v)$
 we have $p= \rho_0 U_0^2 \bar p$
 let us write the viscous part:
 
 $$ \tau =  \rho_0 U_0^2 \left[
 (\dfrac{\tau_y}{ 2 \rho_0 U_0^2 || \bar D||}+\frac{\mu_N }{\rho_0 U_0L} (\frac{U_0}{L})^{N-1}||\bar D||^{N-1}) 2  {\bar D}
 \right]
 $$
 or finally
 $$ \tau =  \rho_0 U_0^2 \left[
 (\dfrac{\bar \tau_y}{ 2  ||\bar D||}+\frac{1 }{Re_N}||\bar D||^{N-1})( 2  {\bar D} )
 \right]
 $$
 with the yield stress without dimension $\bar \tau_y= \dfrac{\tau_y}{  \rho_0 U_0^2}$ and a kind of Reynolds number
 $Re_N=\frac{\rho_0 U_0L} {\mu_N }(\frac{U_0}{L})^{1-N}$.
 
 
 With these two parameters the  Navier Stokes Herschell-Bulkley is:
 $$\frac{d \bar u }{d \bar t } =  -  {\bar \nabla} p  -
 \frac{gL}{U_0^2}
 e_y
 +  {\bar \nabla}   \left[
 (\dfrac{\bar \tau_y}{ 2  ||\bar D||}+\frac{1 }{Re_N}||\bar D||^{N-1})(2  {\bar D} )
 \right]
 $$
 
 In case of Bingham flow $N=1$, the Reynolds is
 $Re_1=\frac{\rho_0 U_0L} {\mu_1 }$
 and we define the Bingham number $Bi=(\dfrac{\tau_y L}{  \mu_1 U_0})$
 with these two parameters the  Navier Stokes Bingham is:
 $$\frac{d \bar u }{d \bar t } =  -   {\bar \nabla} p  -
 \frac{gL}{U_0^2}
 e_y
 +
 \frac{1 }{Re_1}
 {\bar \nabla}   \left[
 (\dfrac{Bi}{ 2  ||\bar D||}+1)(2  {\bar D} )
 \right]
 $$
 
 
 This is a collapse so that it is natural to put one in front of the gravity: $\frac{gL}{U_0^2}=1$. This means that the scales of length and time  $L/T^2=g$, then $T=\sqrt{L/g}$.
 The velocity scale is $U_0=L/T$ which is $\sqrt{gL}$, the gravitational term is 1, the viscous term:
 $$\frac{1}{Re_1}= \frac{\nu}{\sqrt{g L^3}}$$

 
# Code
 
 Includes and definitions
 */
//#include "grid/cartesian.h"
//#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "vof.h"
#include "heights.h"
// Domain extent
#define LDOMAIN 4.000001
//passive fluid
#define RHOF 1e-3
#define mug  1e-3
double mub,Bi;
double tmax;
double etamx;
// Maximum refinement 4./2**8 = 0.015625
// Minimum refinement 2./2**4 = 0.125
#define LEVEL 7 //7 // 9 OK
#define LEVELmin 2
scalar f[];
scalar * interfaces = {f};
scalar eta[];
/**
 Boundary conditions
 */
 u.t[bottom] = dirichlet(0);

/**
 Boundary conditions
*/
int main() {
    system("mkdir SIM");
    L0 = LDOMAIN;
    // number of grid points
    N = 1 << LEVEL;
    // maximum timestep
    DT = 0.025/4 ;
    TOLERANCE = 1e-3;
    
    /**
     Numerical values
     
    The problem is a heap of viscoplastic material of height = length
     released on a flat plate.
     We may define a cone for Abrams slup test (concrete)
     */
    // Initial conditions
#define Ltas 0.9999999999
#define Ls 0.9999999999
#define Htas 0.99999999999

/**
     
The fluid has a density $\rho_0$, time is $T=\sqrt{L/g}$,
We define a quantity, $\mu_b$ which is in fact the inverse of
$Re_1=\frac{\rho_0 U_0L}{\mu_N }$     and wich is here 1.
and we define the Bingham number $Bi=(\frac{\tau_y L}{\mu_N U_0})$, we change it.
 
*/
    
    mub= 1;
    etamx=1000.;
    //tmax = 199.99;
    tmax = 100;
    DT = 0.025;
    Bi = .0;
    run();
    system("cp velo.mp4 veloo.mp4");
  //  system("cp  img.png imgo.png");
   
    mub= 1;
    tmax=50;
    etamx=1000.;
    DT = 0.01 ;
    Bi =  .05*2/sqrt(2);
    run();
    system("cp velo.mp4 veloi.mp4");
    system("cp mu.mp4 mui.mp4");
   // system("cp  img.png imgi.png");


    mub= 1;
    tmax=50;
    etamx=1000.;
    DT = 0.001 ;
    Bi = .14*2/sqrt(2);
    run();
    system("cp velo.mp4 veloii.mp4");
    system("cp mu.mp4 muii.mp4");
   // system("cp  img.png imgii.png");

}
face vector alphav[];
face vector muv[];
scalar rhov[];

event init (t = 0) {
    const face vector g[] = {0.,-1.};
    a = g;
    alpha = alphav;
    mu = muv;
    rho = rhov;
/**
     Defining the initial trapozeoidal shape of the Abrams cone.
*/
    fraction (f, min(Htas - y, Ltas - (fabs(x)+ y * (Ltas-Ls))));
}
event stop (t = tmax);
/**
 total density
 */
#define rho(f) ((f) + RHOF*(1. - (f)))
/**
 Viscosity computing:
 $D_{ij}=(u_{i,j}+u_{j,i})/2$, incompressibility
 $D_{ii}=0$.
 The second invariant $D_2$ written `D2` in Gerris and Basilisk
 $D_2=\sqrt{D_{ij}D_{ij}}$
 in 2D:
 $D_{ij}$ is
 $D_{11}=\frac{\partial u}{\partial x}$,
 $D_{12} =\frac{1}{2}( \frac{\partial u}{\partial y}+ \frac{\partial v}{\partial x})$,
 $D_{21} =D_{12} =\frac{1}{2}( \frac{\partial u}{\partial y}+ \frac{\partial v}{\partial x})$, $D_{22}=\frac{\partial v}{\partial y}$
 Then
 $$D_{ij}D_{ij}=
 D_{11}D_{11}
 +D_{12}D_{21} +
 D_{21}D_{12} + D_{22} D_{22}
 $$
 or:
 $$D_{ij}D_{ij}= ( \frac{\partial u}{\partial x})^2 + (\frac{\partial v}{\partial y})^2 +  2(\frac{1}{2}( \frac{\partial u}{\partial y}+ \frac{\partial v}{\partial x}))^2$$
 or:
 $$D_{ij}D_{ij}=( ( \frac{\partial u}{\partial x})^2 +
 (\frac{1}{2}( \frac{\partial u}{\partial y}+ \frac{\partial v}{\partial x}))^2)
 +(
 (\frac{\partial v}{\partial y})^2 +   (\frac{1}{2}( \frac{\partial u}{\partial y}+ \frac{\partial v}{\partial x}))^2)$$
 as it is coded:
 
 dxx=u.x[1,0] - u.x[-1,0] is $(\frac{\partial u}{\partial x})(2 \Delta x)$
 
 dxy is $\frac{1}{2}(\frac{\partial u}{\partial y} + \frac{\partial v}{\partial y})(2 \Delta x)$
 then add the second dimension.
 
 This gives the equivalent final viscosity
 $$\eta_{eq}=(Bi/(\sqrt{2} D_2) + 1) \mu_b$$
 
 note that we do a regularisation such as for too small $D_2$, the viscosity is the viscosity of a very viscous newtonian fluid
 $$\eta_{eq}= \eta_{max}$$
 
 there are several ways to do that:
 $$\eta_{eq} = min ( (Bi/(\sqrt{2} D_2) + 1) \mu_b ,\eta_{max})$$
 or
 $$\eta_{eq} = ( \frac{Bi}{(\sqrt{2} D_2)  + \mu_1 Bi/\eta_{max} } + 1) \mu_1  $$
 
 */
event properties (i++) {
    trash ({alpha});
    foreach() {
        eta[] = etamx;
        double D2 = 0.;
        foreach_dimension() {
            double dxx = u.x[1,0] - u.x[-1,0];
            double dxy = (u.x[0,1] - u.x[0,-1] + u.y[1,0] - u.y[-1,0])/2.;
            D2 += sq(dxx) + sq(dxy);
        }
        if (D2 > 0.) {
            D2 = sqrt(D2)/(2.*Delta);
            eta[] = (Bi/(sqrt(2)*D2 + mub*Bi/etamx) + 1)*mub ;
           // eta[] = min(eta[],etamx);
        }
    }
    boundary ({eta});
    scalar fa[];
    // filtering density twice
    foreach()
    fa[] =  (4.*f[] +
             2.*(f[-1,0] + f[1,0] + f[0,-1] + f[0,1]) +
             f[1,1] + f[-1,1] + f[1,-1] + f[-1,-1])/16.;
    boundary ({fa});
    
    foreach_face() {
        double fm = (fa[] + fa[-1,0])/2.;
        muv.x[] =  (fm*(eta[] + eta[-1,0])/2. + (1. - fm)*mug);
        //        muv.x[] = 1./(2.*fm/(eta[] + eta[-1,0]) + (1. - fm)/mug);
        alphav.x[] = 1./rho(fm);
    }
    foreach()
    rhov[] = rho(fa[]); 
    boundary ({muv,alphav,rhov});
}
/**
 convergence outputs
 */
void mg_print (mgstats mg)
{
#if 0
    if (mg.i > 0 && mg.resa > 0.)
        fprintf (stderr, "#   %d %g %g %g\n", mg.i, mg.resb, mg.resa,
                 exp (log (mg.resb/mg.resa)/mg.i));
#else
    if (mg.i > 0 && mg.resa > 0.)
        fprintf (stderr, "#   %d %g %g %g\n", mg.i, mg.resb, mg.resa,
                 mg.resb > 0 ? exp (log (mg.resb/mg.resa)/mg.i) : 0);
#endif
}
/**
 convergence outputs
 */
event logfile (i++) {
    stats s = statsf (f);
    fprintf (stderr, "%g %d %g %g %g %g\n", t, i, dt, s.sum, s.min, s.max - 1.);
    mg_print (mgp);
    mg_print (mgpf);
    mg_print (mgu);
    fflush (stderr);
}
/** 
 some outputes for comaprison
*/
event pij0 (t = 0.) {
    char s[80];
    sprintf (s, "out0-%g",Bi);
    FILE * fp = fopen (s, "w");
    output_facets (f, fp);
    fclose (fp);
}
event pij2 (t = 2.5001234) {
    char s[80];
    sprintf (s, "out1-%g",Bi);
    FILE * fp = fopen (s, "w");
    output_facets (f, fp);
    fclose (fp);
}
event pij5 (t = 4.9876) {
    char s[80];
    sprintf (s, "out2-%g",Bi);
    FILE * fp = fopen (s, "w");
      output_facets (f, fp);
    fclose (fp);
}
event pij10 (t = 9.9867212){
    char s[80];
    sprintf (s, "out3-%g",Bi);
    FILE * fp = fopen (s, "w");
    output_facets (f, fp);
    fclose (fp);
}
event pij40 (t = 39.9867212){
    char s[80];
    sprintf (s, "out4-%g",Bi);
    FILE * fp = fopen (s, "w");
    output_facets (f, fp);
    fclose (fp);
}
event pij50 (t = 49.9867212){
    char s[80];
    sprintf (s, "out5-%g",Bi);
    FILE * fp = fopen (s, "w");
    output_facets (f, fp);
    fclose (fp);
}
event pijtmax (t = tmax){
    char s[80];
    sprintf (s, "out6-%g",Bi);
    FILE * fp = fopen (s, "w");
    output_facets (f, fp);
    fclose (fp);
    
}
#if 1
event interface (t +=1 ) {
    char s[80];
    sprintf (s, "SIM/field-%g", t);
    FILE * fp = fopen (s, "w");
    output_field ({f,p,u,uf,pf,eta}, fp, linear = true);
    fclose (fp);
#endif
}

#if gfsv
event gfsview (i += 10){
    static FILE * fp = popen ("gfsview2D -s column_SCC.gfv", "w");
    output_gfs (fp);
}
#endif
/**
 Saving the top position and runout positionfor slump measurements
*/
vector h[];
event timeseries (t += 0.1 ) {
    heights (f, h);
    double maxy = - HUGE,maxx = - HUGE;;
    foreach()
    if ((h.y[] != nodata) && (h.x[] != nodata)) {
        double yi = y + height(h.y[])*Delta;
        double xi = x + height(h.x[])*Delta;
        if (yi > maxy)
            maxy = yi;
        if (xi > maxx)
            maxx = xi;
    }
    char s[80];
    sprintf (s, "hmax-%g",Bi);
    static FILE * fp0 = fopen (s, "w");
    fprintf (fp0, "%g %g %g\n", t, maxx, maxy);
    fflush (fp0);
}
/**
 Movies
 */
#if 1
event pictures (t = {0, 1. , 2., 3., 4.} ) {
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
            output_ppm (f, file = "f3.png", min=0, max = 2,  spread = 2, n = 256, linear = true,
                        box = {{0,0},{4,2}});
    foreach()
    l[] = level;
    boundary ({l});
            output_ppm (f, file = "l3.png", min=0, max = 8,  spread = 2, n = 256, linear = true,
                    box = {{0,0},{4,2}});
}
#endif

#if 0
event movie (t += 0.05) {
    scalar l[];
    static FILE * fp1 = popen ("ppm2mpeg > velo.mpg", "w");
    foreach()
    l[] = f[]*(0.025+sqrt(sq(u.x[]) + sq(u.y[])));
    output_ppm (l, fp1, min = 0, max = 0.1,
                n = 1024, box = {{0,-1},{LDOMAIN-.5,1.5}});
    fflush (fp1);

    static FILE * fp2 = popen ("ppm2mpeg > mu.mpg", "w");
    foreach()
    l[] = f[]*muv.x[];
    boundary ({l});
    output_ppm (l, fp2, min = 0, max = 5, linear = true,
                 n = 1024, box = {{0,-1},{3,1.5}});
    fflush (fp2);
  
    static FILE * fp3 = popen ("ppm2mpeg > f.mpg", "w");
    foreach()
    l[] = f[];
    output_ppm (l, fp3, min = 0, max = 2, linear = true,
                n = 1024, box = {{0,-1},{LDOMAIN,1.5}});
     fflush (fp3);
}
#else
event movie (t += 0.05) {
    scalar l[];
    foreach()
    l[] = level;
    boundary ({l});
    output_ppm (l, file = "level.mp4", min = 0, max = LEVEL,
                n = 1024, box = {{0,-1},{4,3.5}});

    foreach()
    l[] = f[]*(1+sqrt(sq(u.x[]) + sq(u.y[])));
    boundary ({l});
    output_ppm (l, file = "velo.mp4", min = 0, max = 1.5, linear = true,
                n = 1024, box = {{0,-1},{4,1.5}});
   
    foreach()
    l[] = f[]*muv.x[];
    boundary ({l});
    output_ppm (l, file = "muv.mp4", min = 0, max = 5, linear = true,
                n = 1024, box = {{0,-1},{4,1.5}});
    foreach()
    l[] = f[]*p[];
    boundary ({l});
    output_ppm (l, file = "f.mp4", min = 0, max =1, linear = true,
                n = 2048, box = {{0,-1},{4,1.5}});
}
#endif

#ifdef gnuX
event output (t += .1 ) {
    fprintf (stdout, "p[0:4][0:2]  '-' u 1:2    not  w   l \n");
    output_facets (f, stdout);
    fprintf (stdout, "e\n\n");
}
#endif

#if QUADTREE
// if no #include "grid/multigrid.h" then adapt
event adapt (i++) {
scalar K1[],K2[];
foreach()
K1[]=(f[0,1]+f[0,-1]+f[1,0]+f[-1,0])/4*noise();
boundary({K1});

for(int k=1;k<8;k++)
{
    foreach()
    K2[]=(K1[0,1]+K1[0,-1]+K1[1,0]+K1[-1,0])/4;
    boundary({K2});
    foreach()
    K1[]=K2[]*noise();
}
adapt_wavelet({K1,f},(double[]){0.001,0.01}, maxlevel = LEVEL, minlevel = LEVELmin);
//adapt_wavelet ({f,u,muv.x}, (double[]){5e-3,0.02,0.02,0.01}, LEVEL, LEVELmin,list = {p,u,pf,uf,g,f});
}
#endif


/**
# Run
 
 to run
 
~~~bash
 qcc -g -O2 -DTRASH=1 -Wall -DgnuX=1  -o column_SCC column_SCC.c -lm
 qcc -g -O2 -Wall -DgnuX=1  -o column_SCC column_SCC.c -lm
 ./column_SCC | gnuplot
 
 qcc -g -O2 -Wall -Dgfsv=1  -o column_SCC column_SCC.c -lm
  ./column_SCC


make column_SCC.tst;make column_SCC/plots
make column_SCC.c.html ; open column_SCC.c.html
~~~
 
 
# Results
 
## Exemples of  viscous collapse

we first revover the Huppert collapse for newtonian fluid, $Bi=0$, in $h \sim t^{-1/5}$ and $x \sim t^{1/5}$

~~~gnuplot  height function of time
 reset
 set logscale
 set xlabel "t"
 set ylabel "h_m(t),x_m(t)"
 p'hmax-0' u 1:2 w l t 'xm','' u 1:3 w l t 'hm',x**(1./5) t't^{1/5}',x**(-1./5) t'1/t^{1/5}'

~~~
 
 
~~~gnuplot
reset
 p'out2-0' u ($1/5**(1./5)):($2*5**(1./5)) t 't=2.5' w l,\
 'out3-0' u ($1/10**(1./5)):($2*10**(1./5))t 't=10' w l,\
 'out4-0' u ($1/40**(1./5)):($2*40**(1./5)) t 't=40'w l,\
 'out5-0' u ($1/50**(1./5)):($2*50**(1./5)) t 't=50'w l,\
 'out6-0' u ($1/200**(1./5)):($2*100**(1./5)) t 't=100' w l
 ~~~
 
## Exemples of  Bingham collapses
 
We next look at Bingham cases.
 First the height and the run out as function of time
 
~~~gnuplot  height function of time
 reset
 set xlabel "t"
 set ylabel "h_m(t),x_m(t)"
 p'hmax-0.0707107' u 1:2 w l t 'Bi=0.07','' u 1:3 w l t 'Bi=0.07',\
  'hmax-0.19799' u 1:2 w l t 'Bi=0.19799','' u 1:3 w l t 'Bi=0.19799'
 
~~~
 
Second compared to Liu et al, note that the definition of $Bi$ contains a $\sqrt{2}$.
 

 
   
~~~gnuplot comparisons
 reset
 set size ratio -1
 unset tics
 p [1:1494][1:632]'./Img/liu16.png' binary filetype=png   with rgbimage not,\
 'out0-0.0707107' u ($1*(376-121)+121):($2*(400-166)+166) not  w l linec 1,\
 'out1-0.0707107' u ($1*(376-121)+121):($2*(400-166)+166) not w l linec 1,\
 'out2-0.0707107' u ($1*(376-121)+121):($2*(400-166)+166) not w l linec 1,\
 'out3-0.0707107' u ($1*(376-121)+121):($2*(400-166)+166) not w l linec 1,\
 'out4-0.0707107' u ($1*(376-121)+121):($2*(400-166)+166) not w l linec 1,\
 'out3-0.19799' u ($1*(1082-829)+829):($2*(400-166)+166) not w l linec 1,\
 'out4-0.19799' u ($1*(1082-829)+829):($2*(400-166)+166) not w l linec 1,\
 'out5-0.19799' u ($1*(1082-829)+829):($2*(400-166)+166) not w l linec 1

 
~~~
 
 
 
Some  Films
 
 
 [![](./column_SCC/f0.png)](./column_SCC/velo.mp4)
 velocity  (click on image for animation)
 
 [![](./column_SCC/f1.png)](./column_SCC/velo.mp4)
 velocity  (click on image for animation)
 
 [![](./column_SCC/f0.png)](./column_SCC/velo.mp4)
 velocity  (click on image for animation)
 
 [![](./column_SCC/f0.png)](./column_SCC/velo.mp4)
 velocity  (click on image for animation)

 [![](./column_SCC/f0.png)](./column_SCC/muv.mp4)
 viscosity (click on image for animation)
 
 [![](./column_SCC/l3.png)](./column_SCC/level.mp4)
level (click on image for animation)
 
 
# Links
 
 * 1D viscous collapase
 * Multilayer viscous collapse
 * 1D Bingham collapse
 * Multilayer Bingham collapse
 * granular collapse
 
 
 
# Bibliography
 
 * Lagrée, Staron, Popinet
 ["Scaling Laws for the Slumping of a Bingham Plastic Fluid" ](http://www.lmm.jussieu.fr/~lagree/TEXTES/PDF/JOR001265.pdf)
 J. Rheol. 57, 1265 (2013); doi: 10.1122/1.4802052
 
 * Dufour and  Pijaudier-Cabotz "Numerical modelling of concrete flow: homogeneous approach"
 Int. J. Numer. Anal. Meth. Geomech., 2005; 29:395–416
 
 * [gerris](http://gerris.dalembert.upmc.fr/gerris/tests/tests/couette.html).
 
 
 * Y. Liu, N.J. Balmforth, S. Hormozi, D.R. Hewitt, [Two–dimensional viscoplastic dambreaks](https://www.math.ubc.ca/~njb/Research/leoslump.pdf), Journal of Non-Newtonian Fluid Mechanics, Volume 238, December 2016, Pages 65-79,
 
 
 *  Guillaume Vinay Anthony Wachs, Jean-François Agassant, 
 Numerical simulation of non-isothermal viscoplastic waxy crude oil flows
 J. Non-Newtonian Fluid Mech. 128 (2005) 144–162
 
*/
