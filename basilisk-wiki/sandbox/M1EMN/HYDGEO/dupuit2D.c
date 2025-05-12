/** 
# Problem of the flow in a porous media
 

 What is the unsteady flow in a porous media not fully saturated?
 This is the problem of the "aquifère" see Guerin chap 2 without rain.
 
~~~gnuplot
reset
set style fill transparent solid 0.5 noborder
set style function filledcurves y1=0
unset colorbox
set lmargin 6
set xlabel "x"
set label "h(x,t)" at 5,1.25
set arrow 1 from -.1,0 to -.1,1
set label "h1" at -.3,0.75
 
set arrow from 10.1,0 to 10.1,.4
set label "h2" at 10.15,0.45
plot [0:10][0:2] 1-exp(x-10.5) fs solid 1.0 lc rgb "blue" not
 unset arrow 1 

~~~

We have the darcy law with $\beta = 1/B$, the permeability :
$$ 0 = -  \overrightarrow{\nabla} p -  B   \overrightarrow{u}  + \rho  \overrightarrow{g} 
\text{  \; with } \overrightarrow{\nabla} \cdot \overrightarrow{u} = 0$$ 
 
 the water is in a box and goes out in the bottom right corner, the water comes from the left.
 The height of water as function of time and space is unknown,
 but they are imposed left, where $h=h_1$ and right $h=h_2$.

 
## the problem without dimension:

Solution of the "Basic problem" in 2D, potential incompressible flow
$$ \bar u =- \frac{\partial \bar p}{\partial \bar x}, \;\; \bar v = - \frac{\partial \bar p}{\partial \bar y} - 1,\;\;\; 
  \frac{\partial \bar u}{\partial \bar x} +  \frac{\partial \bar v}{\partial \bar y}=0$$
  to be solved as:
  $$\frac{\partial^2 \bar p}{\partial \bar x^2} +  \frac{\partial^2 \bar p}{\partial \bar y^2}=0$$
  in   space $\bar h(\bar x, \bar t ) >\bar y>0$ and $0 < \bar x < L_0$ 
 at the wall $\bar y=0$ velocity slips.
 
 The shape of the aquifer $\bar h(\bar x, \bar t )$ has to be found. At time zero $\bar h(\bar x, \bar t=0 )=1$. 
  
## the problem we solve here:

We introduce a penalization parameter: $\beta_m \gg 1$ 
We have to solve by "projection" the problem without dimension
 $$ 0 = -  \overrightarrow{\nabla} \bar p -   \frac{1}{\beta} \overrightarrow{\bar u} - \bar{\rho} \overrightarrow{e_y}
\text{  with } \overrightarrow{\bar \nabla} \cdot \overrightarrow{\bar u} = 0$$ 
 with $\beta = \beta_m \gg 1$ (100 in practice) and $\bar{\rho}=\bar{\rho_{sec}}\ll 1$ (.005 in practice) outside the water, and $\beta =1$  $\bar{\rho}=1$ in the soil with the water.


in the code we remove the bars

# Code

includes for advection of a field, VOF, running tools and Poisson solver 
*/
#include "advection.h"
#include "vof.h"
#include "run.h"
#include "poisson.h"

#define MAXLEVEL 8 //  11
#define MINLEVEL 5  //

scalar p[], source[];
double Q;   
double H0=2;
double H1=2;
double H2=1;
double rhosec=.005;
double betam=100;
double tmax;
face vector beta[],rhop[];
mgstats mgp;

scalar f[];
scalar * interfaces = {f}, * tracers = NULL;

/** 

## boundary conditions

Domain of size $2L_0$, 
 at the bottom hydrostatic pressure gradient,
 left  $p=h_1-y$ , and right, $ $p=h_2-y$

*/ 
u.n[right]  =  neumann(0);
u.t[right]  =  neumann(0);
u.n[left]   =  neumann(0);
u.t[left]  =   neumann(0);
u.n[bottom]   = dirichlet(0);
u.t[bottom]  =   neumann(0);
u.n[top]   =  neumann(0);

p[left]  =dirichlet(max(0,H1-y));
p[right]  =dirichlet(max(0,H2-y));
p[top]   = dirichlet(0) ;
p[bottom] =  neumann(1);

f[bottom] =  neumann(0);
f[right]  = neumann(0);
f[left]  = neumann(0);
/**
 domain is `L0 x L0`
*/
int main()
{    
    L0=20.;
    Y0=0;
    X0=0;
    init_grid (1 << 4);
    DT = .05;
    tmax = 90;
    run();
}
/**
gravity and initial values
*/
 const face vector g[] = {0.,-1.};

 event init (i = 0) {
     fraction (f, ((x<L0/2?H0:H0) - y -x*0));
   foreach_face() {
        beta.x[] = 1*f[]+betam*(1-f[]); 
        rhop.x[] =   f[]+rhosec*(1-f[]);
    }  
    boundary  ((scalar *){beta}); 
    boundary  ((scalar *){rhop}); 
 
  foreach(){
    p[] = (y<H0? H0-y :0 );
    source[] = 0.;
  }
  boundary ({p});

// We may change the default gradient function (used for advection) to minmod-limited (rather than the centered default).
//  gradient = minmod2;
}
 
event pressurevelocities (i++)
{
/** 
ssmooth the interface
*/
   scalar fa[];
   foreach()
   fa[] =  (4.*f[] +
             2.*(f[-1,0] + f[1,0] + f[0,-1] + f[0,1]) +
             f[1,1] + f[-1,1] + f[1,-1] + f[-1,-1])/16.;
    boundary ({fa});
/**
  evaluate physcical quantities $\beta$ permeability and $\rho$ density as function of concentration
    with (remember $B=1$ by choice of the scales)     
and $\beta_m \gg 1$ the penalization parameter.
$$
\beta = \frac{  1}{    
B }f
+   ( \beta_m )(1 - f);
$$ 
and the density is 1 by adimensionalisation in the fluid and $\rho_{sec}$ in the dry porous solid 
$$
\rho =  f
+   ( \rho_{sec})(1 - f);
$$ 
*/    
  foreach_face() {
      double fm = (fa[] + fa[-1,0])/2.;
        beta.x[] = 1*fm+betam*(1-fm); 
        rhop.x[] = fm+rhosec*(1-fm);
    } 
    boundary  ((scalar *){beta}); 
    boundary  ((scalar *){rhop});  
/**

We have to solve 
$$ 0 = -  \overrightarrow{\nabla}   p -   \frac{1}{\beta} \overrightarrow{  u} - \overrightarrow{\rho e_y}
\text{  with } \overrightarrow{  \nabla} \cdot \overrightarrow{  u} = 0$$ 
 with $\beta = \beta_m$ (with $\beta_m \gg 1$) outside the water, and $\beta =1$ in the water,
 and by choice of scales $\overrightarrow g=- \overrightarrow{  e_y}$.      

this gives : 
$$0 = - \nabla \cdot (\beta \nabla p  ) +  \nabla \cdot (\beta \rho \overrightarrow g ) $$ 

We solve Poisson equation 
$\nabla \cdot (\beta \nabla p  )= s$ with the source term of the Poisson equation
$$ 
s= \nabla \cdot (\beta \rho \overrightarrow g )
$$
*/
scalar divr[];
    foreach(){
        divr[]=0.0;
        foreach_dimension()
           divr[] += (g.x[1]*rhop.x[1]*beta.x[1]-g.x[]*rhop.x[]*beta.x[]) ;
        divr[]/=(Delta);
        } 
  boundary  ({divr});            
  source =  divr;
/**
solve Poisson equation 
$$\nabla \cdot (\beta \nabla p  )= s$$ 
with [http://basilisk.fr/src/poisson.h](http://basilisk.fr/src/poisson.h)
*/

  mgp = poisson (p, source, beta);

/** the velocity is then computed from the gradient (note the `face_gradient_x p=  ((p[i] - p[i-1])/Delta)`) see 
[http://basilisk.fr/src/poisson.h](http://basilisk.fr/src/poisson.h)

*/
   foreach_face()
    u.x[] = - beta.x[]*face_gradient_x (p, 0) + beta.x[]*rhop.x[]*g.x[];
  boundary ((scalar *){u});
}  
/**
error
*/
event logfile (i++)
{
    stats s = statsf (p);
    fprintf (stderr, "%d %g %d %g %g %g\n",
             i, t, mgp.i, s.sum, s.min, s.max);
     fprintf (stderr,"%g \n",dt);
}
/**
Rate of flowing materials across the hole (loss of volume $-dV/dt$)
*/

event debit (t += 0.5 ) {
  static double Vold,V=1;
  Vold=V; 
  V=0;
  foreach()
    V = V + f[]* Delta * Delta;
  Q = -(V-Vold)/.5;   
  if(t>=.1) fprintf (stdout,"%lf %lf %lf \n",t,V/L0/H0,Q); 
  fflush (stdout);
}

/**
 if no `#include "grid/multigrid.h"` then adapt the mesh using the diffusive trick 
*/
 
//#if QUADTREE
event adapt(i+=1){

   refine ( ( f[] > .00051) && level < MAXLEVEL);
   

}
//#endif

/**
Save in a files...
*/
#if 0
event gfsview (i += 100) {
  static FILE * fp = popen ("gfsview2D gfsview.gfv", "w");
  output_gfs (fp, t = t);
}
#endif

event sauve (i+=200)
{
FILE *  fpc = fopen("pressure.txt", "w");
output_field ({p,u.x,u.y,f}, fpc, linear = true);
fclose(fpc);
    
FILE *  fpq = fopen("qh.txt", "w");
    double dx=.1;
    for (double x = 0 ; x < L0; x += dx){
        for (double y = 0 ; y < L0; y += dx){
            fprintf (fpq, "%g %g %g \n",
                     x,    interpolate (p, x, 0), interpolate (u.x, x, 0));}
        fprintf (fpq,"\n");}
    fclose(fpq);
}

event movie (t+=.02;t<tmax)
{
  //static FILE * fp2 = popen ("ppm2mpeg > lf.mpg", "w");
  scalar l[];
  foreach()
    l[] =   f[]*p[]  ;
   output_ppm (l,file="lf.mp4",min = 0.,max = 2,n = 1024,linear = true,box = {{0,-1},{L0,5}});
}

/**
# Results
## Run
To compile and run:

~~~bash
 qcc -O2 -Wall -o dupuit2D dupuit2D.c -lm
./dupuit2D
 
 
make dupuit2D.c.html
cp dupuit2D.c.html to
 sed -i -e 's/\\)//g' to;sed -i -e 's/\\(//g' to;sed -i -e 's/\\\[//g' to;sed -i -e 's/\\\]//g' to;
 mv to dupuit2D.c.html
 open dupuit2D.c.html;

~~~
 
or more clean

~~~bash
 make dupuit2D.tst; make dupuit2D/plots ; make dupuit2D.c.html;
~~~


## Plots

[(A film here)](dupuit2D/lf.mp4)  of pressure  (click for animation)
showing the discharge and the stabilsation on the Dupuit solution.

 
some iso p

~~~gnuplot
set key bottom
set xlabel "x"
set ylabel "p(x,0)" 
plot [:][:] 'pressure.txt' u 1:(abs($2)<.01?($3):NaN) t'num.' 
~~~ 

  
bla bla bla

 
 
Flow rate ($Q=-dV/dt$ as function of time) out of the aquifere.

~~~gnuplot
reset
set xlabel "t"
set ylabel "Q(t)" 
plot [:][0:1] 'out' u 1:3  not
~~~ 

### Dupuit solution :
The flow is in thin layer we have almost an hydrostatic pressure:
 $$p = \rho g h $$
 so that the longitudinal velocity is proportional to the slope (Darcy), and the flux is:
 $$ q=- K \rho g h \frac{\partial h}{\partial x}$$
 and by mass conservation :
 $$\frac{\partial h}{\partial t} =
  K \rho g
 \frac{\partial }{\partial x}( h \frac{\partial h}{\partial x})$$
 case of the Dupuit steady solution, $q$ is constant so that
 $$q = \frac{K  ( h_1^2 -  h_2^2)}{2 L}$$
 The final surface is called Dupuit's parabola:
 $$h=\sqrt{h_1^2- \frac{2 qx}K}.$$
 
~~~gnuplot
 reset
 set xlabel "x"
 h1=2
 h2=1
 qt=(h1*h1-h2*h2)/2./20.
 h(x)=sqrt(h1*h1-2*qt*x)
 p'qh.txt' t'p(x,0)' w p,''u 1:($3*$2) t 'q',qt t 'q theo',h1 not,h2 not ,h(x) t'h theo'
~~~
 

### Early times self similar solution
The Dupuit Boussinesq equation
$$\frac{\partial h}{\partial t} =
 K \rho g
\frac{\partial }{\partial x}( h \frac{\partial h}{\partial x})$$
has a  self similar solution of variable
$$\eta = \frac{ x}{\sqrt{t}}$$
and the mass flow rate is proportional to $t^{-1/2}$.
That   is visible at the begining, after a while, the flow is no more self similar as it goes to the Dupuit parabolic solution.

~~~gnuplot
reset
set logscale x
set logscale y
set xlabel "t"
set ylabel "Q(t)"
set size 1,.5
 plot [:][0.001:] 'out' u 1:3  not , 1.2/sqrt(x) t'self similar'
~~~ 

# Links
 
 * see [http://basilisk.fr/src/hele-shaw.h](http://basilisk.fr/src/hele-shaw.h)
 
 * see [http://basilisk.fr/sandbox/M1EMN/Exemples/darcyLambSneddon.c]()
 
 * see [http://basilisk.fr/sandbox/M1EMN/HYDGEO/dupuit2D.c]()
 
 * see [http://basilisk.fr/sandbox/M1EMN/HYDGEO/toddbear59.c]()
 
# Bibliography

* Adrien Guerin [phd](http://adrienguerin.fr/papers/manuscrit.pdf)
Dynamique de l’écoulement dans un aquifère non confiné, 2015

* Joseph Bear "Dynamic of fluids in porous media"

* [Lagrée P-Y](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/MFEnv.pdf)
 "Equations de Saint Venant et application, Ecoulements en milieux naturels,
écoulements en milieux souterrains" Cours MSF12, M1 UPMC

*/
