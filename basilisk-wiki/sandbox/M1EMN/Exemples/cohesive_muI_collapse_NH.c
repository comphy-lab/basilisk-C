/**
 
# collapse of a rectangular cohesive granular column on a slope (Eugène and Eugène)
 
 
 We test here the hydro solver with a non newtonian viscosity,
 it is  more or less working, with many points and with  small oscillations.
 
## Problem
 
 A heap of fluid following cohesive granular  (with cohesion)  rheology is released along a constant slope.
 
 
![animation of the collapse](cohesive_muI_collapse_NH/a2.gif)
 
 
 
 
## RNSP model equations
 We solve the boundary layer (RNSP) equations corresponding to a thin layer approximation of the Navier Stokes equations:
 $$\left\{
 \begin{aligned}
 \dfrac{\partial  u}{\partial   x}  + \dfrac{\partial  v}{\partial   y} &= 0  \\
 \rho (\dfrac{\partial   u}{\partial   t} + \dfrac{\partial   u^2}{\partial   x}  + \dfrac{\partial  u  v}{\partial   y} )&=  -   \rho g  Z'_b -
 \dfrac{\partial   p}{\partial   x} +
 \dfrac{\partial   \tau_{xy}}{\partial   y} \\
 0 &=   -\rho g -\dfrac{\partial  p}{\partial  y}.
 \end{aligned}
 \right.$$
 following $\mu(I)$ rheology, with $I=\dfrac{d}{\sqrt{p/\rho}}\dfrac{\partial u}{\partial y}$ this gives
 $\tau_{xy}= \tau_Y + (\mu( \dfrac{d}{\sqrt{p/\rho}}\dfrac{\partial u}{\partial y})p),$
 there is cohesion  with additional $\tau_Y$ (Yield) to the classical $(\mu(I)p)$.
 
 
 
 We solve this with the hydro solver   from Popinet.
 We tune the solver  to consider a non newtonian viscosity function of the shear
 $$\nu_{eq}=  \frac{(\tau_Y+ \mu(I)p  )/\rho}{|\frac{\partial u}{\partial y}|}$$
 
 
 the `layered/hydro.h` is changed in [http://basilisk.fr/sandbox/M1EMN/Exemples/hydroNN.h]()
 
 
## Link with 1D model
 
  See discussion in [http://basilisk.fr/sandbox/M1EMN/Exemples/cohesive_muI_collapseML.c]()
 
 


## Link with viscous collapse
 
 The example is linked to the  collapse of a viscous fluid (Huppert 82 “The propagation of two-dimensional and axisymmetric viscous gravity currents over a rigid horizontal surface”), here we have a slope.
 The same RNSP equations are solved 
 [with Multilayer shallow water](http://basilisk.fr/sandbox/M1EMN/Exemples/viscous_collapse_ML.c)
 (using the Multilayer Shallow Water (Saint Venant Multi Couches) strategy of Audusse Sainte-Marie  et al 2011. See De Vita 2020 for details).
 
 See discussion in [http://basilisk.fr/sandbox/M1EMN/Exemples/cohesive_muI_collapseML.c]()
 
 
 
 ![Snapshot of the heap.](cohesive_muI_collapse_NH/snapshot.png)
 
 
 
 
 We compare here  with
 [Popinet (2019)](/Bibliography#popinet2019)
  both resolutions of RNSP system.
  The  hydrostatic
 [http://basilisk.fr/src/layered/hydro.h]()
 it is possible to test the  non hydrostatic
 [http://basilisk.fr/src/layered/nh.h]().
  "hydro.h" has been changed in "hydroNN.h" to include variable Non Newtonian viscosity.
  "difusion.h" has been changed in it and is now  "difusionNN.h".
 
 

 ![animation](cohesive_muI_collapse_NH/animate.gif)
 
 
## Code

 */

#include "grid/multigrid1D.h"
//#include "layered/hydro.h"
#include "hydroNN.h"
//#include "layered/nh.h"


/**
definition of viscosity as a function of shear:
*/
double BBingham;

double nu_eq(double shear,double press){
    double nueq,In,muP;
    In = 1./32*sqrt(sq(shear))/sqrt(press);
    muP = (.4 +.3*In/(In+.4)) *  press;
    nueq = ( BBingham +  muP )/sqrt(sq(shear) + sq(.1e-10));
    return nueq;}
/**
 definition of viscosity from the function
 $\tau = \tau_Y+\mu(I)p$
  we define an equivalent viscosity:
 $$\nu_{eq}= \frac{\tau}{|\frac{\partial u}{\partial z}|}$$
 Note the regularization introduced to avoid division by zero.
 
 `nu=1`is here a flag for diffusion.
 */



double tmax,alpha;

int main() {
    X0 = 0;
    L0 = 8;
    G  = 1.;
    N  = 512;
    nl = 512;
    nu = 1.;
    alpha=0.0;
    tmax = 5.;
    DT = .001;
    BBingham = 0.0 ;
    run();
}

/**
 We impose boundary condition for the positio,n of the surface $\eta$ (as  $h$ not defined in 'hydro.h', but by construction
 $\eta= z_b+h$). */
eta[left] = neumann (0);
eta[right] = neumann (0);

/**
## Initialization  */

event init (i = 0) {
    /**
     We set a zero velocity at the inlet and a free outlet. */
    
        u.n[left] = neumann(0.);
        u.n[right] = neumann(0.);

    /**
     We initialize *h*. */
    foreach() {
        zb[] = -(x)*alpha;
        double h0 = (fabs(x)<2) + 0.0000000001 +  dry;
       foreach_layer()
           h[] = (h0 )/nl;

}
}
/**
## Output
 
 We use gnuplot to visualise the  profile during time in live with `X11`
 or we generate an animation in `gif`
 
![animation of the collapse](cohesive_muI_collapse_NH/a2.gif)
 

*/

void plot_profile (double t, FILE * fp)
{
    fprintf (fp,
             "set title 't = %.2f'\n"
             "p [0:][:1.25]'-' u 1:3:2 w filledcu lc 3 t 'surface','../savagestaron/log' t' Savage Hutter 1 D' w l \n", t);
    foreach()
    fprintf (fp, "%g %g %g\n", x,eta[]-zb[],0.00);
    fprintf (fp, "e\n\n");
    fflush (fp);
}


event profiles (t += .01,t<=tmax) {
    static FILE * fp = popen ("gnuplot --persist 2> /dev/null", "w");
    if(t==0) fprintf (fp,"set term gif animate;set output 'animate.gif';set size ratio .333333\n");
   //  if(t==0) fprintf (fp,"set term x11;set size ratio .333333\n");
    // comment to have direct X11 animation, uncomment to generate the gif
    plot_profile (t, fp);
    
    if(t==tmax){
        fprintf (fp,"! cp animate.gif ../animate.gif\n");
        fprintf (fp,"! cp animate.gif a2.gif \n");
    }
}

/**
 
generate a snapshot
 
*/

event gnuplot (t = end) {
    FILE * fp = popen ("gnuplot", "w");
    fprintf (fp,"! cp animate.gif ../animate.gif\n");
    fprintf (fp,"! cp animate.gif a2.gif \n");
    fprintf (fp,
             "set term png enhanced size 640,200 font \",8\"\n"
             "set output 'snapshot.png'\n");
    plot_profile (t, fp);
}


/**
 
 We print the elevation $h=\eta -z_b$ at last time at $t=t_{max}$.. */

event output (t = {1,2,3,4}) {
    foreach()
    fprintf (stdout, "%g %g %g \n", x, eta[]-zb[],t );
    fprintf (stdout, "\n");
}



/**

## Compilation

Usual make file

~~~bash
 make cohesive_muI_collapse_NH.tst
 make cohesive_muI_collapse_NH/plots;make cohesive_muI_collapse_NH.c.html
 
 
 source ../c2html.sh  cohesive_muI_collapse_NH
~~~
 
 

 
 
# Results

Comparison with 1D Savage-Hutter  (coded in [http://basilisk.fr/sandbox/M1EMN/Exemples/savagestaron.c]()) at ime $t=1$

~~~gnuplot Comparison
 set xlabel "x"
 set ylabel "h(x,t)"
 p [:5][0:1.2]'out' w lp t'NH',\
'../savagestaron/log' w l
~~~

Comparison with 1D Savage-Hutter and
with 2D RNSP Multilayer `multilayer.h`

~~~gnuplot
set xlabel "x"
set ylabel "h(x,t)"
p[:5][:1.2]'out' w p t'2D NH',\
'../savagestaron/log'   t '1D   ' w l,\
'../cohesive_muI_collapse_ML/shapeML-0.00.txt' t'2D ML' w l
~~~
 


 
## Links
 
 * [http://basilisk.fr/sandbox/M1EMN/Exemples/cohesive_muI_collapse_ML.c]() 2D viscous `multilayer.h`
 * [http://basilisk.fr/sandbox/M1EMN/Exemples/cohesive_muI_collapse_NH.c]() 2D viscous with `hydro.h`
 * [http://basilisk.fr/sandbox/M1EMN/Exemples/savagestaron.c]() 1D with `saintvenant.h`
 
## Bibliography
 
 * Lagrée  [M2EMN
 Master 2 Ecoulements en Milieu Naturel](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/mainM2EMN.pdf)
 
 *  [Popinet (2019)](/Bibliography#popinet2019)
 A vertically-Lagrangian, non-hydrostatic, multilayer model
 for multiscale free-surface flows, Journal of Computational Physics

 * Francesco De Vita, Pierre-Yves Lagrée, Sergio Chibbaro, Stéphane Popinet
 [Beyond Shallow Water: appraisal of a numerical approach to hydraulic jumps based upon the Boundary Layer Theory.](http://www.lmm.jussieu.fr/~lagree/TEXTES/PDF/devita20.pdf)
 Volume 79, January–February 2020, Pages 233-246
 European Journal of Mechanics - B/Fluids
 https://doi.org/10.1016/j.euromechflu.2019.09.010
 
 
 */
