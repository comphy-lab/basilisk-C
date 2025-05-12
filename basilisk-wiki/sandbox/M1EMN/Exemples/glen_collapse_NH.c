/**
 
# collapse of a rectangular ice column
 
## Problem of Halfar: collapse of power law fluid (Ostwald de Waele fluid)
 
 A heap of fluid (following power law  rheology) is released along a flat horizontal  plate.
 It spread due to gravity, and it is slowed down due to no Newtonian viscosity.
 This is for example the case of ice on Antarctic, where for ice rheology the Glen's law is applied $\dot\gamma = \tau^3$
 
  See  related Newtonian examples  with [only mass equation and lubrication](http://basilisk.fr/sandbox/M1EMN/Exemples/viscous_collapse_noSV.c) which is of course the sploutch problem.
 
![animation of the collapse](glen_collapse_NH/animate.gif)
 
 
 The 1D kinetic wave case with realistic values for ice is [here](http://basilisk.fr/sandbox/nlemoine/halfar1D.c) and is refered as Halfar's solution.
 
 
 
 
 We do not solve Navier Stokes equations, we solve here the simplified system of RNSP equations with the "vertically-Lagrangian multilayer solver for free-surface flows" (not Multilayer).

 
 
## RNSP Model equations
 We solve the boundary layer Prandtl (RNSP) equations corresponding to a thin layer approximation of the Navier Stokes equations (hydrostatic):
 $$\left\{
 \begin{aligned}
 \dfrac{\partial  u}{\partial   x}  + \dfrac{\partial  v}{\partial   y} &= 0  \\
 \rho (\dfrac{\partial   u}{\partial   t} + \dfrac{\partial   u^2}{\partial   x}  + \dfrac{\partial  u  v}{\partial   y} )&=  -   \rho g  Z'_b -
 \dfrac{\partial   p}{\partial   x} +
 \dfrac{\partial   \tau_{xy}}{\partial   y} \\
 0 &=   -\rho g -\dfrac{\partial  p}{\partial  y}.
 \end{aligned}
 \right.$$
 with Ostwald de Waele  rheology $\tau_{xy}= \rho \nu \left(\dfrac{\partial u}{\partial z}\right)^n$.
 We define an equivalent viscosity
 $$\tau_{xy}=  \rho \nu_{eq} \dfrac{\partial u}{\partial z}
 \text{ with } \nu_{eq}=  \nu \left( \dfrac{\partial u}{\partial z}\right)^{n-1})$$
 the `layered/hydro.h` is changed in [http://basilisk.fr/sandbox/M1EMN/Exemples/hydroNN.h]()

## Link with 1D model
 
 If we integrate over the depth of flow $h$ the incompressiblility equation, we obtain
 $\frac{\partial h}{\partial t} + \frac{\partial Q}{\partial x}=0$,
 where $Q=\int_0^h udy$. If we neglect inertia in the momentum equation, we solve for $u$ in
 $$0 =  -   \rho g  Z'_b -
 \dfrac{\partial   p}{\partial   x} +
 \dfrac{\partial   \tau_{xy}}{\partial y} $$
 were the pressure is hydrostatic $p=\rho g (h-y)$, and with the power law  rheology. This gives then $Q=\int_0^h udy$
 $$Q=  \frac{n\;h^{\frac{n+1}{n}} }{(n+1)(2n+1)}\Bigg(\frac{\rho g}{\mu}(-\frac{\partial h}{\partial x})\Bigg)^{\frac{1}{n}}\; (n + 1 )h
$$
 that we put in the mass conservation. In the case of Herschel Bulkley fluids, we have a threshold stress and hence we obtain the [1D kinetic wave](http://basilisk.fr/sandbox/M1EMN/Exemples/herschel-column-noSV.c) from the paper of Balmforth (in case of Herschel Bulkley fluids).
 
 
 
 
 
 

## Link with the Sploutch
 
 The example is linked to the Sploutch :  collapse of a viscous fluid (Huppert 82 “The propagation of two-dimensional and axisymmetric viscous gravity currents over a rigid horizontal surface”), here we have a slope.
 It is done with only mass equation and lubrication [here](http://basilisk.fr/sandbox/M1EMN/Exemples/viscous_collapse_noSV.c) and
 with shallow water [there](http://basilisk.fr/sandbox/M1EMN/Exemples/viscous_collapse.c).
 The same RNSP equations are solved 
 [with Multilayer shallow water](http://basilisk.fr/sandbox/M1EMN/Exemples/viscous_collapse_ML.c)
 (using the Multilayer Shallow Water (Saint Venant Multi Couches) strategy of Audusse Sainte-Marie  et al 2011. See De Vita 2020 for details). 
 
 
 
 Here we use the Hydro Solver
 [Popinet (2019)](/Bibliography#popinet2019) to see the influence of the inertia neglected in the 1D model.

 
  The  hydrostatic
 [http://basilisk.fr/src/layered/hydro.h]()
 "hydro.h" has been changed in "hydroNN.h" to include variable viscosity.
 "difusion.h" has been chenged in it and is now "difusionNN.h" (Non-Newtonian).
 
 it is possible to test the  non hydrostatic
 [http://basilisk.fr/src/layered/nh.h](), it should be done next.

 
 
 
## Code

### Includes
 */


#include "grid/multigrid1D.h"
double BBingham;
double nGlen;
#include "hydroNN.h"
//#include "layered/nh.h"
/**

### Non Newtonian viscosity
 The definition of viscosity as a function of shear (computed in  "difusionNN.h"):
*/

double nu_eq(double shear,double pipi){
    double nu_eq;
    nu_eq = nu*pow((sqrt(sq(shear) + sq(.1e-8))),1./nGlen-1);
    return nu_eq;}
/**
 definition of viscosity (it is a function of shear: power law flow)
 $$\tau = \nu (\frac{\partial u}{\partial z})^n $$
 where $\nu$ is the consitency and $n$ is the index, we define an equivalent viscosity:
 $$\nu_{eq}=\nu  |\frac{\partial u}{\partial z}|^{n-1}$$
 Note the regularization introduced to avoid division by zero.
 The $\nu$ and the $n$ have no dimension here....

 Of course, as Glen's law is defined by $\dot\gamma=\tau^3$ instead of
 $\tau = (\dot\gamma)^{1/3}$
 the index is $n=1/3$


### Main and BC
 */
double tmax,alpha;


int main() {
    X0 = -1.;
    L0 = 2;
    G  = 1.;
    N  = 1024;
    nl = 256;
    nu = 1.;
    alpha=0;
    tmax = 4;
    nGlen = 3;
    BBingham = 0 ;
    run();
}

/**
 We impose boundary condition for the position of the surface $\eta$ (as  $h$ not defined in 'hydro.h', but by construction
 $\eta= z_b+h$). */
eta[left] = neumann (0);
eta[right] = neumann (0);

/**
### Initialization  */

event init (i = 0) {
    /**
     We set a zero velocity at the inlet and a free outlet. */
        u.n[left] = dirichlet(0);
        u.n[right] = neumann(0.);
    /**
     We initialize *h*. */
    foreach() {
        zb[] = -(x)*alpha;
        double h0 = (fabs(x)<.25) + dry;
        foreach_layer()
        h[] = (h0 )/nl;}
}
/**
### Output
 
 We use gnuplot to visualise the  profile during time in live with `X11`
 or we generate an animatio`in `gif`
 
 ![animation of the collapse](glen_collapse_NH/animate.gif)

*/

void plot_profile (double t, FILE * fp)
{
    fprintf (fp,
             "set title 't = %.2f'\n"
             "p [-1:1. ][:1.25]'-' u 1:3:2 w filledcu lc 3 t 'surface' \n", t);
    foreach()
    fprintf (fp, "%g %g %g\n", x,eta[]-zb[],0.00);
    fprintf (fp, "e\n\n");
    fflush (fp);
}


event profiles (t += .02) {
    static FILE * fp = popen ("gnuplot --persist 2> /dev/null", "w");
     if(t==0) fprintf (fp,"set term gif animate;set output 'animate.gif';set size ratio .333333\n");
    // comment to have direct X11 animation, uncomment to generate the gif
    plot_profile (t, fp);
    
    if(t==tmax){
        fprintf (fp,"! cp animate.gif ../animate.gif\n");
        fprintf (fp,"! cp animate.gif a2.gif \n");
    }
}

/**
 
generate a snapshot at $t=t_{max}$.
 
*/

event gnuplot (t = end) {
    FILE * fp = popen ("gnuplot", "w");
    fprintf (fp,
             "set term png enhanced size 640,200 font \",8\"\n"
             "set output 'snapshot.png'\n");
    plot_profile (t, fp);
}


/**
 
 We print the elevation $h=\eta -z_b$ as function of time. */

event output (t += .05 ; t <= tmax) {
    foreach()
    fprintf (stdout, "%g %g %g \n", x, eta[]-zb[],t );
    fprintf (stdout, "\n");
}



/**

## Compilation

Usual make file

~~~bash
 make glen_collapse_NH.tst
 make glen_collapse_NH/plots;make glen_collapse_NH.c.html
 source ../c2html.sh  glen_collapse_NH
~~~
 
 

 
 
## Results

 
 
![Snapshot of the heap at finale time.](glen_collapse_NH/snapshot.png)
 
 
Plot of the solution at different times
 
~~~gnuplot
 set xlabel "x"
 set ylabel "h(x,t)"
 nG=3
 p [:][0:1.2]'out' u 1:2  w lp t'calc.'
~~~

 The self similar solution of Halfar is (with $n=3$)
 $$h = t^{-1/(3 n + 2)}F(\frac{x}{t^{1/(3 n + 2)}})$$
 so we rescale accordingly and plot for time enough large, say larger than 2:
 
 
 
~~~gnuplot Comparison
 set xlabel "x/t**(1/(3 n+2) )"
 set ylabel "t**(-1/(3 n+2)) h "
 nG=3
 p [-.6:.6][0:1.]'out' u ($1/($3**(1./(nG*3+2)))):($3>2?($2*($3**(1./(3*nG+2)))):NaN) w lp t'calcul'
~~~




 
 
## Links
 
* [http://basilisk.fr/sandbox/nlemoine/halfar1D.c]() Nicolas's Sandbox
 
* [http://basilisk.fr/sandbox/M1EMN/Exemples/herschel-column-noSV.c]()

 
 
## Bibliography
 
 * On the Dynamics of the Ice Sheets P. HALFAR,
JGR VOL.86 ,1981
 
 * On the Dynamics of the Ice Sheets 2 P. HALFAR, JGR VOL. 88,   1983
 
 * Lagrée  [M2EMN
 Master 2 Ecoulements en Milieu Naturel](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/mainM2EMN.pdf)
 see  power law law and Ostwald de Waele fluids problems
 
 *  [Popinet (2019)](/Bibliography#popinet2019)
 A vertically-Lagrangian, non-hydrostatic, multilayer model
 for multiscale free-surface flows, Journal of Computational Physics

 * Francesco De Vita, Pierre-Yves Lagrée, Sergio Chibbaro, Stéphane Popinet
 [Beyond Shallow Water: appraisal of a numerical approach to hydraulic jumps based upon the Boundary Layer Theory.](http://www.lmm.jussieu.fr/~lagree/TEXTES/PDF/devita20.pdf)
 Volume 79, January–February 2020, Pages 233-246
 European Journal of Mechanics - B/Fluids
 https://doi.org/10.1016/j.euromechflu.2019.09.010
 
 
 */
