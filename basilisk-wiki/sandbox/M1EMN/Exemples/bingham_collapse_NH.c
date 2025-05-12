/**
 
# collapse of a rectangular Bingham column on a slope (Eugène and Eugène)
 
## Problem of Liu and Mei/ Balmforth: collapse of Bingham fluid
 
 A heap of fluid following Bingham rheology is released along a constant slope.
 We see the front moving to the right, and the left front going slowly up hill to the left. A real Bingham flow should stop. See Balmforth 1D  related examples  with [only mass equation and lubrication](http://basilisk.fr/sandbox/M1EMN/Exemples/viscous_collapse_noSV.c) 
 
![animation of the collapse](bingham_collapse_NH/a2.gif)
 
 We solve here the simplified system of RNSP equations with the "vertically-Lagrangian multilayer solver for free-surface flows" (not Navier Stokes and not Multilayer).

 
 
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
 with Bingham rheology $\tau_{xy}= \tau_y + \mu \dfrac{\partial u}{\partial z}$. 
 We define an equivalent viscosity
 $$\tau_{xy}=  \rho \nu_{eq} \dfrac{\partial u}{\partial z}
 \text{ with } \nu_{eq}=\rho \nu (1 + \dfrac{\tau_y/{\rho \nu}}{|\dfrac{\partial u}{\partial z}|})$$
 the `layered/hydro.h` is changed in [http://basilisk.fr/sandbox/M1EMN/Exemples/hydroNN.h]()

## Link with 1D model
 
 If we integrate over the depth of flow $h$ the incompressiblility equation, we obtain
 $\frac{\partial h}{\partial t} + \frac{\partial Q}{\partial x}=0$,
 where $Q=\int_0^h udy$. If we neglect inertia in the momentum equation, we solve for $u$ in
 $$0 =  -   \rho g  Z'_b -
 \dfrac{\partial   p}{\partial   x} +
 \dfrac{\partial   \tau_{xy}}{\partial y} $$
 were the pressure is hydrostatic $p=\rho g (h-y)$, and with the Bingham rheology. This gives then $Q=\int_0^h udy$ that we put in the mass conservation, and hence we obtain the [1D kinetic wave](http://basilisk.fr/sandbox/M1EMN/Exemples/bingham_collapse_noSV.c) from the paper of Balmforth.
 Here we use the Hydro Solver to see the influence of the inertia neglected in the 1D model.


## Link with viscous collapse
 
 The example is linked to the  collapse of a viscous fluid (Huppert 82 “The propagation of two-dimensional and axisymmetric viscous gravity currents over a rigid horizontal surface”), here we have a slope.
 It is done with only mass equation and lubrication [here](http://basilisk.fr/sandbox/M1EMN/Exemples/viscous_collapse_noSV.c) and
 with shallow water [there](http://basilisk.fr/sandbox/M1EMN/Exemples/viscous_collapse.c).
 The same RNSP equations are solved 
 [with Multilayer shallow water](http://basilisk.fr/sandbox/M1EMN/Exemples/viscous_collapse_ML.c)
 (using the Multilayer Shallow Water (Saint Venant Multi Couches) strategy of Audusse Sainte-Marie  et al 2011. See De Vita 2020 for details). 
 We compare here  with
 [Popinet (2019)](/Bibliography#popinet2019)
  both resolutions of RNSP system.
  The  hydrostatic
 [http://basilisk.fr/src/layered/hydro.h]()
 "hydro.h" has been changed in "hydroNN.h" to include variable viscosity.
 "difusion.h" has been chenged in it and is now "difusionNN.h" (Non-Newtonian).
 
 it is possible to test the  non hydrostatic
 [http://basilisk.fr/src/layered/nh.h](), it should be done next.

 
 
 
 ![Snapshot of the heap.](bingham_collapse_NH/snapshot.png)
 
## Code

### Includes
 */


#include "grid/multigrid1D.h"
double BBingham;
#include "hydroNN.h"
//#include "layered/nh.h"
/**

### Non Newtonian viscosity
 The definition of viscosity as a function of shear (computed in  "difusion.h"):
*/

double nu_eq(double shear,double pipi){
    double nu_eq;
    nu_eq = nu*(1 + BBingham/sqrt(sq(shear) + sq(.1e-8)));
    return nu_eq;}
/**
 definition of viscosity (it is a function of shear: Bingham flow)
 $$\tau = \nu (\frac{\partial u}{\partial z} + B)$$
 where $\nu$ is the viscosity and $\nu B$ is the yeld stress, we define an equivalent viscosity:
 $$\nu_{eq}=\nu (1 + \frac{B}{|\frac{\partial u}{\partial z}|})$$
 Note the regularization introduced to avoid division by zero.
 The $\nu$ and the $B$ have no dimension....



### Main and BC
 */
double tmax,alpha;


int main() {
    X0 = -0.5;
    L0 = 1.5;
    G  = 1.;
    N  = 2048/2;
    nl = 128*2;
    nu = 1.;
    alpha=1;
    tmax = 1;
    BBingham = 1.25 ;
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
 
 ![animation of the collapse](bingham_collapse_NH/a2.gif)

*/
#if 0
void plot_profile (double t, FILE * fp)
{
    fprintf (fp,
             "set title 't = %.2f'\n"
             " h(x) = (0.9*(1.28338-x*x))**(1/3.) \n"
             "t = %.2f'\n"
             "p [0:2][0:1.5]'-' u ($1/(t**.2)):($2*(t**.2)) w lp lc 3 t 'num' , h(x) t'huppert' \n", t,t);
    foreach()
    fprintf (fp, "%g %g \n", x,eta[] );
    fprintf (fp, "e\n\n");
    fflush (fp);
}
#else
void plot_profile (double t, FILE * fp)
{
    fprintf (fp,
             "set title 't = %.2f'\n"
             "p [-.5:1. ][:1.25]'-' u 1:3:2 w filledcu lc 3 t 'surface', '../bingham_collapse_noSV/shape-1.25.txt' t 'B=1.25 h' w l\n", t);
    foreach()
    fprintf (fp, "%g %g %g\n", x,eta[]-zb[],0.00);
    fprintf (fp, "e\n\n");
    fflush (fp);
}
#endif

event profiles (t += .01) {
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
    //foreach()
    // fprintf (stdout, "%g %g %g \n", x,eta[],t);
    //getchar();
}


/**
 
 We print the elevation $h=\eta -z_b$ at last time. */

event output (t = tmax) {
    foreach()
    fprintf (stdout, "%g %g %g \n", x, eta[]-zb[],t );
    fprintf (stdout, "\n");
}



/**

## Compilation

Usual make file

~~~bash
 make bingham_collapse_NH.tst
 make bingham_collapse_NH/plots;make bingham_collapse_NH.c.html
 source ../c2html.sh  bingham_collapse_NH
~~~
 
 

 
 
## Results

Comparison with 1D Balmforth  with $B=1.25$ at time $t=1$

~~~gnuplot Comparison
 set xlabel "x"
 set ylabel "h(x,t)"
 p [-.4:.5][0:1.2]'out' w lp t'NH',\
'../bingham_collapse_noSV/shape-1.25.txt' t 'B=1.25 h' w l
~~~

Comparison with 1D Balmforth  and 
 with 2D RNSP Multilayer `multilayer.h` with $B=1.25$  at time $t=1$


~~~gnuplot B=1.25
set xlabel "x"
set ylabel "h(x,t)"
p[-0.4:.5][:1.2]'out' w lp t'2D NH',\
'../bingham_collapse_noSV/shape-1.25.txt' t '1D   ' w l,\
'../bingham_collapse_ML/shapeML-1.25.txt' t'2D ML' w l
~~~
 


confine 03/20 (instead of using alcoholic gel, use hydro include file!)
 
## Links

* see the related example  bingham collapse in 1D    [1D kinetic wave](http://basilisk.fr/sandbox/M1EMN/Exemples/bingham_collapse_noSV.c)

* see the related example in 2D Multilayer for comparison
   [with `multilayer.h`](http://basilisk.fr/sandbox/M1EMN/Exemples/bingham_collapse_ML.c)
   
* 
 related examples [only mass equation and lubrication newtonian ](http://basilisk.fr/sandbox/M1EMN/Exemples/viscous_collapse_noSV.c) and
 [with shallow water, newtonian](http://basilisk.fr/sandbox/M1EMN/Exemples/viscous_collapse.c)
 see as well Navier Stokes solution.
 
 * [http://basilisk.fr/sandbox/M1EMN/Exemples/viscous_collapse_ML.c]() viscous `multilayer.h`

 * [http://basilisk.fr/sandbox/M1EMN/Exemples/viscous_collapse_NH.c]() viscous with `hydro.h`
 
## Bibliography
 
 * Lagrée  [M2EMN
 Master 2 Ecoulements en Milieu Naturel](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/mainM2EMN.pdf)
 see Liu & Mei preoblem and Balfmorth problem
 
 *  [Popinet (2019)](/Bibliography#popinet2019)
 A vertically-Lagrangian, non-hydrostatic, multilayer model
 for multiscale free-surface flows, Journal of Computational Physics

 * Francesco De Vita, Pierre-Yves Lagrée, Sergio Chibbaro, Stéphane Popinet
 [Beyond Shallow Water: appraisal of a numerical approach to hydraulic jumps based upon the Boundary Layer Theory.](http://www.lmm.jussieu.fr/~lagree/TEXTES/PDF/devita20.pdf)
 Volume 79, January–February 2020, Pages 233-246
 European Journal of Mechanics - B/Fluids
 https://doi.org/10.1016/j.euromechflu.2019.09.010
 
 
 */
