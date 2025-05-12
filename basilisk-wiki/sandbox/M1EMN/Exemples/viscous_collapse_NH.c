/**
 
# collapse of a rectangular viscous column (Eugène and Eugène)

 
## problem
 
What happens if a viscous fluid is suddenly released from its container?
 
This is the  collapse of a viscous fluid (double viscous dam break),
 from the paper: Huppert 82 “The propagation of two-dimensional and axisymmetric viscous gravity currents over a rigid horizontal surface”
 
## Equations
 it was  done with [only mass equation and lubrication ](http://basilisk.fr/sandbox/M1EMN/Exemples/viscous_collapse_noSV.c) and
 [with shallow water](http://basilisk.fr/sandbox/M1EMN/Exemples/viscous_collapse.c) and
 [with Multilayer shallow water](http://basilisk.fr/sandbox/M1EMN/Exemples/viscous_collapse_ML.c)
(using the Multilayer Shallow Water (Saint Venant Multi Couches See De Vita 2020 for details)
 and now compare with
 [Popinet (2019)](/Bibliography#popinet2019)
 hydrostatic  [http://basilisk.fr/src/layered/hydro.h]()
  and 
 the non hydrostatic
 [http://basilisk.fr/src/layered/layered/nh.h"]() models.

 
 
 
 
 

 ![Snapshot of the heap.](viscous_collapse_NH/snapshot.png)
 
## Code
  
 */
#define NH 0

#if NH
#include "grid/multigrid1D.h"
#include "layered/hydro.h"
//#include "layered/nh.h"
#else
#include "grid/multigrid1D.h"
//#include "grid/cartesian1D.h"
#include "saint-venant.h"
#endif

double tmax;
/**
 
*/
int main() {
    X0 = 0.;
    L0 = 5;
    G  = 1.;
    N  = 128;
    nl = 15;
    nu = 1.;
    tmax = 40;
    run();
}

/**
 We impose boundary condition for$\eta$ (as  $h$ not defined in 'hydro.h' (but by construction
 $\eta= z_b+h$). */
eta[left] = neumann (0);
eta[right] = neumann (0);

/**
## Initialization  */

event init (i = 0) {
    /**
     We set a zero velocity at the inlet and a free outlet. */
    
    for (vector u in ul) {
        u.n[left] = dirichlet(0);
        u.n[right] = neumann(0.);
    }
    /**
     We initialize *h*. */
    foreach() {
        zb[] = 0.;
        double h0 = (x<1) + dry ;
#if NH
        for (scalar h in hl)
            h[] =  ( h0 )/nl ;
#else
        h[] =  h0  ;
        u.x[] = 0;
#endif
}
}
/**
## Output
 
 We use gnuplot to visualise the  profile during time
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
             "p [ ][0:1.5]'-' u 1:3:2 w filledcu lc 3 t ''\n", t);
    foreach()
    fprintf (fp, "%g %g %g\n", x,eta[],zb[]);
    fprintf (fp, "e\n\n");
    fflush (fp);
}
#endif

event profiles (t += .1) {
    static FILE * fp = popen ("gnuplot 2> /dev/null", "w");
    plot_profile (t, fp);
}

/**
 
 runs and to generate a snapshot at $t=t_{max}$.
 
 */

event gnuplot (t = end) {
    FILE * fp = popen ("gnuplot", "w");
    fprintf (fp,
             "set term png enhanced size 640,200 font \",8\"\n"
             "set output 'snapshot.png'\n");
    plot_profile (t, fp);
    foreach()
    fprintf (stdout, "%g %g %g \n", x,eta[],t);
    // getchar();
}


/**
 
 We print the elevation  . */

event output (t = tmax) {
    foreach()
    fprintf (stdout, "%g %g %g \n", x, eta[],t );
    fprintf (stdout, "\n");
}



/**

## Compilation

~~~bash
 
 make viscous_collapse_NH.tst
 make viscous_collapse_NH/plots;make viscous_collapse_NH.c.html
 source ../Exemples/c2html.sh  viscous_collapse_NH
 
~~~
 

## Results

Comparison of theoretical slf similar solution.
The non hydrostatic code creates oscillation...
 
~~~gnuplot Comparison of theoretical and numerical timeseries
 p [0:1.5]'out'u ($1/($3**.2)):($2*($3**.2)) w lp t'NH',\
 '../../Exemples/viscous_collapse_ML/log' u ($1/($4**.2)):($4>80?$2*($4**.2):NaN) t'comp ML.' w l,\
   (9./10*(1.28338-x*x))**(1/3.) t'analytic'
~~~

 
 
 
confine 03/20
 
## Links
 related examples [only mass equation and lubrication ](http://basilisk.fr/sandbox/M1EMN/Exemples/viscous_collapse_noSV.c) and
 [with shallow water](http://basilisk.fr/sandbox/M1EMN/Exemples/viscous_collapse.c)
 see as well Navier Stokes solution.
 

 * [Popinet (2019)](/Bibliography#popinet2019)  from bar.c test case
 
## Bibliography
 
 * Huppert
 "The propagation of two-dimensional and axisymmetric viscous gravity currents over a rigid horizontal surface"
 [JFM 82](http://www.itg.cam.ac.uk/people/heh/Paper47.pdf)
 
 * Lagrée  [M2EMN
 Master 2 Ecoulements en Milieu Naturel](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/mainM2EMN.pdf)
 
 * Francesco De Vita, Pierre-Yves Lagrée, Sergio Chibbaro, Stéphane Popinet
 [Beyond Shallow Water: appraisal of a numerical approach to hydraulic jumps based upon the Boundary Layer Theory.](http://www.lmm.jussieu.fr/~lagree/TEXTES/PDF/devita20.pdf)
 Volume 79, January–February 2020, Pages 233-246
 European Journal of Mechanics - B/Fluids
 https://doi.org/10.1016/j.euromechflu.2019.09.010
 
 * [Popinet (2019)](/Bibliography#popinet2019)
 
 */
