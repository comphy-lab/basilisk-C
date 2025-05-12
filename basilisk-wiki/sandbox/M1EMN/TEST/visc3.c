/**
#    Collapse of a heap with   on inclined plane

 

 ![animation  ](visc3/animate.gif)


We solve with [Multilayer Euler Lagrange solver](hydro.h) in 1D,
the only 1 layers set :
$$\left\{\begin{array}{l}
         \partial_t h_0+\partial_x (hu)_0=0\\
         \partial_t (hu)_0+ \partial_x (hu^2)_0= - g h_0\partial_x \eta   - \partial_x (h q)_0 + \left[q \partial_x z \right]_0
   + \frac{- 3 ( h_0 u_0 )}{
   h_0^2  }\\
  \end{array}\right.$$
 




- where the friction is laminar only in the lower layer :

No upper layer   
 
 
 
- with at  t=0, a given  e heap $h_0(|x|<0.5,t=0)=.5$  so that $\int_{-\infty}^{\infty} h_0(x,0)dx =1$.

 

 We check  the one fluid case (second Huppert problem).
*/

#include "grid/multigrid1D.h"
#include "layered/hydro.h"
// declare parameters

double tmax;
/**
 position of domain `X0`, length `L0`, no dimension `G=1`
run with 2048*8 points (less is not enough)
*/

int main()
{
    X0 = -5.0;
    L0 = 32;    // 15
    N =  512*2*2*2; // 2048*8
    nl = 1;
    G = 1.0;
    DT = .001;
    tmax=5;   // 400
    run();
}
/** start by a rectangular column  of different densities
 */
event init (i = 0)
{
    foreach(){
        zb[] = - x*0;
      
        h[0,0,0]= ((fabs(x-1)<.5));

    }
}

/**
Remember the equilibrium is
$$(hu)_0=-\frac{h_0^2 (h_0 )}{3}(\frac{\partial p_0}{\partial x}) 
 $$

so we can write the split for friction
$$
  \partial_t (hu)_0+=
   - \frac{3 h_0 u_0  }{
   h_0^2 }   
$$
 
*/

 
#if 0
event viscous_term (i++,last)
{
  if (1 > 0.) {
    foreach() {
      foreach_layer()
    foreach_dimension()
      u.x[] += dt*(ha.x[] + ha.x[1])/(hf.x[] + hf.x[1] + dry);
        //   foreach_dimension()
        // vertical_diffusion (point, h, u.x, dt, nu, dut.x[], u_b.x[], lambda_b.x[]);
        foreach_layer()
      foreach_dimension()
        {     double ff = (h[]) < dry ? HUGE : (1. + 3*dt/sq((hf.x[] + hf.x[1])/2));
          u.x[] /= ff;}
 
      foreach_layer()
    foreach_dimension()
      u.x[] -= dt*(ha.x[] + ha.x[1])/(hf.x[] + hf.x[1] + dry);
    }
  }
}

#endif

#if 1
event friction (i++) {
  // Poiseuille friction (implicit scheme)
  foreach()
  {       double ff = h[] < dry ? HUGE : (1. + 3.*dt/h[]/h[]);
    u.x[] /= ff;}
}
#endif
 
 

 
event outputfile (t+=.2;t <= tmax) {
  foreach()
  fprintf (stderr,"%g %g %g %g %g \n", x, h[0,0,0], u.x[0,0,0], t ,zb[]);
//   fprintf (stderr,"%g %g %g %g %g \n", x, h[0,0,0], h[0,0,1], t ,zb[]);
  fprintf (stderr,"\n");
  //foreach()
  // fprintf (stdout,"%g %g %g %g %g \n", x, h[0,0,0], h[0,0,1], t ,zb[]);
  //fprintf (stdout,"\n");
}
/** animation
*/
event animatedplot (t+=.1) {
    static FILE * fp = popen ("gnuplot -persist 2> /dev/null", "w");
    if(t==0) fprintf (fp,"set term gif animate;set output 'animate.gif';set size ratio .333333\n");
    fprintf (fp,"\nset grid\n");
    fprintf (fp,"set title 'Landslide 1D --- t= %.2lf '\n"
      "t= %.2lf  ; "
      "p[%g:%g][-1:1]  '-' u 1:($2+$3) t'free surface' w l lt 3,"
      "'' u 1:4 t'topo/L0' w l lt -1\\\n",
           t,t,X0,X0+L0);
    foreach()
    fprintf (fp,"%g %g %g %g %g\n", x, h[0,0,0], u.x[0,0,0], zb[]/L0, t);
    fprintf (fp,"e\n\n");
    fflush (fp);
    if(t==tmax){
        fprintf (fp,"! cp animate.gif a2.gif \n");}
}
 

/**
end of subroutines

 
#  Run

To compile and run  using `make`
 
~~~bash
make visc2layerML.tst
make visc2layerML/plots
make visc2layerML.c.html
source ../c2html.sh visc2layerML
open isc2layerML.c.html

~~~
 
 
#   Results
 
 
 
 An example of collapse along a slope
 
~~~gnuplot collapse
set xlabel 'x '
set ylabel 'h(x,t)'
p [-1:3]'../visc3_SR/log' u 1:($4==.4?($2) :NaN) t'SR' w l,\
'../visc3/log' u 1:($4==.4?($2) :NaN) t'NH' w l

~~~

  
~~~gnuplot collapse
set xlabel 'x '
set ylabel 'h(x,t)'
p '../visc3_SR/log' u 1:($4==1?($2) :NaN) t'SR  1' w l,\
'log' u 1:($4==1?($2) :NaN) t'NH  1' w l,\
'../visc3_SR/log' u 1:($4==2?($2) :NaN) t'SR 2' w lp,\
'log' u 1:($4==2?($2) :NaN) t'NH 2' w l,\
'../visc3_SR/log' u 1:($4==4?($2) :NaN) t'SR 9' w lp,\
'log' u 1:($4==4?($2) :NaN) t'NH 9' w l
~~~
 

~~~gnuplot collapse
 set xlabel 'x '
 set ylabel 'h(x,t)'
 p '../visc3_SR/log' u 1:($4==1?($2*$3) :NaN) t'SR  1' w l,\
 'log' u 1:($4==1?($2*$3) :NaN) t'NH  1' w l,\
'../visc3_SR/log' u 1:($4==2?($2*$3) :NaN) t'SR 10' w l,\
 'log' u 1:($4==2?($2*$3) :NaN) t'NH 10' w l,\
 '../visc3_SR/log' u 1:($4==4?($2*$3) :NaN) t'SR 99' w l,\
 'log' u 1:($4==4?($2*$3) :NaN) t'NH 99' w l
~~~
 

 
 
 
 
# Links


 * [http://basilisk.fr/sandbox/M1EMN/Exemples/viscolsqrt.c]()
 * [http://basilisk.fr/sandbox/M1EMN/Exemples/viscous_collapsesqrt_ML.c]()

 with Bingham
 
 * [Bingham 1D collapse on a incline](bingham_collapse_noSV.c)

 
# Bibliography

 * [Huppert H.](http://www.itg.cam.ac.uk/people/heh/Paper47.pdf)
 
 
 
 */

 