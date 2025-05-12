/**
# Enforcing the ratio $\frac{R}{\Delta}$ 

Some like to hint at the fidelity of their two-phase simulations in
terms of the dimensionless ratio; $\frac{R}{\Delta} < \text{Const}$,
with $R$ a radius of curfavure of the interface and $\Delta$ the grid
size. Here we implement a method to enforce this ratio along an
interface.

![It seems to work](rc/plot.mp4)

The "challenge" is to make in cooperate with `adapt_wavelet()`,
eventough it is not a wavelet-based criterion.
 */
#include "fractions.h"
#include "curvature.h"
#include "view.h"

double DoC = 1./20.; //Use 20 to 40 cells per radius
/**
`prolongate_ratio()` is defined for this purpose only.
 */
void prolongate_ratio (Point point, scalar s) {
  foreach_child() {
    if (s[] != nodata)
      s[] += s[]*Delta;
  }
}
/**
We test the method with a circle of radius $R$.
 */

#define CIRC(R) (sq(x - R + 0.3) + sq(y) - sq(R))

scalar f[], curv[];

int main() {
  f.prolongation = f.refine = fraction_refine; //f is the fraction field...
  X0 = Y0 = -L0/2;
  init_grid (1 << 7);
  view (width = 700, height = 700);
  //Loop over radii
  for (double R = 10; R > 0.1; R /= 1.01) {
    fraction (f, CIRC(R));                 //Compute drop
    boundary ({f});                        //is this needed?
    curvature (f, curv);                   //Compute curvature
    curv.prolongation = prolongate_ratio;  //_then_ overload the prolongation attribute
    /**
For testing, we do not set an effective maximum level. Note that
eventough we have not used proper boundary conditions for `f`, the
resolution remains limited when the interface touches the domains
boundaries..
     */
    adapt_wavelet ({curv}, (double[]){DoC}, maxlevel = 99);
    //Make a movie
    draw_vof ("f", lc = {1., 0., 1.}, lw = 2);
    cells();
    save ("plot.mp4");
  }

  /**
It is interesting to compare against the "default" refinement based on the vof field.
  */
  for (double R = 10; R > 0.1; R /= 1.01) {
    fraction (f, CIRC(R));              
    boundary ({f});                      
    adapt_wavelet ({f}, (double[]){0.01}, maxlevel = 99);
    
    draw_vof ("f", lc = {1., 0., 1.}, lw = 2);
    cells();
    save ("plot2.mp4");
  }
  /**
     Its appears similar, except we do not know how the
     refinement-citerion value (`0.01`) compares to
     $\frac{R}{\Delta}$. Looking closely, we can see that the
     curvature-based method uses a more constant grid along the
     interface, wheareas the fov-based method is more robust in time.
     
     ![Compare](rc/plot2.mp4)
  */
}
  

