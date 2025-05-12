/**
# Prolongation/Refinement in 1D

We test the accuracy of refinement/prolongation methods in 1D. 1st,
2nd and 3rd order-accurate methods are defined under
[`src/grid/`](/src/grid/multigrid-common.h). Aditionally a 4th and
5th-order method are defined below.
*/

#include "grid/bitree.h" 

static inline void refine_order4 (Point point, scalar s) {
  foreach_child()
    s[] = 1./64.*(-3*coarse(s, 2*child.x) + 17*coarse(s, child.x) +
		  55*coarse(s, 0) + -5*coarse(s, -child.x));
}

static inline void refine_order5 (Point point, scalar s) {
  double a =  (-3*s[-2] + 22.*s[-1] + 128.*s[] + -22*s[1] + 3.*s[2])/128.;
  double b = 2*s[] - a;
  foreach_child() {
    if (child.x == -1)
      s[] = a;
    else
      s[] = b;
  }
}

/**
  The test function is the compact Gaussian-'bell' shape ($s =
  e^{-x^2}$), `#define`d via the grid-cell averaged value (`GCA`).
*/

#define GCA (sqrt(pi)/2.*(erf(x + (Delta/2)) - erf(x - (Delta/2)))/Delta)

scalar s1[], s2[], s3[], s4[], s5[],
  * list = {s1, s2, s3, s4, s5};

int main() {
  s1.refine = refine_injection;
  s2.refine = refine_bilinear;
  s3.refine = refine_linear;
  s4.refine = refine_order4;
  s5.refine = refine_order5;
  L0 = 20;
  X0 = -L0/2;
  init_grid (16);
  foreach()
    for (scalar s in list)
      s[] = GCA;
  /**
     The convergence test is started:
  */
  for (int lev = depth() + 1; lev <= 19; lev++) {
    refine (level < lev); //Initialize new vualues using `si.refine`
    double e[5]  = {0., 0., 0., 0., 0.};
    double em[5] = {0., 0., 0., 0., 0.};
    foreach(){
      int j = 0;
      for (scalar s in list) {
	double el = fabs(s[] - GCA); 
	e[j] += el*Delta;
	if (el > em[j])
	  em[j] = el;   
	s[] = GCA;
	j++;
      }
    }
    /**
       printing the results; 
    */
    printf ("%d", 1 << depth());
    int j = 0;
    for (scalar s in list) {
      printf ("\t%g\t%g", e[j], em[j]);
      j++;
    }
    printf ("\n");
  }
}

/**
# Results

~~~gnuplot The results look good
set xr [16 : 1e6]
set xlabel 'N'
set ylabel 'L1 error'
set grid
set logscale xy 2
set key bottom left box on
plot 'out' u 1:2 t '1 point', 'out' u 1:4 t '2 point', 'out' u 1:6 t '3 point',	\
     'out' u 1:8 t '4 point', 'out' u 1:10 t '5 point'
~~~

Note that there exists an *optimimal* grid size the balances rounding
errors versus resolution.
 */
