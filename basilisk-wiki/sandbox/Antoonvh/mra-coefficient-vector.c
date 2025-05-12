/**
![](mra2D.png)

# One-dimensional multi-resolution coefficients in two dimensions

One-dimensional prolongation to obtain multi-resolution
analysis(MRA)-coefficients is easy and [well
understood](refineing1d.c). We may extend it to more dimensions by
introducing a multi-resolution-coeficient vector to dilineate between
the dimensions. After "full" restriction, first half of a parent
siblings are restricted to obtain a coarse-cell alliged estimate. Then
the one-dimensional techniques can be applied.
*/

#include "utils.h"

vector w[];
scalar s[];
/**
First, we define functions that compute the average of the child cells
at the lhs and rhs of a parent cell. This is extended to the other
dimensions using the `foreach_dimension()`
[rotation](/src/rotate.lex). such that: `left_average_x ->
bottom_average_y -> front_average_z` and, `right_average_x ->
top_average_y -> back_average_z`.
*/
foreach_dimension() {
  double left_average_x (Point point, scalar s) {
#if dimension == 1 //Dont use in 1D
    return fine(s,0,0,0);
#elif dimension == 2
    return (fine(s,0,0,0) + fine(s,0,1,0))/2.;
#else //dimension == 3
    return (fine(s,0,0,0) + fine(s,0,1,0) +
	    fine(s,0,0,1) + fine(s,0,1,1))/4.;
#endif  
  }
  
  double right_average_x (Point point, scalar s) {
#if dimension == 1 //Dont use in 1D
    return fine(s,1,0,0);
#elif dimension == 2
    return (fine(s,1,0,0) + fine(s,1,1,0))/2.;
#else //dimension == 3
    return (fine(s,1,0,0) + fine(s,1,1,0) +
	    fine(s,1,0,1) + fine(s,1,1,1))/4.;
#endif
  }
  /**
Next, the 2,3,4 and 5 point prolongations are implemented. Note that
these are also rotated to the higher dimensions. The 1-point
(injection) method is not subject to rotation. 
   */
  static inline double bilinear_x (Point point, scalar s) {
    return (3.*coarse(s,0,0,0) + coarse(s,child.x,0,0))/4.;
  }

  static inline double linear_x (Point point, scalar s) {
    return (coarse(s,0,0,0) +
	    (double)child.x*(coarse(s,1,0,0) - coarse(s,-1,0,0))/8.);
  }

  static inline double fourth_order_x (Point point, scalar s) {
    return (-3*coarse(s, 2*child.x,0,0) + 17*coarse(s, child.x,0,0) +
	    55*coarse(s, 0,0,0) + -5*coarse(s, -child.x,0,0))/64.;
  }

  static inline double fifth_order_x (Point point, scalar s) {
    return (-3*coarse(s,2*child.x,0,0) + 22.*coarse(s, child.x,0,0) +
	    128.*coarse(s, 0,0,0) + -22*coarse(s, -child.x,0,0)
	    + 3.*coarse(s, -2*child.x,0,0))/128.;
  }
}
/**
# Computing the wavelet vector

The following function fills the $Z$-components of the MRA-coefficient
vector for the $Z$ dimensional scalar field. Using the prolongation
method with `order`-th accuracy.
 */

trace
void wavelet_vector (scalar s, vector w, int order) {
  restriction ({s});
  /**
     It is computed in each cell in the grid, including (grand) parent
     cells.
   */
  for (int l = 0; l < depth(); l++) {
    foreach_level (l) {
      /**
We write the method for the $x$ direction and rotate to the other
dimensions using the aforementioned `foreach_dimension()` code
wrapper.
      */
      foreach_dimension() {
	// Store the lhs and rhs values for the moment.
	double l = left_average_x  (point, s);
	double r = right_average_x (point, s);
	foreach_child() {
	  if (order < 2)      // 1-point: Injection
	    w.x[] = child.x*(child.x < 0 ?          //lhs?
			     l - coarse(s,0,0,0) :  //Yes.
			     r - coarse(s,0,0,0));  //no.
	  else if (order < 3) // 2-point: bilinear
	    w.x[] = (child.x < 0 ?
		     l - bilinear_x (point, s) :
		     r - bilinear_x (point, s));
	  else if (order < 4) // 3-point: linear
	    w.x[] = child.x*(child.x < 0 ?
			     l - linear_x (point, s) :
			     r - linear_x (point, s));
	  else if (order < 5) // 4-point: 4th order
	    w.x[] = (child.x < 0 ?
		     l - fourth_order_x (point, s) :
		     r - fourth_order_x (point, s));
	  else                // 5-point: 5th order
	    w.x[] = child.x*(child.x < 0 ?
			     l - fifth_order_x (point, s) :
			     r - fifth_order_x (point, s));
	}
      }
    }
  }
}

int main() {
  L0 = 10.;
  X0 = Y0 = Z0 = -L0/2;
  /**
     Start of the convergence experiment:
   */
  for (N = 32; N <= 1024 ; N*=2) {
    init_grid (N);
    foreach()
      s[] = exp (-sq(x) - sq(2.*y) - sq(z/2.));
    /**
    we need ghost-cell values at the boundaries for `order` > 1.
    */
    boundary ({s});
    printf("%d", N);
    /**
       Testing the five different prolongation operators;
    */
    for (int i = 1; i < 6; i++) {
      wavelet_vector (s, w, i);
      foreach_dimension()
	printf("\t%g", normf(w.x).avg);
      if (N == 512) {// We output images
	char fname[99];
	sprintf (fname, "wx_%d.png", i);
	output_ppm (w.x, file = fname);
	sprintf (fname, "wy_%d.png", i);
	output_ppm (w.y, file = fname);
      }
    }
    fputs ("\n", stdout);
  }
}

/**
## Results

It "works" as intended:

~~~gnuplot Convergence for the various methods in the x direction
set xr [23:1500]
set logscale xy 2
set key left bottom box
set grid
set xlabel 'N'
set ylabel 'avg. wavelet coef.'
plot 'out' u 1:2 t 'w.x injection',\
     'out' u 1:4 t 'w.x bilinear',\
     'out' u 1:6 t 'w.x linear',\
     'out' u 1:8 t 'w.x 4th',\
     'out' u 1:10 t 'w.x 5th',\
     0.2*x**(-1) w l lw 3 t 'first order',\
     5000*x**-5 w l lw 3 t 'fifth order'
~~~


~~~gnuplot Convergence for the various methods in the y direction
set xr [23:1500]
set logscale xy 2
set key left bottom box
set grid
set xlabel 'N'
set ylabel 'avg. wavelet coef.'
plot 'out' u 1:3 t 'w.y injection',\
     'out' u 1:5 t 'w.y bilinear',\
     'out' u 1:7 t 'w.y linear',\
     'out' u 1:9 t 'w.y 4th',\
     'out' u 1:11 t 'w.y 5th',\
     0.3*x**(-1) w l lw 3 t 'First order',\
     50000*x**(-5) w l lw 3 t 'Fifth order'
~~~

Furtheremore, we may visually inspect the 10 structures of $^z\chi_i$:

![$^1\chi_x$](mra-coefficient-vector/wx_1.png)

![$^2\chi_x$](mra-coefficient-vector/wx_2.png)

![$^3\chi_x$](mra-coefficient-vector/wx_3.png)

![$^4\chi_x$](mra-coefficient-vector/wx_4.png)

![$^5\chi_x$](mra-coefficient-vector/wx_5.png)

![$^1\chi_y$](mra-coefficient-vector/wy_1.png)

![$^2\chi_y$](mra-coefficient-vector/wy_2.png)

![$^3\chi_y$](mra-coefficient-vector/wy_3.png)

![$^4\chi_y$](mra-coefficient-vector/wy_4.png)

![$^5\chi_y$](mra-coefficient-vector/wy_5.png)

Looks good!
 */
