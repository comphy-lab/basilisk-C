 /**
# The fonction test of "axi.h"

Here we use a simple axisymetric poiseuille flow to show that axi.h is changing coordinates from cartesienne to cylindrical automatically. 
So in order to compare the axisymetric resultes of 
changing variable (which can done by us with math and solve it in basilisk 2D planar case) , 
and using axi.h (which is totally done by basilisk) ,
we need to write the equation differentely. SYM is a control of using axi.h, the code need to be run 2 times with SYM = 1 & 0.
*/
#define SYM 1
#if SYM  
#include "axi.h"
#endif 
#include "navier-stokes/centered.h"
#include "view.h"

/**
##  main fonction
we periodic the inlet/outlet surface */
int main() { 
  periodic (right);
  TOLERANCE = 1e-4;
  N=64;
  run(); 
}
/**
No slip walls as BC. we set a field un to restore the velocity field of  */
u.n[top]  = dirichlet(0);
u.t[top]  = dirichlet(0);
scalar un[];

/**
## Initial event*/

event init (t = 0) {

/**
### With axi.h
So here comes the differente equations.
after adimensinalization and apply the conditions, we slove the equation 
$$
\frac{\partial ^2 \bar{u}}{\partial \bar{x}^2} + 
\frac{\partial ^2 \bar{u}}{\partial \bar{y}^2} =
 - \frac{\partial \bar{p}}{\partial \bar{x}} = - 1
$$
Writw in basilisk form we simply apply $a =(1,0)$
*/
#if SYM 
 const face vector g[] = {1.,0.};
/**
### Without axi.h
we apply the changement of variable, $x^2 +y^2 =r^2$:
$$
\frac{\partial ^2 \bar{u}}{\partial \bar{r}^2} +
\frac{1}{\bar{r}} \frac{\partial \bar{u}}{\partial \bar{r}} = -1
$$
which can be write in planer form (so we can solve it by basilisk)
$$
\frac{\partial ^2 \bar{u}}{\partial \bar{r}^2} = - 0.5
$$
Write in basilisk form we simply apply $a =(0.5, 0)$
*/
#else
  const face vector g[] = {0.5,0.};
#endif
  a = g;
  mu = fm;
  foreach()
    un[] = u.x[];
}

event logfile (t += 0.1; i <= 100) {
  double du = change (u.x, un);
  if (i > 0 && du < 1e-6)
    return 0; /* stop */
}

event profile (t = end) {
 view (width = 640, height = 640, tx=-0.5,ty= -0.5);
	box();
 squares("u.x", linear = true, min =0., max = 0.25);
#if SYM
  output_field ((scalar *){u}, fopen ("velo_axi", "w"), linear = true);
 //	save("view_axi.png");
   output_ppm (u.x, file = "view_axi.png", min = 0,   max= 0.25 );
#else
  output_field ((scalar *){u}, fopen ("velo_noaxi", "w"),  linear = true);
	save("view_noaxi.png");
#endif
}

/**
# results

![AxiVelocityFieldView](2Dtube/view_axi.png)
![NoaxiVelocityFieldView](2Dtube/view_noaxi.png)

~~~gnuplot compare
plot 'velo_axi' u 2:($3<1? $3:NaN) w lp t'With axi.h' ,\
'velo_noaxi' u 2:($3<1? $3:NaN) w l t'Without axi.h'
~~~
*/