/** 
# Testing `reference_height()`

In this example use the sample results to test the reference_height function
which is useful to highlight density contrasts. 

<center>
<table>
<tr>
<td>![](test_available_potential/init_c.png){ width="50%" }</td>
<td>![](test_available_potential/init_yref.png){ width="50%" }</td>
</tr>
<tr>
<td><center>c</center></td> 
<td><center>yref</center></td> 
</tr>
</table>
</center>

The corresponding profile is shown below
~~~gnuplot (Spatial) cumulative probability function
set ylabel 'y'
set xlabel 'c'
p "reference_state_gsl_0.asc" u 2:3
~~~

*/

#define MAXLEVEL 9  
#include "view.h"
#include "available_potential.h"

double rho1 = 3; 
double rho2 = 1; 


int main()
{
  L0 = 1.0;
  X0 = Y0 = Z0 = -L0 / 2;
  N = 1 << MAXLEVEL;
  init_grid(N);

  /** We load the fields from */
  scalar f[], c[];
  vector u[];
  foreach(){
    f[] = (0 > y);
    c[] = 0.5*( 1. + tanh( (0.25 + 0.05*cos(4*pi*x) + y)/(0.01 + 0.005*cos(6*pi*x))) ) * f[];  
  }
  boundary({c,f});
  
  stats s = statsf(c);

  scalar yref[];
  reference_height(yref, f, c, s.min, s.max, true, L0/2.);

  /** Then, visualize the results  */

  draw_vof ("f", filled = -1, fc = {0.5,0.5,0.5});
  squares("c", linear = false, spread=-1);
  save ("init_c.png");

  draw_vof ("f", filled = -1, fc = {0.5,0.5,0.5});
  squares("yref", linear = false, spread=-1);
  save ("init_yref.png");

}