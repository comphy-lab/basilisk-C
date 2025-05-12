/** 
# Testing `profile_foreach_region()`

In this example use the sample results from
[test_init_2Dto3D.c](../input_fields/test_init_2Dto3D.c) to test the
horizontally averaged profiles. An example of the 2D volume fraction field is
shown below

<center>
<table>
<tr>
<td>![](test_profiles_foreach_region/init_f.png){ width="100%" }</td>
<td>![](test_profiles_foreach_region/init_p.png){ width="100%" }</td>
<td>![](test_profiles_foreach_region/init_u.png){ width="100%" }</td>
<td>![](test_profiles_foreach_region/init_v.png){ width="100%" }</td>
</tr>
<tr>
<td><center>f</center></td> 
<td><center>u</center></td> 
<td><center>v</center></td> 
<td><center>p</center></td> 
</tr>
</table>
</center>

The corresponding profile is shown below
~~~gnuplot Horizontally averaged profiles
set multiplot layout 2,2
set ylabel 'y'
set xlabel 'f'
p "profiles.asc" u 2:1 title "Average", "profiles_favre.asc" u 2:1 title "Weighted Average",
set xlabel 'p'
p "profiles.asc" u 3:1 title "Average", "profiles_favre.asc" u 3:1 title "Weighted Average",
set xlabel 'u'
p "profiles.asc" u 4:1 title "Average", "profiles_favre.asc" u 4:1 title "Weighted Average",
set xlabel 'v'
p "profiles.asc" u 5:1 title "Average", "profiles_favre.asc" u 5:1 title "Weighted Average",
~~~

*/

#define MAXLEVEL 9  
#define D0 L0
#include "profiles_foreach_region.h"
#include "view.h"

double rho1 = 3; 
double rho2 = 1; 

// Function to generate a triangle wave
double triangle_wave(double t, double amplitude, double period) {  
  return (4 * amplitude / period) * fabs(t - (period / 2) * floor(2 * t / period + 0.5)) - amplitude/2;
}

int main()
{
  L0 = 1.0;
  X0 = Y0 = Z0 = -L0 / 2;
  N = 1 << MAXLEVEL;
  init_grid(N);

  /** We load the fields from */
  scalar f[], p[];
  vector u[];
  foreach(){
    f[] = (triangle_wave(x,0.33,0.25) > y);
    p[] = (y + 0.1*noise())*f[];    
    u.x[] = (  sq(y) + sin(2*pi*x)*cos(2*pi*y) + 0.1*noise())*f[] + (sin(3*pi*x)*cos(3*pi*y) + 0.1*noise())*(1-f[]);
    u.y[] = (cube(y) + sin(4*pi*x)*cos(4*pi*y) + 0.1*noise())*f[] + (sin(5*pi*x)*cos(5*pi*y) + 0.1*noise())*(1-f[]);
  }
  
  /** Then, visualize the results just to make sure  */

  draw_vof("f");
  squares("f", linear = false);
  save ("init_f.png");

  squares("u.x", linear = false);
  save ("init_u.png");

  squares("u.y", linear = false);
  save ("init_v.png");

  squares("p", linear = false);
  save ("init_p.png");

  /** We define a density */
  scalar rho[];
  foreach(){		    
    rho[] = rho1*f[] + rho2*(1-f[]);
  }
  
  /** Then, we extract the profiles with and without weights */
  scalar * list = {f, p, u};
  profile_foreach_region(list, filename="profiles.asc", ymin=-D0/2., ymax=D0/2., m2=N/8, mode="w");  
  profile_foreach_region(list, rho, filename="profiles_favre.asc", ymin=-D0/2., ymax=D0/2., m2=N/8, mode="w");

}