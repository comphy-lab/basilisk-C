/**
# Culvert flushing

   In this example, we placed a culvert structure in the center of the field closed 
   by a dam in the computation domaine. */


#include "saint-venant.h"
#include "culvert_struct.h"
#include "view.h"
/**
   The domain is 10 metres squared, centered on the origin. Time is in
   seconds. */


/**
   ![Bed elevation](culvert/topo.png)
*/

scalar etamax[], hmax[], zbmax[], umax[];

struct culvert c1;

double q_culvert;
double tmax,XMAX;

int LEVEL;

int main()
{
  L0 = 10.0; 
  X0 = 0.;
  Y0 = 0. ;
  G = 10.;
  tmax = 80.;

  /**
     The structure "culvert" needs 13 parameters in order to compute 
     the culvert dicharge "q_culvert" at each time-steps: \n
     the width, inlet and outlet bed elevations (CBE1,CBE2), the coordinates of the inlet and outlets, the culvert LENGTH, the contraction coefficient $\mu$ and inlet, outlet and linear friction coefficients.    */

  c1.x_inlet=   2.5;
  c1.y_inlet =  5.0;
  c1.x_outlet = 7.5;
  c1.y_outlet = 5.0;
  c1.LENGTH = pow(pow(c1.x_inlet-c1.x_outlet,2)+pow(c1.y_inlet-c1.y_outlet,2),1/2);
  c1.WIDTH = 1.;
  c1.CBE1 = 0;
  c1.CBE2 = 0;
  c1.C56 = 1.;
  c1.C1 = 1.3;
  c1.C2 =1.2;
  c1.C3 =1.5;
  c1.W = 1.;
  c1.mu = pow((c1.C1+c1.C2+c1.C3),(-1./2.));

 
  LEVEL = 6;
  N = 1 << LEVEL;
  run();
  LEVEL = 7;
  N = 1 << LEVEL;
  run();
  LEVEL = 8;
  N = 1 << LEVEL;
  run();
}

double etas;
u.n[right] = neumann(0);

event init (t = 0)
{
  foreach() {
    zb[] = (x>5.5) ? -0.5*x+3.25:10;
    zb[] = (x<4.5) ? 0:zb[];
    h[] = (x<4.5) ? 5.:0;
  }

}
/**
   The culvert flux is then computed : */
event inflow (i++) {
  scalar l[];
  foreach()
    l[] = level;
  q_culvert = culvert_flux(c1);
  /**
## Implementation


     At each time-step we can then estimate the discharge as a function of the culvert geometry and the upstream and downstream water depth: 
     $$Q_{c} = f_{c}(h_{1},h_{2}).$$

     To implement culvert discharge calculations equation in the source-term coupling approach (see figure 1), we need to define a splitting as follows:
	
     $$h_{1,k}=h_{1,k}-\dfrac{Q_{c}}{p}\;\;\;\;\; k \in \{1,...,p\},$$
     $$h_{2,k}=h_{2,k}+\dfrac{Q_{c}}{p}\;\;\;\;\; k \in \{1,...,p\},$$

     with $h_{1,k}$, $k \in \{1,...,p\}$ the $p$ inlet cells (well) and $h_{2,k}$, $k \in \{1,...,p\}$ the $p$ outlet cells (source).

  */
  c1.Dx_inlet =  L0/pow(2,interpolate(l,c1.x_inlet,c1.y_inlet));
  c1.Dx_outlet =  L0/pow(2,interpolate(l,c1.x_outlet,c1.y_outlet));
  foreach(){
    h[]-=(((fabs(x-c1.x_inlet)<=c1.Dx_inlet/2) && (fabs(y-c1.y_inlet)<=c1.Dx_inlet/2) && (x-c1.x_inlet<c1.Dx_inlet/2) && (y-c1.y_inlet<c1.Dx_inlet/2))  ? q_culvert*dt/c1.Dx_inlet/c1.Dx_inlet :0 );
    h[]+=(((fabs(x-c1.x_outlet)<=c1.Dx_outlet/2) && (fabs(y-c1.y_outlet)<=c1.Dx_outlet/2) && (x-c1.x_outlet<c1.Dx_outlet/2) && (y-c1.y_outlet<c1.Dx_outlet/2)) ? q_culvert*dt/c1.Dx_outlet/c1.Dx_outlet :0 );  }
}

/**
## Outputs

   We compute the water flow in the culvert and make GIF movies. */

event movie(t+=2.) { 
  view(width=800, height=800,tx=-0.5,ty=-0.5); 
 clear(); 
 squares("h", min=0, max=3, linear=1); 
 isoline("h", 0.01, lw=2.0); 
 cells();
 char timestring[15]; 
 sprintf(timestring, "t=%g", t); 
 draw_string(timestring, lw = 2);
 char name[80];
 sprintf (name, "h-%d.gif",LEVEL);
 save(name); 
}




event logfile (i++) {
  char name[80];
  sprintf (name, "log-%d", LEVEL);
  static FILE * fp = fopen (name, "w");
  fprintf(fp, "%g %g %g %g %g\n", t,q_culvert,culvert_FT(c1),interpolate(h,c1.x_inlet,c1.y_inlet),interpolate(h,c1.x_outlet,c1.y_outlet));
}


event outputtopo (t=1) {
  output_ppm (zb, min = -2, max = 1, file = "topo.png");}



event end (t = tmax) {
  printf ("i = %d t = %g\n", i, t);
}
/**
## Results

### Analytic solution FT = 5

As the flowtype is initially 5, an analytic solution would be solving the equation (see [culvert_struct.h](/sandbox/zied/culvert_flushing/culvert_struct.h)): 

$S_0\frac{dh}{dt} + mD^2\sqrt{2gh} = 0$ 

with $S_0$ the basin area on the left side of the topographic dam, $m$ the contraction coefficient and $D$ the culvert width.

Let's define the following variables:


$\alpha = S_0$

$\beta = mD^2\sqrt{2g}$

The solution h is then:

$h(t) = \frac{\beta^2 t^2}{4\alpha^2} - \frac{C\beta^2}{2\alpha}t + \frac{C^2\beta^2}{4\alpha^2}$

with  $C = \frac{2S_0\sqrt{h_0}}{mD^2\sqrt{2g}}$ in order to have $h(0) = h_0$


So:

$h(t) = \frac{m^2D^4g}{2S_0^2}t^2 - \frac{mD^2\sqrt{2gh_0}}{S_0}t + h_0$

and 

$Q_c(t) = mD^2\sqrt{2g h(t)}$

Using the previous analytic solution for flowtype 5, we can establish: 



### Analytic solution FT = 2

Those equations are valid until a time $T_{5,2}$ when the water depth reaches $1.5D$, the flowtype switches then to 2 and the associated equation is $Q_c(t) = mh_c D \sqrt{2g(S_1 - z_2 + h_c)} = \frac{2}{3} m D \sqrt{\frac{10}{3}g}   h^{\frac{3}{2}}$ in the present test case.

We can obtain $T_{5,2}$ using the analytic solution for the first flowtype,The equation provided is a quadratic equation in terms of $h(t)$. The general form of a quadratic equation is $ax^2 + bx + c = 0$, where a, b, and c are constants and x is the variable. In this case, the variable is $h(t)$ and the constants are:

$a = \frac{m^2D^4g}{2S_0^2}$
$b = -\frac{mD^2\sqrt{2gh_0}}{S_0}$
$c = h_0-1.5D$

To solve the equation, we can use the quadratic formula:
$x = \frac{-b \pm \sqrt{b^2-4ac}}{2a}$

Therefore, the solution of the equation is:

$T_{5,2} = \frac{\frac{mD^2\sqrt{2gh_0}}{S_0} - \sqrt{\frac{m^2D^4g}{S_0^2}h_0 - 1.5mD^3\sqrt{2g}}}{\frac{m^2D^4g}{2S_0^2}}$



~~~gnuplot Evolution of the flow rate in the culvert
set size ratio 0.5
set key top left
set xlabel 'Time (s)'
set ylabel 'Flow rate (m^3/s)'
set grid xtics ytics mxtics mytics
set xrange [0:100]
set yrange [0:6]
set style line 1 lc rgb '#0060ad' lt 1 lw 1
set style line 2 lc rgb '#dd181f' lt 1 lw 1
set style line 3 lc rgb '#00aa00' lt 1 lw 1
set style line 4 lc rgb 'black' lt 1 lw 2 dt 2
q(x)= (x <40.7) ? 1./2.*sqrt(2*10*((x**2)/1620.-(201.25/(810.*sqrt(5.)))*x+5.)):1./6.*sqrt(200./3*(4374./(x**2+26.6*x+176.9))**3)
plot './log-8' u 1:2 w l ls 1 t 'lvl 8',\
'./log-7' u 1:2 w l ls 2 t 'lvl 7',\
'./log-6' u 1:2 w l ls 3 t 'lvl 6',q(x) w l ls 4 t 'Analytic Solution'
~~~

~~~gnuplot Evolution of the flow rate
set key top right
set x2label 'Time (s)'
set y2label 'Flow rate(m^3/s)'
set ylabel 'Water depth (m)'
set yrange [0:5]
set y2tics
set style line 1 lc rgb '#0060ad' lt 1 lw 1.5
set style line 2 lc rgb '#dd181f' lt 1 lw 1.5
set style line 3  lt 21 dt 2
set style line 4  lt 20 dt 3
set style line 5 lc rgb '#666666' lt 1 lw 1.5
set style line 6 lc rgb 'black' lt 1 lw 2 dt 2
set cblabel 'Flow type'
set cbtics offset character -1, 0, 0
set cbrange [1:6]
set palette defined (1 "red", 2 "orange", 3 "yellow", 4 "green", 5 "blue", 6 "purple")
q(x)= (x <40.7) ? 1./2.*sqrt(2*10*((x**2)/1620.-(201.25/(810.*sqrt(5.)))*x+5.)):1./6.*sqrt(200./3*(4374./(x**2+26.6*x+176.9))**3)

plot './log-6' every 100 using 1:4 with points pointtype 9 ps 0.5 lc rgb 'blue' t 'Upstream water depth',\
'./log-6' every 100 using 1:(2./3.*$4) w l lc rgb 'blue' dt 2 t 'Critical water depth h_c',\
'./log-6' every 100 using 1:5 with points pointtype 10 ps 0.5 lc rgb 'black' t 'Downstream water depth',\
'./log-6' u 1:2:3 w l ls 1 lw 2 lc palette axis x1y2 t 'Q_c(t)',\
q(x) w l ls 6 axis x1y2 t 'Flowrate analytic Solution'

~~~


Water depth evolution for multiple resolution:

![LEVEL = 6](culvert/h-6.gif)

![LEVEL = 7](culvert/h-7.gif)

![LEVEL = 8](culvert/h-8.gif)



## See also

* [Flow rates for multiple rivers](multiriverinflow.c)
*/
