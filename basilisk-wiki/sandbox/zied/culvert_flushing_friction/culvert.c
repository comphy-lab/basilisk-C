/**
   # Culvert in a river

   In this example, we placed a culvert structure in the center of a river closed 
   by a bridge in the computation domaine. */


#include "saint-venant.h"
#include "culvert_struct.h"
#include "view.h"

/**
   The domain is 10 metres squared, centered on the origin. Time is in
   seconds. */


scalar etamax[], hmax[], zbmax[], umax[];

struct culvert c1;

double q_culvert;

double tmax,XMAX;

int LEVEL;

double n = 0;
double manning_culvert;
// fixme: nmanning[] is allocated even n is a constant
// this should use constant fields instead
scalar nmanning[];



/**
   ![Bed elevation](culvert/topo.png)
*/

/**
   The structure "culvert" needs 13 parameters in order to compute 
   the culvert dicharge "q_culvert" at each time-steps: \n
   the width, inlet and outlet bed elevations (CBE1,CBE2), the coordinates of the inlet and outlets, the culvert LENGTH, the contraction coefficient $\mu$ and inlet, outlet and linear friction coefficients.    */


int main()
{
  dry = 1e-4;
  L0 = 10.0; 
  X0 = 0.;
  Y0 = 0. ;
  G = 10;
  tmax = 200.;
  c1.x_inlet=   2.5;
  c1.y_inlet =  5.0;

  c1.x_outlet = 7.5;
  c1.y_outlet = 5.0;

  c1.LENGTH = pow(pow(c1.x_inlet-c1.x_outlet,2)+pow(c1.y_inlet-c1.y_outlet,2),1/2);
  c1.WIDTH = 1.;
  c1.CBE1 = 0;
  c1.CBE2 = 0;

  c1.LHL = 0.;
  c1.FRIC = 0.015;
  c1.C56 = 1.;
  c1.C1 = 1.5;
  c1.C2 =1.2;
  c1.C3 =1.3;
  c1.W = 1.;
  c1.mu = pow((c1.C1+c1.C2+c1.C3),(-1./2.));
  LEVEL = 6;
  N = 1 << LEVEL;

  manning_culvert=0.;

  run();
  
  manning_culvert=0.1;
  
  run();
  
  manning_culvert=0.8;

  run();
  
  manning_culvert=0.3;

  run();
  
  manning_culvert=10.;

  run();



}

double etas;
u.n[right] = neumann(0);

event init (t = 0)
{
  foreach() {
    nmanning[]=manning_culvert;
    zb[] = (x>4.5&&x<5.5) ? 10:0;
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


event manning_friction (i++)
{
  foreach()
    if (h[] > dry) {      
      if (nmanning[] != 0)
	n = nmanning[];      
      double s = 1. + dt*G*sq(n)*norm(u)/pow(h[],4/3.);
      foreach_dimension()
	u.x[] /= s;
    }
  boundary ((scalar *){u});
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
  sprintf (name, "h-%g.gif",manning_culvert);
  save(name); 
}



event logfile (i++) {
  char name[80];
  sprintf (name, "log-%g",manning_culvert);
  static FILE * fp = fopen (name, "w");
  fprintf(fp, "%g %g\n", t,q_culvert);
}

event outputtopo (t=1) {
  output_ppm (zb, min = 0, max = 10, file = "topo.png");}



event end (t = tmax) {
  printf ("i = %d t = %g\n", i, t);
}
/**
   ## Results

   ~~~gnuplot Evolution of the flow rate in the culvert
   set xlabel 'Time (s)'
   set ylabel 'Flow rate(m^3/s)'
   q(x)= (x <90) ? 1./2.*sqrt(2*10*((x**2)/1620-(201.25/(810.*sqrt(5.)))*x+5.)):0
   plot  './log-0' u 1:2 w l t 'n=0',\
   './log-0.3' u 1:2 w l t 'n=0.3',\
   './log-0.1' u 1:2 w l t 'n=0.1',\
   './log-0.8' u 1:2 w l t 'n=0.8',\
   './log-10' u 1:2 w l t 'n=10',q(x)
   ~~~


   Water depth evolution for multiple manning coefficients:

   ![n=0](culvert/h-0.gif)

   ![n=0.1](culvert/h-0.1.gif)

   ![n=0.3](culvert/h-0.3.gif)

   ![n=0.8](culvert/h-0.8.gif)

   ![n=10](culvert/h-10.gif)


   ## See also

   * [Flow rates for multiple rivers](multiriverinflow.c)
   */
