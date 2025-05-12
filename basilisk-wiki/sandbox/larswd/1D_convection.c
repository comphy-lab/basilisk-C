/**
This test case solves the 1D unsteady convection problem given by an initial profile defined by 

$$ c(x,0) = e^{-(x-x_0)^2}.$$
Where $U$ is the velocity in the $x$-direction and $x_0$ is the initial peak concentration. The analytical solution is given by

$$c_a(x,t) = e^{-(x-x_0)^2},$$
and we enforce dirichlet boundary conditions given by $c(0,t) = c_a(0,t)$ and $c(L_0,t) = c_a(L_0,t)$. 

*/
#include "grid/multigrid1D.h"
#include "layered/hydro.h"
#include "mltracer.h"

double U       = 1;
double kappa_c = 0.0;
double rho_c   = 1;
double L       = 4;
double t_max   = 6;
double t_step  = 1.5;
double x0      = 1;


double c_analytical(double x, double t){
  return exp(-sq(x - U*t - x0));
}

double bf(Point point, double t, Ml_Tracer t_idx){
  if ( x < L0/2.){
    return c_analytical(0,t);
  } else {
    return c_analytical(L0,t);
  }
}

/**
Now we initialize the domain, and the scalar and index we will use to create our tracer. We enforce periodic boundary conditions. Note, the boundary conditions on the tracer are entirely separated from the boundary conditions on the multilayer solver. Hence, the domain is not periodic for the tracers unless manually enforced by the user. 
*/
scalar c;
int t_idx;
int k = 0;
int main(){
  nl = 1;
  L0 = L;
  N = 64;
  periodic(right);
  run();
}

/**
Here, we initialize the tracer by setting it equal to the analytical solution at $t=0$. Then, we create a tracer with the chosen coefficients and scalar field and a tolerance of 1e-4, which we then give the boundary function as BC. 
*/

event init(i=0){
  c = new scalar[nl];
  foreach(){
    zb[] = -1;
    foreach_layer(){
      h[]   = -zb[]/nl;
      u.x[] = U;
      c[]   = c_analytical(x, 0);
    }
  }
  t_idx = init_tracer(c, kappa_c, rho_c, 1e-4);
  t_idx = set_tracer_BC(t_idx, &bf);
  

}

/**
This event is all the user need to call for tracer transport. 
*/

event transport(i++){
  if (i > 0)
    t_idx = transport_tracer(t_idx);
}

event log_concentration (t <= t_max; t += t_step) {
  static int j = 0;
  char name[10]; 
  sprintf(name, "log%d",j++);
  FILE * fp = fopen(name, "w");
  foreach()   
    fprintf (fp,"%g %g %g\n", x, c[], c_analytical(x,t));
  fclose(fp);
}

/**
The solver replicates the analytical solution

~~~gnuplot
set xlabel 'x'
set ylabel 'Temperature'
set grid 
plot 'log0' w l lw 2 t 'Initial concentration' ,\
     'log0' u 1:3 t 'Scalar (t = 0)' ,\
     'log1' w l lw 2 t 'Middle concentration' ,\
     'log1' u 1:3 t 'Scalar (t = 1.5)' ,\
     'log4' w l lw 2 t 'Final Temperature' ,\
     'log4' u 1:3 t 'Scalar (t = 6)'
~~~
*/


