/**
This test case solves the unsteady diffusion problem described in Example 8.1 in * Introduction to Computational Fluid Dynamics, The Finite Volume Method: 2nd Edition* by Versteeg and Malalasekara.

*/
#include "grid/multigrid1D.h"
#include "layered/hydro.h"
#include "mltracer.h"


/**
This example considers diffusive heat transport in a one-dimensional rod.
The analytical solution is given by

$$T(x,t) = \frac{4T_0}{\pi} \sum_{n=1}^{\infty} \frac{(-1)^{n+1}}{2n-1}e^{-\alpha \lambda_n^2 t}\cos(\lambda_n x),$$

with $\lambda_n = (2n-1)\pi/(2L)$.

We start by defining the constants used by Versteeg and Malalasekara. 
*/
double kappa_c = 10;
double rho_c   = 1e7;
double L       = 0.02;
double t_max   = 160;
double t_step  = 5;
double c0      = 200;


double c_analytical(double x, double t){
  double ans  = 0; 
  double alph = kappa_c/rho_c; 
  for (int n = 1; n < 500; n++){
    double lambda_n = (2*n - 1)*pi/(2*L0);
    if (n % 2 == 0){
      ans -= exp(-alph*sq(lambda_n)*t )*cos(lambda_n*x)/(2*n - 1);
    } else {
      ans += exp(-alph*sq(lambda_n)*t )*cos(lambda_n*x)/(2*n - 1);
    }
  }
  return 4*c0/pi * ans;
}

/**
The boundary conditions are homogenous Neumann on the left hand side and zero Dirichlet on the right hand side. We define an appropriate boundary condition function here. The default boundary condition is a zero dirichlet condition on all horizontal boundaries. Note, for now, the top and bottom boundaries are hard coded to be impermeable boundaries, i.e no flux through top or bottom. 

TODO: 
Simplify the boundary condition functions. Make them more basilisk-like in how to enforce them. 
*/
double bf(Point point, double t, Ml_Tracer t_idx){
  ml_tracer tracer = tl[t_idx];
  if ( x < L0/2.){
    return tracer.c[];
  } else {
    return 0;
  }
}

/**
Now we initialize the domain, and the scalar and index we will use to create our tracer. 
*/
scalar c;
int t_idx;
int main(){
  nl = 1;
  L0 = L;
  N = 8;
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
      u.x[] = 0.0;
      c[]   = c_analytical(x, 0);
    }
  }
  t_idx = init_tracer(c, kappa_c, rho_c, 1e-4);
  t_idx = set_tracer_BC(t_idx, &bf);
  
  FILE * fp = fopen("beginlog", "w");
  foreach()   
    fprintf (fp,"%g %g %g\n", x, c[], c_analytical(x,0));
  fclose(fp);
}

/**
This event is all the user need to call for tracer transport. 
*/

event transport(i++){
  if (i > 0)
    t_idx = transport_tracer(t_idx);
}

event stop (t = t_max) {
  FILE * fp = fopen("endlog", "w");
  foreach()   
    fprintf (fp,"%g %g %g\n", x, c[], c_analytical(x,t));
  fclose(fp);
  return 1;
}

/**
The solver replicates the analytical solution

~~~gnuplot
set xlabel 'x'
set ylabel 'Temperature'
set grid 
plot 'beginlog' w l lw 2 t 'Initial Temperature' ,\
     'beginlog' u 1:3 t 'Scalar (t = 0)' ,\
     'endlog' w l lw 2 t 'Final Temperature' ,\
     'endlog' u 1:3 t 'Scalar (t = t_{end})'
~~~
*/


