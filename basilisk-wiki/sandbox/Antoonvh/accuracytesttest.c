/**
# Checking implementations for their intended convergence rate.
   
One may believe it is sufficient to test a solver implementation
for it's convergence properties. One this page we provide a counter
example. The aim is to numerically compute paths of particles in a
(2D) rotational flow field (i.e. $\nabla \times \mathbf{u} = \omega
\neq 0$). Here we test a few time integration schemes for the
equations,

$$\frac{\partial x}{\partial t} = u_x(x,y),$$
$$\frac{\partial y}{\partial t} = u_y(x,y).$$

The test case concerns a particle in a flow field that is described
by a solid-body rotation ($\omega = 1$):

$$\frac{\partial x}{\partial t} = -y,$$
$$\frac{\partial y}{\partial t} = x,$$

subject to initial conditions $\{x_0, y_0\} = \{1, 0\}$. The
analytical solution is a circular path $x = \mathrm{cos}(t), y =
\mathrm{sin}(t)$.
 */
#include <stdio.h>
#include <math.h>
#define sq(x) ((x)*(x))
#define distance(xp) (pow(sq(xp[0] - 1.) + sq(xp[1]), 0.5))

void dxdt (double x[2], double v[2]){ // Velocities from position
  v[0] = -x[1];
  v[1] = x[0];
}

void advance (double x[2], double v[2], double dt){
  x[0] += dt*v[0];
  x[1] += dt*v[1];
}

int main(){
  FILE * fp = fopen("data", "w"); // File for the convergence data
  /**
     We conduct a convergence study where we increase the number of
     time steps per revolution. (`N`).
  */
  for (int N = 20; N <= 1280; N *= 2){
    /**
       We initialize the particle positions for the two different
       approaches and start the time loop.
    */
    double xfe[2] = {1, 0}; // Forward Euler
    double xmp[2] = {1, 0}; // Mid-point method
    double dt = 2.*M_PI/(double)N;
    for (int i = 1; i <= N; i++){
      /** 
##Forward Euler
	  
Forward Euler is straight forward:
      */
      double v[2];
      dxdt(xfe, v); //Compute velocities
      advance(xfe, v, dt);
      /**
##The mid-point method
	  
We choose not to implement the mid-point method correctly,
i.e. we introduce a bug here. We "forget" to compute the
velocity at $t_n$ for the mid-point-method's
location. Rather, we use the one that is obtained from the
forward-Euler method. Since this path is diverging, one may
expect the result to be *totally* wrong.
      */
      double xs[2] = {xmp[0], xmp[1]}; //A scratch
      //dxdt(xs, v);                   // <- This should be uncommented.
      advance(xs, v, dt/2.);           //Advance to the mid point
      dxdt(xs, v);                     //Compute velocity there
      advance(xmp, v, dt);             //Update original position
    }
    /**
       At the end of the run we compute the distance from the position
       corresponding ot the analytical solution. 
    */
    fprintf(fp, "%d\t%g\t%g\n", N, distance(xfe), distance(xmp));
  }
  /**
     We can view the result of our convergence study:
     
     ~~~gnuplot Error convergence
         set logscale xy 2
	 set xr [10 : 2400] 
	 set xlabel 'timesteps'
	 set ylabel 'Distance Error'
	 set grid
	 set key outside
	 set size square
	 plot 'data' u 1:2 t 'Forward Euler',		\
	 'data' u 1:3 t 'Mid-point method' ,		\
	 10*x**(-1) lw 2 t 'First order',		\
	 100*x**(-2) lw 2 t 'Second order'
     ~~~
      
     Eventough the implementaiton of the mid-point method makes no
     sense, it still converges with second order. As such, one should
     use a more critical test for checking implementations. 
     
## See also

The proper implementation:

* [More tests for particle advection](timetest.c)
  */
}
