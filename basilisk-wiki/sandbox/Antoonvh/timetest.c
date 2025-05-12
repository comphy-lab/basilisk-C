/**
# A test of four time-integration schemes for particle advection
   
   We aim to numerically compute paths of particles in a (2D)
   rotational flow field (i.e. $\nabla \times \mathbf{u} = \omega 
   \neq 0$). Here we test a few time integration schemes for the
   equations,
   
   $$\frac{\partial x}{\partial t} = u_x(x,y),$$
   $$\frac{\partial y}{\partial t} = u_y(x,y).$$

   The test case concerns a particle in a flow field that is described
   by a solid body rotation ($\omega = 1$):
   
   $$\frac{\partial x}{\partial t} = -y,$$
   $$\frac{\partial y}{\partial t} = x,$$
   
   subject to initial conditions $\{x_0, y_0\} = \{1, 0\}$. The
   analytical solution is a circular path $x = \mathrm{cos}(t), y =
   \mathrm{sin}(t)$. The system appears to be characterized by a
   Hamiltonian: $H = \frac{1}{2} \left( x^2 + y^2 \right)$, such that,

   $$\frac{\partial x}{\partial t} = -\frac{\partial H}{\partial y},$$
   $$\frac{\partial y}{\partial t} = \frac{\partial H}{\partial x}.$$

   Furthermore, $\frac{\partial H}{\partial t}= 0$. We test the
   forward Euler method, the mid-point method, the backward Euler
   method and a x/y-split semi-implicit method. We compute the error
   after one revolution ($t_e = 2\pi$) and trace the path for 20
   revolutions to check if the errors organize in a particular
   pattern.
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
  FILE * fpp = fopen("paths", "w");
  /**
     We conduct a convergence study where we increase the number of
     time steps per revolution. (`N`).
  */
  for (int N = 20; N <= 1280; N *= 2){
    /**
       We initialize the particle positions for the four different
       approaches and start the time loop.
    */
    double xfe[2] = {1, 0}; // Forward Euler
    double xmp[2] = {1, 0}; // Mid-point method
    double xbe[2] = {1, 0}; // Backward Euler
    double xsi[2] = {1, 0}; // x/y-split. Semi-implicit. method
    double dt = 2.*M_PI/(double)N;
    for (int i = 1; i <= 20*N; i++){
      /** 
      ##Forward Euler
	  
      Forward Euler is straight forward:
      */
      double v[2];
      dxdt(xfe, v); //Compute velocities
      advance(xfe, v, dt);
      /** 
      ##The mid-point method
	  
      The mid-piont method requires an additional stratch for the
      location.
      */
      double xs[2] = {xmp[0], xmp[1]};
      dxdt(xs, v);
      advance(xs, v, dt/2.);
      dxdt(xs, v);
      advance(xmp, v, dt);
      /** 
      ##Backward Euler
	  
      Our backward Euler method uses an iterative Newton-Rapson
      procedure to find the velocities corresponding to the next time
      step within a certain `tolerance`. This warrant an additional
      scratch for the velocities. We use the forward Euler solution as
      an initial guess.
      */
      double tolerance = 1e-13; // Very small tolerance
      double error = 1;
      double vs[2];
      int nmax = 40;
      int it = 0;
      dxdt (xfe, vs); 
      while (error > tolerance && it < nmax){
	xs[0] = xbe[0]; xs[1] = xbe[1];
	advance(xs, vs, dt);
	dxdt(xs, v);
	error = sq(vs[0]-v[0]) + sq(vs[1]-v[1]);
	vs[0] = v[0]; vs[1] = v[1];
	it++;
      }
      if (it >= nmax)
	fprintf(stderr, "N-R procedure did not convergence\n");
      advance(xbe, v, dt);
      /**
      ## x/y-split Semi-implicit method
	 
      The equation for the x component is solved in the
      forward-in-time direction wheareas we use a backward formulation
      for the y component. It is rather easy to implement for this
      perticular test case:
      */
      if (i % 2 ==0){
	xsi[0] -= xsi[1]*dt; // x Forward
	xsi[1] += xsi[0]*dt; // y Backward
      }else{
	xsi[1] += xsi[0]*dt; // y Forward
	xsi[0] -= xsi[1]*dt; // x Backward
      }
      /**
      ##Output
	 
      We compute the absolute distance from the analitycal solution
      after a single revolution and write the result to a file

      */
      if (i == N) // t = 2*3.14...
	fprintf(fp, "%d\t%g\t%g\t%g\t%g\n", N, distance(xfe), distance(xmp), distance(xbe), distance(xsi));
      /**
	 We can view the result:

	 ~~~gnuplot Error convergence
	 set logscale xy
	 set xr [ 10:2400] 
	 set xlabel 'timesteps'
	 set ylabel 'Distance Error'
	 set key outside
	 set size square
	 plot 'data' u 1:2 t 'Forward Euler',		\
	 'data' u 1:3 t 'Mid-point method' ,		\
	 'data' u 1:4 t 'Backward Euler' ,		\
	 'data' u 1:5 t 'x/y S-imp. method',		\
	 10*x**(-1) lw 2 t 'First order',		\
	 50*x**(-2) lw 2 t 'Second order'
	 ~~~

	 It Appears the Euler methods are first order accurate and the
	 other two convergerce with $N^{-2}$.
	 
	 It will prove insightfull to see how the errors organize over
	 time and hence we output the particle paths. We study the
	 case where $N = 80$, which lays comfortably inside the
	 convergence domain for all methods.
       */
      if (N == 80) 
	fprintf(fpp, "%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", i,
		xfe[0], xfe[1], sq(xfe[0]) + sq(xfe[1]),
		xmp[0], xmp[1], sq(xmp[0]) + sq(xmp[1]),
		xbe[0], xbe[1], sq(xbe[0]) + sq(xbe[1]),
		xsi[0], xsi[1], sq(xsi[0]) + sq(xsi[1]));
    }
    /**
       We plot the paths

       ~~~gnuplot Particle Paths
       reset
       set xr [-2:2]
       set yr [-2:2]
       set xlabel 'x'
       set ylabel 'y'
	 set size square
	 set key outside
	 plot 'paths' u 2:3 w l lw 2 t 'Forward Euler',		\
	 'paths' u 5:6 w l lw 4 t 'Mid-point method',		\
	 'paths' u 8:9 w l lw 2 t 'Backward Euler',		\
	 'paths' u 11:12 w l lw 2 t 'x/y S-imp. method'
       ~~~

       This inspired to plot the evolution of the Hamiltonians over
       time. 

       ~~~gnuplot Notice the log scale
       set logscale y
       unset xr 
       unset yr 
       set xr[0:1600]
       set yr[0.001 : 1000]
       set xlabel 'iteration'
       set ylabel '2 H'
       plot 'paths' u 1:4 w l lw 2 t 'Forward Euler' ,			\
       'paths' u 1:7 w l lw 4 t 'Mid-point method' ,			\
       'paths' u 1:10 w l lw 2 t 'Backward Euler' ,			\
       'paths' u 1:13 w l lw 2 t 'x/y S-imp. method'
       ~~~
	 
       The forward Euler method displays a monotomically increasing
       Hamiltonian, and the backward-Euler method produces a
       decreasing one. Dr. Delandmeter was right, and the mid-point
       method is neither here nor there. Only the semi-explicit method
       enherites the conservation of the Hamiltonian from the physical
       system, meaning that it is symplectic. This can be a desirable
       property for long-term time integration projects. Noting that
       is should not be confused with accuracy.
    */
  }
}
