/**
# 1D advection-diffusion problem
Consider the scalar field $C$, defined on a 1D domain $x \in [0,1]$. 

It's evolution is given by:

$$\frac{\partial C(x,Pe)}{\partial x} = \frac{1}{Pe}\frac{\partial^2 C(x,Pe)}{\partial x^2}$$

and we will consider boundary conditions $C(0,Pe) = 0., \ C(1,Pe) = 1.$.

The analytical solution reads as:

$$
C(x,Pe) = \frac{e^{Pe \cdot x} - 1.}{e^{Pe}-1.}
$$

We will compare this analycal solution with the following 
basilisk implementation.
*/
#include "grid/multigrid1D.h"
#include "run.h"
#include "diffusion.h"
#include "advection.h"
scalar C[];
scalar * tracers = {C};
face vector D[];

double Pe = 5;

FILE * output_file; // global output file pointer

scalar CN[]; // RESIDUAL FIELD, required to build up a stopping criterion

/**
Define a function to compute the analytical solution...*/

double Analytical(double x, double Pe) {

	return (exp(Pe*x)-1.0)/(exp(Pe)-1.0);

}


int main() {
	L0 = 1.; // domain size
	X0 = 0.; // smallest X
	init_grid(256); // grid resolution
	DT = 1.0;
	output_file = fopen ("Output.dat", "w");
	run();
	fclose(output_file);
}
/**
Boundary conditions on C field */
C[left] = dirichlet(0.);
C[right]= dirichlet(1.);

event init (t=0) {
// Diffusion coefficient, containing the Peclet number
	foreach_face()
		D.x[]=1./Pe;
// Setting initial conditions
	foreach()
		C[] = 0.;
// Providing the velocity field
	foreach_face()
		uf.x[] = 1.;
// Setting up all of the BCs
	boundary({C});
// Set up the initial RESIDUAL field
  foreach()
		CN[]=1.0;// a dummy value...
}

event integration(i++){
// Advection is performed within the "advection.h" call
// on the _tracer_ field.
// Diffussion instead is done here...
	diffusion(C,dt,D);
}

// This one is for plotting purposes only.
// Please not that DT and final computation time
// are set down here.
event printdata (t += DT; t <= 10.) {
  foreach()
    fprintf (stderr, "%+6.5e %+6.5e %6.5e \n",t, x, C[]);
  fprintf (stderr, "\n\n");
// Add a stopping criterion, if convergence is achieved earlier...
  double RES = change(C,CN);
	if(i>1 && RES < 1e-3){// i.e. if I'm converging towards a steady solution...
		foreach()
			fprintf(output_file,"%+6.5e %+6.5e %+6.5e \n",x,C[],Analytical(x,Pe));
		fflush(output_file);
		return 1; /*stop*/
	}
}

/**
### 1D advection-diffusion, general view
~~~gnuplot 1D advection diffusion, general plot
set grid
set xlabel "x"
set ylabel "C(x,Pe)"
set title "1D advection diffusion"
set key top left

# Plot data from file, only lines with column 6 equal to 0 (i.e. yShift = 0)
plot 'Output.dat' using ($1):($2) with points pt 7 ps 2 lc rgb "black" title "numerical",\
     'Output.dat' using ($1):($3) with lines lw 5 lc rgb "red" title "analytical"
~~~
Nice match! But can we be more _quantitative_?

### 1D advection-diffusion, residual
~~~gnuplot 1D advection diffusion, general residual
set grid
set xlabel "x"
set ylabel "(C(x,Pe)-A(x,Pe))/A(x,Pe)"
set title "1D advection diffusion"
set key top left
set logscale y

# Plot data from file, only lines with column 6 equal to 0 (i.e. yShift = 0)
plot 'Output.dat' using ($1):(($2-$3)/$3) with points pt 7 ps 2 lc rgb "black" title "Residual"
~~~
Residual here is high for x=0...
...because analytical value must be zero ;)
*/
