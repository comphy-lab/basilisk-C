/**
# Differentially heated cavity in 2D

This is similar to [cav2d.c]() but uses `bview` to generate the outputs, for instance:

![Temperature](cav2d_view/cav2d_T.png "Temperature") ![Stream-function](cav2d_view/cav2d_psi.png "Stream-function")

## Model equations
For this example we use the dimensionless
Boussinesq equations written in the following
form:
$$
\nabla \cdot \vec{u} = 0
$$
$$
\partial_t \vec{u} + \nabla \cdot \left(\vec{u} \otimes \vec{u}
\right) = -\nabla p + \nabla \cdot \left( PrRa^{0.5} \nabla \vec{u} \right) +
Pr\theta\vec{e}_y
$$
$$
\partial_t \theta + \nabla \cdot \left( \vec{u} \theta \right) = \nabla \cdot
\left( Ra^{-0.5} \nabla \theta \right)
$$
Under the Boussinesq approximation our system
is fully defined in terms of two dimensionless
parameters: the *Rayleigh* and *Prandtl* numbers,
in addition to the geometry and boundary conditions.

Instead of writing new code, we combine the centered
the (compressible) Navier-Stokes equations with an
advection-diffusion problem for the temperature field.
*/

	#include "convection_boussinesq.h"
	#include "view.h"

/**
## Dimensionless parameters
For this example for $Ra=10^3$ and $Pr=0.71$ which leads to a steady-state
solution. We use a $16^2$ grid and refine progressively up to $256^2$.
*/

	#define MINLEVEL 4
	#define MAXLEVEL 7

	double EndTime= 300.;
	int main() {
		L0 = 1.;
		X0 = Y0 = -0.5;
		DT = 0.1;
		TOLERANCE = 1e-5;
		Ra = 1e3; Pr = 0.71; N = 1<<MINLEVEL ; run();
	}

/**
## Boundary conditions
The left and right walls are iso-thermal
$$
\theta = \mp 1/2
$$
at $x = \pm 1/2$,
while all other boundaries are adiabatic
$$
\partial_n \theta = 0
$$
on all other walls.
*/
	T[left] 	= dirichlet(-0.5);
	T[right] 	= dirichlet(0.5);
/**
In addition to no-slip walls
$$
\vec{u} = 0
$$
on all boundaries
*/
	u.n[top] 		= dirichlet(0.);
	u.n[bottom] = dirichlet(0.);
	u.n[right] 	= dirichlet(0.);
	u.n[left] 	= dirichlet(0.);

	u.t[top] 		= dirichlet(0.);
	u.t[bottom] = dirichlet(0.);
	u.t[right] 	= dirichlet(0.);
	u.t[left] 	= dirichlet(0.);
/**
## Initial conditions
Initial conditions correspond to a linear temperature
profile and no motion.
*/

	event init (t=0) {
		foreach(){
			T[] = x;
			foreach_dimension()
				u.x[] = 0.;
		}
		boundary ({T,u});
	}
/**
## Outputs

For this problem we extract the temperature field
$T(x,y)$ and the stream-function $\psi(x,y)$
$$
\delta\psi(x,y) = u \delta y - v \delta x
$$
This function is obtained by solving a poisson problem for the vorticity.
*/

void streamfunction (vector u, scalar psi) {
	scalar omega[];
	vorticity (u, omega);
	boundary ({psi, omega});
	poisson (psi, omega);
}

/**## Results
To identify if the steady-state is reached we follow the maximum change in the
horizontal velocity over two consecutive time units. If the change is smaller
than $10^{-8}$, we write the results and stop the simulation.
*/
	scalar un[];
	event init_un (i = 0) {
		foreach()
			un[] = u.x[];
	}

	event save_and_refine (t += 1.0; t <= EndTime) {
		static int nf = MINLEVEL;
		double deltau = change (u.x, un);
		if (deltau < 1e-8 && t > 5.0){

			char name[80];
			sprintf(name, "cav2d_l%.2d.dump", nf);
			dump (name);

			if (nf >= MAXLEVEL) return 1;
			refine (level < nf+1);
			nf++; N*=2;
		}
	}

/**
Finally, we output the temperature and streamfunction fields at the steady
state
*/

	event pictures (t = end) {
		view(width = 400, height=400);
		squares ("T", linear = true);
		isoline ("T", n = 21, min=-0.5, max=0.5);
		box();
		save ("cav2d_T.png");

		scalar psi[];
		streamfunction (u, psi);
		boundary ({psi});

		squares ("psi", linear = true);
		isoline ("psi", n = 21);
		box();
		save ("cav2d_psi.png");
	}
