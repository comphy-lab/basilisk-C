/**
# Differentially heated cavity in 3D

This is similar to [cav2d.c]() but uses `bview` to generate the outputs, for instance:

![Temperature](cav3d_view/cav3d_T.png "Temperature") ![$\lambda_2$-criterion](cav3d_view/cav3d_l2.png "$\lambda_2$-criterion")

*/

#include "grid/octree.h"
#include "convection_boussinesq.h"
#include "view.h"

/**
## Dimensionless parameters
For this example for $Ra=10^3$ and $Pr=0.71$ which leads to a steady-state
solution. We use a $16^3$ grid and refine progressively up to $64^3$.
*/


#define MINLEVEL 4
#define MAXLEVEL 6

double EndTime= 50.;

int main() {
	L0 = 1.;
	X0 = Y0 = Z0 = -0.5;
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

T[left] = dirichlet(-0.5);
T[right] = dirichlet(0.5);

/**
In addition to no-slip walls
$$
\vec{u} = 0
$$
on all boundaries.
*/

u.n[top] 		= dirichlet(0.);
u.n[bottom] = dirichlet(0.);
u.n[right] 	= dirichlet(0.);
u.n[left] 	= dirichlet(0.);
u.n[front] 	= dirichlet(0.);
u.n[back] 	= dirichlet(0.);

u.t[top] 		= dirichlet(0.);
u.t[bottom] = dirichlet(0.);
u.t[right] 	= dirichlet(0.);
u.t[left] 	= dirichlet(0.);
u.t[front] 	= dirichlet(0.);
u.t[back] 	= dirichlet(0.);

u.r[top] 		= dirichlet(0.);
u.r[bottom] = dirichlet(0.);
u.r[right] 	= dirichlet(0.);
u.r[left] 	= dirichlet(0.);
u.r[front] 	= dirichlet(0.);
u.r[back] 	= dirichlet(0.);

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
## Results

To identify if the steady-state is reached we follow the maximum
change in the horizontal velocity over two consecutive time units. If the
change is smaller than $10^{-6}$, we write the results and stop the simulation.
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
		sprintf(name, "cav3d_l%.2d.dump", nf);
		dump (name);

		if (nf >= MAXLEVEL) return 1;
		refine (level < nf+1);
		nf++; N*=2;
	}
	fprintf (stderr, "%f %.9g \n", t, deltau);
}

/**
Finally, we output the temperature and streamfunction fields at the steady
state
*/
#include "lambda2.h"
event pictures (t = end) {

	view(theta = pi/6, phi = pi/6, width = 400, height=400, fov=35);
	for (float val = -0.5; val < 0.5; val=val+0.1) {
		isosurface ("T", val, color = "T", min = -0.5, max = 0.5,  linear = true, map = cool_warm);
	}
	box(notics=1);
	save ("cav3d_T.png");

	scalar l2[];
	lambda2 (u, l2);
	isosurface ("l2", 0);
	box(notics=1);
	save ("cav3d_l2.png");
}
