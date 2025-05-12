/**
# Differentially heated cavity in 2D
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

/**
We define the following function to write our fields using binary files in a
gnuplot-compatible format as well as in *VTK format - UnstructuredGrid (.vtu)*
files. Fields are stored at the end of the simulation (t=EndTime) or when
a steady-state is reached.
*/

	#include "../output_fields/output_vtu_foreach.h"
	void backup_fields (scalar T, vector u, int nf)
	{
		char name[80];
		FILE * fp ;
		nf > 0 ? sprintf(name, "cav2d_T_%6.6d.bin", nf) : sprintf(name, "cav2d_T.bin");
		fp = fopen(name, "w");
		output_matrix (T, fp, N, linear = true);
		fclose (fp);

		scalar psi[];
		streamfunction (u, psi);
		boundary ({psi});

		nf > 0 ? sprintf(name, "cav2d_Psi_%6.6d.bin", nf) : sprintf(name, "cav2d_Psi.bin");
		fp = fopen(name, "w");
		output_matrix (psi, fp, N, linear = true);
		fclose (fp);

		nf > 0 ? sprintf(name, "cav2d_%6.6d_n%3.3d.vtu", nf,pid()) : sprintf(name, "cav2d_n%3.3d.vtu",pid());
		fp = fopen(name, "w");
		output_vtu_ascii_foreach ((scalar *) {T,psi}, (vector *) {u}, N, fp, false);
		fclose (fp);
	}

	event logfile (t <= EndTime) {
		backup_fields(T,u,0);
	}
/**

Finally, to identify if the steady-state is reached we follow the maximum
change in the horizontal velocity over two consecutive time units. If the
change is smaller than $10^{-8}$, we write the results and stop the simulation.
For this particular case, we are interested on comparing the heat-flux evaluated
at the left and right walls against reference results from
[G. DE VAHL. DAVIS (1983)].

*/

	scalar un[];
	event init_un (i = 0) {
		foreach()
			un[] = u.x[];
	}
	#include "global_nusselt.h"
	event logfile (t += 1.0; t <= EndTime) {
		static int nf = MINLEVEL;
		double deltau = change (u.x, un);
		if (deltau < 1e-8 && t > 5.0){
			backup_fields(T,u,N);

			FILE * fp ;
			if (nf > MINLEVEL) {
				fp = fopen("cav2d.asc", "a");
			} else {
				fp = fopen("cav2d.asc", "w");
				fprintf (fp, "[1]Ra [2]Pr [3]N [4]nuleft [5]nuright [6]nuvol [7]umin [8]umax [9]vmin [10]vmax\n");
			}

			double nu_l = 0., nu_r = 0, nu_vol=0. ;
			nu_r=nusselt_right(T); nu_l=nusselt_left(T); nu_vol=nusselt_vol(T,u);

	  		stats velx = statsf (u.x); stats vely = statsf (u.y);

			fprintf (fp, "%.9g %.9g %d %.9g %.9g %.9g %.9g %.9g %.9g %.9g \n",Ra,Pr,N,nu_l,nu_r,nu_vol, velx.min, velx.max, vely.min, vely.max );
			fclose (fp);

			refine (level < nf+1);
			if (nf >= MAXLEVEL) return 1;
			nf++; N*=2;
		}
	  scalar div[];
	  foreach() {
	    div[] = 0.;
	    foreach_dimension()
	      div[] += u.x[1] - u.x[];
	    div[] /= Delta;
	  }
	  stats s0 = statsf (div);
		fprintf (stderr, "%f %.9g %.9g %.9g \n", t, deltau, s0.sum/s0.volume, s0.max);
	}
/**
## Results

After processing by gnuplot with [cav2d.plot]() we get

![](cav2d/cav2d.png)

*/
