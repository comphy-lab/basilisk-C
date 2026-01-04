/**
# Differentially heated cavity in 3D
This example is a simple extension of the 2D case presented here [cav2d.c](). Like before the model equations
correspond to the Boussinesq equations.

*/

#include "grid/octree.h"
#include "convection_boussinesq.h"

/**
## Dimensionless parameters
For this example for $Ra=10^3$ and $Pr=0.71$ which leads to a steady-state
solution. We use a $16^3$ grid and refine progressively up to $64^3$.
*/


#define MINLEVEL 4
#define MAXLEVEL 6

double EndTime= 10.;

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
## Outputs

For this problem we extract the temperature field
$T(x,y)$ and the velocity field $\vec{u}(x,y)$.
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
	nf > 0 ? sprintf(name, "cav3d_%4.4d_n%3.3d.vtu", nf,pid()) : sprintf(name, "cav3d_n%3.3d.vtu",pid());
	fp = fopen(name, "w"); output_vtu_bin_foreach ((scalar *) {T}, (vector *) {u}, N, fp, false); fclose (fp);

@if _MPI
	char subname[80];
	if (pid()==0){
		nf > 0 ? sprintf(name, "cav3d_%4.4d.pvtu", nf) : sprintf(name, "cav3d.pvtu");
		nf > 0 ? sprintf(subname, "cav3d_%4.4d", nf) : sprintf(subname, "cav3d");
		fp = fopen(name, "w"); output_pvtu_bin ((scalar *) {T}, (vector *) {u}, N, fp, subname); fclose (fp);
	}
	MPI_Barrier(MPI_COMM_WORLD);
@endif
}

event logfile (t += EndTime/4.;t <= EndTime) {
	backup_fields(T,u,0);
}

/**

Finally, to identify if the steady-state is reached we follow the maximum
change in the horizontal velocity over two consecutive time units. If the
change is smaller than $10^{-6}$, we write the results and stop the simulation.
For this particular case, we are interested on comparing the heat-flux evaluated
at the left and right walls against reference results from
[E. Tric et al. / Int. J. Heat Mass Transfer 43 (2000)].

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
	if (deltau < 1e-6 && t > 5.0){
		backup_fields(T,u,N);
		FILE * fp ;
		if (nf > MINLEVEL) {
			fp = fopen("cav3d.asc", "a");
		} else {
			fp = fopen("cav3d.asc", "w");
			fprintf (fp, "[1]Ra [2]Pr [3]N [4]nuleft [5]nuright [6]nuvol [7]umin [8]umax [9]vmin [10]vmax [11]wmin [12]wmax\n");
		}

		double nu_l = 0., nu_r = 0, nu_vol=0. ;
		nu_r=nusselt_right(T); nu_l=nusselt_left(T); nu_vol=nusselt_vol(T,u);
		stats velx = statsf (u.x); stats vely = statsf (u.y); stats velz = statsf (u.z);

		fprintf (fp, "%.9g %.9g %d %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g\n",Ra,Pr,N,nu_l,nu_r,nu_vol, velx.min, velx.max, vely.min, vely.max, velz.min, velz.max);
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

We process the resulting file cav3d.ptu using paraview to obtain iso-thermal surfaces (left) and velocity vector fields (right).

![](cav3d/cav3d.png)

*/
