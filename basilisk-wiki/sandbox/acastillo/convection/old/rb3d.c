/**
# Rayleigh-Bénard convection in 3D
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
Under the Boussinesq approximation our system is fully defined in terms of two
dimensionless parameters: the *Rayleigh* and *Prandtl* numbers, in addition to
the geometry and boundary conditions.

Instead of writing new code, we combine the centered the (compressible)
Navier-Stokes equations with an advection-diffusion problem for the temperature
field.
*/


#include "grid/octree.h"
#include "convection_boussinesq.h"

#define minlevel 4
#define maxlevel 7

double EndTime= 10.;

int main() {
	L0 = 1.;
	X0 = Y0 = Z0 = -0.5;
	DT = 0.1;
	TOLERANCE = 1e-5;
	Ra = 1e6; Pr = 0.71; N = 1<<minlevel ; run();
}

T[top] = dirichlet(-0.5);
T[bottom] = dirichlet(0.5);

u.n[top] = dirichlet(0.);
u.n[bottom] = dirichlet(0.);
u.n[right] = dirichlet(0.);
u.n[left] = dirichlet(0.);
u.n[front] = dirichlet(0.);
u.n[back] = dirichlet(0.);

u.t[top] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);
u.t[right] = dirichlet(0.);
u.t[left] = dirichlet(0.);
u.t[front] = dirichlet(0.);
u.t[back] = dirichlet(0.);

event init (t=0) {
	foreach(){
		T[] = -y + 0.03*noise();
		foreach_dimension()
			u.x[] = 0.;
	}
	boundary ({T,u});
}


event adapt (i+=1) {
	adapt_wavelet((scalar *){T,u}, (double[]){1e-4,1e-4,1e-4}, maxlevel);
}

#define H5FILE_NAME "fields.h5"
#include "output_fields/output_xmf_foreach.h"
#include "available_potential.h"
void backup_fields (scalar T, vector u, int nf)
{
	FILE * fp ;
	char name[80], stamp[1024];
	sprintf(stamp, "%6.6i", nf);

	scalar yref[];
	reference_height(yref,T,-0.5,0.5,128);

	scalar pid[];
	foreach()
		pid[] = pid();

	nf > 0 ? sprintf(name, "fields_%6.6d.xmf", nf) : sprintf(name, "fields.xmf");
	fp = fopen(name, "w"); output_xmf_h5_foreach ((scalar *) {T,yref,pid}, (vector *) {u}, N, fp, stamp); fclose (fp);
}

scalar un[];
event init_un (i = 0) {
	foreach()
		un[] = u.x[];
}
event logfile (t += 1.0; t <= EndTime) {

	double deltau = change (u.x, un);
	if (deltau < 1e-6 && t > 5.0){
		backup_fields(T,u,0);
		return 1;
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

event wrapup (t += 1.0; t <= EndTime) {
	static int nf=0;
	backup_fields(T,u,nf);
	nf++;
}
