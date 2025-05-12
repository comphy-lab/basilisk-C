/**
# Rayleigh-Bénard convection in 2D
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

#include "embed.h"
#include "convection_boussinesq.h"
#include "view.h"

#define minlevel 9
#define maxlevel 9

double EndTime= 40.;
scalar yref[];

int main() {
	L0 = 1.;
	X0 = Y0 = Z0 = -0.5;
	DT = 0.01;
	TOLERANCE = 1e-6;
	Ra = 1e7; Pr = 0.71; N = 1<<minlevel ; run();
}

T[top] = dirichlet(-0.5);
T[bottom] = dirichlet(0.5);

yref[top] = dirichlet(-0.5);
yref[bottom] = dirichlet(0.5);

u.n[top] = dirichlet(0.);
u.n[bottom] = dirichlet(0.);
u.n[right] = dirichlet(0.);
u.n[left] = dirichlet(0.);

u.t[top] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);
u.t[right] = dirichlet(0.);
u.t[left] = dirichlet(0.);

T[embed]  = dirichlet(-sign(y)*0.5);
u.n[embed]  = dirichlet(0.);
u.t[embed]  = dirichlet(0.);


event init (t=0) {

	vertex scalar phi[];
  foreach_vertex()
		phi[] = ((y > x - 0.75) && (y > -x - 0.75)) && ((y < x + 0.75) && (y < -x + 0.75)) ;
	boundary ({phi});
  fractions (phi, cs, fs);
	fractions_cleanup (cs, fs);

	foreach(){
		T[] = ((y > x - 0.75) && (y > -x - 0.75)) && ((y < x + 0.75) && (y < -x + 0.75)) ?  -y*cs[] : -sign(y)*0.5;

		foreach_dimension()
			u.x[] = 0.;
	}
	boundary ({T,u});
}


// event adapt (i+=1) {
// 	adapt_wavelet((scalar *){T,u}, (double[]){1e-4,1e-4,1e-4}, maxlevel);
// }

#include "output_fields/output_vtu_foreach.h"
void backup_fields (scalar T, vector u, int nf)
{
	FILE * fp ;
	char name[80], stamp[1024];
	sprintf(stamp, "%6.6i", nf);

	nf > 0 ? sprintf(name, "rb2d_%6.6d_n%3.3d.vtu", nf,pid()) : sprintf(name, "rb2d_n%3.3d.vtu",pid());
	fp = fopen(name, "w");
	output_vtu_bin_foreach ((scalar *) {T}, (vector *) {u}, N, fp, false);
	fclose (fp);
}
/**
We also define a function to plot the temperature and streamfunction using
bview.
*/
	void plot_fields(scalar T)
	{
		draw_vof ("cs", filled = -1, fc = {1,1,1});
		squares ("T", linear = true);
		isoline ("T", n = 21, min=-0.5, max=0.5);
		box();
		save ("rb2d_T.png");
	}

scalar un[];
event init_un (i = 0) {
	foreach()
		un[] = u.x[];
}

event logfile (t += 1.0; t <= EndTime) {
	double deltau = change (u.x, un);
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

event wrapup (t += 5.0; t <= EndTime) {
	static int nf=0;
	backup_fields(T,u,nf);
	plot_fields(T);
	nf++;
}

event movie (t += 0.1; t <= EndTime)
{
	draw_vof ("cs", filled = -1, fc = {1,1,1});
  squares ("T", linear = true, min=-0.5, max=0.5);
  isoline ("T", n = 21, min=-0.5, max=0.5);
  box();
  save ("movie.mp4");
}
