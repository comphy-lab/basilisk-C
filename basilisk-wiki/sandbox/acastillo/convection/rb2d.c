/**
# Rayleigh-Bénard convection in 2D

We consider a square (2D) Rayleigh-Bénard cell and perform the
available mechanical energy balance as in (Castillo-Castellanos et al. (2016)).

![](rb2d/rb2d.png "Time-series")

We also show how to enhance the visualization of thermal plumes using the
reference height required to compute the available potential energy
(Sutherland (2010), Tseng & Ferziger (2001))

![Temperature](rb2d/rb2d_T.png "Temperature") ![Reference height](rb2d/rb2d_yref.png "Reference height")

## Model equations
For this example we use the dimensionless Boussinesq equations written in the
following form:
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

#include "convection_boussinesq.h"
#include "view.h"

double EndTime= 100;

int main() {
	L0 = 1.;
	X0 = Y0 = -0.5;
	DT = 0.1;
	TOLERANCE = 1e-5;
	Ra = 5e5; Pr = 0.71; N = 128 ; run();
}

T[top		] = dirichlet(-0.5);
T[bottom] = dirichlet(0.5);

u.n[top		] = dirichlet(0.);
u.n[bottom] = dirichlet(0.);
u.n[right	] = dirichlet(0.);
u.n[left	] = dirichlet(0.);

u.t[top		] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);
u.t[right	] = dirichlet(0.);
u.t[left	] = dirichlet(0.);

event init (t=0) {
	foreach(){
		T[] = -y + 0.03*noise();
		foreach_dimension()
			u.x[] = 0.;
	}
	boundary ({T,u});

	FILE * fp = fopen("rb2d.asc", "w");
	fprintf (fp, "t ekin epot bpot nu_top nu_bot nu_vol nu_eps nu_tmp nu_mix \n");
	fclose (fp);
}

/**
## Results
For this problem we are interested on the the stream-function
$\psi(x,y)$
$$
\delta\psi(x,y) = u \delta y - v \delta x
$$
*/
void streamfunction (vector u, scalar psi) {
	scalar omega[];
	vorticity (u, omega);
	boundary ({psi, omega});
	poisson (psi, omega);
}
/** and the reference height required to compute the available potential energy
(Sutherland (2010), Tseng & Ferziger (2001))
$$
y_{ref} - y_{bot} = \int_{T_{min}}^{T} P(T) dT
$$
*/
scalar yref[];
yref[top	 ] = dirichlet(-0.5);
yref[bottom] = dirichlet(0.5);
#include "available_potential.h"
event reference (i++,last){
	reference_height(yref,T,-0.5,0.5,N);
	boundary ({yref});
}

/** Using the temperature, velocity and reference height fields we may compute
the mechanical energy balance as in (Castillo-Castellanos et al. (2016))
$$
\dfrac{dE_{kin}}{dt}  = Pr Ra^{−0.5}[Nu_{vol} − (\epsilon + 1)]
$$
$$
\dfrac{dE_{pot}}{dt}  = Pr Ra^{−0.5}[-Nu_{vol} + \frac{1}{2}(Nu_{top} + Nu_{bot})]
$$
$$
\dfrac{dE_{bpot}}{dt} = Pr Ra^{−0.5}[\Phi_d − \frac{1}{2}(Nu_{top} + Nu_{bot})]
$$
 */
#define dA Delta
event time_series (t += 0.05; t <= EndTime){

	vector dT[], dyref[];
	gradients ({T}, {dT});
	gradients ({yref}, {dyref});

	tensor du[];
	foreach(){
		du.x.y[] = (u.x[0,1] - u.x[0,-1])/2./Delta;
		du.y.x[] = (u.y[1,0] - u.y[-1,0])/2./Delta;
		du.x.x[] = (u.x[1,0] - u.x[-1,0])/2./Delta;
		du.y.y[] = (u.y[0,1] - u.y[0,-1])/2./Delta;
	}

	double nu_top = 0., nu_bot = 0.;
	foreach_boundary (top,reduction(+:nu_top))
		nu_top += dA*(T[] - T[ghost])/Delta;

	foreach_boundary (bottom,reduction(+:nu_bot))
		nu_bot += dA*(T[ghost] - T[])/Delta;

	double ekin = 0., epot = 0., bpot = 0.;
	foreach(reduction(+:ekin), reduction(+:epot), reduction(+:bpot)){
		epot += -dv() * Pr * T[] * y;
		bpot += -dv() * Pr * T[] * yref[];
		foreach_dimension()
			ekin += dv() * 0.5 * sq(u.x[]);
	}

	double nu_vol=0., nu_eps=0., nu_tmp=0., nu_mix=0.;
	foreach(reduction(+:nu_vol), reduction(+:nu_eps), reduction(+:nu_tmp), reduction(+:nu_mix)){
		nu_vol += dv() * (sqrt(Ra)*u.y[]*T[] - dT.y[]);
		nu_eps += dv() * (sq(du.x.x[]) + sq(du.x.y[]) +  sq(du.y.x[]) + sq(du.y.y[]));
		foreach_dimension(){
			nu_tmp += dv() * (dT.x[]*dT.x[]);
			nu_mix += dv() * (dT.x[]*dyref.x[]);
		}
	}

	FILE * fp = fopen("rb2d.asc", "a");
	fprintf (fp, "%f %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g \n", t, ekin, epot, bpot, nu_top, nu_bot, nu_vol, nu_eps, nu_tmp, nu_mix);
	fclose (fp);
}

event movie (t += 0.1; t <= EndTime)
{
  squares ("yref", linear = true, min=-0.45, max=0.45);
  isoline ("T", n = 11, min=-0.45, max=0.45);
  box();
  save ("movie.mp4");
}

event pictures (t = end) {
	scalar psi[];
	streamfunction (u, psi);
	boundary ({psi});

	view(width = 400, height=400);
	squares ("T", linear = true, min=-0.45, max=0.45);
	isoline ("psi", n = 11);
	box();
	save ("rb2d_T.png");

	squares ("yref", linear = true, min=-0.45, max=0.45);
	isoline ("psi", n = 11);
	box();
	save ("rb2d_yref.png");
}
