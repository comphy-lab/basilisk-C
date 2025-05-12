/**
# Discharge of a silo through a lateral orifice with a bottom inclination
We propose an implementation of the Jop Pouliquen Forterre µ(I) rheology for the flow in a hourglass.
We propose to identify the relative role of the acceleration and the friction on the discharge flow by varying the silo geometry (bottom inclination angle), adding the influence of a moderate friction on the front and back wall of the silo, and non-slip boundary condition is chosen for the silo bottom. 


# Code
Includes and definitions
*/

#include "grid/quadtree.h"
#include "navier-stokes/centered.h"
#include "vof.h"
// Domain extent
#define LDOMAIN 4.
// heap definition
double  H0, R0, D, W, tmax, Q, Wmin, muwall, WDOMAIN, Hmask, Lz, angle;

/**
passive fluid small density to preserve 0 pressure
and small viscocity
*/
#define RHOF 1e-4
#define mug  1e-5
// Maximum refinement
#define LEVEL 8
//#define MINLEVEL 7
char s[80];
FILE * fpf, *fwq, *fqout;
scalar f[];
scalar * interfaces = { f };
face vector alphav[];
face vector muv[];
scalar rhov[];
/**
Boundary conditions for granular flow, pressure must be zero at the surface.
The pressure is zero in the hole $x=0$ and $0 < y< W$, but the lithostatic gradient is given elsewhere
on the right wall.
No slip boundary conditions on the other walls.
*/
p[top] = dirichlet(0);
u.n[top] = neumann(0);
u.t[bottom] = dirichlet(0);
u.n[bottom] = dirichlet(0);
u.n[right] = dirichlet(0);
u.t[right] = dirichlet(0);
f[right] = neumann(0);
u.t[left] = (fabs(y) >= 0. &&
	fabs(y) <= (W + 0.)) ? neumann(0) : dirichlet(0);
u.n[left] = (fabs(y) >= 0.&&
	fabs(y) <= (W + 0.)) ? neumann(0) : dirichlet(0);
p[left] = (fabs(y) >= 0. &&
	fabs(y) <= (W + 0.)) ? dirichlet(0) : neumann(0);


int main(int argc, char ** argv) {
	L0 = LDOMAIN;
	// number of grid points
	N = 1 << LEVEL;
	// maximum timestep
	DT = 0.005;
	// coefficient of friction of wall
	muwall = 0.1;
	TOLERANCE = 1e-3;
	//silo height
	H0 = 3.84;
	R0 = 20.000;
	// Grain size
	D = 1. / 90.;
	// size of the hole
	W = 0.5;
	fwq = fopen("outWQ", "w");
	fclose(fwq);
	//film size
	Lz = LDOMAIN;
	//Silo thickness
	WDOMAIN = 2.0;
	const face vector g[] = { 0., -1., 0 };
	a = g;
	alpha = alphav;
	mu = muv;
	rho = rhov;
	//inclinaison angle
	angle = 40.;
	Q = 0.;
	//maximum simulation time
	tmax = 20.*0.625 / W;
	fpf = fopen("interface.txt", "w");
	fqout = fopen("qout.txt", "w");
	run();
	fclose(fqout);
	fclose(fpf);
	fprintf(stdout, "\n");
	fwq = fopen("outWQ", "a");
	fprintf(fwq, " %lf %lf \n", W, Q);
	fclose(fwq);
}


/**
initial heap, a rectangle
*/
// note the v
event init(t = 0) {
	mask(x > 1.*L0 / 4. ? right : none);

	scalar phi[];
	foreach_vertex()
		phi[] = min(H0 - y, R0 - x);
	fractions(phi, f);
	/**
	lithostatic pressure, with a zero pressure near the hole
	*/
	foreach()
		p[] = (fabs(y - (W / 2. + 0.)) <= W / 2. &&
		fabs(x) <= .1) ? 0 : max(H0 - y, 0);
}
/**
total density
*/

#define rho(f) ((f) + RHOF*(1. - (f)))
/**
Viscosity computing $D_2=D_{ij}D_{ji}$;

In the pure shear flow
$D_{11}=D_{22}=0$ et $D_{12}=D_{21}=(1/2)(\partial u/\partial y)$,
so that
$D_2=\sqrt{D_{ij}D_{ij}} =\sqrt{ 2 D_{12}^2} = \frac{\partial u}{ \sqrt{2}\partial y}$.
In a pure shear flow, $\partial u/\partial y= \sqrt{2} D_2$.
The inertial number $I$ is $D \sqrt{2} D_2/\sqrt(p)$
and $\mu = \mu_s+ \frac{\Delta \mu}{1+I/I_0}$
the viscosity is $\eta = \mu(I)p/D_2$:

note that if $\eta$ is too small an artificial small viscosity $\rho D \sqrt{gD}$
is taken see Lagrée et al. 11 § 2.3.1
*/
event properties(i++) {
	trash({ alphav });
	scalar eta[];
	foreach() {
		eta[] = mug;
		if (p[] > 0.) {
			double D2 = 0.;
			foreach_dimension() {
				double dxx = u.x[1, 0] - u.x[-1, 0];
				double dxy = (u.x[0, 1] - u.x[0, -1] + u.y[1, 0] - u.y[-1, 0]) / 2.;
				D2 += sq(dxx) + sq(dxy);
			}
			if (D2 > 0.) {
				D2 = sqrt(2.*D2) / (2.*Delta);
				double In = D2*D / sqrt(p[]);
				double muI = .4 + .28*In / (.4 + In);
				double etamin = sqrt(D*D*D);
				eta[] = max((muI*p[]) / D2, etamin);
				eta[] = min(eta[], 100);
			}
		}
	}
	boundary({ eta });
	scalar fa[];
	foreach()
		fa[] = (4.*f[] +
		2.*(f[-1, 0] + f[1, 0] + f[0, -1] + f[0, 1]) +
		f[1, 1] + f[-1, 1] + f[1, -1] + f[-1, -1]) / 16.;
	boundary({ fa });
	foreach_face() {
		double fm = (fa[] + fa[-1]) / 2.;
		muv.x[] = (fm*(eta[] + eta[-1]) / 2. + (1. - fm)*mug);
		alphav.x[] = 1. / rho(fm);
	}
	foreach()
		rhov[] = rho(fa[]);
	boundary({ muv, alphav, rhov });
}
/**
convergence outputs
*/
void mg_print(mgstats mg)
{
	if (mg.i > 0 && mg.resa > 0.)
		fprintf(stderr, "#   %d %g %g %g\n", mg.i, mg.resb, mg.resa,
		exp(log(mg.resb / mg.resa) / mg.i));
}
/**
convergence stats
*/
event logfile(i += 1) {
	stats s = statsf(f);
	fprintf(stderr, "%g %d %g %g %g %g\n", t, i, dt, s.sum, s.min, s.max - 1.);
	mg_print(mgp);
	mg_print(mgpf);
	mg_print(mgu);
	fflush(stderr);
}


/**
wall friction
$$\frac{du}{dt} = \frac{-2 \mu_w p u}{W ||u||}$$
equation with discretization:
$$\frac{u^{n+1}-u^n}{\delta t} = \frac{-2 \mu_w p u^{n+1}}{W ||u^n||} $$
*/


event friction(i++) {
	vector uold[], um[];
	int it = 0;
	double errmuw;

	foreach() {
		foreach_dimension(){
			uold.x[] = u.x[];
		}
	}

	do{
		it++;
		foreach() {
			double a = norm(u) == 0 ? HUGE : 1. + 2.*muwall*dt*p[] / (norm(u)*WDOMAIN);
			foreach_dimension(){
				um.x[] = u.x[];
				u.x[] = uold.x[] / a;
			}
		}
		boundary((scalar *){ u });
		errmuw = 0;
		foreach() {
			errmuw += sqrt(sq(u.x[] - um.x[]) + sq(u.y[] - um.y[]));
		}
		errmuw /= pow(2, 2 * LEVEL);
	}
	/**
	convergence if $|u^{m+1}-u^{m}| < 1e-6$ or more 25 iterations
	*/
	while (it < 25 && errmuw>1e-6);
	fprintf(stderr, "errmuw = %g ----  %d   \n", errmuw, it);
}


/**
non-slip boundary condition in the inclination bottom area
*/
event unul(i++) {
	foreach() {
		foreach_dimension(){
			u.x[] = u.x[] * (y >tan(angle / 180.0*pi)*(x)+0.);
		}
	}

}


/**
save some interfaces
*/
event interface (t = 0; t += 1.; t <= tmax) {
#if dimension == 2
	output_facets(f, fpf);
#endif
	char s[80];
	sprintf(s, "field-%g.txt", t);
	FILE * fp = fopen(s, "w");
	output_field({ f, p, u, uf, pf }, fp, linear = true);
	fclose(fp);
}



/**
Rate of flowing materials across the hole
*/
event debit(t += 0.1) {
	static double V = 1;
	V = 0;
	foreach()
		V = V + f[] * Delta * Delta;
	if (t >= 0.) fprintf(stdout, "%lf %lf \n", t, V);
	fflush(stdout);
	fprintf(fqout, "%lf %lf \n", t, V);
}
/**
film output
*/
#if 1
event movie(t += 0.05) {
	static FILE * fp1 = popen("ppm2mpeg > level.mpg", "w");
	scalar l[];
	foreach()
		l[] = level;
	output_ppm(l, fp1, min = 0, max = LEVEL,
		n = 512, box = { { 0, 0 }, { Lz, Lz } });

	foreach()
		l[] = f[] * (1 + sqrt(sq(u.x[]) + sq(u.y[])))*(y >tan(angle / 180.0*pi)*(x)+0.);
	boundary({ l });
	static FILE * fp2 = popen("ppm2mpeg > velo.mpg", "w");
	output_ppm(l, fp2, min = 0, max = 2., linear = true,
		n = 512, box = { { 0, 0 }, { Lz, Lz } });

	static FILE * fp3 = popen("ppm2mpeg > f.mpg", "w");
	foreach()
		l[] = f[] * p[];
	output_ppm(l, fp3, min = 0, linear = true,
		n = 512, box = { { 0, 0 }, { Lz, Lz } });
}
event pictures(t == 3) {
	output_ppm(f, file = "f.png", min = 0, max = 2, spread = 2, n = 512, linear = true,
		box = { { 0, 0 }, { 2, 2 } });
}
#endif


/**
# Run

to run

~~~bash
qcc -g -O2 -Wall -o Lateral LateralSandglassAngle40.c -lm
./Lateral > out
~~~



/**
#Bibliography
L. Staron, P.-Y. Lagrée, & S. Popinet (2014)
"Continuum simulation of the discharge of the granular silo, A validation test for the μ(I) visco-plastic flow law"
Eur. Phys. J. E (2014) 37: 5 DOI 10.1140/epje/i2014-14005-6

L. Staron, P.-Y. Lagrée & S. Popinet (2012)
"The granular silo as a continuum plastic flow: the hour-glass vs the clepsydra"
Phys. Fluids 24, 103301 (2012); doi: 10.1063/1.4757390

Y. Zhou, P. Ruyer, and P. Aussillous, Discharge flow of a bidisperse granular media from a silo: Discrete particle simulations Phys. Rev. E 92,
062204 (2015).

Zhou, Y., Lagree, P. Y., Popinet, S., Ruyer, P., and Aussillous, P. (2017). Experiments on, and ´
discrete and continuum simulations of, the discharge of granular media from silos with a lateral
orifice. Journal of Fluid Mechanics, 829:459–485.
*/

