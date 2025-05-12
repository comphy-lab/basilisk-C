/** ## Hole expansion 3D.
*/
/**
## Model equations
$$
\rho(\partial_{t}\vec{u} + \vec{u}.\nabla \vec{u}) = -\nabla p + \mu \nabla^{2}\vec{u} + \gamma \kappa \delta_{s}n
$$
$$
\nabla . \vec{u} = 0
$$
Non-dimensionalizing the equations by :
$$
\vec{u} = \vec{u}^{*}u_{c} \quad u_{c} = \sqrt{\frac{2\gamma}{\rho H}}, \quad \nabla = \nabla^{*}H, \quad p = p^{*}\frac{2\gamma}{H}, \quad t = t^{*}\tau \quad \tau = \frac{H}{u_{c}} = \sqrt{\frac{\rho H^{3}}{2\gamma}}
$$
With 2H the thickness's liquid sheet, $\tau$  the time scale for low Oh otherwise $\tau_{vis} = \frac{Oh H}{u_{c}}$(High Oh).
\newline
Dimensionless equations :
$$
\tilde{\nabla} \tilde{ . \vec{u}} = 0
$$
$$
\tilde{\rho}(\tilde{\partial_{t}}\tilde{\vec{u}} + \tilde{\vec{u}}.\tilde{\nabla} \tilde{\vec{u}}) = -\tilde{\nabla} \tilde{p} + Oh \tilde{\mu} \tilde{\nabla^{2}}\tilde{\vec{u}} + \tilde{\kappa} \delta_{s}n
$$
With the Ohnesorge number oh, which compares viscous and capillary effects
$$
Oh = \frac{\mu}{\sqrt{2\gamma \rho H}}
$$
The dimensionless density $\tilde{\rho}$ and viscosity $\tilde{\mu}$ are both 1 in the liquid phase (water) and $\frac{\rho_{g}}{\rho_{L}}$ $\frac{\mu_{g}}{\mu_{L}}$, respectively, in the gas phase (air).
*/

#include "grid/octree.h"// for 3D dimensions
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "view.h"
#include "navier-stokes/perfs.h"
#include "output_crown.h"// to get the evolution of the rim (x,y) at z=0

#define MINLEVEL 5
//#define MAXLEVEL 8
/**
Geometrical parameter 
*/
#define L0 150.// Domain size
//#define H 1. // H is Half thickness of the liquid layer
#define R0 H // radius's Sphere
#define line z-H // Liquid Layer's equation

#define corrugation 1 // set to 1 to get an initial perturbation on the fraction field
#define smooth 0 // set to 1 to get a smooth initial circular shape

#if corrugation
	#define R 10. // radius of the circle
	#define A 4. // Amplitude of the perturbation
	#define frq 2. // frequency
	#define circle (sq(x)+sq(y)-sq(R + A*fabs(sin(atan2(y,x)*frq)))) // Corrugate circle
	#define Co (R + A*fabs(sin(atan2(y,x)*frq))) // Corrugate radius
	#define Sphere sq(x-Co*cos(atan2(y,x)))+sq(y-Co*sin(atan2(y,x)))+sq(z)-sq(R0) // Sphere used to smooth the edge of the rim
#endif

#if smooth
	#define Co 11 // origin's Sphere
	#define Sphere (sq(x-Co)+sq(y-Co))+sq(z)-sq(R0) // Sphere used to smooth the edge of the rim
#endif

/**
The ratios are the one for a water sheet which is retracting in the air. */

#define rhoeau 1000.
#define rhoair 1.2
#define rhoLiquid 1000.
#define muRatio 10.

double tEnd = 50.;
double Oh =0.1;
double H =1.;
int MAXLEVEL = 8;
double TaylorCulick;

int main(int argc, const char* argv[]) { 
  if (argc > 1)
    Oh = atof (argv[1]);
  if (argc > 2)
    H = atof (argv[2]);
  if (argc > 3)
    MAXLEVEL = atoi (argv[3]);
	size (L0); 
	X0=Y0=-L0/2.;
	init_grid (1 << 8);
	double mu = Oh;
	rho2 = rhoeau/rhoLiquid, mu2 = mu;
	rho1 = rhoair/rhoLiquid, mu1 = mu/muRatio;
	f.sigma = 1.;
	TaylorCulick = sqrt((f.sigma)/(rho2*2*H));
	run();
}

/**
## Boundary conditions
*/
/**
## Initial conditions
*/
event init (t = 0) {
	int iteration = 0;
	do 
	{
#if corrugation
		fraction(f, (circle <=0.) ? Sphere : line);
#endif
#if smooth
		fraction(f, sqrt(sq(x)+sq(y))<=Co ? Sphere : line);
#endif
		iteration++;
	 }while (adapt_wavelet ({f}, (double[]){1e-4}, MAXLEVEL).nf !=0 && 
	iteration<=10);
	fprintf(stdout,"%d\n", iteration);
	fflush(stdout);
	boundary ({f,u});
	dump();
}
/**
## Mesh adaptation
We use an adaptive mesh with a maximum level of refinement of 10. We adapt 
the mesh with respect to the interface the velocity.*/

event adapt (i++) {
	double uemax = 0.05;
	adapt_wavelet ({f,u}, (double[]){0.005,uemax,uemax,uemax}, MAXLEVEL, MINLEVEL);
}
event tol (i++; t<=0.5){
	TOLERANCE=1e-4;
}
/**
## Outputs
*/

#if 1
event dumping (t+=2.){
  char snap[80];
  sprintf(snap,"snap%g",t*sqrt((2*f.sigma)/(rho2*pow(2*H,3))));
  dump(snap);
}
#endif

event extractPosition (i++) {
	vector h[];
	heights (f, h);
	double yMin = +HUGE;
	foreach(reduction(min:yMin)){
		if (h.y[] != nodata) {
			double yi = y + height(h.y[])*Delta;
			if (fabs(yi) < yMin){
				yMin = yi;
			}
		}
	}
	double xMin = +HUGE;
	foreach(reduction(min:xMin)){
		if (h.x[] != nodata) {
			double xi = x + height(h.x[])*Delta;
			if (fabs(xi) < xMin){
				xMin = xi;
			}
		}
	}
	char name[80];
	sprintf(name, "data%g",Oh);
	static FILE * fout = fopen (name, "w");
	fprintf(fout,"%g %g %g\n", t, t*sqrt((2*f.sigma)/(rho2*pow(2*H,3))), sqrt(sq(xMin) + sq(yMin)));
	fflush(fout);
}
/**
We output the fraction field to have an idea of the evolution of the instability.*/
event interface (i++) {
	clear();
	view (fov = 17.2793, quat = {-0.55487,-0.0688854,-0.124054,-0.819737}, tx = 0.0262771, ty = -0.0825254, bg = {0.3,0.4,0.6}, width = 820, height = 700, samples = 4);
	draw_vof("f");
	char time[80];
	sprintf(time, "t=%.2gtau",t*sqrt((2*f.sigma)/(rho2*pow(2*H,3))));
	draw_string(time,pos=0);
#if corrugation
	char ratio[80];
	sprintf(ratio, "Oh=%g", Oh);
	draw_string(ratio,pos=3);
#endif
#if smooth
	char ratio[80];
	sprintf(ratio, "Oh=%g, R0/H=%g", Oh,(Co-R0)/H);
	draw_string(ratio,pos=3);
#endif
	char video[80];
	sprintf(video, "Hole%g.mp4",Oh);
	save(video);

	scalar l[];
	char gridfile[80];
	sprintf(gridfile, "grid%g.mp4",Oh);
	foreach()
		l[] = level;
	output_ppm (l, file = gridfile, n=512, min=MINLEVEL, max=MAXLEVEL);
}

/**
We output the interface of the fluid to track the evolution of the
liquid bulge. */

event plotInterface (t += 0.1; t <= tEnd) {
  	char names[80];
	sprintf(names, "file%d", pid());
	FILE * fp = fopen (names, "w");
	output_facets (f,fp);
	fclose(fp);
	char command1[80];
	sprintf(command1, "LC_ALL=C  cat fil* > interface%g%g",t*sqrt((2*f.sigma)/(rho2*pow(2*H,3))),Oh);
        system(command1);

	char crown[80];
	sprintf(crown, "profil%d", pid());
	FILE * fcr = fopen (crown, "w");
	output_crown(f,fcr);
	fclose(fcr);
	char command2[80];
	sprintf(command2, "LC_ALL=C  cat profil* > crown%g%g",t*sqrt((2*f.sigma)/(rho2*pow(2*H,3))),Oh);
        system(command2);
}
/** python command to get the evolution of the rim and lambda_max(t) */
event python(t=tEnd){
	char command3[80];
#if corrugation
	sprintf(command3, "python RP.py\t%g\t%g",Oh,R);
#endif
#if smooth
	sprintf(command3, "python RP.py\t%g\t%g",Oh,Co-R0);
#endif
        system(command3);
}
/**
## Results
![fraction field Oh=0.1.](Hole/Hole0.1.mp4)
![refinement Oh=0.1.](Hole/grid0.1.mp4)

![fraction field Oh=0.01.](Hole/Hole0.01.mp4)
![refinement Oh=0.01.](Hole/grid0.01.mp4)
