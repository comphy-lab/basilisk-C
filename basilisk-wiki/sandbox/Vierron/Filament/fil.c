/** ## Axisymmetric liquid filament retraction 2D.
*/
/**
## Model equations
$$
\frac{\partial u_r}{\partial t}  + u_r \frac{\partial u_r}{\partial r} + u_z \frac{\partial u_r}{\partial z} =  - \frac{1}{\rho} \frac{\partial p}{\partial r} + \frac{\mu}{\rho} \left\{-\frac{u_r}{r^2} + \frac{1}{r} \frac{\partial}{\partial r} \left( r \frac{\partial u_r}{\partial r} \right) + \frac{\partial^2 u_r}{\partial z^2} \right\} + \gamma \kappa \delta_{s}n\\
\frac{\partial u_z}{\partial t}  + u_r \frac{\partial u_z}{\partial r} + u_z \frac{\partial u_z}{\partial z} =  - \frac{1}{\rho} \frac{\partial p}{\partial z} + \frac{\mu}{\rho} \left\{ \frac{1}{r} \frac{\partial}{\partial r} \left( r \frac{\partial u_z}{\partial r} \right) + \frac{\partial^2 u_z}{\partial z^2} \right\} + \gamma \kappa \delta_{s}n
$$
$$
\frac{1}{r}\frac{\partial}{\partial r} (r u_r) +\frac{\partial u_z}{\partial z} = 0.
$$
Non-dimensionalizing the equations by :
$$
\vec{u} = \vec{u}^{*}u_{c} \quad u_{c} = \sqrt{\frac{2\gamma}{\rho R}}, \quad r = r^{*}R, \quad p = p^{*}\frac{2\gamma}{R}, \quad t = t^{*}\tau \quad \tau = \frac{R}{u_{c}} = \sqrt(\frac{\rho R^{3}}{2\gamma})
$$
With R the radius's filament , $\tau$  the time scale for low Oh otherwise $\tau_{vis} = \frac{Oh R}{u_{c}}$(High Oh).
\newline
Dimensionless equations :
$$
\tilde{\nabla} \tilde{ . \vec{u}} = 0
$$
$$
\tilde{\rho}(\tilde{\partial_{t}}\tilde{\vec{u}} + \tilde{\vec{u}}.\tilde{\nabla} \tilde{\vec{u}}) = -\tilde{\nabla} \tilde{p} + Oh \tilde{\mu} \tilde{\nabla^{2}}\tilde{\vec{u}} + \frac{R}{2}\tilde{\kappa} \delta_{s}n
$$
With the Ohnesorge number oh, which compares viscous and capillary effects
$$
Oh = \frac{\mu}{\sqrt{2\gamma \rho R}}
$$
The dimensionless density $\tilde{\rho}$ and viscosity $\tilde{\mu}$ are both 1 in the liquid phase (water) and $\frac{\rho_{g}}{\rho_{L}}$ $\frac{\mu_{g}}{\mu_{L}}$, respectively, in the gas phase (air).
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "view.h"
#include "navier-stokes/perfs.h"

#define MINLEVEL 5
/**
Geometrical parameter
*/
#define L0 100
#define R 1. //H is Half thickness
#define R0 R
#define xo 31
#define Circle sq(x-xo)+sq(y)-sq(R0)
#define line y-R

/**
The ratios are the one for a water filament which is retracting in the air. */

#define rhoeau 1000.
#define rhoair 1.2
#define rhoLiquid 1000.
#define muRatio 10.

int MAXLEVEL;
double tEnd = 30.;
double Oh;
double TaylorCulick;
double Ohs[2] = {0.5, 0.1};
int main(int argc, char *argv[]) { 
	MAXLEVEL = 8;
	size (L0); 
	init_grid (1 << (MAXLEVEL-1));
	for (int i=0; i<=1; i++){
		Oh = Ohs[i];
		double mu = Oh;
		rho2 = rhoeau/rhoLiquid, mu2 = mu;
		rho1 = rhoair/rhoLiquid, mu1 = mu/muRatio;
		f.sigma = 1/2.;
		TaylorCulick = sqrt((f.sigma)/(rho2*R));
		run();
	}
}

/**
## Boundary conditions
*/
/**
## Initial conditions
*/
event init (t = 0) {
	fraction(f, x<=xo ? Circle : line);
	boundary ({f,u});
}
/**
## Mesh adaptation
We use an adaptive mesh with a maximum level of refinement of 10. We adapt 
the mesh with respect to the interface the velocity.*/

event adapt (i++) {
	double uemax = 0.05;
	adapt_wavelet ({f,u}, (double[]){0.005,uemax,uemax}, MAXLEVEL, MINLEVEL);
}
/**
## Outputs
*/

/**
We output the interface of the fluid to track the evolution of the
liquid bulge. */
double Xo;
event extractPosition (i++) {
	vector h[];
	heights (f, h);
	double xMin = +HUGE;
	foreach(reduction(min:xMin)){
		if (h.x[] != nodata) {
			double xi = x + height(h.x[])*Delta;
			if (xi < xMin){
				xMin = xi;
			}
		}
	}
	if(t==0){
		Xo = xMin;			
	}
	char name[80];
	sprintf(name, "datafil%g",Oh);
	static FILE * fout = fopen (name, "w");
	fprintf(fout,"%g %g %g %g %g\n", t, t*sqrt((2*f.sigma)/(rho2*pow(R,3))),(xMin-Xo)/(R), interpolate(u.x,0,xMin)/(TaylorCulick*sqrt(2)), TaylorCulick*sqrt(2));
	fflush(fout);
}


/**
We output the fraction field to have an idea of the evolution of the instability. */
event interface (i++) {
	clear();
	view (fov = 7.93684, quat = {0,0,0,1}, tx = -0.593915, ty = -0.109611, bg = {1,1,1}, width = 1028, height = 500, samples = 1);
	box();
	draw_vof("f",lw = 2);
 	scalar omega[];
	scalar omegabis[];
	vorticity(u, omega);
	double y1 = L0/10.;
	foreach() {
	    if (y > y1)
	      omegabis[] = nodata;
	    else
	      omegabis[] = omega[];
 	}
	squares("omegabis");
	double y2=L0/10.;
	scalar ubis[];
	foreach() {
	    if (y > y2)
	      ubis[] = nodata;
	    else
	      ubis[] = u.x[]/(TaylorCulick*sqrt(2));
 	}
	begin_translate (y = L0/10.);
		squares("ubis", min=0, max=1.);
		draw_vof("f",lw = 2);
	end_translate ();
	begin_translate (y = L0/5.);
		cells();
		draw_vof("f",lw = 2);
	end_translate ();
	char time[80];
	sprintf(time, "t=%.2gtau",t*sqrt((2*f.sigma)/(rho2*pow(R,3))));
	draw_string(time,pos=0);
	char ratio[80];
	sprintf(ratio, "Oh=%g", Oh);
	draw_string(ratio,pos=3);
	char video[80];
	sprintf(video, "vof%g.mp4",Oh);
	save(video);

	scalar l[];
	char gridfile[80];
	sprintf(gridfile, "grid%g.mp4",Oh);
	foreach()
		l[] = level;
	output_ppm (l, file = gridfile, n=512, min=MINLEVEL, max=MAXLEVEL, box = {{X0,Y0},{L0,L0}});
}
/**
We output the interface of the fluid to track the evolution of the
liquid bulge. */

event plotInterface (t += 0.1; t <= tEnd) {
  	char names[80];
	sprintf(names, "interface%d", pid());
	FILE * fp = fopen (names, "w");
	output_facets (f,fp);
	fclose(fp);
	char command[80];
	sprintf(command, "LC_ALL=C  cat interfa* > allaxi%g%g.dat",t,Oh);
        system(command);
}

/**
### Moderate Oh
![fraction field Oh=0.5.](fil/vof0.5.mp4)
![refinement Oh=0.5.](fil/grid0.5.mp4)
![fraction field Oh=0.1.](fil/vof0.1.mp4)
![refinement Oh=0.1.](fil/grid0.1.mp4)

### fit TaylorCulick's speed
~~~gnuplot speed
set xlabel 'times t/tau'
set ylabel 'Ux edge/U_{TC}'
set title "Evolution of the rim's velocity"
plot 'datafil0.5' u 2:4 w l title "Oh=0.5",\
      'datafil0.1' u 2:4 w l title "Oh=0.1"
~~~
*/
