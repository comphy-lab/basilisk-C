/** ## Axisymetric Planar sheet retraction 2D.
![Hole expansion](axiss2/movie.mp4)
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
\vec{u} = \vec{u}^{*}u_{c} \quad u_{c} = \sqrt{\frac{2\gamma}{\rho H}}, \quad r = r^{*}H, \quad p = p^{*}\frac{2\gamma}{H}, \quad t = t^{*}\tau \quad \tau = \frac{H}{u_{c}} = \sqrt(\frac{\rho H^{3}}{2\gamma})
$$
With H the thickness's liquid sheet, $\tau$  the time scale for low Oh otherwise $\tau_{vis} = \frac{Oh H}{u_{c}}$(High Oh).
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
Oh = \frac{\mu}{\sqrt{2 \gamma \rho H}}
$$
The dimensionless density $\tilde{\rho}$ and viscosity $\tilde{\mu}$ are both 1 in the liquid phase (water) and $\frac{\rho_{g}}{\rho_{L}}$ $\frac{\mu_{g}}{\mu_{L}}$, respectively, in the gas phase (air).
*/

#include "axi.h" // y axis => r and x axis => z. (r,z)
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "view.h"
#include "navier-stokes/perfs.h"

#define MINLEVEL 5
/**
Geometrical parameters
*/
#define L0 60 //domain size
#define H 1. //H is Half thickness in the code
#define R0 H // radius bulge
#define yo 11 // r position center radius
#define Circle sq(x)+sq(y-yo)-sq(R0)
#define line x-H

/**
The ratios are the one for a water sheet which is retracting in the air. */

#define rhoeau 1000.
#define rhoair 1.2
#define rhoLiquid 1000.
#define muRatio 10.

int MAXLEVEL;
double tEnd = 30.;
double Oh;
double TaylorCulick;
double Ohs[3] = {0.1, 0.01, 0.005};
int main(int argc, char *argv[]) { 
	MAXLEVEL = 8;
	size (L0); 
	init_grid (1 << (MAXLEVEL-2));
	for (int i=0; i<=1; i++){
		Oh = Ohs[i];
          #if 0
		if(Oh==0.1){
			for(int lvl=8; i<=9; lvl++){
				MAXLEVEL =  lvl;
				double mu = Oh;
				rho2 = rhoeau/rhoLiquid, mu2 = mu;
				rho1 = rhoair/rhoLiquid, mu1 = mu/muRatio;
				f.sigma = 1.;
				TaylorCulick = sqrt((f.sigma)/(rho2*2*H));
				run();
			}
		}
          #endif
          #if 1
		double mu = Oh;
		rho2 = rhoeau/rhoLiquid, mu2 = mu;
		rho1 = rhoair/rhoLiquid, mu1 = mu/muRatio;
		f.sigma = 1.;
		TaylorCulick = sqrt((f.sigma)/(rho2*2*H));
		run();
          #endif
	}
}

/**
## Boundary conditions
*/
/**
## Initial conditions
*/
event init (t = 0) {
	fraction(f, y<=yo ? Circle : line);
        boundary ({f,u});
}
/**
## Mesh adaptation
We use an adaptive mesh ([Adaptation Algorithm](http://basilisk.fr/sandbox/Antoonvh/the_adaptive_wavelet_algorithm)) with a maximum level of refinement of MAXLEVEL. We adapt the mesh with respect to the interface and the velocity (u.r ; u.z).*/

event adapt (i++) {
	double uemax = 0.05;
	adapt_wavelet ({f,u}, (double[]){0.005,uemax,uemax}, MAXLEVEL, MINLEVEL);
}
/**
## Outputs
*/

/**
We output the interface of the fluid to track the evolution of the
liquid bulge. We calculate the minimun position to interpolate the speed of retraction at the edge.*/
double Yo;
event extractPosition (i++) {
	vector h[];
	heights (f, h);
	double yMin = +HUGE;
	foreach(reduction(min:yMin)){
		if (h.y[] != nodata) {
			double yi = y + height(h.y[])*Delta;
			if (yi < yMin){
				yMin = yi;
			}
		}
	}
	if(t==0){
		Yo = yMin;			
	}
	char name[80];
	sprintf(name, "dataaxi%g",Oh);
	static FILE * fout = fopen (name, "w");
	fprintf(fout,"%g %g %g %g %g\n", t, t*sqrt((2*f.sigma)/(rho2*pow(2*H,3))),(yMin-Yo)/(2*H), interpolate(u.y,0,yMin)/(TaylorCulick*sqrt(2)), TaylorCulick*sqrt(2));
	fflush(fout);
	if(Oh==0.1){
		char nam[80];
		sprintf(nam, "logs%d",MAXLEVEL);
		static FILE * flog = fopen (nam, "w");
	  	fprintf (flog, "%d %f %d %d %d %g %ld\n",i,
		   t*sqrt((2*f.sigma)/(rho2*pow(2*H,3))), N, mgp.i, mgu.i, interpolate(u.y,0,yMin)/(TaylorCulick*sqrt(2)), grid->n);
                fflush(flog);
	}
}

/**
We output the evolution of the fraction field with the vorticity using [view.h](http://basilisk.dalembert.upmc.fr/src/view.h) and the [draw.h command](http://basilisk.dalembert.upmc.fr/src/draw.h).
[How to find view arguments](https://groups.google.com/forum/#!searchin/basilisk-fr/view$20parameters|sort:date/basilisk-fr/Z8goUFJPivA/c5N9ALcpBAAJ)*/
event interface (i++) {
	clear();
	view (fov = 19.298, quat = {0,0,0,1}, tx = -0.364702, ty = -0.504662, bg = {1,1,1}, width = 860, height = 760, samples = 4);
	box();
	draw_vof("f",lw = 2);
 	scalar omega[];
  	vorticity(u, omega);
	squares("omega");
	begin_translate (x = L0/6.);
	 	cells();
	end_translate ();
	char time[80];
	sprintf(time, "t=%.2gtau",t*sqrt((2*f.sigma)/(rho2*pow(2*H,3))));
	draw_string(time,pos=0);
	char ratio[80];
	sprintf(ratio, "Oh=%g, R0/H=%g", Oh,(yo-R0)/2*H);
	draw_string(ratio,pos=3);
	char video[80];
	sprintf(video, "vof%g.mp4",Oh);
	save(video);

	scalar l[];
	char gridfile[80];
	sprintf(gridfile, "grid%g.mp4",Oh);
	foreach()
		l[] = level;
	output_ppm (l, file = gridfile, n=512, min=MINLEVEL, max=MAXLEVEL, box = {{X0,Y0},{L0/6,L0}});
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
We run a gnuplot script to plot the evolution of the axisymetric hole in a liquid sheet. It's works with gnuplot 5.2 version
*/
#if 0
event movie(t=tEnd){
    char gnu[200];
    sprintf(gnu, "gnuplot -e 'Oh=%g; pas=%g; tend=%g' 3d3.gnu | ffmpeg -f png_pipe -s:v 1024x768 -i pipe: -y movie.mp4",Oh,step,tEnd*sqrt((2*f.sigma)/(rho2*pow(2*H,3))));
    //sprintf(gnu, "gnuplot -e 'Oh=%g; pas=%g; tend=%g' 3d4.gnu | ffmpeg -f png_pipe -s:v 1024x768 -i pipe: -y movie.mp4",Oh,step,tEnd*sqrt((2*f.sigma)/(rho2*pow(2*H,3)))); save pictures on the hard drive
	system(gnu);
}
#endif


/**
### Moderate Oh

![fraction field Oh=0.1.](axiss2/vof0.1.mp4)
![refinement Oh=0.1.](axiss2/grid0.1.mp4)

![fraction field Oh=0.005](axiss2/vof0.005.mp4)
![refinement Oh=0.005](axiss2/grid0.005.mp4)

### fit TaylorCulick's speed
~~~gnuplot speed
set tics font ",12"
set xlabel "t/tau" font ",12"
set ylabel "u/U_{TC}" font ",12"
set title font ",12"
set style line 1 linewidth 2
set key font ",12"
set yrange [0:1.2]
set title "Evolution of the rim's velocity"
plot 'dataaxi0.1' u 2:4 w l lw 2 title "Oh=0.1",\
      'dataaxi0.01' u 2:4 w l  lw 2 title "Oh=0.01",\
      'dataaxi0.005' u 2:4 w l lw 2 title "Oh=0.005"
~~~
### fit N.Savva et al 2009
[NIKOS SAVVA and JOHN W. M. BUSH (2009)](https://doi.org/10.1017/S0022112009005795). Viscous sheet
retraction.
~~~gnuplot Savva
reset
set key bottom right
set xrange [0:15.]
set yrange [0:1.]
set xlabel "times t/tau"
set ylabel 'Uy/U_{TC}'
set title "0h_{1}=0.2, 0h_{2}=0.04"
plot 'dataaxi0.2' every 6 u 2:4 t 'Oh_{1}=0.2',\
 'savvaA.csv' u 1:2 w l lw 2 t 'Savva_{1}',\
 'dataaxi0.04' every 6 u 2:4 t 'Oh_{2}=0.04',\
 'savvaB.csv' u 1:2 w l lw 2 t 'Savva_{2}'
~~~
### Convergence test Oh=0.1
~~~gnuplot convergence
reset
set tics font ",12"
set xlabel "t/tau" font ",12"
set ylabel "u/U_{TC}" font ",12"
set title font ",12"
set style line 1 linewidth 2
set key font ",12"
plot for [i=8:10:1] 'logs'.i using 2:6 w l lw 2t sprintf("N=2^{%i}",i)
~~~
*/
