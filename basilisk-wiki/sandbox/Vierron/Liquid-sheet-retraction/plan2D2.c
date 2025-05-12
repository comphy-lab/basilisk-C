/** ## Planar sheet retraction 2D.
We simulate the retraction of a water planar sheet under the capillary force effect surrounded of air.
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
\vec{u} = \vec{u}^{*}u_{c} \quad u_{c} = \sqrt{\frac{2\gamma}{\rho H}}, \quad X = X^{*}H, \quad p = p^{*}\frac{2\gamma}{H}, \quad t = t^{*}\tau \quad \tau = \frac{H}{u_{c}} = \sqrt(\frac{\rho H^{3}}{2\gamma})
$$
With H the liquid sheet's thickness, $\tau$  the time scale for low Oh otherwise $\tau_{vis} = \frac{Oh H}{u_{c}} =  \sqrt(\frac{\rho H^{3}}{2\gamma})$(High Oh).
\newline
Dimensionless equations :
$$
\tilde{\nabla} \tilde{ . \vec{u}} = 0
$$
$$
\tilde{\rho}(\tilde{\partial_{t}}\tilde{\vec{u}} + \tilde{\vec{u}}.\tilde{\nabla} \tilde{\vec{u}}) = -\tilde{\nabla} \tilde{p} + Oh \tilde{\mu} \tilde{\nabla^{2}}\tilde{\vec{u}} + \tilde{\kappa} \delta_{s}n
$$
With the Ohnesorge number Oh, which compares viscous and capillary effects
$$
Oh = \frac{\mu}{\sqrt{2\gamma\rho H}}
$$
The dimensionless density $\tilde{\rho}$ and viscosity $\tilde{\mu}$ are both 1 in the liquid phase (water) and $\frac{\rho_{g}}{\rho_{L}}$ $\frac{\mu_{g}}{\mu_{L}}$, respectively, in the gas phase (air).
*/

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "view.h"
#include "navier-stokes/perfs.h"


#define MINLEVEL 5
/**
Geometrical parameters
*/
#define L0 80 //Domain size
#define H 1. //H is Half thickness (we simulate half thickness)
#define R0 H // radius circle
#define xo 21 // position center circle
#define Circle sq(x-xo)+sq(y)-sq(R0) 
#define line y-H

/**
The ratios are the one for a water sheet which is retracting in the air. */

#define rhoeau 1000.
#define rhoair 1.2
#define rhoLiquid 1000.
#define muRatio 10.

int MAXLEVEL;
double tEnd = 50;
double Oh;
double TaylorCulick;
double Ohs[3] = {0.1, 0.01, 0.005};
int main(int argc, char *argv[]) { 
	MAXLEVEL = 8;
	size (L0); 
	init_grid (1 << (MAXLEVEL));
	for (int i=0; i<=2; i++){
		Oh = Ohs[i];
          #if 0
		if(Oh==0.1){
			for(int lvl=8; i<=10; lvl++){
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
	fraction(f, x<=xo ? Circle : line);
	boundary ({f,u});
}
/**
## Mesh adaptation
We use an adaptive mesh ([Adaptation Algorithm](http://basilisk.fr/sandbox/Antoonvh/the_adaptive_wavelet_algorithm)) with a maximum level of refinement of MAXLEVEL. We adapt the mesh with respect to the interface and the velocity (u.x ; u.y).*/

event adapt (i++) {
	double uemax = 0.05;
	adapt_wavelet ({f,u}, (double[]){0.005,uemax,uemax}, MAXLEVEL, MINLEVEL);
}
/**
## Outputs
*/
/**
We output the interface of the fluid to track the evolution of the
liquid bulge. We calculate the minimun position to interpolate the speed of retraction. Then the evolution of the maximum thickness. */
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
  	double yMax = -HUGE;
	foreach(reduction(max:yMax)){
		if (h.y[] != nodata) {
			double yi = y + height(h.y[])*Delta;
			if (yi > yMax){
				yMax = yi;
			}
		}
	}
	char name[80];
	sprintf(name, "data%g",Oh);
	static FILE * fout = fopen (name, "w");
	fprintf(fout,"%g %g %g %g %g %g\n", t,(xMin-Xo)/(2*H),t*sqrt((2*f.sigma)/(rho2*pow(2*H,3))), yMax/H, interpolate(u.x,xMin,0)/(TaylorCulick*sqrt(2)),TaylorCulick*sqrt(2));
	fflush(fout);
	if(Oh==0.1){
		char nam[80];
		sprintf(nam, "logs%d",MAXLEVEL);
		static FILE * flog = fopen (nam, "w");
	  	fprintf (flog, "%d %g %d %d %d %g %ld\n",i,
		   t*sqrt((2*f.sigma)/(rho2*pow(2*H,3))), N, mgp.i, mgu.i, interpolate(u.x,xMin,0)/(TaylorCulick*sqrt(2)), grid->n);
                fflush(flog);
	}

}

/**
We output the evolution of the fraction field with the vorticity using [view.h](http://basilisk.dalembert.upmc.fr/src/view.h) and the [draw.h command](http://basilisk.dalembert.upmc.fr/src/draw.h).
[How to find view arguments](https://groups.google.com/forum/#!searchin/basilisk-fr/view$20parameters|sort:date/basilisk-fr/Z8goUFJPivA/c5N9ALcpBAAJ)*/
event interface (i++) {
	clear();
	view (fov = 12.5982, quat = {0,0,0,1}, tx = -0.498782, ty = -0.221298, bg = {1,1,1}, width = 1080, height = 630, samples = 4);
	box();
	draw_vof("f",lw = 2);
 	scalar omega[];
  	vorticity(u, omega);
	squares("omega");
	begin_translate (y = L0/4.);
	 	cells();
	end_translate ();
	char time[80];
	sprintf(time, "t=%.2gtau",t*sqrt((2*f.sigma)/(rho2*pow(2*H,3))));
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
	output_ppm (l, file = gridfile, n=512, min=MINLEVEL, max=MAXLEVEL, box = {{X0,Y0},{L0,L0/4}});
}

/**
We output the interface of the fluid to track the evolution of the
liquid bulge. */

event plotInterface (t+= 0.5; t <= tEnd) {
  	char names[80];
	sprintf(names, "interface%d", pid());
	FILE * fp = fopen (names, "w");
	output_facets (f,fp);
	fclose(fp);
	char command[80];
	sprintf(command, "LC_ALL=C  cat interfa* > all2D%g%g.dat",t,Oh);
        system(command);// allow to use linux command in the c code to concatenate our files
}
/**
## Results

![fraction field Oh=0.1](plan2D2/vof0.1.mp4)
![refinement Oh=0.1](plan2D2/grid0.1.mp4)

### low Oh
![fraction field Oh=0.01](plan2D2/vof0.01.mp4)
![refinement Oh=0.01](plan2D2/grid0.01.mp4)

![fraction field Oh=0.005](plan2D2/vof0.005.mp4)
![refinement Oh=0.005](plan2D2/grid0.005.mp4)

### fit TaylorCulick's speed
~~~gnuplot speed
set tics font ",12"
set xlabel "t/tau" font ",12"
set ylabel "u/U_{TC}" font ",12"
set title font ",12"
set style line 1 linewidth 2
set key font ",12"
set key bottom right
set yrange [0:1]
plot 'data0.1' u 3:5 t "Oh=0.1",\
 'data0.01' u 3:5 t "Oh=0.01",\
 'data0.005' u 3:5 t "Oh=0.005"
~~~
### fit G.AGBAGLAH et al 2009.
We plot the evolution of the maximun thickness of the bulge and superpose with the results of [G.AGBAGLAH et al 2009](https://doi.org/10.1088/0169-5983/41/6/065001)
~~~gnuplot bourrelet
reset
set title ' Oh_{1}=0.1, Oh_{2}=0.01, Oh_{3}=0.005'
set ylabel '2yMax/2H'
set xlabel 't/tau'
set yrange [0:5.5]
set xrange [0:16]
set key left top
set key font ",10"
plot 'data0.1' u 3:4 every 4 ps 0.5 t 'Oh=0.1',\
'data0.01' u 3:4 every 4 ps 0.5  t 'Oh=0.01',\
'data0.005' u 3:4 every 4 ps 0.5  t 'Oh=0.005',\
'010.csv' u 1:($2/5) w l lw 2 t 'AGBAGLAH_{1}',\
'001.csv' u 1:($2/5) w l lw 2 t 'AGBAGLAH_{2}',\
'0005.csv' u 1:($2/5) w l lw 2 t 'AGBAGLAH_{3}'
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
set key bottom right
set title 'Oh=0.1, R0/H=5'
plot for [i=8:10:1] 'logs'.i using 2:6 w l t sprintf("N=2^{%i}",i)
~~~
*/
