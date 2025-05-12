/** ## Rayleigh Benard instability with adaptative mesh */
#include "convection_boussinesq_adapt.h"
#include "navier-stokes/perfs.h"
#include "global_nusselt.h"
#include "profil5c.h"

#define MINLEVEL 5
#define MAXLEVEL 8

double EndTime= 200.;
double ytop=1.;

int main() {
  size(4.);
  Ra = 1e4; 
  Pr = 0.72;
  init_grid(1<<(MAXLEVEL));
  run(); 
}

/** 
## Boundaries Conditions 
*/

T[left] = neumann(0.);
T[right] = neumann(0.);

T[embed] = dirichlet(-0.5);
T[bottom] = dirichlet(0.5);

u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);

u.n[right] = dirichlet(0.);
u.n[left] = dirichlet(0.);
u.t[right] = dirichlet(0.);
u.t[left] = dirichlet(0.);

p[left] = neumann(0.);
pf[left] = neumann(0.);

p[right] = neumann(0.);
pf[right] = neumann(0.);


/**
## Initial conditions
Initial conditions correspond to a linear temperature
profile and no motion.
*/
event init (t=0) {
        vertex scalar phi[];
        foreach_vertex() {
           phi[] = ytop-y;
        }
        boundary ({phi});
        fractions (phi, cs, fs);
 	foreach(){
		T[] = -y;
		foreach_dimension()
	 		u.x[] = 0.;
 	}
        boundary ({T,u});
}

/** 
## Adaptative mesh 
*/
event adapt (i++) {
	adapt_wavelet ({T,u}, (double[]){1e-2,5e-2,5e-2}, MAXLEVEL, MINLEVEL);
}
#if 1
/** 
## Output Data 
*/
event output (t += 1.0; t <= EndTime) {
  #if 1
	double nu_vol=0., nu_t=0., nu_b=0.;
	nu_vol = nusselt_vol(T,u);
	nu_t = nusselt_top(T);
	nu_b = nusselt_bot(T);
	static FILE * fout = fopen ("data", "w");
	if (t<1.){
		fprintf (fout, "# [1]t [2]Ra [3]Pr [4]N [5]umin [6]umax [7]vmin [8]vmax [9]nutop [10]nubot [11]nuvol \n");
	}
	stats velx = statsf (u.x), vely = statsf (u.y);
	fprintf (fout, " %g %.9g %.9g %d %.9g %.9g %.9g %.9g %.9g %.9g %.9g \n",
	t, Ra,Pr,N, velx.min, velx.max, vely.min, vely.max, nu_t, nu_b, nu_vol);
	fflush (fout);
  #endif
/** 
## Output video Temperature field + Grid 
*/
	output_ppm (T, file="temperature.mp4", n = 1024, box = {{0.,0},{L0,ytop}});
	scalar l[];
	foreach(){
		l[] = level;
	}
	output_ppm (l, file="grid.mp4", n=512, min=MINLEVEL, max=MAXLEVEL, box = {{0.,0},{L0,ytop}});
}

/** 
## Calcul profil de temperature + velocity
*/
event pro (t=EndTime) {
	profile({T,u},"profils",0, ytop);
}
#endif
/**

## Results

~~~gnuplot
set output "Nu=f(t).png"
set xlabel 'times t'
set ylabel 'Nusselts'
titre=system("awk 'NR==2 {print $2}' data")
set title sprintf("Ra=%s",titre)
plot 'rbi2dadapt/data' u 1:9 w l title 'Nu-top',\
      '' u 1:10 w l title 'Nu-bot'
~~~

![Temperature field.](rbi2dadapt/temperature.mp4)
![Animation of the level of refinement.](rbi2dadapt/grid.mp4)
*/