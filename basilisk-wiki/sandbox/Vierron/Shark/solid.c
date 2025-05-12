#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "utils.h"
#include "distance.h"
#include "fractions.h"
#include "view.h"
#include "navier-stokes/perfs.h"
#include "lambda2.h"

/**
## Importing the geometry 
*/
double maxl = -HUGE;

void fraction_from_stl (scalar f, FILE * fp, double eps, int maxlevel)
{

  /**
  We read the STL file and compute the bounding box of the model. */
  
  coord * p = input_stl (fp);
  coord min, max;
  bounding_box (p, &min, &max);
  foreach_dimension()
  //double length;

    if (max.x - min.x > maxl){
      maxl = max.x - min.x;
    }
  
  /**
  We initialize the distance field on the coarse initial mesh and
  refine it adaptively until the threshold error (on distance) is
  reached. */

  scalar d[];
  distance (d, p);
  while (adapt_wavelet ({d}, (double[]){eps*L0}, maxlevel, 5).nf);

  /**
  We also compute the volume fraction from the distance field. We
  first construct a vertex field interpolated from the centered field
  and then call the appropriate VOF functions. */

  vertex scalar phi[];
  foreach_vertex()
    phi[] = (d[] + d[-1] + d[0,-1] + d[-1,-1] +
	     d[0,0,-1] + d[-1,0,-1] + d[0,-1,-1] + d[-1,-1,-1])/8.;
  boundary ({phi});
  fractions (phi, f);
}
/**
 Parameters
*/
int LEVEL = 9;
double Re = 1.e5;
double speed = 1.;
double tEnd=20.;
scalar Shark[];
face vector muv[];

/**
# Main function
*/
int main(){
  size (20);
  origin (- L0/2,- L0/2,- L0/2);
  init_grid (32);
  mu = muv;
  //TOLERANCE = 1e-5;

  Shark.refine = Shark.prolongation = fraction_refine;
  run();
}

/**
We set a constant viscosity corresponding to a Reynolds number of 300, based on the Shark's length and the inflow velocity (1).
*/
event properties (i++){
  foreach_face()
    muv.x[] = fm.x[]*maxl*speed/Re;
}

/**
# Boundary conditions
*/
u.n[top] = dirichlet(-speed);
u.n[bottom] = neumann(0.);

p[top] = neumann(0);
p[bottom] = dirichlet(0);
pf[top] = neumann(0);
pf[bottom] = dirichlet(0);
Shark[top] = 0.;
/**
# Initial conditions
*/
event init(t=0){
  if (!restore (file = "restart")) {
    FILE * fp = fopen ("shark.stl", "r");
    fraction_from_stl (Shark, fp, 5e-4, LEVEL);
    fclose (fp);
  }
   dump();
	foreach(){
		u.y[]=-speed;
	}
	boundary ({Shark,u.y});	
}

event velocity (i++){
# if 1
  foreach()
    foreach_dimension()
      u.x[] = (speed - Shark[])*u.x[];
#endif
}

/**
# Mesh adaptation
*/
event adapt (i++) {
  adapt_wavelet ({u,Shark}, (double[]){5e-2,5e-2,5e-2,1e-3}, LEVEL,5);
}

/**
# Outputs
*/
event movie(t+= 0.01; t<=tEnd){
        view (fov = 15.1858, quat = {-0.00958367,0.65255,0.757675,0.00382529}, tx = 0, ty = 0, bg = {0.3,0.4,0.6}, width = 722, height = 598, samples = 1);
        travelling(0.1,tEnd/2,quat = {-0.632565,0.191283,0.219854,-0.717591},tx = 0, ty = 0, fov = 15.1858);
        travelling(tEnd/2+0.1,tEnd-0.1,quat = {-0.211687,-0.592443,-0.725576,-0.278822},tx = 0, ty = 0, fov = 15.1858);
	clear();
	box();
	draw_vof ("Shark", filled=-1, fc = {0.2,0.2,0.8});
  	translate (z = -L0/2){
	  cells();
	}
	scalar l2[];
	lambda2 (u, l2);
	isosurface ("l2", -1.);
	save ("test.mp4");
}
/**
# Movie
![Shark.](solid/test.mp4)
([Shark movie](http://basilisk.dalembert.upmc.fr/sandbox/Vierron/Shark/solid/test.mp4))
*/
