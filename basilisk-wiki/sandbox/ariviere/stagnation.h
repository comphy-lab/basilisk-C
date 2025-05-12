/** #Bubble in a stagnation point flow.

We first create the axi-symmetric stagnation point flow. This is done by imposing the analytical boundary conditions (inflow on the top and outflow on the right). After a transient, the flow converges to the analytical solution. 
The stationary flow is then used to study the deformations and breakup of an initially spherical bubble which lies at the stagnation point. 

The velocity field of a stagnation point flow is (using Basilisk notation, y is the radial direction and x the axial direction.):
$$  \bm{u} = -Ex \bm{e_x} +\frac{1}{2}E y \bm{e_y}$$ 

In what follows, the flow enters radially and outs along the axis of symmetry, so that $E <0$. 
 */



//#include "grid/quadtree.h"
#include "axi.h"
#include "navier-stokes/centered.h"
#include "view.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "navier-stokes/perfs.h"
#include "maxruntime.h"

//Physical ratios
#define RHOR 850.0 //density ratio
#define MUR 25.0 //kinematic viscosity ratio; real value 55; to be consistent with turbubble simulations

#define HALFBOX 5.

double amp=-1., R0=1., SIGMA=1.;//velocity amplitude; bubble radius; surface tension
int MINLEVEL=5, MAXLEVEL=8;
double MAXTIME = 5.;

/** Boundary conditions are chosen so that they match the extract solution for the flow at a stagnation point.
 */

//boundary conditions
u.n[left] = dirichlet(0.);
u.n[right] = neumann(-amp);
u.n[top]    = dirichlet(1./2.*amp*2.*HALFBOX);
u.n[bottom] = dirichlet(0.);

p[top] = neumann(-1./2.*sq(amp)*HALFBOX);
pf[top] = neumann(-1./2.*sq(amp)*HALFBOX);
p[bottom] = neumann(0.);
pf[bottom] = neumann(0.);

p[left] = neumann(0.);
pf[left] = neumann(0.);
p[right] = dirichlet(-1./8.*sq(amp*y));
pf[right] = dirichlet(-1./8.*sq(amp*y));


int main(int argc, char * argv[]){
    if (argc>1){
	amp = atof(argv[1]);//typical shear
	R0 = atof(argv[2]);//in practice, R0 is kept constant at 1. It is the length reference
	SIGMA = atof(argv[3]);
	MINLEVEL = atoi(argv[4]);
	MAXLEVEL = atoi(argv[5]);
	MAXTIME = atof(argv[6]);
    }
    size(2.*HALFBOX);
    init_grid(1<<MINLEVEL);

    //Physical properties
    rho1 = 1.0;//liquid is the reference
    rho2 = rho1/RHOR;

    mu1 = 0.01;
    mu2 = mu1/MUR;
    f.sigma = SIGMA;

    run();
}

event init(i = 0){
    /**
     The same simulation can be used to create the stagnation point flow and to
     inject a bubble in a previously computed flow. To inject a bubble in a
     flow, be sure to include in the folder a dump file called 'restart'. 
     */
    if (!restore("restart") && !restore("dump")){
	foreach(){
	    f[] = 1.0;
	}
	boundary({f});
    }
    else if (restore(file="restart")){//bubble injection
	fraction(f, sq(x) + sq(y) - sq(R0));
    }
    else{}
    p.nodump=false;//just added
}
/** 
 We adapt on both the velocity field and the position of the interface.*/
event adapt(i++){//grid adaptation
    double uemax = 1.e-3;
    double femax = 1.e-2;
    adapt_wavelet((scalar *){f, u}, (double[]){femax, uemax, uemax}, MAXLEVEL, MINLEVEL);
}


/**
   To compute the spherical harmonics decomposition, we extract the coordinates
 of the interface. This is very low level and should be improved in the future.*/

event getInterface(i++){
    char fixname[100];//x
    char fiyname[100];//y
    char fikname[100];//curvature
    sprintf(fixname, "fix_%d.dat", pid() );
    sprintf(fiyname, "fiy_%d.dat", pid() );
    sprintf(fikname, "fik_%d.dat", pid() );
    FILE * fix = fopen(fixname, "a");
    FILE * fiy = fopen(fiyname, "a");
    FILE * fik = fopen(fikname, "a");
    fprintf (fix, "%g ", t);
    fprintf (fiy, "%g ", t);
    fprintf (fik, "%g ", t);

    scalar curve[]; //curvature
    curvature(f, curve);
    boundary((scalar*){curve});

    //position with the heigh function
    scalar xpos[];
    scalar ypos[];

    position(f, xpos, {1,0});
    position(f, ypos, {0,1});
  
    foreach(){
      if (interfacial(point,f)){
	fprintf (fix, "%g ", xpos[]);
        fprintf (fiy, "%g ", ypos[]);
	fprintf (fik, "%g ", curve[]);
      }
    }
    fprintf (fix, "\n");
    fprintf (fiy, "\n");
    fprintf (fik, "\n");
    fclose(fix);
    fclose(fiy);
    fclose(fik);
    return 0;
}

/** We generate a movie of the simulation. The bubble will be in white, the
 * surrounding flow in black.*/

event movie(t+=0.05){
  view(fov=20.,tx=-0.5);
  clear();
  draw_vof("f", filled = 1);
  box(notics=true);
  mirror({0, 1}){
    draw_vof("f", filled = 1);
    box(notics=true);
  }
  save("movie.mp4");
}

/** We output a dump file at the end of the simulation. */
event end(t = MAXTIME){
    printf("%g", t);
   dump("end");
} 
