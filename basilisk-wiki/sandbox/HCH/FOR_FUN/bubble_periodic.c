/**
A test case for bubble rising in the periodic BC
*/
#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h"

#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
#include "two-phase.h"

#include "tension.h"
#if dimension == 3
# include "lambda2.h"
#endif
#include "view.h"

#include "maxruntime.h"

#define RHOR 10.
#define MUR 10.

# define Ga 29.9
# define Bo 2.

# define MAXTIME 20
#define WIDTH 1.61290
#define R0 0.5

int LEVEL = 4;

int main (int argc, char * argv[]) {
  maxruntime (&argc, argv);
  if (argc > 1)
    LEVEL = atoi (argv[1]);
  
  size (WIDTH);
  origin (-L0/2.,-L0/2.,-L0/2.);
  init_grid (64);
  periodic(right);
  periodic(front);
  periodic(top);
  
  rho1 = 1.;
  rho2 = 1./RHOR;
  mu1 = 1./Ga;
  mu2 = 1./(MUR*Ga);
  f.sigma = 1./Bo;
  
  TOLERANCE = 1e-4;
  run();
}


event init (t = 0) {
  if (!restore (file = "restart")) {
  fraction (f, sq(x) + sq(y) + sq(z) - sq(R0));
  }
}

/**
The acceleration is set to ensure the conservation of momentum flux across the boundary [Buner & Tryggvason](https://doi.org/10.1017/S0022112002001179). Since $\mathbf{g}$ is unity in this case the acceleration reads 
$$
(1 - \rho_0/\rho) = (1 - ((1-\alpha)\rho_1 + \alpha\rho_2)/\rho)
$$
Moreover, according to [two-phase-generic.h](/src/two-phase-generic.h#properties), the $rho$ on each phase is altered according to vertically adjacent cells.
*/
event acceleration (i++) {
face vector av = a;
foreach_face(y){
double ff1 = (f[0,0,0]+f[0,-1,0])/2.;

double rho_m = fm.y[]/rho(ff1);

av.y[] -= (1. - (0.875212*rho1+0.124788*rho2)*rho_m);
}
}

event movie (t += 1e-1){

  clear();
  draw_vof ("f");
  save ("movie.mp4");
}

event logfile (i = 0 ; t <= MAXTIME; i += 10) {
  double sb = 0.;
  double vby = 0., v_y = 0.;
  double Re = 0.;
  foreach(
	  reduction(+:vby) reduction(+:v_y)
	  reduction(+:sb)) {
    double dv = (1. - f[])*dv();
    vby += u.y[]*dv;
    v_y += u.y[]*dv();
    sb += dv;
  }
  
  Re = (vby/sb - v_y/(pow(WIDTH,3)))*Ga;
  
  fprintf (stderr,
	   "%.8f %.8f %.8f %.8f %.8f\n", 
	   t, sb,
	   vby/sb, v_y/(pow(WIDTH,3)), Re);
  fflush (stderr);
}

event dumplog(t = 0; t <= MAXTIME; t++)
{
  char name[80];
  sprintf(name, "snapshot/dump-%03d", (int) t);
  dump(file = name);
}