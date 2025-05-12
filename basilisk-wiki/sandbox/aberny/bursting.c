/**
# Bursting bubble

With this code, we want to simulate a symplified version of a bursting
bubble in a water bath. We based this simulation on the work of
Duchemin et al. published in 2002. */

#include "axi.h"
#include "navier-stokes/centered.h"
// #define mu(f) ((mu1*mu2)/(clamp(f,0,1)*(mu2-mu1)+mu1))
#include "two-phase.h"
#include "tension.h"
#include "tag.h"

/**
Definition of the geometrical parameters. The 2 radii correspond
respectively to the radius of the bubble and to the radius of a smooth
interface. */

#define LEVEL 14
#define L0 10
#define R1 1.

/**
The parameters L1, L2 and L3 are used to define the hole between the upper 
interface and the spherical droplet. Once the hole is done, we smooth the 
interface by adding a circle between the upper interface and the drop. This 
circle as to be continuous with the upper interface but also with the bubble 
interface. $2\times L2$ defines the distance between the top of the bubble and 
the interface. L1 defines the radius of the smoothing circle, but also the 
distance between the center of this cylinder with the upper interface.*/

#define L2 0.05
#define L3 0.15
#define L1 (-(sq(R1+2.*L2)-sq(R1)+sq(L3))/(4.*(L2+R1)))

#define Sigma 1.

uf.n[left] = 0.;
uf.n[right] = 0.;
uf.n[top] = 0.;
uf.n[bottom] = 0.;

double uemax = 0.3;

int main(int argc, char *argv[]) {
  size (L0);
  origin (-L0/2., 0.);
  init_grid (1 << (9));

  /**
  By default, the value of the Laplace number is 10 000. But, we can
  change that with an initial input argument. */
  double La = 10000.;
  if (argc >= 2){
      La = atof (argv[1]);
    }

  /**
  The Ohnesorge number is defined such that: $Oh = \frac{1}{\sqrt{La}} */
  
  double Oh = sqrt(1./La);

  /**
  We then define the properties of the fluids. In this case, the fluid
  1 correspond to water, and the fluid 2 correspond to air.*/

  rho1 = 1., mu1 = Oh;
  rho2 = 1./998., mu2 = Oh/55;
  
  /**
  The surface tension is defined with the Ohnesorge number*/

  f.sigma = Sigma;

  /**
  We print the characteristics of the fluid in the logfile. */
  fprintf (stderr, "props %f %f %f %f %f\n", rho1, mu1, La, R1, Sigma);
  
  run();
}

double geometry(double x, double y) {

  /**
  The geometry need the definition of 2 circles (C*) and 3 lines
  (D*). Then we need to proceed correctly to the intersection and the
  union of those parameters.*/
  
  double C1 = sq(x + R1 + 2*L2) + sq(y) - sq(R1);
  double C2 = sq(x - L1) + sq(y - L3) - sq(L1);
  double D1 = - x - 1e-8;
  double D2 = y - L3;
  double D3 = - x + (2*L1);

  /**
  Our geometry is:
  $$geom = C1 \cup (D3 \cap (D1 \cup D2))$$ */

  double D1D2 = min(D1, D2);
  double D1D2D3 = max(D1D2, D3);
  double D1D2D3C1 = min(D1D2D3, C1);
  double D1D2D3C1C2 = min(-D1D2D3C1, C2);
  return -D1D2D3C1C2;
} 

/**
We initialise the geometry of the interface.

![Initial interface of the simulation. The length of the rim is half of the 
size of the hole](bursting/initial.png)
*/
  
double intemax = 0.005;

event init (t = 0) {
  double iteration = 0;
  do {
    iteration++;
    fraction(f, geometry(x,y))
  } while( adapt_wavelet({f,u}, (double []){intemax,uemax,uemax,uemax},
     maxlevel = LEVEL, 5).nf != 0 && iteration <= 10);

  static FILE * fp = fopen("initial.png", "w");
  output_ppm(f, fp, 400, min = 0, max = 1, box = {{-2.5,0},{0.5,3}});
  static FILE * fp2 = fopen("initial.gfs", "w");
  output_gfs(fp2);
}


/**
We use an adaptive mesh. Our refinement criteria are based on the
interface and the velocity. */

event adapt (i++) {
  adapt_wavelet ({f,u}, (double []){intemax,uemax,uemax,uemax},
		 maxlevel = LEVEL, 5);
}

/**
We output the shape of the interface at regular intervals. */


/**
This gives the following image

~~~gnuplot Evolution of the interface
unset key
set size ratio -1
plot [-2:2]'interface' u 2:1 w l, 'interface' u (-$2):1 w l lt 1
~~~

We output the position of the interface of the fluid allong the axis of
symmetry at every timestep. We track this position with the height
function. We also compute the corresponding velocity. */

static double tprev = -1, xprev = -1,xDropPrev = -1;
static double compt = 0, toto = 0;

event droplets (i ++; t<=2) {
  scalar m[];
  foreach()
    m[] = f[] > 1.5e-1;
  int n = tag(m);

  /**
  We want to observe the formation of the first drop. We detect when the tag 
  function return more than 1 phase of the internal fluid. We then output 15 
  steps in gfs file to observe the behaviour of the fluid at this moment.*/

  if (n>1 && t>0.15) {
    compt++;
    if (compt>=1 && compt<=15) {
      char rupture[80];
      sprintf(rupture, "rupture-%07ld.gfs",i);
      FILE* fprup = fopen (rupture, "w");
      output_gfs(fprup);
    }
  }  

  /**
  We tag the droplet, when there is some.*/

  double v[n];
  coord b[n];
  double dropVelo[n];

  for (int j = 0; j < n; j++)
    dropVelo[j] = v[j] = b[j].x = b[j].y = b[j].z = 0.;

  foreach_leaf()
    if (m[] > 0) {
      int j = m[] - 1;
      v[j] += dv()*f[];
      dropVelo[j] += dv()*f[]*u.x[];
      coord p = {x,y,z};
      foreach_dimension()
      b[j].x += dv()*f[]*p.x;
    }

  /**
  We compute the droplets velocity. We need at least 2 steps with a droplet to 
  have a correct velocity. We add for that an extra variable, which is 
  counting the number of steps with the presence of a drop. We also compute 
  the volume of the first drop.*/

  double xDrop, veloDrop, volume, radius;

  for (int j = 1; j < n; j++) {

    /**
    Computation of the drop velocity, and the drop position.*/

    xDrop = b[j].x/v[j];
    veloDrop = dropVelo[j]/v[j];

    /**
    Computation of the drop radius.*/ 

    volume = v[j];
    radius = sqrt(2*volume/M_PI);

    fprintf (fout, "%d %g %d %g %g %g %g %g %d\n", i, t,
       j, v[j], xDrop, veloDrop, volume, radius, n);
  }
  fflush (fout); 

/**
We also compute the velocity of the top of the jet. We use for that the height 
function.*/

  vector h[];
  double xMax = -HUGE;;
  heights (f, h);
  foreach()
    if (y <= Delta && h.x[] != nodata && m[]==0) {
      double xi = x + height(h.x[])*Delta;
      if (xi > xMax)
        xMax = xi;
    }
  double velo = tprev >= 0 ? (xMax - xprev)/(t - tprev) : 0.;
  double deltaT = t - tprev;
  printf ("%d %g %g %g %g\n", i, t, deltaT, xMax, velo);
  fflush (stdout);

  if (xMax > 0. && xprev < 0.) {
    fprintf (stderr, "end %g %g\n", t, velo);
    char name[80];
    sprintf(name, "jetX-0.gfs");
    FILE* fp = fopen (name, "w");
    output_gfs (fp);
  }  

  /**
  The stop condition of the simulation is based on 2 conditions. First, we 
  need to overpass the axi x=0. Secondly, we need to observe at least 2 steps 
  consecutive steps with a drops.*/

  if (xMax>=0. && n>1) {
    fprintf(stderr, "firstDrop %g %g %g %g\n", t, veloDrop, volume, radius);
    return 1;
  }

  tprev = t, xprev = xMax, xDropPrev = xDrop;
}

/**
This gives the following figure

~~~gnuplot Evolution of the velocity of the centerline
set xlabel 'Time'
set ylabel 'Velocity'
set grid
plot [0.01:]'out' u 1:4 w l t ''
~~~
*/

event interface(i+= 10) {
  static FILE * fp = popen ("ppm2mpeg > tracer.mpeg", "w");
  static FILE * fp2 = popen ("ppm2mpeg > grid.mpeg", "w");
  output_ppm(f, fp, 512, min=0, max=1, box = {{-2.5,0},{0.5,3}});
  scalar l[];
  foreach()
    l[] = level;
  output_ppm(l,fp2, min=5, max=LEVEL, box = {{-2.5,0},{0.5,3}});
}

/**
We can follow the evolution of the simulation with a video

![The bursting bubble](bursting/tracer.mpeg)
*/

#if 0
event gfsview (i += 10) {
  static FILE * fp = popen ("gfsview2D bursting.gfv", "w");
  output_gfs (fp);
}
#endif

event logfile (i += 10) {
  fprintf (stderr, "%d %g %d %d %d\n", i, t, mgp.i, mgpf.i, mgu.i);
}

/**
We output the final state of the simulation. In the final state, we include 
the tag function, to see where basilisk identify the different phase.*/

event finalState (t = end) {
  scalar m[];
  foreach()
    m[] = f[] > 1.5e-1;
  int n = tag(m);
  output_gfs (file = "final.gfs");
  static FILE * fp = fopen("final.ppm", "w");
  output_ppm(f, fp, 512, min=0, max=1, box = {{-2.5,0},{0.5,3}});
  dump (file = "dump");
}

/**
We can control the refinement parameter by loocking at the evolution 
of the mesh during the simulation.

![Evolution of the mesh](bursting/grid.gif) */
