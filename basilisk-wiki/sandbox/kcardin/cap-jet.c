/**
This simulation is based on ``Bursting.c`` from the aberny sanbox.

The Ohnesorge number and cavity aspect ratio arguments are provided from the terminal: ``./bursting "10000." "5.0" 2> >(tee 10000_5p0_11)``*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "tag.h"
#include "view.h"

/**
Definition of the geometrical parameters. The 2 radii correspond
respectively to the radius of the cavity and to the radius of the corner between the cavity wall and free surface. */

#define LEVEL 11
#define R1 1.
#define R2 0.1

double La;
double H;
int count;

#define Sigma 1.

uf.n[left] = 0.;
uf.n[right] = 0.;
uf.n[top] = 0.;
uf.n[bottom] = 0.;

double uemax = 0.3;

int main(int argc, char *argv[]) {
  size (L0);
  origin (-L0/2., 0.);
  init_grid (1 << (10));

  
  if (argc >= 2){
      La = atof (argv[1]);
    }
  
  /**
  The Ohnesorge number is defined such that: $Oh = \frac{1}{\sqrt{La}}$ */
  
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
  We set the height of the channel*/
  H = atof (argv[2]);
  /**
  We print the characteristics of the fluid and some other info in the logfile. */
  fprintf (stderr, "props %f %f %f %f %f\n", rho1, mu1, La, R1, Sigma);
  
  
  fprintf (stderr, "Aspect Ratio %f \n", H);
  fprintf (stderr, "La %f\n", La);
  fprintf (stderr, "Level %d\n", LEVEL);
  
  /**
  We set the droplet count to 0*/
  count = 0;
  run();
}

/**
  We define the initial geometry of the channel*/
double geometry(double x, double y) {
  double F1 = - x - 1e-8;
  double F2 = y-(R1+R2);
  double F1F2 = min(F1, F2);
  double F3 = -x-R2;
  double F1F2F3 = max(F1F2, F3);
  double F4 = y-R1;
  double F1F2F3F4  = min(F1F2F3, F4);
  double F5 = -x-H*R1;
  double F1F2F3F4F5 = max(F1F2F3F4, F5);
  double C1 = sq(x + H*R1) + sq(y) - sq(R1);
  double F1F2F3F4F5C1 = min(F1F2F3F4F5,C1 );
  double C2 = sq(x + R2) + sq(y - R1 -R2) - sq(R2);
  double F1F2F3F4F5C1C2 = min(-F1F2F3F4F5C1, C2);
  return -F1F2F3F4F5C1C2;
} 



/**
We initialise the geometry of the interface.*/
  
double intemax = 0.005;

event init (t = 0) {
  double iteration = 0;
  do {
    iteration++;
    fraction(f, geometry(x,y));
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


event movie (t += 0.01; t <= 3) {
  char str[99];
  char str2[99];
  view (psi = -pi/2.);
  /*cells();*/
  draw_vof ("f");
  draw_vof ("f", filled = 1, fc = {0.2, 0.4, 0.9});
  sprintf (str, "%0.2f %g %g %d", t, H, La, LEVEL);
  sprintf (str2, "%g_%g_%d.mp4", H,La,LEVEL);
  draw_string (str, size = 30);
  mirror({0,-1}){
    draw_vof ("f");
    draw_vof ("f", filled = 1, fc = {0.2, 0.4, 0.9});}
  /*cells();*/
  save (str2);
  /*save ("movie.mp4");*/
}


event droplets (t += 0.01)
{
  /**
  We determine the maximum on axis velocity*/
  double xMax = -HUGE;
  scalar bound[];
  foreach()
  	bound[]=0.;
  foreach_boundary (bottom){
        bound[] = 1.; 
  	if ((u.x[]*f[]) > xMax)
        	xMax = u.x[]*f[];
        }
  fprintf (stderr, "%g %g\n",t, xMax );
  
  /**
  We tag connected regions*/
  scalar m[];
  foreach()
    m[] = f[] > 1e-3;
  int n = tag (m);

  /**
  Once each cell is tagged with a unique droplet index, we can easily
  compute the volume *v* and position *b* of each droplet.*/

  double v[n];
  coord b[n];
  double dropVelo[n]; /*Sum of color function weighted velocity in tagged region*/

  
  for (int j = 0; j < n; j++)
    v[j] = b[j].x = b[j].y = b[j].z =dropVelo[j]= 0. ;
  foreach_leaf()
    if (m[] > 0) {
      int j = m[] - 1;
      v[j] += dv()*f[];
      dropVelo[j] += dv()*f[]*u.x[];
      coord p = {x,y,z};
      foreach_dimension()
	b[j].x += dv()*f[]*p.x; /*Center of mass of tagged region*/
    }



  /**
  We output the volume and position of each droplet to
  standard output. */
  if (n>1 && count<=4) {
  double veloDrop;
  count +=1;
  for (int j = 0; j < n; j++){
    veloDrop = dropVelo[j]/v[j]; /*Color function weighted average velocity in tagged region*/
    fprintf (stderr, "drop %g %d %g %d %g %g %g %g\n", H, i, t,
	     j, v[j], b[j].x/v[j], b[j].y/v[j],veloDrop);
    }     
  /**
  We stop the simulation once the 4th droplet is formed*/
  }
  else if (count>4){
    return 1;
  }
}


