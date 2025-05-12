/** External flow driven coalescence with constant strain rate */
/**
![Flow driven coalescence](flow_driven/movie.mp4)
*/
#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
//#include "reduced.h"
#include "tag.h"
#include "view.h"
#include "navier-stokes/perfs.h"


double stop_sim = 20.;
int LEVEL = 7; //maxlevel
int initlevel = 7;
double stension = 0.01;
#define R 1. //Radius
#define RHOR 1.1 //suspending/drop 
#define MUR 5.26 //suspending/drop
#define GRA 0. //gravity
#define HX 2.*R //Initial distance from droplets center to the
		//collision plane
#define HY 0.
#define LENGTH_DOMAIN 7.*R
#define UEMAX 1e-4
#define STEPS1 (int)(0.2*stop_sim*pow(3, (LEVEL - 5))) 

double OH = 0.02; //Ohnesorge number = mud/sqrt(rhod*sigma*R)

//Uinf = non-dimensional far field velocity =
		     //G*sqrt(rhod*R^3/sigma), which is equal to G, as
		     //sqrt(rhod*R^3/sigma) is 1(characteristic time
		     //scale)
double Uinf = 0.095; 


u.t[top] = dirichlet(-Uinf*x);
u.n[top] = dirichlet(0.5*Uinf*y);
//p[top] = neumann(0);
//pf[top] = neumnn(0);

u.n[left] = dirichlet(-Uinf*x);
u.t[left] = dirichlet(0.5*Uinf*y);
//p[left] = neumann(0);
//pf[left] = neumann(0);

int main(int argc, char * argv[])
{ 
  if(argc > 1)
    LEVEL = atoi(argv[1]);
  if(argc > 2)
    initlevel = atoi(argv[2]);
  if(argc > 3)
    Uinf = atof(argv[3]);
  if(argc > 4)
    OH = atof(argv[4]);
  if(argc > 5)
    stension = atof(argv[5]);
  if(argc > 6)
    stop_sim = atof(argv[6]);

  size(LENGTH_DOMAIN);
  origin(-L0,0);
  init_grid(1 << initlevel);
  
  //G.y = -GRA;
  rho1 = stension, mu1 = OH*stension;
  rho2 = rho1*RHOR, mu2 = mu1*MUR, f.sigma = stension;
  // static face vector velf[];
  TOLERANCE = 1e-4;// chage this for finer solution 
  
  run();
}

// init function

event init(i = 0)
{
  if (!restore (file = "dump")) {
    fraction (f, - (sq(x+HX) + sq(y) - sq(R)));
  }
  foreach()
    {
      if(!f[])
	{
	  u.x[] = -Uinf*x;
	  u.y[] = 0.5*Uinf*y;
	}
    }
  boundary({u.x, u.y});
}

// adaptive function

event adapt(i++)
{
  double uemax = UEMAX;
  adapt_wavelet ({f,u}, (double[]){1e-2,uemax,uemax}, LEVEL);
}
 
// minimum and center film thickness

event filmthickness(i++) {
  double min = R, cent = -R;
  scalar pos[];
  position (f, pos, {1,0});
  foreach_face(y, reduction(max:cent))
  {
    if(f[] && f[]!=1 && !y && (cent < pos[]))
      cent = pos[];
  }
  min = -statsf(pos).max;
  char name[80];
  sprintf (name, "height");
  static FILE * fp = fopen (name, "w");
  fprintf (fp, "%g %g %g\n", t, min, -cent);
  fflush (fp);
}

event movies(i+=STEPS1)
{
  view(fov = 7, tx = -0.003073, ty = 0.00850717);
  draw_vof("f", lw = 3);
  mirror({0,-1,0})
    draw_vof("f", lw = 3);
  mirror({1,0,0}) {
    draw_vof("f", lw = 3);
    mirror({0,-1,0})
      draw_vof("f", lw = 3);
  }
  cells();
  save("movie.mp4");
}

event stop(t=stop_sim)
{
 dump("dump");
}                 
                                                                                                                                                                 
/**
   Results: 
   Film thickness is naturally scaled to radius of the droplet and time is naturally scaled to the characteristic time:$$\sqrt{\frac{\rho R^3}{\sigma}}$$,  which is equal to 1. 
   ~~~gnuplot Minimum & center film thickness as a function of time 
   set grid
   set logscale y
   set xlabel 't'
   set ylabel 'z'
   set xrange[0:15]
   plot 'height' u ($1<10.5?$1:1/0):($2>0.056?$2:1/0) w l title "minimum thickness", "../data.txt" u 1:2 w l title "Sambath et al.", 'height' u ($1<10.5?$1:1/0):($3>0.056?$3:1/0) w l title "center thickness"
   ~~~

*/
