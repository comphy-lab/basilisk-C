/**
# Taylor-Culick retraction of a liquid drop
This file is largely inspired by the work of Alexis Berny. Thanks to him !*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"


#define LEVEL 9
#define INLEVEL 5
#define MINLEVEL 5

/**
We are simulating only one half of the liquid drop*/
#define L0 5.
#define L 2.


/**
Parameters */

#define R 0.5
//Ohnesorge number
#define Oh 1.
#define muRatio 1000.
#define rhoRatio 100.

double tC = R*R*Oh;
double tEnd = 30*R*R*Oh;
//Taylor-Culick velocity (viscosity / density inside the drop = 1)
double utc = 1./(1.*R*Oh);


/**
The boundary conditions are slip walls*/


int main() {

  /**
  The domain */
  size (L0);
  init_grid (1 << INLEVEL);

  //mass conservation tolerance
  TOLERANCE = 1e-5;
  NITERMAX = 200;


  /**
  We define the parameters of the fluids. The internal fluid (the liquid) is 
  the fluid 1. The external fluid is the fluid 2.*/
  double rhoInt = 1;
  double muInt = 1.;

  rho1 = rhoInt, mu1 = muInt;
  rho2 = rhoInt/rhoRatio, mu2 = muInt/muRatio;
  f.sigma = (mu1 /Oh)*(mu1 /Oh)*1./rho1/R;
    
  run();
}


/**
We define the geometry with boolean operations. We define a
rectangular area (the intersection of Line and Line_vert). We also
define a half circle.  Then, our geometry will be the union of the
rectangle with the half circle. */

double geometry (double x, double y) {
  double Line = y-R;
  double Line_vert = x-L;
  double Circle = sq(x-L)+sq(y)-sq(R);
  double right_part = max(-Circle,-Line_vert);
  return min(-Line, right_part);
}

/**
We do not use AMR, because previous tests have shown that it can generate
spurious oscillations of the tip velocity. */

event init (t = 0) {
    if (!restore (file = "dump-xxxx")){

       /**
       We initialise the geometry of the interface.*/
       fraction(f, geometry(x,y));

       //Width of the region with maximal refinement
       double Wm=3*R;
       //Width between two adjacent layers
       double Wa=0.5*R;

       //Main region 
       refine (x < (L + Wm - R ) && y < (Wm) && level < LEVEL);

       for(int k=0; k < (LEVEL-MINLEVEL); k +=1){

         //East region
         refine (x > (L + Wm - R + k*Wa) && x < (L + Wm - R + (k+1)*Wa) && y > 0. && y < (Wm) && level < LEVEL-k-1);
         //Top region
         refine (x < (L + Wm - R + (k+1)*Wa) && y > (Wm) && y < (Wm+(k+1)*Wa) && level < LEVEL-k-1);

       }
    }

  scalar le[];
  foreach(){
    le[] = level;
  }
  char name_grid[100];
  sprintf (name_grid, "grid.png");
  FILE* fgrid = fopen (name_grid, "w");
  output_ppm (le, file = name_grid,min=MINLEVEL,max=LEVEL,n = 1000,box = {{0.,0.},{L0,L0}});
  fclose (fgrid);

}


/**
We output : the interface position and scalar / vector fields . */

event saveDatas (t += tC; t <= tEnd) {

  /**
  We only output the interface in this function. */
  
  char name[80];
  sprintf (name, "interface-%f.txt", t);
  FILE* fp = fopen (name, "w");
  output_facets (f, fp);
  fclose (fp);

  scalar vort[];
  vorticity (u, vort);
  char namev[80];
  sprintf (namev, "vort-%f.txt", t);
  FILE * fpv = fopen (namev, "w");

  /**
  Txt format. */
  output_field ({vort}, fpv, n=1000);
  fclose (fpv);

  char named[80];
  sprintf (named, "dump-%d", i);
  dump (file = named);

}

/**
We output the maximum x position of the liquid finger. Indeed this
position should linearly evolve in time, with a few variations coming
from the capillary wave. */

double xPrev = -1, tPrev = -1;

event extractPosition (i ++) {

  /**
  We define a new vector field, h. */
  
  vector h[];
  
  /**
  We reconstruct the height function field and take the corresponding 
  maximum along the x axis. */
  
  heights (f, h);
  double xMax = -HUGE;;
  foreach()
    if (h.x[] != nodata) {
      double xi = x + height(h.x[])*Delta;
      if (xi > xMax)
  xMax = xi;
  }

  /**
  We also output the velocity of the end of the bulge.*/

  double veloTip = tPrev >= 0 ? (xMax - xPrev)/(t - tPrev) : 0.;
  char name[80];
  sprintf (name, "velocity_pos.dat");
  static FILE * fp = fopen (name, "w");
  fprintf (fp," %g %g %g\n", t, xMax, veloTip);
  fflush (fp);
  tPrev = t, xPrev = xMax;

}

/**
We output, in the standard output file, the step with the corresponding
time. */

event logfile(i++) {
  printf ("i = %d t = %g\n", i,t);
  fflush(stdout);
}


/* Images for post processing */

event movies (t += tC*1e-1; t <= tEnd){
  scalar vort[];
  vorticity (u, vort);
  scalar m[];
  foreach() {
    if (f[]<0.95 && f[] >0.05)
      m[] = -10;
    else
      m[] = 0.;
  }
  //we output both the interface and the normalized vorticity.
  char name_vort[100];
  sprintf (name_vort, "vort-%g.png", t);
  FILE* fvort = fopen (name_vort, "w");
  double om = (mu1 /Oh)*1./rho1/R/R;
  output_ppm (vort, fvort, mask =m,min=-om,max=om, map=cool_warm, n = 1000,box = {{0.,0.},{L0,3*R}});
  fclose (fvort);


}

event end (t = tEnd) {}

![Tip velocity as function of time](data/velocity_oh1.png)