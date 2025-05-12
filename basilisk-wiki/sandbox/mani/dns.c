/** This file must do a DNS of coalescence of two drops approaching
    each other*/
// Some things are always fixed: Boundaries (symmetric)

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "reduced.h"
#include "tag.h"
#include "view.h"
#include "navier-stokes/perfs.h"

// vary only BO, AR and rho1
#define STOP 15.
#define D 1. //Diameter
#define RHOR 0.9 //drop/suspending 
#define MUR 1.2 //drop/suspending
#define GRA 1.
#define HX 3.5*D
#define HY D
#define LEVEL 7
#define DELTA ((HX + 0.5*D)/(1<<LEVEL))

int maxlevel = 7;
double uemax = 1e-6;
double AR = 40.; //Archimedes number = delta(rho)*rhos*g*d^3/mus^2
double BO = 1.; //Bond number = delta(rho)*g*dÂ²/sigma
face vector velf[];

int main(int argc, char * argv[])
{
  if(argc > 1)
    AR = atof (argv[1]);
  if(argc > 2)
    BO = atof (argv[2]);
  if(argc > 3)
    rho1 = atof (argv[3]);
  
  size(HX + 0.5*D);
  origin(0,-L0);
  init_grid( 1 << LEVEL);
  
  G.y = -GRA;
  rho1 = 9., mu1 = mu2*MUR;
  rho2 = rho1/RHOR, mu2 = sqrt(rho2*(rho2 - rho1)/AR), f.sigma = (rho2 - rho1)/BO;
  
  TOLERANCE = 1e-4;
  
  run();
}

// init function

event init(i = 0)
{
  if (!restore (file = "restart")) {
    fraction (f, max (- (sq(x) + sq(y) - sq(0.5*D)),
		      - (sq(x) + sq(y + (D+HY)) - sq(0.5*D))));
  }
}

// adaptive function
event adapt(i++)
{
  adapt_wavelet ({f,u}, (double[]){1e-4,uemax,uemax}, maxlevel);
}

/** 

## Minimum distance between drops

This event tracks the minimum distance between two droplets and
velocity of the moving droplet. */

double bn1, bn2, ext = 1.;
double hmin = 2*D;
event thickness_and_velocity(i++)
{
  scalar b[];

  if(hmin > 2*DELTA)
    {
      foreach()
	b[] = f[]>1e-6;
      tag(b);

      if(!i)
	foreach()
	  {
	    if(y <= -(D+HY) && f[] && b[])//change this when you restart from a dump file
	      bn2 = b[];
	    else if(y >= -0.25*D && f[] && b[])	
	      bn1 = b[];// top most bubble
	  }
    }
  else
    {
      foreach()
	b[] = (f[]==1);

      tag(b);

      foreach()
	if(f[]>1e-6 && f[]<1. - 1e-6)
	  for(int i = -1; i <=1; i++)
	    for(int j = -1; j <= 1; j++)
#if dimension>2
	      for(int k = -1; k<=1; k++)
#endif
		if(f[i,j,k] == 1)
		  b[] = b[i,j,k];// remember that you assume still the
				 // tag values are same as before
    }
  vector n[];
  scalar alpha[];
  reconstruction (f, n, alpha);
  coord segment[2];

  for(int l = 0; l < (2./DELTA); l++)
    {
      double pt[4] = {0, -HX, l*0.5*DELTA, l*0.5*DELTA};
      foreach()
	{
	  if(x > l*Delta && x < (l+1)*Delta && f[] > 1e-6 && f[] < 1. - 1e-6) 
	    {	
	      coord m = {n.x[],n.y[]};
	      if (facets (m, alpha[], segment) == 2)
		{
		  if(b[] == bn1)
		    {
		      if(!(pt[0] < (y + segment[0].y*Delta) && pt[0] < (y + segment[1].y*Delta)))
			{
			  if((y + segment[0].y*Delta) < (y + segment[1].y*Delta))
			    {
			      pt[0] = (y + segment[0].y*Delta);
			      pt[2] = (x + segment[0].x*Delta);
			    }
			  else
			    {
			      pt[0] = (y + segment[1].y*Delta);
			      pt[2] = (x + segment[1].x*Delta);
			    }
			}
		    }
		
		  else if(b[] == bn2)
		    {
		      if(!(pt[1] > (y + segment[0].y*Delta) && pt[1] > (y + segment[1].y*Delta)))
			{
			  if((y + segment[0].y*Delta) > (y + segment[1].y*Delta))
			    {
			      pt[1] = (y + segment[0].y*Delta);
			      pt[3] = (x + segment[0].x*Delta);
			    }
			  else
			    {
			      pt[1] = (y + segment[1].y*Delta);
			      pt[3] = (x + segment[1].x*Delta);
			    }
			}
		    }  

		  printf ("%g %g %g %g\n%g %g %g %g\n\n", 
			     x + segment[0].x*Delta, y + segment[0].y*Delta, t, hmin, 
			     x + segment[1].x*Delta, y + segment[1].y*Delta, t, hmin);
		}
	    }		    
	}

      //pt[2], pt[0] = x, y of topmost droplet and pt[3], pt[1] = x, y
      //of bottom one
      static FILE * fp1 = fopen("film_data.txt", "w");
      fprintf(fp1, "%g %g %g %g %g %g %d\n", t, pt[2], pt[0], pt[3], pt[1], pt[0] - pt[1], i);
      fflush(fp1);
      if(hmin >= (pt[0] - pt[1]))
	hmin = (pt[0] - pt[1]);
    }	  
  if(hmin > DELTA)
    {
      static FILE * fp = fopen("hmin_data.txt", "w");
      fprintf(fp, "%g %g\n", t, hmin);
      fflush(fp);
      //Rise velocity tracking
      double yb = 0., vb = 0., sb = 0.;
      foreach(reduction(+:yb) reduction(+:vb) reduction(+:sb)) {
	if(b[] == bn2)
	  {     
	    double dv = (f[])*dv();      
	    yb += y*dv;
	    vb += u.y[]*dv;
	    sb += dv;
	  }
      } 
      static FILE * fp2 = fopen("vel_data.txt", "w");
      fprintf (fp2, "%g %g %g %g\n", t, yb/sb, vb/sb, dt);                    
      fflush (fp2);
    }
 
  //Film velocity field
  foreach()
    {
      if( x < 0.5*D && y > -(hmin + D) && f[] == 0 && hmin <= HY)
	{
	  velf.x[] = u.x[];
	  velf.y[] = u.y[];
	  fprintf(stderr, "%g %g %g %g %g %g %d\n", t, x, y, u.x[], u.y[], hmin, i);
	  fflush(stderr);
	}
      else
	{
	  velf.x[] = 0.;
	  velf.y[] = 0.;
	}
    }
}

event movies(i+=50)
{
  box();
  view (tx = -0.484867, ty = 0.523272);
  draw_vof("f");
  squares("u.y");
  save("movie.mp4");
}
/*event dump_files(i+=300)
{
  char name[80];
  sprintf(name,"snap%g", t);
  dump(file = name);
}*/

event stop(t=STOP)
{
  dump("restart");
}                                                                                                                                                                                                   
event non_dim_nos(i=0)
{
  static FILE * fp = fopen("non_dimensional_numbers.txt", "w"); 
  fprintf (stderr," RHOR = %g\n MUR = %g\n AR = %g\n BO = %g\n D = %g\n g = %g\n rho_drop = %g\n", RHOR, MUR, AR, BO, D, GRA, rho1);          

}


/**

## Results

~~~gnuplot Minimum distance between drops as a function of time
set grid
set xlabel 'Time'
set ylabel 'Min-distance'
plot [0:14][0:]'hmin_data.txt' u 1:2 w l notitle
~~~


~~~gnuplot Rise velocity as a function of time
set grid
set xlabel 'Time'
set ylabel 'Velocity'
plot [0:14][0:] 'vel_data.txt' u 1:(2.25*$3) w l notitle
~~~
~~~gnuplot Evolution of film between drops over time

set xlabel 'x'
set ylabel 'y'
set term gif animate
set output 'film.gif' 
do for [i = 1:292] {

plot "film_data.txt" u 2:($7==(5*i)?($3<-0.1?$3:1/0):1/0) notitle w p, "film_data.txt" u 4:($7==(5*i)?($5>-3.5?$5:1/0):1/0) notitle w p 

}

~~~

~~~gnuplot Evolution of velocity inside film between drops over time

set xlabel 'x'
set ylabel 'y'
set term gif animate
set output 'velfilm.gif' 
do for [i = 160:292] {

plot "log" u ($7==5*i?$2:1/0):($7==5*i?$3:1/0):($7==5*i?$4:1/0):($7==5*i?$5:1/0) w vec notitle

}

~~~


![DNS of two approaching drops. The color field 
is the horizontal component of the velocity field.](tes/movie.mp4)

*/