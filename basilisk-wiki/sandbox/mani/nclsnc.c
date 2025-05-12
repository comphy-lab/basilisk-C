/**
![Avoiding coalescence using multiple VOF tracers. The color field 
 is the horizontal component of the velocity field.](nclsnc/movie_coalescence.mp4)
 
 The other interesting example in contrast to this is [control coalescence.](cntrl_clsnc/movie_coalescence.mp4)
*/

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
//double drainage_time = 2.;
//#include "control_no-coalescence.h"
#include "no-coalescence.h"
#include "view.h"
#include "navier-stokes/perfs.h"

double STOP=20.;
#define D 1.
int LEVEL=7;
double SIZE=6.;
#define height (SIZE*D)
#define L (SIZE*D)
#define I (0.5*(height-L))

double geometry(double x, double y) 
{

  double a=0., b=0., neg=-1., circ=1., i = 1., j = 1.;

  for (a=(I+(0.75*D)) ; a<=(I+L-0.75*D); a+=(1.5*D))
    {
      for (b=(height-(0.75*D)); b>=(0.75*D); b-=(1.5*D))
	{ 
	  circ = -sq(x-a) - sq(y-b) + sq(0.45*i*D);
	  if (circ >= 0)
	    {
	      neg = circ;  
	    }
	  if(j==1.)
	    j = -1.;
	  else
	    j = 1.;
	  i-=j*0.3;    	
	}
      if(j==1.)
	j = -1.;
      else
	j = 1.;
      i-=j*0.3;
    }
      
  return(neg);
}

/** no slip on all boundaries*/
u.t[top] = dirichlet(0.);
u.t[bottom]  = dirichlet(0.);
u.t[left] = dirichlet(0.);
u.t[right] = dirichlet(0.);

/**
   We make sure there is no flow through the top and bottom boundary,
   otherwise the compatibility condition for the Poisson equation can be
   violated. */

uf.n[left] = 0.;
uf.n[right] = 0.;
uf.n[top] = 0.;
uf.n[bottom] = 0.;


int main()
{

  size (height);
  
  init_grid (1 << LEVEL);

  rho1 = 1., mu1 = 0.01;

  rho2 = 1., mu2 = 0.01, f.sigma = 1.;

  
  run();
}

event init (t = 0)
{

  fraction(f, geometry(x,y));

  foreach()
    {
      if(x<=0.5*height)
	{
	  u.x[] = 2*f[];
	}
      else
	{
	  u.x[] = -2*f[];
	}
      if(y<=0.5*height)
	{
	  u.y[] = 2*f[];
	}
      else
	{
	  u.y[] = -2*f[];
	}
    }

}

/*scalar Ti[];
event tagging(i++)
{
  scalar A[];
  foreach()
    A[] = f[]>EPS;
  tag(A);
  foreach()
    Ti[] = A[];

    }*/

event movies(t=0; t<=STOP; i+=20)
{
  scalar ab[];
  box();
  view (tx = -0.484923, ty = -0.425958);
  for(scalar c in interfaces)
    {foreach()
	ab[] = c[];
      draw_vof("ab", lw = 2.);
    }
  squares("u.x");
  save("movie_coalescence.mp4");
}






