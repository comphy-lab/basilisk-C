/** 
This code can reproduce the results in:

Jian,Z.,Deng P.,&Thoraval M.(2020). Air sheet contraction.Journal of Fluid Mechanics,899,A7.
*/
#include "navier-stokes/centered.h"
#define FILTERED 1
#include "two-phase.h"
#include "tension.h"
#define film_thickness_simulation 0.1
#define film_length_simulation 79.95
#define domain_length 100.00
#define Maxlevel 13
#define Minlevel 9
/** The geometrical shape of the air sheet is defined by this function. geometry(x,y) is positive outside the air sheet and negative inside the air sheet.*/
double geometry(double x, double y)
{
  double Line_up = domain_length/2.0 + film_thickness_simulation/2.0 - y;
  double Line_bottom = y - (domain_length/2.0 - film_thickness_simulation/2.0);
  double rectangle1 = min(Line_up,Line_bottom);
  double Line_right = film_length_simulation - x;
  double rectangle = min(rectangle1,Line_right);
  double circle = sq(film_thickness_simulation/2.)-sq(x-film_length_simulation)-sq(y-domain_length/2.);
  double geometry = max(rectangle,circle);
  return -geometry;
}
int main()
{
  size(domain_length);
  origin(0,0);
  init_grid(2048);
  rho1 = 1.;
  rho2 = 1./814.;
 /** this case is for air sheet with $Oh = 0.05$.*/
  mu1 = 0.005;
  mu2 = mu1/50.;
  f.sigma = 0.1;
  TOLERANCE = 1e-4;
  run();
}
/**
area around the interface is refined at the initial step.
*/
event init(i=0)
{
  refine(y-50.05<0.03 && y-50.05>-0.03 && x-79.95<0 && level <13);
  refine(y-49.95<0.03 && y-49.95>-0.03 && x-79.95<0 && level <13);
  refine(sq(x-79.95)+sq(y-50)<sq(0.08)&&sq(x-79.95)+sq(y-50)>sq(0.02)&&level <13);
  fraction(f,geometry(x,y));
}
/** 
As the adaptive mesh refinement will make the mesh coarse because the air sheet far from the rim is quiescent. We enforce the refinement after the adaptive mesh refinement for the area around the interface.
*/
event adapt(i++)
{
  adapt_wavelet({f,u.x,u.y},(double []){1e-6,1e-3,1e-3},13,9);
  refine(y-50.05<0.03&& y-50.05>-0.03 && level <13);
  refine(y-49.95<0.03&&y-49.95>-0.03&&level<13);
}

char filename[100];
/** Simulation output is stored.*/
event extract_position(t+=0.01)
{
  vector h[];
  heights(f,h);
  double xMax = -HUGE;
  double y3 = 0;
  foreach()
  {
    if(h.x[]!=nodata)
    {
      if((x+Delta*height(h.x[]))>xMax)
      {
        xMax = x+Delta*height(h.x[]);
        y3 = y;
      }
    }
  }
  fprintf(stderr,"%f %f %f\n",t,xMax,y3);
}
event end(t=100.)
{
  printf("END");

}
/**
After simulation is finished, post-processing of the dumpfiles is done to extract information.*/



