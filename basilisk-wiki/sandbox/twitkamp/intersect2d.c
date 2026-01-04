// #include "grid/multigrid.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "tavares/contact-embed.h"
#include "view.h"

#define R0 0.5
#define xc 0.
#define yc 0.575 
#define T 20

#define rl 1

#define MAXLEVEL 7
#define MINLEVEL 4

//Chongsen defines these maybe so the embed 
//does not cross the cell exactly
#define ed        L0/pow(2,MAXLEVEL)
#define oxy       ed/4.

double theta0, volume_vof_init;

bool segment_facet_intersect (coord a[2], scalar f, Point point, coord * n, coord isect[2]) {
  coord cc = {x, y, z};
  n[0]  = mycs (point, f); // interface normal
  double alpha = plane_alpha (f[], n[0]);

  double d0 = 0, d1 = 0;

  foreach_dimension() {
    d0 += n[0].x * (a[0].x - cc.x)/Delta;
    d1 += n[0].x * (a[1].x - cc.x)/Delta;
  }

  // segment endpoints distance to the plane
  double sd0 = d0 - alpha;
  double sd1 = d1 - alpha;

  // test if segment crosses
  if (sd0 * sd1 > 0.0) // same sign: no crossing
    return false;

  // find intersection point by linear interpolation
  double w = fabs(sd0)/(fabs(sd0) + fabs(sd1));

  foreach_dimension() {
    isect[0].x = a[0].x + w*(a[1].x - a[0].x);
  }

  // If you want to output the segment lying *on* the facet
  // you need to recover its full extent, which is the
  // intersection of the interface with the cell faces
  // (i.e. the polygonal VOF facet). That is more complex,
  // needing plane tracing across the cell faces.
  // Here, we only return the intersection point:
  isect[1] = isect[0]; // degenerate segment

  return true;
}


int main() {

  L0 = 4 * rl;

  origin(-L0/2, -L0/2);
    
  init_grid (1 << MAXLEVEL);
  /**
  We use a constant viscosity. */
  mu1 = mu2 = 0.1;
  
  /**
  We set the surface tension coefficient. */
  
  f.sigma = 1.;

  /**
  We vary the contact_angle. */

  for (theta0 = 45; theta0 <= 45; theta0 += 15) {
  	const scalar c[] = theta0*pi/180.;
  	contact_angle = c;
  	run();
  }
}

event init (t = 0)
{
  //embedded_solid
  solid(cs, fs , -x + 2. * oxy + y);
  // solid(cs, fs , y - oxy);

  // solid(cs, fs, y - x );
  //droplet
  fraction(f, - (sq(y + oxy) + sq(x + oxy) - sq(rl)));
}


event intersection_comp(t += 0.1, t <= T){


  char intersect_file[80];
  sprintf (intersect_file, "intersection-%d-%g-t_%.4f.dat", N, theta0, t);

  FILE * fp_inter = fopen (intersect_file,"w");

  coord isect[2];
  coord seg_dom[2];
  foreach(serial){
    if ((cs[] > 0.01 && cs[] < 0.99) && (f[] > 0.01 && f[] < 0.99) ) {
    // if (cs[] > 0.01 && cs[] < 0.99) {
      coord n = facet_normal (point, cs, fs);
      double alpha = plane_alpha (cs[], n);
      coord segment[2];
      // facets just outputs the number of intersections (0, 1, 2), 1 is the degenerate case. 
      // We are only interested in intersection so facets = 2
      // it also fills segment with the segment value
      if (facets (n, alpha, segment) == 2){
        seg_dom[0].x = x + segment[0].x*Delta; seg_dom[0].y = y + segment[0].y*Delta; 
        seg_dom[1].x = x + segment[1].x*Delta; seg_dom[1].y = y + segment[1].y*Delta; 
      }
      bool inters = segment_facet_intersect(seg_dom, f, point, &n, isect);
      if (inters){
        fprintf(fp_inter, "%g, %d, %g, %g, %g, %g\n", t, inters, x, y, isect[0].x, isect[0].y);
      } 
    }
  }
}

event facets(t += 0.1, t < T){

  char name_facets[80];
  sprintf (name_facets, "facets-%d-%g-t_%.4f.dat", N, theta0, t);

  FILE * fp_facet = fopen (name_facets,"w");
  output_facets (f, fp_facet);

}

event end (t = end)
{
  /**
  At the end, we output the equilibrium shape. */
  char name_facets[80];
  sprintf (name_facets, "facets-%d-%g.dat", N, theta0);
  FILE * fp_facet = fopen (name_facets,"w");
  output_facets (f, fp_facet);
  output_facets(cs, fout);

}
