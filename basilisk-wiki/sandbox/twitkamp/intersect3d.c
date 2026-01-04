// #include "grid/octree.h"
// #include "embed.h"
// #include "navier-stokes/centered.h"
// #include "two-phase.h"
// #include "tension.h"
// #include "contact-embed.h"
// #include "view.h"

// #define R0 0.5
// #define xc 0.
// #define yc 0.575 
// #define T .1

// #define rl 1

// #define MAXLEVEL 5
// #define MINLEVEL 4

// //Chongsen defines these maybe so the embed 
// //does not cross the cell exactly
// #define ed        L0/pow(2,MAXLEVEL)
// #define oxy       ed/4.

// double theta0, volume_vof_init;


// bool segment_facet_intersect (coord a[3], scalar f, Point point, coord * n, coord isect[0]) {
//   coord cc = {x, y, z};
//   n[0]  = mycs (point, f); // interface normal
//   double alpha = plane_alpha (f[], n[0]);

//   double d0 = 0., d1 = 0.;

//   foreach_dimension() {
//     d0 += n[0].x * (a[0].x - cc.x)/Delta;
//     d1 += n[0].x * (a[1].x - cc.x)/Delta;
//   }

//   // segment endpoints distance to the plane
//   double sd0 = d0 - alpha;
//   double sd1 = d1 - alpha;

//   // test if segment crosses
//   if (sd0 * sd1 > 0.0) // same sign: no crossing
//     return false;

//   // find intersection point by linear interpolation
//   double w = fabs(sd0)/(fabs(sd0) + fabs(sd1));

//   foreach_dimension() {
//     isect[0].x = a[0].x + w*(a[1].x - a[0].x);
//   }
//   fprintf(ferr, "isect = %g %g %g\n", isect[0].x, isect[0].y, isect[0].z);
//   // If you want to output the segment lying *on* the facet
//   // you need to recover its full extent, which is the
//   // intersection of the interface with the cell faces
//   // (i.e. the polygonal VOF facet). That is more complex,
//   // needing plane tracing across the cell faces.
//   // Here, we only return the intersection point:
//   return true;
// }



// int main() {

//   L0 = 4 * rl;

//   origin(-L0/2., -L0/2., -L0/2.);
    
//   init_grid (1 << MAXLEVEL);
//   /**
//   We use a constant viscosity. */
//   mu1 = mu2 = 0.1;
  
//   /**
//   We set the surface tension coefficient. */
  
//   f.sigma = 1.;

//   /**
//   We vary the contact_angle. */

//   for (theta0 = 45; theta0 <= 45; theta0 += 15) {
//   	const scalar c[] = theta0*pi/180.;
//   	contact_angle = c;
//   	run();
//   }
// }

// event init (t = 0)
// {
//   vertex scalar phi[];

//   foreach_vertex(){
//     phi[] = -x + 2. * oxy + y;
//   }

//   fractions (phi, cs, fs);
//   boundary ({cs, fs});
//   fractions_cleanup (cs, fs);
//   boundary ({cs, fs});
//   restriction ({cs, fs});

//   //droplet
//   fraction(f, - (sq(y + oxy) + sq(x + oxy) + sq(z + oxy) - sq(rl)));
// }

// event intersection3D(t += .1, t <= T){
//   char intersect_file[80];
//   sprintf (intersect_file, "intersection-%d-%g-t_%.4f.dat", N, theta0, t);

//   FILE * fp_inter = fopen (intersect_file,"w");

//   coord isect[0];
//   foreach(serial){
//     if ((cs[] > 0.01 && cs[] < 0.99) && (f[] > 0.01 && f[] < 0.99) ) {
//     // if (cs[] > 0.01 && cs[] < 0.99) {
//       coord n = facet_normal (point, cs, fs);
//       coord v[12];
//       double alpha = plane_alpha (cs[], n);
//       int m = facets (n, alpha, v, 1.);
//       fprintf(ferr, "int m = %d \n", m);
//       coord segments[m];
//       // I think the logic makes sense but I need to update segment_facet_intersect to fill in the next isect 
//       for (int i = 0; i < m; i++){
//         segments[i].x = x + v[i].x * Delta;
//         segments[i].y = y + v[i].y * Delta;
//         segments[i].z = z + v[i].z * Delta;

//         bool inters = segment_facet_intersect(&segments[i], f, point, &n, isect[0]);
//         if (inters){
//           fprintf(fp_inter, "%g, %g, %g\n", isect[0].x, isect[0].y, isect[0].z);
//         } 
//       }
//       fprintf(fp_inter, "\n");
//     }
//   }
// }

// event contact_line_output (t += .1; t <= T) {
//   char fname[80];
//   sprintf(fname, "contact-line-theta%g-t%.2f.dat", theta0, t);
//   FILE * fp = fopen(fname, "w");

//   foreach() {
//     // Identify cells near the interface and near the embedded boundary
//     if (cs[] < 1. && f[] > 1e-3 && f[] < 1. - 1e-3) {
//       fprintf(fp, "%g %g %g\n", x, y, z);
//     }
//   }

//   fclose(fp);
// }


// event end (t = end)
// {

//   /**
//   At the end, we output the equilibrium shape. */
//   char name_facets[80];
//   sprintf (name_facets, "facets-%d-%g-t-%.2f.dat", N, theta0, t);
//   FILE * fp_facet = fopen (name_facets,"w");
//   output_facets (f, fp_facet);
// }




#include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "contact-embed.h"
#include "view.h"

#define R0 0.5
#define xc 0.
#define yc 0.575 
#define T 10

#define rl 1

#define MAXLEVEL 5
#define MINLEVEL 4

//Chongsen defines these maybe so the embed 
//does not cross the cell exactly
#define ed        L0/pow(2,MAXLEVEL)
#define oxy       ed/4.

double theta0, volume_vof_init;


bool segment_facet_intersect (coord a[2], scalar f, Point point, coord * n, coord isect[0]) {
  coord cc = {x, y, z};
  n[0]  = mycs (point, f); // interface normal
  double alpha = plane_alpha (f[], n[0]);

  double d0 = 0., d1 = 0.;

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

  fprintf(stderr, "isect %g %g %g\n", isect[0].x, isect[0].y, isect[0].z);
  // If you want to output the segent lying *on* the facet
  // you need to recover its full extent, which is the
  // intersection of the interface with the cell faces
  // (i.e. the polygonal VOF facet). That is more complex,
  // needing plane tracing across the cell faces.
  // Here, we only return the intersection point:

  return true;
}



int main() {

  L0 = 4 * rl;

  origin(-L0/2., -L0/2., -L0/2.);
    
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
  vertex scalar phi[];

  foreach_vertex(){
    phi[] = -x + 2. * oxy + y;
  }

  fractions (phi, cs, fs);
  boundary ({cs, fs});
  fractions_cleanup (cs, fs);
  boundary ({cs, fs});
  restriction ({cs, fs});

  //droplet
  fraction(f, - (sq(y + oxy) + sq(x + oxy) + sq(z + oxy) - sq(rl)));
}

event intersection3D(t += .1, t <= T){
  char intersect_file[80];
  sprintf (intersect_file, "intersection-%d-%g-t_%.4f.dat", N, theta0, t);

  FILE * fp_inter = fopen (intersect_file,"w");

  coord isect[0];
  foreach(serial){
    if ((cs[] > 0.01 && cs[] < 0.99) && (f[] > 0.01 && f[] < 0.99) ) {
    // if (cs[] > 0.01 && cs[] < 0.99) {
      coord n = facet_normal (point, cs, fs);
      coord v[12];
      double alpha = plane_alpha (cs[], n);
      int m = facets (n, alpha, v, 1.);
      coord corners[m];
      coord segment[2];
      int one_end[6] = {0, 0, 0, 1, 1, 2};
      int the_other_end[6] = {1, 2, 3, 2, 3, 3};
      // I think the logic makes sense but I need to update segment_facet_intersect to fill in the next isect 

      // this part is clearly not correct. I should create the segments then send them one by one in the gunction
      // Currently im just sending a single point? 
      for (int i = 0; i < m; i++){
        corners[i].x = x + v[i].x * Delta;
        corners[i].y = y + v[i].y * Delta;
        corners[i].z = z + v[i].z * Delta;
        fprintf(ferr, "seg %d %g %g %g\n", i, corners[i].x, corners[i].y, corners[i].z);
      }

      for (int seg = 0; seg < 6; seg++){
        
        int one_end_seg = one_end[seg];
        int the_other_end_seg = the_other_end[seg];
        segment[0].x = corners[one_end_seg].x; segment[0].y = corners[one_end_seg].y; segment[0].z = corners[one_end_seg].z;
        segment[1].x = corners[the_other_end_seg].x; segment[1].y = corners[the_other_end_seg].y; segment[1].z = corners[the_other_end_seg].z;
        bool inters = segment_facet_intersect(&segment, f, point, &n, isect);
        if (inters){
          fprintf(fp_inter, "%d %g, %g, %g\n", inters, isect[0].x, isect[0].y, isect[0].z);
        } 
      }

      fprintf(fp_inter, "\n");
    }
  }
}

event contact_line_output (t += .1; t <= T) {
  char fname[80];
  sprintf(fname, "contact-line-theta%g-t%.2f.dat", theta0, t);
  FILE * fp = fopen(fname, "w");

  foreach() {
    // Identify cells near the interface and near the embedded boundary
    if (cs[] < 1. && f[] > 1e-3 && f[] < 1. - 1e-3) {
      fprintf(fp, "%g %g %g\n", x, y, z);
    }
  }

  fclose(fp);
}


event end (t = end)
{

  /**
  At the end, we output the equilibrium shape. */
  char name_facets[80];
  sprintf (name_facets, "facets-%d-%g-t-%.2f.dat", N, theta0, t);
  FILE * fp_facet = fopen (name_facets,"w");
  output_facets (f, fp_facet);
}


