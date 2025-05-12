/**
# Geometry functions
We define several geometric functions that our code can utilize, such as sphere
and rectangle. Initially, we need to define operations for intersection, union,
and difference between two geometric objects.
*/

#include "common.h"

#define intersection(a,b)   min(a,b)
#define union(a,b)          max(a,b)
#define difference(a,b)     min(a,-(b))

/**
## Rectangle
Parameters:

 * coord p: reference point
 * coord p1: rectangle bottom-left(-back) point
 * coord p2: rectangle top-right(-front) point
*/
double rectangle (coord p, coord p1, coord p2) {      
  #if dimension == 2
    double rec = intersection(intersection(intersection(p.x-p1.x, -p.x+p2.x),
          p.y-p1.y), -p.y+p2.y);
  #else // dimension == 3
    double rec = intersection(intersection(intersection(intersection(
          intersection(p.x-p1.x, -p.x+p2.x), p.y-p1.y), -p.y+p2.y),
          p.z-p1.z), -p.z+p2.z);
  #endif
  return rec;
}


/**
## Cylinder
Parameters:

 * coord p: reference coord
 * coord p1: first coord at the end of the cylinder
 * double r: radius of the cylinder
 * double h: length of the cylinder (only in 3D)
 * char axis: principal axis of the cylinder. 'x', 'y' and 'z' (only in 3D)
*/
#if dimension == 2
double cylinder (coord p, coord p1, double r) {
        double cyl = sq(r) - sq(p.x-p1.x) - sq(p.y-p1.y);
        return cyl;
}
#else // dimension == 3
double cylinder (coord p, coord p1, double r, double h, char axis) {
        double cyl, cond1, cond2;
        if (axis == 'x') {
                cyl = sq(r) - sq(p.y-p1.y) - sq(p.z-p1.z);
                cond1 = p.x - p1.x+h/2.;
                cond2 = -p.x + p1.x+h/2.;
        }
        else if (axis == 'y') {
                cyl = sq(r) - sq(p.x-p1.x) - sq(p.z-p1.z);
                cond1 = p.y - p1.y+h/2.;
                cond2 = -p.y + p1.y+h/2.;
        }
        else if (axis == 'z') {
                cyl = sq(r) - sq(p.x-p1.x) - sq(p.y-p1.y);
                cond1 = p.z - p1.z+h/2.;
                cond2 = -p.z + p1.z+h/2.;
        }
        else {
                fprintf(stderr, "Error on axis value in cylinder function.\
                    Accepted values are 'x', 'y' and 'z'.\n");
                return -1.;
        }
        double lcyl = intersection(intersection(cyl, cond1), cond2);;
        return lcyl;
}
#endif


/**
## Sphere
Parameters:

 * coord p: reference coord
 * coord p1: coord at the center of the sphere
 * double r: sphere radius
*/
double sphere (coord p, coord p1, double r) {
        #if dimension == 2
        double sph = sq(r) - sq(p.x-p1.x) - sq(p.y-p1.y);
        #else // dimension == 3
        double sph = sq(r) - sq(p.x-p1.x) - sq(p.y-p1.y) - sq(p.z-p1.z);
        #endif
        return sph;
}

/**
## Curved rectangle
Parameters:

 * coord p: reference coord
 * coord p1: rectangle bottom-left(-back) point
 * coord p2: rectangle top-right(-front) point
 * double r: corner radius
 */
double curvedRectangle (coord p, coord p1, coord p2, double r) {     
  assert (dimension == 2); // fixme: 2D only for the moment

  coord prec1_1 = {p1.x, p1.y+r};
  coord prec1_2 = {p2.x, p2.y-r};
  coord prec2_1 = {p1.x+r, p1.y};
  coord prec2_2 = {p2.x-r, p2.y};
  double recs = union(rectangle(p, prec1_1, prec1_2),
      rectangle(p, prec2_1, prec2_2));

  coord pcyl1 = {p1.x+r, p2.y-r};
  coord pcyl2 = {p2.x-r, p2.y-r};
  coord pcyl3 = {p1.x+r, p1.y+r};
  coord pcyl4 = {p2.x-r, p1.y+r};
  double cyls = union(union(union(cylinder(p, pcyl1, r),
      cylinder(p, pcyl2, r)), cylinder(p, pcyl3, r)), cylinder(p, pcyl4, r));
  double crec = union(recs, cyls);
        
  return crec;
}