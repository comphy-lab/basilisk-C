/**
# Distance to a fraction field

We use a slightly modified [`distance()`](/src/examples/distance.c) function
to find the distance field to a surface that is represented by volume
fractions.
 */
#include "fractions.h"
#include "distance2.h"
/**
Therefore, we need triangles; Using [Alexis Berny's
method](/sandbox/aberny/output_stl.h) for triangularization of a
fraction field.
 */
struct cff {
  scalar c;
  face vector f; // optional
};

coord * coords_from_fractions (struct cff input) {
  coord * p;
  scalar c = input.c;
  if (!input.f.x.i)
    input.f.x.i = -1;
  //We first count how many triangles there will be:
  long int nr_tri = 0, nr_p = 0;
  foreach(reduction(+:nr_tri)) {
    if (c[] > 1e-6 && c[] < 1. - 1e-6) {
      coord n = facet_normal (point, c, input.f);
      double alpha = plane_alpha (c[], n);
      coord v[12];
      nr_tri += facets (n, alpha, v, 1.) - 2;
    }
  }
  //Each triangle has 3 coordinates and we need a `nodata` coordinate:
  p = (coord*) malloc ((nr_tri*3 + 1)*sizeof(coord));
  // Fill `p` with triplets
  foreach(){
    if (c[] > 1e-6 && c[] < 1. - 1e-6) {
      coord n = facet_normal (point, c, input.f);
      double alpha = plane_alpha (c[], n);
      coord v[12];
      int m = facets (n, alpha, v, 1.); //m is the number of coordinates
      for (int j = 0; j < m - 2; j++) {
	coord temp1 = {x + v[0].x*Delta  , y + v[0].y*Delta  , z + v[0].z*Delta}; 
	coord temp2 = {x + v[j+1].x*Delta, y + v[j+1].y*Delta, z + v[j+1].z*Delta};
	coord temp3 = {x + v[j+2].x*Delta, y + v[j+2].y*Delta, z + v[j+2].z*Delta};
	p[nr_p++] = temp1;
	p[nr_p++] = temp2;
	p[nr_p++] = temp3;
      }
    }
  }
  coord w = {nodata, nodata, nodata};
  p[nr_p] = w;
  return p;
}
/**
# distance to surface

The user interface is provided via `distance_to_surface()`
 */


struct dts {
  scalar c;           // Volume fraction field
  face vector f;      // Optional face fraction field
  scalar d;           // Optional output: centered distance field
  vertex scalar phi;  // Optional output: distance at vertices
};

void distance_to_surface (struct dts dff) {
  assert (dimension == 3); //fixme: there should be a 2D analog
  struct cff nja = (struct cff){dff.c, dff.f};
  coord * p = coords_from_fractions (nja);
  
  scalar d = automatic(dff.d);

  distance (d, p);

  if (dff.phi.i) { //Distance at vertices
    vertex scalar phi = dff.phi;
    boundary ({d});
    foreach_vertex()
      phi[] = (d[] + d[-1] + d[0,-1] + d[-1,-1] +
	       d[0,0,-1] + d[-1,0,-1] + d[0,-1,-1] + d[-1,-1,-1])/8.;
    boundary ({phi});
  }
}
