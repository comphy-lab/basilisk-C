/**
# Computation of a distance field from a discretised curve or surface

The goal is to compute the (signed) minimal distance from a set of
oriented segments (2d) or triangles (3D) to any point on the
grid. This signed distance function is also dynamically updated when
the mesh is refined or coarsened. 

The scheme is robust and will give results even for inconsistent
surface representations (i.e. surfaces with holes, manifold egdes,
incompatible orientations etc...). Faces do not need to be properly
connected i.e. the scheme works also for "triangle soups".

We need a few utility function such as computation of the minimal
distance between a point and a segment or a point and a triangle.
*/

#include <stdint.h>
#include "PointTriangle.h"

/**
This function computes the coordinates of the bounding box of a set of
segments or triangles. */

void bounding_box (coord * p, coord * min, coord * max)
{
  foreach_dimension()
    (*min).x = HUGE, (*max).x = - HUGE;
  while (p->x != nodata) {
    foreach_dimension() {
      if ((*p).x < (*min).x)
	(*min).x = (*p).x;
      if ((*p).x > (*max).x)
	(*max).x = (*p).x;
    }
    p++;
  }
}
  
/**
An extra field, holding a pointer to the elements (segments or
triangles) intersecting the neighborhood of the cell, is associated
with the distance function. The neighborhood is a sphere centered on
the cell center and with a diameter $3\Delta$. */

attribute {
  scalar surface;
}

#define double_to_pointer(d) (*((void **) &(d)))
#define BSIZE 3. // if larger than 1, cells overlap

/**
To compute the minimal distance and its sign, we need to store
information on the closest (2 in 2D, 12 in 3D) elements. */

typedef struct {
  double d2; // minimal distance
  coord * v; // the element
  int type;  // 0,1,2: vertex, 3: edge, 4: face
} closest_t;

#if dimension == 2
#  define ND 2
#else // dimension == 3
#  define ND 12
static double vertex_angle (coord * p, int type)
{
  if (type >= 3) // edge or face
    return pi;
  coord * v = p + type, * v1 = p + (type + 1)%3, * v2 = p + (type + 2)%3;
  coord vv1 = vecdiff (*v1, *v), vv2 = vecdiff(*v2,*v);
  return acos(vecdot(vv1,vv2)/sqrt(vecdot(vv1,vv1)*vecdot(vv2,vv2)));
}

static coord face_normal (coord * q, int type)
{
  coord ab = vecdiff(*(q+1),*q), ac = vecdiff(*(q+2),*q);
  coord n = vecdotproduct(ab,ac);
  double nn = sqrt(vecdot(n,n));
  assert (nn > 0.);
  nn = vertex_angle(q, type)/nn;
  foreach_dimension()
    n.x *= nn;
  return n;
}
#endif // dimension == 3

static void update_distance (Point point, coord ** i, scalar d)
{
  scalar surface = d.surface;
  Array * a = array_new();
  coord c = {x,y,z}, closest = {0};
  closest_t q[ND];
  for (int i = 0; i < ND; i++)
    q[i].d2 = HUGE;
  int nd = 0;
  double r2 = sq(BSIZE*Delta/2.);
  bool first = (level == 0);
  while (*i) {
    coord * p = *i;
#if dimension == 2
    coord r;
    double s, d2 = PointSegmentDistance (&c, p, p + 1, &r, &s);
#elif dimension == 3
    double s, t, d2 = PointTriangleDistance (&c, p, p + 1, p + 2, &s, &t);
#endif
    // keep pointers/distances/types of up to ND closest elements
    for (int i = 0; i < ND; i++)
      if (d2 < q[i].d2) {
	for (int j = ND - 1; j > i; j--)
	  q[j] = q[j-1];
	q[i].d2 = d2, q[i].v = p;
#if dimension == 2
	// vertices
	if (s == 0.)
	  q[i].type = 0;
	else if (s == 1.)
	  q[i].type = 1;
	else
	  // edge
	  q[i].type = 3;
	if (i == 0)
	  closest = r;
#elif dimension == 3
	// vertices
	if (s == 0. && t == 0.)
	  q[i].type = 0;
	else if (s == 1. && t == 0.)
	  q[i].type = 1;
	else if (s == 0. && t == 1.)
	  q[i].type = 2;
	else if (s == 0. || t == 0. || s + t == 1.)
	  // edge
	  q[i].type = 3;
	else
	  // face
	  q[i].type = 4;
	if (i == 0)
	  foreach_dimension()
	    closest.x = ((*q[0].v).x*(1. - s - t) + s*(*(q[0].v+1)).x +
			 t*(*(q[0].v+2)).x);
#endif // dimension == 3
	if (i >= nd)
	  nd = i + 1;
	break;
      }
    // add elements which are close enough to the local list
    if (d2 < r2 || first)
      array_append (a, &p, sizeof(coord *));
    first = false, i++;
  }
  if (a->len) {
    // set surface[] to list, ended with NULL
    coord * p = NULL;
    array_append (a, &p, sizeof(coord *));
    p = (coord *) array_shrink (a);
    assert (sizeof(double) >= sizeof(void *));
    memcpy (&surface[], &p, sizeof(void *));

    int orient;
#if dimension == 2
    if (q[0].type == 3)
      // edge
      orient = PointSegmentOrientation (&c, q[0].v, q[0].v + 1);
    else {
      // vertex
      if (nd == 1) // a single vertex, cannot find sign
	// get sign from parent
	orient = sign(bilinear (point, d));
      else { // two vertices
	orient = PointSegmentOrientation (&c, q[0].v, q[0].v + 1);
	if (orient != PointSegmentOrientation (&c, q[1].v, q[1].v + 1)) {
	  coord n = {0};
	  for (int i = 0; i < 2; i++) {
#if DEBUG	    
	    fprintf (stderr, "q %g %g %g %g\n",
		     x, y, q[i].v->x - x, q[i].v->y - y);
#endif
	    coord ab = vecdiff(*(q[i].v + 1),*q[i].v);
	    double nn = sqrt(vecdot(ab,ab));
	    assert (nn > 0.);
	    n.x -= ab.y/nn, n.y += ab.x/nn;
	  }
#if DEBUG
	  fprintf (stderr, "vertex %g %g %g %g %g %g %g %g\n",
		   x, y, closest.x - x, closest.y - y,
		   closest.x, closest.y, n.x, n.y);
#endif
	  coord diff = vecdiff(closest,c);
	  orient = sign(vecdot(n,diff));
	}
      }
    }
#elif dimension == 3
    if (q[0].type < 3) {
      coord n = face_normal (q[0].v, q[0].type), diff = vecdiff(closest,c);
      int nv = 1;
      orient = sign(vecdot(n,diff));
      for (int i = 1; i < nd && q[i].type < 3; i++)
	if (vecdist2 (closest, *(q[i].v + q[i].type)) < sq(1e-6)) {
	  coord n1 = face_normal (q[i].v, q[i].type);
	  foreach_dimension()
	    n.x += n1.x;
	  nv++;
	  if (orient > -2 && sign(vecdot(n1,diff)) != orient)
	    orient = -2;
	}
      if (nv < 3 ) {// less than 3 vertices, cannot find sign
	if (level != 0)
	  orient = sign(bilinear (point, d));
	else orient = 1;
      }else if (orient == -2) {
	// the vertices do not have the same orientation
	// get the proper orientation from the pseudo-normal n
	orient = sign(vecdot(n,diff));
#if DEBUG
	fprintf (stderr, "vertex %g %g %g %g %g %g %d %g %g %g %g %g %g %d\n",
		 x, y, z, closest.x - x, closest.y - y, closest.z - z, orient,
		 closest.x, closest.y, closest.z, n.x, n.y, n.z, nv);
#endif
      }
    }
    else if (q[0].type == 3) {
      // edge
      if (nd == 1 || q[1].type != 3) // a single edge, cannot find sign
	// get sign from parent
	orient = sign(bilinear (point, d));
      else { // two edges
	orient = PointTriangleOrientation (&c, q[0].v, q[0].v+1, q[0].v+2);
	if (orient !=
	    PointTriangleOrientation (&c, q[1].v, q[1].v+1, q[1].v+2)) {
	  coord n1 = face_normal (q[0].v, 3), n2 = face_normal (q[1].v, 3), n;
	  foreach_dimension()
	    n.x = n1.x + n2.x;
	  coord diff = vecdiff(closest,c);
	  orient = sign(vecdot(n,diff));
#if DEBUG
	  fprintf (stderr, "edge %g %g %g %g %g %g %d %g %g %g %g %g %g\n",
		   x, y, z, closest.x - x, closest.y - y, closest.z - z, orient,
		   closest.x, closest.y, closest.z, n.x, n.y, n.z);
#endif
	}
#if DEBUG
	else
	  fprintf (stderr, "edge %g %g %g %g %g %g %d\n",
		   x, y, z, closest.x - x, closest.y - y, closest.z - z, orient);
#endif
      }
    }
    else { // face
#if DEBUG
      fprintf (stderr, "face %g %g %g %g %g %g\n",
	       x, y, z, closest.x - x, closest.y - y, closest.z - z);
#endif
      orient = PointTriangleOrientation (&c, q[0].v, q[0].v+1, q[0].v+2);
    }
#endif // dimension == 3
    d[] = sqrt (q[0].d2)*orient; 
  }
  else { // !a->len
    free (a);
    surface[] = 0.;
    if (level > 0)
      d[] = bilinear (point, d);
    else
      d[] = 0.;
  }
}

#undef ND

static void refine_distance (Point point, scalar d)
{
  scalar surface = d.surface;
  if (surface[] == 0.)
    foreach_child() {
      surface[] = 0.;
      d[] = bilinear (point, d);
    }
  else {
    coord ** ap = (coord **) double_to_pointer (surface[]);
    int s = 0;
    foreach_child() {
      update_distance (point, ap, d);
      s += sign(d[]);
    }

    /**
    To increase robustness to inconsistent input, we check whether all
    children are included within the minimum distance sphere. If this
    is the case then the children and parent must have the same
    orientation. We enforce this, using the "average" orientation. */

    if (fabs(d[]) > sqrt(dimension)*Delta/4.) {
      if (abs(s) != 1 << dimension) {
	s = sign(s);
	foreach_child()
	  d[] = s*fabs(d[]);
      }
      if (sign(d[]) != sign(s))
	d[] = - d[];
    }
  }
}

static void restriction_distance (Point point, scalar d) {}

static void coarsen_distance (Point point, scalar d) {
  scalar surface = d.surface;
  foreach_child()
    free (double_to_pointer (surface[]));
}

static void delete_distance (scalar d) {
  scalar surface = d.surface;
  foreach_level (0)
    free (*((void **)double_to_pointer (surface[])));
  for (int l = 0; l <= depth(); l++)
    foreach_level (l)
      free (double_to_pointer (surface[]));
  delete ({surface});
}

trace
void distance (scalar d, coord * p)
{
  scalar surface = d.surface;
  if (surface.i)
    delete_distance (d);
  surface = new scalar;
  surface.restriction = no_restriction;
#if TREE
  surface.prolongation = no_restriction;
  surface.refine = no_restriction; // handled by refine_distance()
  d.prolongation = refine_bilinear;
  d.refine = refine_distance;  
  d.coarsen = coarsen_distance;
#endif
  d.surface = surface;
  d.delete = delete_distance;
  d.restriction = restriction_distance;

  Array * a = array_new();
  while (p->x != nodata) {
#if dimension == 3
    // filter degenerate triangles
    coord ab = vecdiff(*(p+1),*p), ac = vecdiff(*(p+2),*p);
    coord n = vecdotproduct(ab,ac);
    if (vecdot(n,n) > 0.)
#endif
      array_append (a, &p, sizeof (coord *));
    p += dimension;
  }
  p = NULL;
  array_append (a, &p, sizeof (coord *));
  p = (coord *) array_shrink (a);

  foreach_level(0)
    update_distance (point, (coord **) p, d);
  free (p);
  
  boundary_level ({d}, 0);
  for (int l = 0; l < depth(); l++) {
    foreach_coarse_level (l)
      refine_distance (point, d);
    boundary_level ({d}, l + 1);
  }
}
