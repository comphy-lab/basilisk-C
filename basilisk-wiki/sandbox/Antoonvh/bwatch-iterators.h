/**
# Iterators bwatch


## Concept 

Iterate over

* All camera rays  
* Each cell that the ray passes  
* Each cell where a ray intersects a facet  

 */
// A helper function
bool interfacial2 (scalar f, Point point) {
  if (f[] > 1e-6 && f[] < 1 - 1e-6)
    return true;
  return false;
}

void set_cr (coord * cr, coord * proj, coord * vert) {
  foreach_dimension() 
    cr->x = cam.poi.x - cam.O.x;
  normalize (cr);
  foreach_dimension() {
    proj->x = cam.up.x*cr->x;
    vert->x = cam.up.x - proj->x;
  }
}

void set_apoint (coord * apoint, coord hori, coord vert, int ii, int jj) {
  foreach_dimension() {
    apoint->x = (cam.poi.x
		 + cam.fov*((ii + 0.5)/cam.nx - 0.5)*hori.x
		 + ((float)cam.ny/cam.nx)*cam.fov*((jj + 0.5)/cam.ny - 0.5)*vert.x);  
  }
}

macro foreach_ray() {
  {
    coord cr, proj, vert, hori;
    set_cr (&cr, &proj, &vert);
    hori = cross (cr, vert);
    normalize (&vert);
    normalize (&hori);
    for (int jj = cam.ny; jj > 0; jj--) {    
      for (int ii = 0; ii < cam.nx; ii++) {
	coord apoint;
	set_apoint (&apoint, hori, vert, ii, jj);
	ray _r;
	_r.depth = 0;
	_r.O = cam.O;
	_r.dir.x = apoint.x - cam.O.x;
	_r.dir.y = apoint.y - cam.O.y;
	_r.dir.z = apoint.z - cam.O.z;
    
	normalize (&_r.dir);
	{...}
      }
    }
  }
}

macro foreach_ray_cell_intersection (ray r) {
  foreach_cell() { //MG acceleration
    coord _a[2] = {0};
    double _dist = ray_cell_intersect (r, point, _a);
    if (_dist < HUGE) {
      if (is_leaf(cell))
	{...}
    } else // Ray missed, skip children
      continue;
  }
}

macro foreach_ray_facet_intersection (ray r, scalar f) {
  {
    double _distg = HUGE;
    foreach_cell() { //MG acceleration
      coord _a[2] = {{0,0,0}, {0,0,0}}; //Pedantic...
      double _dist = ray_cell_intersect (r, point, _a);
      if (_dist < 0 || _dist < _distg) {
	if (interfacial2 (f, point)) {
	  if (is_leaf(cell)) {
	    coord _n;
	    if (segment_facet_intersect (_a, f, point, &_n)) {
	      normalize (&_n);
	      _distg = _dist;
	      {...}
	    }
	  }
	} else  // No fragment, skip children
	  continue;
      } else   // Ray missed or blocked, skip children
	continue; 
    }
  }
}

/**
## Workhorses

sketch should use optimized iterators,


 */
#include "PointTriangle.h"

void find_nearby_distances (scalar s) {
  assert (s.vofd.i);
  boundary ({s});
  scalar d = s.vofd;
  foreach() {
    coord cc = {x, y, z};
    double dist = 3*Delta*sign (s[] - 0.5);
    foreach_neighbor (1) {
      if (s[] > 1e-6 && s[] < 1. - 1e-6) {
	coord n = mycs (point, s);
	double alpha = plane_alpha (s[], n);
	coord v[12];
	int m = facets (n, alpha, v, 1.);
	// Using Alexis Berny's triangulation
	for (int j = 0; j < m - 2; j++) {
	  coord t1 = {x + v[0].x*Delta  , y + v[0].y*Delta  , z + v[0].z*Delta}; 
	  coord t2 = {x + v[j+1].x*Delta, y + v[j+1].y*Delta, z + v[j+1].z*Delta};
	  coord t3 = {x + v[j+2].x*Delta, y + v[j+2].y*Delta, z + v[j+2].z*Delta};
	  double s, t, dtr = PointTriangleDistance (&cc, &t1, &t2, &t3, &s, &t);
	  if (dtr > 0)
	    if (fabs(sqrt(dtr)) < fabs(dist))
	      dist = sqrt(dtr)*PointTriangleOrientation (&cc, &t1, &t2, &t3);
	}
      }
    }
    d[] = dist;
  }
  boundary ({d});
}

long int find_possible (scalar s, double val) {
  scalar posi = s.possible;
  int n = 0;
  foreach(reduction (+:n)) {//leaf
    posi[] = 0;
    bool pos = false, neg = false; 
    foreach_neighbor(1) {
      if (s[] < val)
	neg = true;
      if (s[] > val)
	pos = true;
    }
    if (pos && neg) {
      n++;
      posi[] = 1;
    }
  }
  for (int l = depth() - 1; l >=0 ; l--) {
    foreach_coarse_level(l) {
      double p = 0;
      foreach_child()
	if (posi[])
	  p = 1;
      posi[] = p;
    }
  }
  return n;
}

bool maybe (scalar posi, Point point) {
  if (posi[])
    return true;
  return false;
}

macro foreach_possible_ray_equi_intersection (ray r, scalar s, scalar posi, double DG) {
  foreach_cell() { //MG acceleration
    if (maybe(posi, point)) {
      coord _a[2] = {0};
      double _dist = ray_cell_intersect (r, point, _a);
      if (_dist < DG) {
	if (is_leaf(cell))
	  {...}
      } else // Ray missed or blocked, skip children
	continue;
    } else // Does not contain the isosurface
      continue;
  }
}

static inline bool precomp_segment_facet_intersect (coord a[2], double alpha, Point point, coord n) {
  coord cc = {x, y, z};
  double ALP = 0, ALP2 = 0;
  foreach_dimension() {
    ALP  += n.x*(a[0].x - cc.x)/Delta;
    ALP2 += n.x*(a[1].x - cc.x)/Delta;
  }
  if ((ALP2 - alpha)/(ALP - alpha) > 0.05) //~/wiki/sandbox/Antoonvh 5% gap filling
    return false;
  double w = fabs((ALP2 - alpha)) / (fabs(ALP - alpha) + fabs(ALP2 - alpha));
  foreach_dimension() {
    a[1].x = w*a[0].x + (1 - w)*a[1].x;
  }
  {;}
  return true;
}

static inline coord get_normal (scalar f, Point point) {
  coord n1 = {0};
  vector v = f.normals;
  foreach_dimension()
    n1.x = v.x[];
  return n1;
}

static inline double get_alpha (scalar f, Point point) {
  scalar s = f.possible;
  return s[];
}
  
macro foreach_precomp_ray_facet_intersection (ray r, scalar f, double DG) {
  {
    double _distg = HUGE;
    foreach_cell() { //MG acceleration
      if (interfacial2 (f, point)) {
	coord _a[2] = {{0,0,0}, {0,0,0}}; //Pedantic...
	double _dist = ray_cell_intersect (r, point, _a);
	if (_dist >= 0 && _dist < _distg && _dist < DG) {
	  if (is_leaf(cell)) {
	    coord _n = get_normal (f, point);
	    double alpha = get_alpha (f, point);
	    if (precomp_segment_facet_intersect (_a, alpha, point, _n)) {
	      normalize (&_n);
	      _distg = _dist;
	      {...}
	    }
	  }
	} else  // Ray missed or blocked, skip children
	  continue;
      } else   // No fragment, skip children
	continue; 
    }
  }
}

/**
   Iterator for volumetric rendering
 */
long int find_possible_vol (scalar s, double mval) {
  scalar posi = s.possible;
  int n = 0;
  foreach(reduction(+:n)) {//leaf
    posi[] = 0;
    bool nv = false;
    foreach_neighbor(1)
      if (s[] > mval)
	nv = true;
    if (nv)  {
      posi[] = 1;
      n++;
    } 
  }
  // restriction
  for (int l = depth() - 1; l >= 0; l--) {
    foreach_coarse_level(l) {
      double p = 0;
      foreach_child()
	if (posi[])
	  p = 1;
      posi[] = p;
    }
  }
  return n;
}

macro foreach_ray_cell_intersection_volume (ray r, double dm, scalar posi) {
  foreach_cell() { //MG acceleration
    if (maybe(posi, point)) {
      coord _a[2] = {0};
      double _dist = ray_cell_intersect (r, point, _a);
      if (_dist < dm) {
	if (is_leaf(cell)) //do something
	  {...}
      } else // Ray missed or blocked, skip children
	continue;
    } else // Nothing to see, skip children 
      continue;
  }
}
