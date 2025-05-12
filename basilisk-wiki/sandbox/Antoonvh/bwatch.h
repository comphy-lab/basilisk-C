/**
# Bwatch

`Bwatch` is a ray-casting toolbox written in Basilisk C. It should
work with Octrees or multigrid3D. 

## "Installation"

<img style = float:right src="minbwatch/min.png" width="340px"/>

Apart from the default Basilisk code, the implementation is scattered
among a few header files.

Required files are:  

* `Bwatch.h` (this file) contains a bunch of generic definitions
* [`Sketch.h`](sketch.h) implements high-level rendering functions and
most of the user interface
* [`Bwatch-iterators.h`](bwatch-iterators.h) defines low-level
functions and iterators to help find ray-object intersections
efficiently.
* [`Radix.h`](radix.h) provides a sorting algorithm that maybe used for volumetric rendering.  

Copy and place them in your current folder, a
`$BASILISK_INCLUDE_PATH` or against my advice in `$BASILISK`.

Optional is:  

* [installing an adaptive raycaster](raycasterv.c) or ([this
one](raycaster.c)) if you want to use
[`store_adaptive()`](sketch.h#store-adaptive) 
* Having [`convert`](https://linux.die.net/man/1/convert) installed on
your system for `image` and `text` rendering.
* [Watching along](sketch.h#watch-the-sketch-process) the
rendering process

## How does it work?

Rays are casted from the camera into the scene. The nearest
intersection is computed with the obects that the user called. The
standard lightning/shading model for coloured materials follows that of
[Bui Tuong
Phong](https://en.wikipedia.org/wiki/Phong_reflection_model). 
 
## tips

1. Compile your bwatch code with the `-fopenmp` (gcc) or `-qopenmp`
(icc) compiler flag (before the `*.c` file) to accelerate rendering.  
2. When your rendering invokes secondary rays (e.g. transmission, reflection etc.), make sure the scene has a background. For example, a large `dull` sphere. 
3. Compile with `-disable-dimensions` to avoid qcc crashing in `get_dimension()` routine. 

## Implementation
*/
#include "fractions.h"
/**
A list of sketch functions is stored. Rendering takes place when an
image is `store()`d.
 */

#define N_SKETCH 256          // Maximum sketch functions
typedef double (*sketch)();   // a Sketch function prototype   
sketch sketch_list[N_SKETCH]; // List sketching functions
int functot;                  // number of sketching function calls
scalar * shading = NULL;      // VOF objects that cast shade
scalar smoke;                 // Concentration field that cast shade
//Prim * Obj  = NULL;         // Primitives should cast a shade as well
int max_ray_depth = 5;        // Number of ray recasts
attribute {                   // attenuation coefficient of smoke fields  
  double att;                 
}

// A camera class and "the" cam
struct Camera {
  coord O;         // camera pos
  coord up;        // upward pointing vector (related to image vertical)
  coord rhs;       // In case up-vector is parallel to the viewing direction
  coord poi;       // Pointing towards
  double fov;      // horizontal width fustrum through poi.
  int nx, ny;      // Resolution
  double min, max; // Sensitivities
  double (*f)();   // Transfer function
};

struct Camera cam = {.O   = {0., 0., 3},   // view from back to front (-z)
		     .up  = {0., 1., 0.},
		     .rhs = {1., 0., 0.},  
		     .poi = {0., 0., 0.},  // view from back to front (-z)
		     .fov = 1.1,
		     .nx = 450,
		     .ny = 400,
		     .min = 0.1,
		     .max = 100,
		     .f = log};
/**
## User function for the Camera setup
 */

bool watch (struct Camera inp) {
  foreach_dimension() {
    if (inp.O.x)
      cam.O.x = inp.O.x;
    if (inp.up.x)
      cam.up.x = inp.up.x;
    if (inp.rhs.x)
      cam.rhs.x = inp.rhs.x;
    if (inp.poi.x)
      cam.poi.x = inp.poi.x;
  }
  if (inp.fov)
    cam.fov = inp.fov;
  if (inp.nx)
    cam.nx = inp.nx;
  if (inp.ny)
    cam.ny = inp.ny;
  if (inp.min)
    cam.min = inp.min;
  if (inp.max)
    cam.max = inp.max;
  if (inp.f)
    cam.f = inp.f;
  return true;
}

// A ray class
typedef struct ray {
  coord O;
  coord dir;
  int   depth;   // for recasted rays
} ray;

// Lights
typedef struct light {
  int type;             // ambient (1), Sun (2), point (3), ...
  unsigned char col[3]; // Color
  double I;             // Intensity (at cam.poi)
  coord O;              // Location of point light
  coord dir;            // Direction of sun light
} light;

light lights[N_SKETCH]; //Reuse N_SKETCH
/**
## Material properties

Properties for matierial objects
*/
#include "utils.h"

struct Material {
  unsigned char col[3];   // prescribed color
  scalar s;               // color from field data
  Colormap map;           // Corresponding colorbar ...
  vector v;               // RGB from vector
  double min, max;        // ... cbar range
  bool linear;            // Cbar interpolation;
  double SPR;             // Specular albedo   (Blinn)
  double SPexp;           // Specular exponent (Blinn)
  double ind;             // refraction index (normal points outwards)
  double T;               // Transparancy  (0 - 1 intended)
  double R;               // Reflectivity  (0 - 1 intended)
  bool dull;              // Dull object: no light interaction just `col`
  unsigned char col2[3];  // Or a gradient to col2 (not so dull!)
  coord n1;               // normal for `col`, col2 at n = -n1;
};

typedef struct Material material;
/**
## A few (vector) helper functions
*/
void bias (ray * r) {
  double dr = L0*1e-4;
  foreach_point (r->O.x, r->O.y, r->O.z, serial)
    dr = Delta;
  foreach_dimension()
    r->O.x += r->dir.x*dr;
}

void bias1 (ray * r) {
  foreach_dimension()
    r->O.x += r->dir.x*L0*1e-4;
}

bool any_col (unsigned char px[3]) {
  if (px[0] || px[1] || px[2])
    return true;
  return false;
}

bool any_comp (coord a) {
  if (a.x || a.y || a.z)
    return true;
  return false;
}

coord cross (coord a, coord b) {
  coord new;
  foreach_dimension()
    new.x = a.y*b.z - a.z*b.y;
  return new;
}

double dot (coord a, coord b) {
  double d = 0;
  foreach_dimension()
    d += a.x*b.x;
  return d;
}

coord reflect (coord in, coord n) {
  coord out;
  foreach_dimension()
    out.x = in.x - 2*dot (in, n)*n.x;
  return out;
}
/**
   Normalization with the approximate inverse square root. Quake3
Arena method of Greg
Welsh, [see](https://en.wikipedia.org/wiki/Fast_inverse_square_root).

Premature optimzation is the fast inverse square root of all evil...
- *Donald Knuth*
 */
static inline double Q_isqrt (double l) {
  double x1 = l*0.5;
  union {
    double f;
    long long i;
  } conv  = { .f = l };
  conv.i = 0x5fe6eb50c7b537a9 - (conv.i >> 1);
  for (int i = 0; i < 2; i++)
    conv.f *= 1.5 - (x1 * conv.f * conv.f);
  return conv.f;
}

void normalize_(coord * p) {
  double l = 0;
  foreach_dimension()
    l += sq(p->x);
  l = Q_isqrt (l);
  foreach_dimension()
    p->x *= l;
}
#define normalize normalize_

// Thanks to scratch pixel:
// https://www.scratchapixel.com/lessons/3d-basic-rendering/
// introduction-to-shading/reflection-refraction-fresnel
coord refract(coord in, coord normal, double ind) { 
  normalize (&in);
  normalize (&normal);
  double cosi = dot (in, normal); 
  double etai = 1, etat = ind; 
  coord n;
  foreach_dimension()
    n.x = normal.x; 
  if (cosi < 0) 
    cosi = -cosi;
  else {
    etai = etat;
    etat = 1.;
    foreach_dimension()
      n.x = -normal.x;
  }
  double eta = etai/etat; 
  double k = 1. - sq(eta) * (1. - sq(cosi));
  assert (k >= 0);
  coord out = {0.};
  foreach_dimension()
    out.x = k < 0 ? 0 : eta*in.x + (eta*cosi - sqrt(k))*n.x;
  return out; 
}

double vec_length (coord v) {
  double l = 0;
  foreach_dimension()
    l += sq(v.x);
  return l > 0 ? sqrt(l): 0;
}

// Thanks you scratchpixel:
// https://www.scratchapixel.com/lessons/3d-basic-rendering/
// introduction-to-shading/reflection-refraction-fresnel
double fresnel(coord in, coord normal, double ind) { 
  normalize (&in);
  normalize (&normal);
  double cosi = dot(in, normal); 
  double etai = 1, etat = ind;
  if (cosi > 0) {
    etai = ind;
    etat = 1.;
  } 
  // Compute sini using Snell's law
  double sint = etai/etat * sqrt(max(0., 1 - sq(cosi))); 
  // Total internal reflection
  if (sint >= 1)
    return 1.;
  else { 
    double cost = sqrt(max(0., 1 - sq(sint))); 
    cosi = fabs(cosi); 
    double Rs = ((etat*cosi) - (etai*cost)) / ((etat*cosi) + (etai*cost)); 
    double Rp = ((etai*cosi) - (etat*cost)) / ((etai*cosi) + (etat*cost)); 
    return (sq(Rs) + sq(Rp))/2; 
  } 
}


/**
## Intersections

One may compute

* ray-cell intersection 
* ray-cuboid interections
* ray-plane intersection
* segment-facet intersection
* ray-triangle intersection
 */
static inline double ray_cell_intersect (ray r, Point point, coord a[2]) {
  coord cc = {x,y,z};
  double in = -HUGE, out = HUGE;
  foreach_dimension() {
    in  = max(in,  r.dir.x != 0 ? ((cc.x - (sign(r.dir.x)*Delta/2.)) - r.O.x)/r.dir.x : -HUGE);
    out = min(out, r.dir.x != 0 ? ((cc.x + (sign(r.dir.x)*Delta/2.)) - r.O.x)/r.dir.x :  HUGE);
  }
  if (in >= out || out < 0) 
    return HUGE;
  in = in < 0 ? 0 : in; //The origin is in the cell.
  foreach_dimension() {
    a[0].x = r.O.x +  in*r.dir.x;
    a[1].x = r.O.x + out*r.dir.x;
  }
  return in;
}

static inline double ray_box (ray r, coord bb[2]) {
  double in, out;
  in = -HUGE; out = HUGE;
  foreach_dimension() {
    if (r.dir.x > 0) {
      in  = max(in,  (bb[0].x - r.O.x)/r.dir.x);
      out = min(out, (bb[1].x - r.O.x)/r.dir.x);
    } else if (r.dir.x < 0) {
      in  = max(in,  (bb[1].x - r.O.x)/r.dir.x);
      out = min(out, (bb[0].x - r.O.x)/r.dir.x);
    }
  }
  if (in >= out || out < 0) 
    return HUGE;
  in = in < 0 ? 0 : in; //The origin is in the box.
  return in;
}

static inline double ray_plane_intersect (ray r, coord n, double alpha, coord * a) {
  double ldotn = 0;
  foreach_dimension()
    ldotn += n.x*r.dir.x;
  if (ldotn == 0) //parallel
    return HUGE;
  double d = alpha;
  foreach_dimension()
    d -= n.x*r.O.x;
  d /= ldotn;
  if (d <= 0) //plane certainly not in view.
    return -1;
  d = d > 1e4*L0 ? 1e4*L0 : d;
  foreach_dimension() {
    a[0].x = r.O.x + d*r.dir.x;
  }
  return d;
}

static inline double ray_sphere_intersect (ray r, coord C, double R, coord * a, coord * n) {
  normalize (&r.dir);
  coord oc;
  foreach_dimension()
    oc.x = r.O.x - C.x;
  double det = (sq(dot(r.dir, oc)) - (sq(vec_length (oc)) - sq(R)));
  if (det < 0)
    return HUGE;
  det = sqrt(det);
  double dist  = -dot(r.dir, oc) - det;
  double dist2 = -dot(r.dir, oc) + det;
  if (dist < 0 && dist2 < 0)
    return HUGE;
  else if (dist < 0)
    dist = dist2;
  else if (dist > dist2)
    dist = dist2;
  foreach_dimension() {
    a->x = r.O.x + r.dir.x*dist;
    n->x = a->x - C.x;
  }
  normalize (n);
  return dist;
}

// use ray_plane_intersect?
static inline bool segment_facet_intersect (coord a[2], scalar f, Point point, coord * n) {
  coord cc = {x, y, z};
  n[0]  = mycs (point, f); //pointing outwards.
  double alpha = plane_alpha (f[], n[0]);
  double ALP = 0, ALP2 = 0;
  foreach_dimension() {
    ALP  += n[0].x*(a[0].x - cc.x)/Delta;
    ALP2 += n[0].x*(a[1].x - cc.x)/Delta;
  }
  if ((ALP2 - alpha)/(ALP - alpha) > 0.05) // 3% gap filling
    return false;
  double w = fabs((ALP2 - alpha)) / (fabs(ALP - alpha) + fabs(ALP2 - alpha));
  foreach_dimension() {
    a[1].x = w*a[0].x + (1 - w)*a[1].x;
  }
  return true;
}
// https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
double ray_triangle (ray r, coord t[3]) {
  double eps = 1e-6;
  coord e1, e2, h, s ,q;
  double a, f, u, v;
  foreach_dimension() {
    e1.x = t[1].x - t[0].x;
    e2.x = t[2].x - t[0].x;
  }
  h = cross (r.dir, e2);
  a = dot (e1, h);
  if (a > -eps && a < eps)
    return HUGE;
  f = 1./a;
  foreach_dimension()
    s.x = r.O.x - t[0].x;
  u = f*dot (s, h);
  if (u < 0 || u > 1.0)
    return HUGE;
  q = cross (s, e1);
  v = f*dot (r.dir, q);
  if (v < 0 || (v + u) > 1)
    return HUGE;
  double d = f*dot (e2, q);
  if (d > eps)
    return d;
  return HUGE;
}
/**
## Iterators

see,
 */
attribute {
  scalar possible; // Marker for possible interesting child
  vector center;   // Procomputed center of vof facets? (check me!)
  vector normals;  // precomputed normal of vof facets/isosurfaces
  scalar vofd;     // Nearby vof distance field
}
#include "bwatch-iterators.h"
/**
## Color from colorbar 

Quantative data can be mapped to an RGB-color code via a colorbar. The
function below is a version of `colormap_color()` in
[output.h](/src/output.h).
 */

bool colormap_pigmentation (unsigned char c[3], double cmap[NCMAP][3],
			    double val, double min, double max) {
  if (val == nodata) {
    c[0] = c[1] = c[2] = 0; // nodata is black
    return false;
  }
  int i;
  double coef;
  if (max != min)
    val = (val - min)/(max - min);
  else
    val = 0.;
  if (val <= 0.) i = 0, coef = 0.;
  else if (val >= 1.) i = NCMAP - 2, coef = 1.;
  else {
    i = val*(NCMAP - 1);
    coef = val*(NCMAP - 1) - i;
  }
  assert (i < NCMAP - 1);
  for (int j = 0; j < 3; j++)
    c[j] = 255*(cmap[i][j]*(1. - coef) + cmap[i + 1][j]*coef);
  return true;
 }
/**
## Lightning

Lights can be shaded by vof facets (in `shading`) and should be
reflected by reflectors (in `rlist`)
 */
bool shaded (coord start, light source) {
  ray l;
  l.O = start;
  foreach_dimension()
    l.dir.x = -source.dir.x;
  for (scalar c in shading) 
    foreach_ray_facet_intersection(l, c)
      return true;
  return false;
}

double Intensity_by_reflection (coord start, light source) {
  return 0;
}

void default_lights (void) {
  // Ambient and sun;
  lights[0].type = 1;
  lights[0].I = 2;
  lights[0].col[0] = 255;
  lights[0].col[1] = 255;
  lights[0].col[2] = 255;
  
  lights[1].type = 2;
  lights[1].I = 80;
  lights[1].dir = (coord){-1, -1.8, -2};
  lights[1].col[0] = 255;
  lights[1].col[1] = 255;
  lights[1].col[2] = 250;
}

double absorbsion (coord a, int l) {
  ray r = {.O = a, .dir = lights[l].dir};
  foreach_dimension()
    r.dir.x *= -1;
  normalize (&r.dir);
  double sl = 0;
  foreach_ray_cell_intersection_volume (r, HUGE, smoke.possible) {
    coord mean = {0, 0, 0};
    double len = 0;
    foreach_dimension() {
      mean.x = (_a[0].x + _a[1].x)/2.;
      len += sq(_a[0].x - _a[1].x);
    }
    double val = interpolate (smoke, mean.x, mean.y, mean.z);
    len = sqrt(len);
    if (len > 0)
      sl += fabs(val)*len;
  }
  if (sl > 0)
    return exp(-sl/smoke.att);
  else
    return 1;
}

double diffuse_Lambert (coord n, coord a) {
  double I = 0;
  int l = 0; //light
  while (lights[l].type > 0) { //Fixme: No mirrored light sources yet
    if (lights[l].type > 1) {
      double Il = 0;
      coord ldir = {0,0,0};
      if (lights[l].type == 2)
	ldir = lights[l].dir;
      if (!shaded (a, lights[l])) {
	normalize (&ldir);
	foreach_dimension()
	  Il -= n.x * ldir.x;
	double IT = lights[l].I;
	if (smoke.i) {
	  IT *= absorbsion (a, l);
	}
	I += max (0, IT*Il);
      }
    }
    l++;
  }
  return I;
}

bool specular (ray r, coord n, coord a, double gloss, double alR, unsigned char c[3]) {
  int l = 0;
  normalize (&n);
  while (lights[l].type > 0) { //Fixme: No mirrored light sources yet
    if (lights[l].type > 1) {
      if (lights[l].type == 2) { //Fixme: Only sun
	if (shaded (a, lights[l]))
	  return false;
	coord R = reflect (lights[l].dir, n);
	normalize (&R);
	double RdN = -dot(r.dir, R);
	if (RdN < 0)
	  return false;
	double I = max(alR*lights[l].I*pow(RdN, gloss), 1e-6);
	for (int i = 0; i < 3; i++) {
	  c[i] = lights[l].col[i]* min(max(log(I) - log(cam.min),0)/log(cam.max/cam.min), 1);
	  c[i] = min (c[i], 255);
	}
      }
    }
    l++;
  }
  return true;
}

double light_intensity (coord n, coord a) {
  double I = diffuse_Lambert(n, a);
  int l = 0;
  //ambients
  while (lights[l].type > 0) {
    if (lights[l].type == 1) 
      I += lights[l].I;
    l++;
  }
  I = min (max(cam.f(I) - cam.f(cam.min),0)/cam.f(cam.max/cam.min), 1);
  return I;  
}

struct _get_pixel{
  unsigned char c[3];
  coord n;
  coord a;
  ray r;
  double SPexp;
  double SPR;
};
  
void get_pixel (struct _get_pixel * inp) {
  double I = light_intensity (inp->n, inp->a);
  for (int i = 0; i < 3; i++)
    inp->c[i] *= I;
  if (inp->SPexp && inp->SPR) {
    unsigned char tmp[3];
    if (specular (inp->r, inp->n , inp->a, inp->SPexp, inp->SPR, tmp))
      for (int i = 0; i < 3; i++) 
	inp->c[i] = min(inp->c[i] + tmp[i], 255);
  }
}
/**
# Pixel color from a material hit (maybe recursive)

Mind the presidence certain properties take over others:

1. Dull
2. Refractive
3. Fully reflective (R = 1)
4. Prescribed color 
5. Color from scalar field data and color bar
6. RGB Color code from vector field 

4, 5 and 6 maybe supplemented with *either* 7 or 8,

7. Reflection (R < 1)
8. Transmission (transparancy, without refraction) 
 */
double get_color_ray(); //Prototype

trace
bool get_color_material (ray r, coord a, coord n, material mat, unsigned char px[3]) {
  // Dull object?
  if (mat.dull) {
    if (any_comp(mat.n1) && any_comp(n)) {
      normalize (&n); normalize (&mat.n1);
      double l = (dot (n, mat.n1) + 1.)/2.;
      for (int i = 0; i < 3; i++)
	px[i] = l*mat.col[i] + (1. - l)*mat.col2[i];
      return true;
    }
    for (int i = 0; i < 3; i++)
      px[i] = mat.col[i];
    return true;
  }
  // Refractive object?
  else if (mat.ind) {
    if (r.depth <= max_ray_depth) {
      normalize (&n);  normalize (&r.dir);
      ray refl, refr;
      refl.depth = refr.depth = r.depth + 1;
      refl.O = refr.O = a;
      refl.dir = reflect (r.dir, n);
      normalize (&refl.dir);
      double w = fresnel (r.dir, n, mat.ind);
      if (w < 1) {
	refr.dir = refract (r.dir, n, mat.ind);
	normalize (&refr.dir);
      }
      bias (&refl); bias (&refr);
      unsigned char ptrefl[3], ptrefr[3];
      get_color_ray (refl, ptrefl);
      if (w < 1) 
	get_color_ray (refr, ptrefr);
      else
	w = 1;
      for (int i = 0; i < 3; i++)
	px[i] = w*ptrefl[i] + (1 - w)*ptrefr[i];
      if (mat.SPexp && mat.SPR) {
	unsigned char tmp[3];
	if (specular (r, n , a, mat.SPexp, mat.SPR, tmp))
	  for (int i = 0; i < 3; i++) 
	    px[i] = min(px[i] + tmp[i], 255);
      }
      return true;
    }
    return false;
  }
  // Full Reflector?
  if (mat.R >= 1) {
    if (r.depth <= max_ray_depth) {
      ray new = {.O = a, .dir = reflect (r.dir, n), .depth = r.depth + 1};
      bias (&new);
      if ((get_color_ray (new, px) == HUGE))
	for (int i = 0; i < 3; i++)
	  px[i] = mat.col[i];
      return true;
    }
    return false;
  }
  // It must be a "tradional" object:
  assert (any_col (mat.col) || mat.s.i || mat.v.y.i || mat.R >= 1 );
  struct _get_pixel out;
  out.r = r;
  out.SPexp = mat.SPexp;
  out.SPR   = mat.SPR;
  out.n = n;
  out.a = a;
  if (any_col (mat.col))
    for (int i = 0; i < 3; i++)
      out.c[i] = mat.col[i];
  else if (mat.s.i) {
    double cmap[NCMAP][3];
    scalar s = mat.s;
    mat.map (cmap);
    foreach_dimension()
      a.x += L0/(1<<depth())*1e-5*r.dir.x;
    double val = interpolate (s, a.x, a.y, a.z, mat.linear);
    if (val == nodata)
      return false;
    colormap_pigmentation (out.c, cmap, val, mat.min, mat.max);
  }
  else if (mat.v.y.i) {
    vector v = mat.v;
    foreach_dimension()
      a.x += L0/(1<<depth())*1e-5*r.dir.x;
    double val = interpolate (v.x, a.x, a.y, a.z, mat.linear);
    if (val == nodata)
      return false;
    out.c[0] =  255*min(max(val - mat.min, 0.)/(mat.max - mat.min), 1.);
    val = interpolate (v.y, a.x, a.y, a.z, mat.linear);
    if (val == nodata)
      return false;
    out.c[1] =  255*min(max(val - mat.min, 0.)/(mat.max - mat.min), 1.);
    val = interpolate (v.z, a.x, a.y, a.z, mat.linear);
    if (val == nodata)
      return false;
    out.c[2] =  255*min(max(val - mat.min, 0.)/(mat.max - mat.min), 1.);
  }
  get_pixel (&out);
  for (int i = 0; i < 3; i++)
    px[i] = out.c[i];
  if ((!mat.R && !mat.T) || r.depth > max_ray_depth)
    return true;
  if (mat.R) {
    unsigned char temp[3];
    ray new = {.O = a, .dir = reflect (r.dir, n), .depth = r.depth + 1};
    bias (&new);
    if ((get_color_ray (new, temp) != HUGE))
      for (int i = 0; i < 3; i++)
	px[i] = (1 - mat.R)*px[i] + mat.R*temp[i];
    return true;
  }
  if (mat.T) {
    unsigned char temp[3];
    ray new = {.O = a, .dir = r.dir, .depth = r.depth + 1};
    bias1 (&new);
    if ((get_color_ray (new, temp) != HUGE))
      for (int i = 0; i < 3; i++)
	px[i] = (1 - mat.T)*px[i] + mat.T*temp[i];
    return true;
  }
  //this should not happen:  
  assert (0);
  for (int i = 0; i < 3; i++)
    px[i] = 0; // Darkness
  return false;
}
/**
## User interface

The sketch and store functions are elsewhere:
*/
#include "sketch.h"
#undef normalize
/**
## Tests

* [Test camera angles](bwatch.c)
* [Test every bwatch function](minbwatch.c)

## Usage

* [Visualize a dream](balls.c)
* [Visualize a dump file](visl2.c)
* [Vortex ring colission with 4th order solver](ring4.c)
* [Turbulence at a cloud-atmosphere interface](subsiding-shell.c)
* [Untying a vortex knot](trefoil4.c)

## All pages using `bwatch`

* [Overview of all pages using `bwatch`](http://www.basilisk.fr/_search?patterns=bwatch.h%22)


## to do

  1. ~~~A flexible camera~~~  
  2. ~~~MG-accelerated `foreach_segment_3D()` ray-casting~~~
  3. ~~~Do something ray-ish~~~
  4. ~~~implement and combine more than one sketch function~~~
  5. ~~~Implement a fun sketch function~~~ 
  6. ~~~Consider lights and colours~~~
  7. Reflection of light sources
  8. ~~~volumetric renderingplume.c~~~
  9. ~~~Use bwatch for something~~~
  10. Think about optimizations to mitigate the long rendering times
  11. ~~~Fill gaps in vof facets~~~  
  12. ~~~Refractive objects~~~
  13. Implement point-light sources 
  14. Primitives cast shades
  15. ~~~Partially reflective and transparant (alpha buffering) objects~~~
  16. ~~~Gouraud shading for facets~~~
  17. ~~~Isosurfaces~~~
  18. ~~~RGB-code from vector field~~~ 
  19. ~~~Parallel acceleration with OpenMP~~~
  20. MPI compatiblity (and acceleration)
  */
