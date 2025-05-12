/**
# Sketching functions for Bwatch

Sketching utilities consist of a few functions and definitions.
*/
#include "bwatch-iterators.h"
/**
## Sketch vof facets  
`sketch_vof()` mimics bview's `draw_vof()` a little
*/
struct _sketch_vof {
  scalar f;
  bool precomp; //precompute normals
  material mat;
};

struct _sketch_vof svlist[N_SKETCH];   // Arguments to sketch vof

trace
double get_pixel_sketch_vof (struct _sketch_vof inp, ray r, double DG, unsigned char * px) {
  double distance = HUGE;
  coord n, a;
  if (inp.precomp) {
    foreach_precomp_ray_facet_intersection(r, inp.f, DG) {
      distance = _distg;
      n = _n;
      a = _a[0];
    }
  }
  else {
    foreach_ray_facet_intersection(r, inp.f) {
      distance = _distg;
      n = _n;
      a = _a[0];
    }
  }
  if (distance < HUGE && distance < DG)
    if (get_color_material (r, a, n, inp.mat, px))
      return distance;
  return HUGE;
}

// User function:
bool sketch_vof (struct _sketch_vof inp) {
  sketch_list[functot] = &get_pixel_sketch_vof;
  if (!inp.mat.dull && !inp.mat.ind && !any_col (inp.mat.col) &&
      !inp.mat.s.i && !inp.mat.v.x.i && !(inp.mat.R >= 1.))
    inp.mat.col[0] = inp.mat.col[1] = inp.mat.col[2] = 250;
  if (!inp.mat.map)
    inp.mat.map = jet;
  if (!inp.mat.min && !inp.mat.max) {
    inp.mat.min = -1;
    inp.mat.max =  1;
  }
  if (!inp.mat.SPexp && !inp.mat.SPR) {
    inp.mat.SPexp = 40;
    inp.mat.SPR = 0.1;
  }
  svlist[functot] = inp;
  functot++;
  shading = list_add (shading, inp.f);
  return true;
}
/**
## A Slice

`quadriangles` may mimic bview's `squares()` a bit. 
 */
struct _quadriangles {
  scalar s;
  double alpha;
  coord n;
  material mat;
};

struct _quadriangles qlist[N_SKETCH]; // Arguments to quadriangles

trace
double get_pixel_quadriangles (struct _quadriangles inp, ray r, double DG, unsigned char * px) {
  coord a;
  double d = ray_plane_intersect (r, inp.n, inp.alpha, &a);
  if (d <= 0 || d >= DG)
    return HUGE;
  if (locate (a.x, a.y, a.z).level < 0) // Only slices in the domain
    return HUGE;
  normalize (&inp.n);
  if (get_color_material (r, a, inp.n, inp.mat, px))
    return d;
  return HUGE;
}

// User interface
void quadriangles (struct _quadriangles inp) {
  sketch_list[functot] = &get_pixel_quadriangles;
  if (!any_comp (inp.n))
    inp.n.z = 1;
  if (!inp.mat.dull && !inp.mat.ind && !any_col (inp.mat.col) &&
      !inp.mat.s.i && !inp.mat.v.y.i && !(inp.mat.R >= 1)) {
    if (inp.s.i)
      inp.mat.s = inp.s;
    else 
      inp.mat.col[0] = 250, inp.mat.col[1] = 150, inp.mat.col[2] = 50;
  }
  if (!inp.mat.map)
    inp.mat.map = jet;
  if (!inp.mat.min && !inp.mat.max) {
    inp.mat.min = -L0;
    inp.mat.max =  L0;
  }
  qlist[functot] = inp;
  functot++;
}
/**
## The grid

`lattice()`  may mimic bview's `cells()` a bit.
*/
struct _lattice {
  double alpha;
  coord n;
  double width;
  material mat;
  bool HQ; // High quality, some background object required
};

struct _lattice llist[N_SKETCH]; // Arguments to lattice

double get_pixel_lattice (struct _lattice inp, ray r, double DG, unsigned char * px) {
  coord a;
  double d = ray_plane_intersect (r, inp.n, inp.alpha, &a);
  if (d < 0 || d > 1e4*L0 || d > DG)
    return HUGE;
  double val = HUGE;
  foreach_point (a.x, a.y, a.z, serial) {
    coord cc = {x, y, z};
    //grid-alligned slices only
    foreach_dimension() {
      if (inp.n.z) {
	int sx = sign (a.x - cc.x);
	int sy = sign (a.y - cc.y);
	double dx = cc.x + sx*Delta/2 - a.x;
	double dy = cc.y + sy*Delta/2 - a.y;
	double dist = min (fabs(dx), fabs(dy));
	if (dist < (inp.HQ ? 3 : 1)*inp.width && d > Delta*.1) {
	  normalize (&inp.n);
	  if (inp.HQ) {
	    double b = 0, c = 0;
	    foreach_dimension() {
	      b += sq(cam.O.x - a.x);
	      c += sq(cam.O.x - cam.poi.x);
	    }
	    b = sqrt(b);
	    c = sqrt(c);
	    double ps = cam.fov*b/(cam.nx*c);
	    coord np = {0,0,0}, rp = {0, 0, 0};
	    np.z = 1.;
	    rp.z = r.dir.z;
	    if(dist == fabs(dx)) 
	      rp.y = r.dir.y;
	    else
	      rp.x = r.dir.x;
	    normalize (&rp);
	    double cc = fabs(dot(rp, np));
	    double ss = cc > 0. ? ps/cc : ps;
	    //printf ("%g %g %g\n", ss, ps, cc);
	    double temp = dist + ss/4. < inp.width ? 1 :
	      dist + ss/4. < 2*inp.width ?
	      1 - (dist + ss/4. - inp.width)/(inp.width) : 0;
	    temp += fabs(dist - ss/4.) < inp.width ? 1 :
	      fabs(dist - ss/4.) < 2*inp.width ?
	      1 - (fabs(dist - ss/4.) - inp.width)/(inp.width) : 0;
	    temp = (2. - temp)/2.;
	    inp.mat.T = max(min(temp, 1), 0);
	    foreach_dimension()
	      a.x -= 0.04*Delta*r.dir.x;
	  }
	  if (get_color_material (r, a, inp.n, inp.mat, px)) 
	    val = d - Delta/20.;
	}
      }
    }
  }
  return val;
}
// User interface
void lattice (struct _lattice inp) {
  sketch_list[functot] = &get_pixel_lattice;
  if (!inp.mat.dull && !inp.mat.ind && !any_col (inp.mat.col) &&
      !inp.mat.s.i && !inp.mat.v.x.i && !(inp.mat.R >= 1)) 
    for (int i = 0; i < 3; i++)
      inp.mat.col[i] = 1;
  if (!any_comp(inp.n))
    inp.n.z = 1.;
  if (!inp.width)
    inp.width = cam.fov/(40*N);
  if (inp.HQ)
    inp.width /= 1.5;
  llist[functot] = inp;
  functot++;
}
/**
## A disk

A round planar disk
*/
struct _disk {
  double R;
  coord P;
  coord n;
  material mat;
};

struct _disk dlist[N_SKETCH]; // Arguments to disk

double get_pixel_disk (struct _disk inp, ray r,  double DG, unsigned char * px) {
  coord a;
  double alpha = 0;
  foreach_dimension()
    alpha += inp.n.x*inp.P.x;
  double d = ray_plane_intersect (r, inp.n, alpha, &a);
  if (d >= DG || d < 0.0001) //No further computations
    return HUGE;
  double D = 0;
  foreach_dimension()
    D += sq(a.x - inp.P.x);
  if (D > sq(inp.R))
    return HUGE;
  normalize (&inp.n);
  normalize (&r.dir);
  if (get_color_material (r, a, inp.n, inp.mat, px))
    return d;
  return HUGE;
}

bool disk (struct _disk inp) {
  sketch_list[functot] = &get_pixel_disk;
  if (!inp.R)
    inp.R = L0/2.;
  if (!any_comp (inp.n))
    inp.n.z = 1;
  if (!inp.mat.dull || !inp.mat.ind || !any_col (inp.mat.col) ||
      !inp.mat.s.i || !inp.mat.v.x.i || !(inp.mat.R >= 1.)) 
    inp.mat.R = 1;
  dlist[functot] = inp;
  functot++;
  return true;
}
/**
## Sketch a sphere
*/
struct _sphere {
  double R;      //Radius
  coord C;       //centre
  struct Material mat;
};

struct _sphere slist[N_SKETCH];    // Arguments to sphere

trace
double get_pixel_sphere (struct _sphere inp, ray r, double DG, unsigned char * px) {
  coord a, n;
  double d = ray_sphere_intersect (r, inp.C, inp.R, &a, &n); 
  if (d > 0 && d < DG) { // possibly in sight
    if (get_color_material (r, a, n, inp.mat, px))
      return d;
  }
  return HUGE;
}

// User interface
void sphere (struct _sphere inp) {
  sketch_list[functot] = &get_pixel_sphere;
  if (!inp.R)
    inp.R = 0.5;
  if (!inp.mat.dull && !inp.mat.ind && !any_col (inp.mat.col) &&
      !inp.mat.s.i && !inp.mat.v.y.i && !(inp.mat.R >= 1.)) 
    inp.mat.col[0]= 200, inp.mat.col[1] =  100, inp.mat.col[2] = 150;
  if (inp.mat.dull && !any_col(inp.mat.col)) {
    inp.mat.col[0]  = 50, inp.mat.col[1]  = 100, inp.mat.col[2]  = 150;
    inp.mat.col2[0] = 50, inp.mat.col2[1] = 150, inp.mat.col2[2] = 100;
    inp.mat.n1 = (coord){0.1, 0.9, 0};
  }
  if (!inp.mat.SPexp && !inp.mat.SPR) {
    inp.mat.SPexp = 40;
    inp.mat.SPR = 0.1;
  }
  if (!inp.mat.map)
    inp.mat.map = jet;
  if (!inp.mat.min && !inp.mat.max) {
    inp.mat.min = -L0;
    inp.mat.max = L0;
  }
  slist[functot] = inp;
  functot++;
}
/**
## Equiplane

An (Curved) Plane where a scalar field `s` has the value `val`.  It
can also be used to sketch a smooth vof-interface reconstruction.
 */
struct _equiplane {
  scalar s;
  double val;
  bool insideout;
  bool vof;
  material mat;
} _equiplane;

struct _equiplane eqlist[N_SKETCH];

trace
double get_pixel_equiplane (struct _equiplane inp, ray r, double DG, unsigned char * px) {
  scalar s = inp.s;
  // Draw the distance equiplane for vof
  if (inp.vof) {
    s = s.vofd;
    if (inp.insideout) //Invert 
      inp.insideout = false;
    else
      inp.insideout = true;
  }
  scalar posi = s.possible;
  vector v = s.normals;
  double d = HUGE;
  coord n, a;
  foreach_possible_ray_equi_intersection(r, s, posi, DG) {
    if (_dist < d) {
      double vals[2];
      for (int i = 0; i < 2; i++) 
	vals[i] = interpolate(s, _a[i].x, _a[i].y, _a[i].z);
      if ((vals[0] - inp.val) * (vals[1] - inp.val) < 0.) {
	double w = fabs(vals[0] - inp.val)/fabs((vals[0] - inp.val) - (vals[1] - inp.val));
	if (w < 0 || w > 1)
	  continue;
	double mrdir = 0;
	normalize (&r.dir);
	foreach_dimension() {
	  a.x = w*_a[1].x + (1. - w)*_a[0].x;
	  if (fabs(r.dir.x) > mrdir) {
	    mrdir = fabs(r.dir.x);
	    d = (a.x - r.O.x)/r.dir.x;
	  }
	}
	n.x = interpolate (v.x, a.x, a.y, a.z); 
	n.y = interpolate (v.y, a.x, a.y, a.z); 
	n.z = interpolate (v.z, a.x, a.y, a.z);
      }
    }
  }
  if (any_comp(n) && d < DG && d > 0) {
    if (inp.insideout) {
      foreach_dimension()
	n.x *= -1;
    }
    normalize (&n);
    if (get_color_material (r, a, n, inp.mat, px))
      return d;
  }
  return HUGE;
}

void equiplane (struct _equiplane inp) {
  sketch_list[functot] = &get_pixel_equiplane;
  if (!inp.mat.dull && !inp.mat.ind && !any_col (inp.mat.col) &&
       !inp.mat.s.i && !inp.mat.v.x.i && !(inp.mat.R >= 1.)) 
     inp.mat.col[0]= 100, inp.mat.col[1] =  200, inp.mat.col[2] = 150;
  if (!inp.mat.SPexp && !inp.mat.SPR) {
    inp.mat.SPexp = 40;
    inp.mat.SPR = 0.1;
  }
  if (!inp.mat.map)
    inp.mat.map = jet;
  if (!inp.mat.min && !inp.mat.max) {
    inp.mat.min = -L0;
    inp.mat.max = L0;
  }
  eqlist[functot] = inp;
  functot++;
}
/**
## Image in the scene

Place an image in the scene
 */
#define PIXEL_INDEX (((ix + inp.res*jy)*3))

struct _image {
  char * fname;
  int res;
  unsigned char * imgrgb; //[res*res*3];
  coord n;
  double alpha;
  coord up;
  unsigned char greenscreen[3]; // Transparent color code (non black green screen)
  material mat;
};

struct _image ilist[N_SKETCH]; 

double get_pixel_image (struct _image inp, ray r, double DG, unsigned char * px) {
  coord a;
  double d = ray_plane_intersect (r, inp.n, inp.alpha, &a);
  if (d <= 0 || d >= DG)
    return HUGE;
  if (locate (a.x, a.y, a.z).level < 0) // Only slices in the domain
    return HUGE;
  normalize (&inp.n);
  {;}
  //Only grid aligned projections
  coord Ori = {X0, Y0, Z0};
  foreach_dimension() {
    if (inp.n.z) {
      int ix = max(min(inp.res*(a.x - Ori.x)/L0, inp.res - 1), 0);
      int jy = max(min(inp.res*(Ori.y + L0 - a.y)/L0, inp.res - 1), 0);
      for (int i = 0; i < 3; i++) 
	inp.mat.col[i] = max(inp.imgrgb[PIXEL_INDEX + i], 1);
    }
  }
  if (any_col (inp.greenscreen)) {
    int count = 0;
    for (int i = 0; i < 3; i++) 
      if (inp.greenscreen[i] == inp.mat.col[i])
	count++;
    if (count == 3) //its the green-screen val
      return HUGE;
  }
  if (get_color_material (r, a, inp.n, inp.mat, px))
    return d;
  return HUGE;
}

bool image (struct _image inp) {
  if (!any_comp(inp.n))
    inp.n.z = 1;
  if (!inp.alpha)
    inp.alpha = 0;
  if (!inp.res)
    inp.res = (1 << grid->maxdepth);
  if (inp.fname && inp.imgrgb == NULL) {
    inp.imgrgb = malloc (3*sizeof(unsigned char)*sq(inp.res));
    char cmd[999];
    sprintf (cmd, "convert -quiet %s -interpolate Nearest -filter point -resize %dx%d^  \
                 -gravity center -extent %dx%d -depth 8 2not8exist5.ppm",
	     inp.fname, inp.res, inp.res, inp.res, inp.res);
    assert(!system (cmd));
    FILE * fpi = fopen ("2not8exist5.ppm", "rb");
    char line_one[3]; //three to hold "P6?"
    int height, width, max_color;
    fscanf(fpi, "%s\n%d %d\n%d\n", line_one, &height, &width, &max_color);
    //printf ("%s %d %d %d\n", line_one, height, width, max_color);
    fread (inp.imgrgb, 3*sizeof(unsigned char), sq(inp.res), fpi);
    fclose (fpi);
  } else
    assert (0);
  ilist[functot] = inp;
  sketch_list[functot] = &get_pixel_image;
  functot++;
  return true;
}
/**
### Display text

It is possible to add text using `convert` and `image`
 */
struct _sketch_text {
  char * str;
  int fs;              // Font size (compare against .res)
  unsigned char tc[3]; // text color
  char * pos;          // Imagemagick Position (center, north, south west etc..)
  int res;             // Text image resolution
  char * ops;          // Options to pass to convert
  coord n;
  double alpha;
  coord up;
  material mat;
};
  
bool sketch_text (struct _sketch_text inp) {
  char fname[99] = "3not9exist6.ppm";
  char cmd[999];
  if (!inp.fs) 
    inp.fs = 100;
  if (!inp.res)
    inp.res = 1000;
  unsigned char tc[3] = {0,0,0}; //black
  char d_ops[2] = " ";
  if (!inp.ops)
    inp.ops = d_ops;
  if (any_col(inp.tc)) {
    tc[0] = inp.tc[0];
    tc[1] = inp.tc[1];
    tc[2] = inp.tc[2];
  }
  unsigned char gs[3]; //Set Green screen color
  for (int i = 0; i < 3; i++) {
    int tci = tc[i];
    gs[i] = tci - 2 <= 0 ? tc[i] + 2 : tc[i] - 2;
  }
  char tmp[99] = "northwest";
  if (!inp.pos)
    inp.pos = tmp;
  sprintf (cmd, "convert -quiet %s -size %dx%d -background \"rgb(%d, %d, %d)\"\
	   -fill \"rgb(%d,%d,%d)\" -pointsize %d -gravity %s caption:'%s' -depth 8 %s",
	   inp.ops, inp.res, inp.res, gs[0], gs[1], gs[2], tc[0], tc[1], tc[2],
	   inp.fs, inp.pos, inp.str, fname);
  //printf ("%s \n", cmd);
  system (cmd);
  struct _image inpn;
  inpn.fname = fname; 
  inpn.greenscreen[0] = gs[0];
  inpn.greenscreen[1] = gs[1];
  inpn.greenscreen[2] = gs[2];
  inpn.res = inp.res;
  inpn.n = inp.n;
  inpn.alpha = inp.alpha;
  inpn.up = inp.up;
  inpn.mat = inp.mat;
  inpn.imgrgb = NULL;
  return image (inpn);
}
/**
## A volumetric object

It is possible to sketch volumetric objects as a smoke concentration
field. (`0 ~ mval < s[]`) 
*/
#include "radix.h" //A sorting algorithm
bool sketch_a_volume = false;
struct _volume {
  scalar s;
  double sc;            // 1/e scale for s times length
  double mval;          // Minimal value (positive and close to zero)
  unsigned char col[3]; // prescribed Color
  bool cols;            // Color by scalar value (color map)
  vector colorv;        // Color map 3D vector direction xyz -> rgb
  Colormap map;         // color Map (default is cool_warm)
  double min;           // min cbar
  double max;           // max cbar
  double shading;       // shading effect for default lights
};

struct _volume vollist; //max 1

double get_attenuation (coord a, struct _volume inp) {
  ray r = {.O = a, .dir = lights[1].dir};
  foreach_dimension()
    r.dir.x *= -1;
  normalize (&r.dir);
  scalar s = inp.s;
  double sl = 0;
  foreach_ray_cell_intersection_volume (r, HUGE, s.possible) {
    coord mean = {0, 0, 0};
    double len = 0;
    foreach_dimension() {
      mean.x = (_a[0].x + _a[1].x)/2.;
      len += sq(_a[0].x - _a[1].x);
    }
    double val = interpolate (inp.s, mean.x, mean.y, mean.z);
    if (val > inp.mval) {
      len = sqrt(len);
      if (len > 0)
	sl += fabs(val)*len;
    }
  }
  return sl;
}

trace
void mod_pixel_volume (struct _volume inp, ray r, double DG, unsigned char * px) {
  scalar s = inp.s;
  int np = 50, incf = 2; // Array size, increase factor;
  double * dvl = NULL;
  int elm = 0;
  if (inp.cols) {
    elm = 3 + (inp.shading > 0) + 3*(inp.colorv.x.i > 0); 
    dvl = malloc (elm*sizeof(double)*np); //Distance, values and length;
  }
  int celld = 0; // Number if cells
  double sl = 0; //integral of s times length;
  foreach_ray_cell_intersection_volume(r, DG, s.possible) {
    coord mean = {0, 0, 0};
    double len = 0;
    foreach_dimension() {
      mean.x = (_a[0].x + _a[1].x)/2.;
      len += sq(_a[0].x - _a[1].x);
    }
    double val = interpolate (inp.s, mean.x, mean.y, mean.z);
    if (val > inp.mval) {
      len = sqrt(len);
      if (len > 0)
	sl += fabs(val)*len;
      if (inp.cols) { // Record data 
	celld++;
	if (celld >= np) { // increase array size
	  np *= incf;
	  dvl = realloc (dvl, elm*sizeof(double)*np);
	}
	double dirm = 0;
	foreach_dimension() { //distance
	  if (fabs(r.dir.x) > dirm)
	    dvl[elm*(celld - 1)] = (mean.x - r.O.x)/r.dir.x;
	}
	assert (dvl[elm*(celld - 1)] > 0);
	dvl[elm*(celld - 1) + 1] = val;
	dvl[elm*(celld - 1) + 2] = len;
	if (inp.shading)
	  dvl[elm*(celld - 1) + 3] = get_attenuation (mean, inp);
	if (inp.colorv.x.i) {
	  double xp = mean.x, yp = mean.y, zp = mean.z;
	  int dim = 0;
	  foreach_dimension()
	    dvl[elm*(celld - 1) + 4 + dim++] = interpolate (inp.colorv.x, xp, yp, zp);
	}
      }
    }
  }
  if (inp.cols && sl > 0) { // Compute effective absorbed color
    int ind[celld], dista[celld];
    for (int it = 0; it < celld; it++) 
      dista[it] = (int)(1 << (depth() + 1))*(dvl[elm*it]/L0); // distance in integer units
    radixsort (celld, dista, ind);
    double col[3] = {0, 0, 0};
    double cmap[NCMAP][3];
    inp.map (cmap);
    double sl2 = 0;
    double TW = 0;
    for (int it = 0; it < celld; it++) {
      unsigned char coll[3];
      if (inp.colorv.x.i) {
	int mx = 0;
	for (int i = 0; i < 3; i++) {
	  coll[i] = 255*min(1, fabs(dvl[elm*ind[it] + 4 + i])/inp.max); 
	  if (coll[i] > mx)
	    mx = coll[i];
	}
	for (int i = 0; i < 3; i++) 
	  coll[i] = mx > 0 ? min (255, (255*coll[i])/mx) : 1; 
      } else 
	colormap_pigmentation (coll, cmap, dvl[elm*ind[it] + 1], inp.min, inp.max);
      double sle = dvl[elm*ind[it] + 2]*fabs(dvl[elm*ind[it] + 1]);
      double Weight = exp(-sl2/inp.sc) - exp(-(sl2 + sle)/inp.sc);
      double fac = 1;
      if (inp.shading) { 
	double I = lights[0].I;
	I += lights[1].I*exp(-dvl[elm*ind[it] + 3]/(smoke.att));
	fac =  min (max(cam.f(I) - cam.f(cam.min),0)/cam.f(cam.max/cam.min), 1);
      }
      TW += Weight;
      for (int i = 0; i < 3; i++) 
	col[i] += (double)coll[i]*Weight*fac;
      sl2 += sle;
    }
    if (TW > 0) {
      for (int i = 0; i < 3; i++) 
	inp.col[i] = (unsigned char)(col[i]/TW);
    }
  }
  if (sl > 0) {
    double W = exp(-sl/inp.sc);
    px[0] = (1 - W)*inp.col[0] + W*px[0];
    px[1] = (1 - W)*inp.col[1] + W*px[1];
    px[2] = (1 - W)*inp.col[2] + W*px[2];
  }
  free (dvl); dvl = NULL;
}

bool volume (struct _volume inp) {
  if (!inp.sc)
    inp.sc = 1.;
  if (!any_col(inp.col)) {
    inp.col[0] = 255;
    inp.col[1] = 255;
    inp.col[2] = 255;
  }
  if (!inp.max)
    inp.max = 1;
  if (!inp.mval)
    inp.mval = 0;
  if (!inp.min)
    inp.min = inp.mval;
  if (inp.shading) {
    scalar s = inp.s;
    smoke = s;
    s.att = inp.sc/inp.shading;
  }
  if (!inp.map)
    inp.map = cool_warm;
  vollist = inp;
  sketch_a_volume = true;
  return sketch_a_volume;
}

/**
## Triangles

Draw a list of triangles, stored as a `nodata` terminated list of
vertex triplets (e.g. as read by `inport_stl`).
 */
struct _triangles {
  coord * c;
  coord bb[2];  //Bounding box (automatically computed)
  material mat;
};

struct _triangles tlist[N_SKETCH];

trace
double get_pixel_triangles (struct _triangles inp, ray r, double DG, unsigned char * px) {
  coord * p = inp.c;
  double dist = DG;
  coord a, n;
  bool any = false;
  if (ray_box (r, inp.bb) < dist) {
    while (p->x != nodata) {
      double d = ray_triangle (r, p);
      if (d < dist) {
	dist = d;
	any = true;
	foreach_dimension()
	  a.x = r.O.x + r.dir.x*d;
	coord a1, a2;
	foreach_dimension() {
	  a1.x = p[1].x - p[0].x;
	  a2.x = p[2].x - p[0].x;
	}
	n = cross (a1, a2);
	normalize (&n);
      }
      p += dimension; //3?
    }
    if (any) 
      if (get_color_material (r, a, n, inp.mat, px))
	return dist;
  }
  return HUGE;
}

struct _triangles get_bb (struct _triangles in) {
  struct _triangles inp = in;
  foreach_dimension() {
    inp.bb[0].x = HUGE;
    inp.bb[1].x = -HUGE;
  }
  while (in.c->x < nodata) {
    foreach_dimension() {
      if (in.c->x < inp.bb[0].x)
	inp.bb[0].x = in.c->x;
      if (in.c->x > inp.bb[1].x)
	inp.bb[1].x = in.c->x;
    }
    in.c++;
  }  
  return inp; 
}

bool triangles (struct _triangles inp) {
  sketch_list[functot] = &get_pixel_triangles;
  if (!inp.mat.dull && !inp.mat.ind && !any_col (inp.mat.col) &&
      !inp.mat.s.i && !inp.mat.v.x.i && !(inp.mat.R >= 1.)) 
    inp.mat.col[0]= 200, inp.mat.col[1] =  150, inp.mat.col[2] = 100;
  if (!inp.mat.SPexp && !inp.mat.SPR) {
    inp.mat.SPexp = 40;
    inp.mat.SPR = 0.1;
  }
  if (!inp.mat.map)
    inp.mat.map = jet;
  if (!inp.mat.min && !inp.mat.max) {
    inp.mat.min = -L0;
    inp.mat.max = L0;
  }
  inp = get_bb (inp);
  tlist[functot] = inp;
  functot++;
  return true;
}
/**
   ## A plain sheet for sketching...
   
`plain()` mimics bview's `clear()` a bit.
*/
void plain (void) {
  for (int func = 0; func < functot; func++) 
    if (sketch_list[func] == &get_pixel_image) { 
      free (ilist[func].imgrgb);
      ilist[func].imgrgb = NULL;
    }
  functot = 0;
  free (shading);
  shading = NULL;
  sketch_a_volume = false;
}

/**
## *The* ray-caster function

The implementation of a helper function that gets the `rgb` values of
a ray by cycling though all sketch function calls...*/
double get_color_ray (ray r, unsigned char * px) {
  double Distg = HUGE;
  for (int func = 0; func < functot; func++) {
    double dist = HUGE;
    unsigned char pt[3];
    if (sketch_list[func] == &get_pixel_quadriangles) 
      dist = sketch_list[func](qlist[func],  r, Distg, pt);
    else if (sketch_list[func] == &get_pixel_sketch_vof) 
      dist = sketch_list[func](svlist[func], r, Distg, pt);
    else if (sketch_list[func] == &get_pixel_disk) 
      dist = sketch_list[func](dlist[func],  r, Distg, pt);
    else if (sketch_list[func] == &get_pixel_lattice) 
      dist = sketch_list[func](llist[func],  r, Distg, pt);
    else if (sketch_list[func] == &get_pixel_sphere) 
      dist = sketch_list[func](slist[func],  r, Distg, pt);
    else if (sketch_list[func] == &get_pixel_equiplane) 
      dist = sketch_list[func](eqlist[func], r, Distg, pt);
    else if (sketch_list[func] == &get_pixel_image) 
      dist = sketch_list[func](ilist[func],  r, Distg, pt);
    else if (sketch_list[func] == &get_pixel_triangles) 
      dist = sketch_list[func](tlist[func],  r, Distg, pt);
    if (dist < Distg && dist > 0) {
      Distg = dist;
      for (int i = 0; i < 3; i++)
	px[i] = pt[i];
    }
  }
  if (sketch_a_volume) // Modify this pixel color 
    mod_pixel_volume (vollist, r, Distg, px);
  return Distg;
}

static inline double point_value (scalar s, Point point) {
  return s[];
}
/**
## Watch the sketch process

You may watch the sketching process *live* in a (minimalistic) X11
window. You should have the [`libX11`](https://www.x.org/wiki/)
library [installed](https://packages.debian.org/nl/sid/libx11-dev)
and compile using something like:

~~~literatec
$ qcc -O2 -DWATCH_ALONG bwatch.c -lm -lX11
$ ./a.out
~~~

Note that it may slow the rendering, especially with openmp.

 */ 
#if WATCH_ALONG
#include <X11/Xlib.h> 
GC gc_2;
Window win = (Window)NULL;
Display *dsp;
int screen;
Colormap cmap;

void XSetup (void) {
  if (win == (Window)NULL) {
#if _OPENMP
    XInitThreads();
#pragma omp barrier
#endif
    dsp    = XOpenDisplay(NULL);
    screen = DefaultScreen(dsp);
    win = XCreateWindow(dsp, DefaultRootWindow(dsp),
			0, 0, cam.nx, cam.ny,
			0, 0, 0, 0, 0, 0);
    XGCValues gcvalues_2;
    gcvalues_2.function = GXcopy;
    gcvalues_2.plane_mask = AllPlanes;
    gcvalues_2.foreground = 0x00FF00;
    gcvalues_2.background = 0xFFFFFF;
    gc_2 = XCreateGC(dsp, win,
		     GCFunction|GCPlaneMask|GCForeground|GCBackground,
		     &gcvalues_2);
    XStoreName (dsp, win, "Bwatch: Watch Along\0");
    XMapWindow(dsp, win);
    cmap= XDefaultColormap(dsp, screen);
   }
}

void XSketchPx (int x, int y, unsigned char px[3]) {
  XLockDisplay (dsp);
  XColor   pX11;
  pX11.red = 256*px[0];
  pX11.green = 256*px[1];;
  pX11.blue = 256*px[2];;        
  if (XAllocColor(dsp,cmap,&pX11)==0)
    assert (0);
  XSetForeground(dsp,gc_2,pX11.pixel);
  XDrawPoint(dsp, win, gc_2, x, cam.ny - y );
  if (tid() == 0)
    XFlush (dsp);
  XUnlockDisplay (dsp);
}

void XStop (void) {
  XDestroyWindow(dsp, win);
  XCloseDisplay(dsp);
  win = (Window)NULL;
}

event stop_X11 (t = end) {
  XStop();
}
#endif

/**
## Preparation 

Before rendering, it maybe usefull to compute some helper fields to
store data for rendering or help guide optimized grid iterators.
 */

void preparation_steps (void){
#if WATCH_ALONG
  XSetup();
#endif
  for (int func = 0; func < functot; func++) {
    //Gradients and possibles for isosurfaces,
    // maybe compute the nearby distances to vof fields.
    if (sketch_list[func] == &get_pixel_equiplane) {
      scalar s = eqlist[func].s;
      if (eqlist[func].vof) {
	char dname[99];
      	sprintf (dname, "dfor%s%d", s.name, func); 
      	scalar d = new_scalar (dname);
	s.vofd = d;
      	find_nearby_distances (s);
	s = d; // we will draw the isosurface of d not s
      }
      char vname[99], pname[99];
      sprintf (vname, "vfor%s%d", s.name, func);
      sprintf (pname, "pfor%s%d", s.name, func);
      vector v = new_vector (vname);
      scalar p = new_scalar (pname);
      s.normals = v;
      s.possible = p;
      find_possible (s, eqlist[func].val);
      gradients ({s}, {v});
      foreach() {
	coord n = {0};
	foreach_dimension()
	  n.x = v.x[];
	if (any_comp(n)) {
	  normalize (&n);
	  foreach_dimension()
	    v.x[] = n.x;
	}
      }
      boundary ((scalar*){v});
    }
    // Precompute normals and plane alphas for vof facets
    else if (sketch_list[func] == &get_pixel_sketch_vof) {
      if (svlist[func].precomp) {
	scalar f = svlist[func].f;
	char vname[99], aname[99];
	sprintf (vname, "vfor%s%d", f.name, func);
	sprintf (aname, "afor%s%d", f.name, func);
	vector v = new_vector (vname);
	scalar alpha = new_scalar (aname);
	f.normals = v;
	f.possible = alpha; //alpha
	foreach() {
	  alpha[] = nodata;
	  foreach_dimension()
	    v.x[] = nodata;
	  if (interfacial2 (f, point)) {
	    coord n = mycs (point, f); //pointing outwards.
	    double l = vec_length(n);
	    alpha[] = plane_alpha (point_value(f, point), n)/l;
	    normalize (&n);
	    foreach_dimension()
	      v.x[] = n.x;
	  }
	}
	boundary ({alpha, v});
      }
    }
  }
  if (sketch_a_volume) { //Pre compute cells with volume
    scalar s = vollist.s;
    char pname [99];
    sprintf (pname, "p_for_%s_volume", s.name);
    scalar p = new_scalar (pname);
    s.possible = p;
    find_possible_vol (s, vollist.mval);
  }
}

void aftercare (void) {
  for (int func = 0; func < functot; func++) {
    if (sketch_list[func] == &get_pixel_equiplane) {
      scalar s = eqlist[func].s;
      if (eqlist[func].vof)
	s = s.vofd;
      vector v = s.normals;
      scalar posi = s.possible;
      delete ({posi, v});
      if (eqlist[func].vof) 
	delete ({s});
    }
    else if (sketch_list[func] == &get_pixel_sketch_vof) 
      if (svlist[func].precomp) {
	scalar s = svlist[func].f;
	vector v = s.normals;
	scalar alpha = s.possible;
	delete ((scalar*){alpha, v});
      }
  }
  if (sketch_a_volume) {
    smoke.i = 0;
    scalar s = vollist.s;
    s = s.possible;
    delete ({s});
  }
}
/**
## Store

Write a `.ppm` to the `FILE * fp`. 
*/

// Default background
#define BGR (0)
#define BGG (25.*(ii + jj)/cam.nx + 120)
#define BGB (50.*jj/cam.nx + 120)

trace
bool store (FILE * fp) {
  if (lights[0].type == 0) // No lights
    default_lights();
  fprintf (fp, "P6\n%d %d\n%d\n", cam.nx, cam.ny, 255); // .ppm header
  preparation_steps();
  
  foreach_ray() {
    unsigned char px[3] = {BGR, BGG, BGB};
    get_color_ray (_r, px);
    fwrite (px, 1, 3, fp);
#if WATCH_ALONG
    XSketchPx (ii, jj, px);
#endif
  }
  fflush (fp);
  aftercare();
  return true;
}

#if _OPENMP
trace
bool store_OPENMP (FILE * fp) {
  if (lights[0].type == 0) // No lights
    default_lights();
  fprintf (fp, "P6\n%d %d\n%d\n", cam.nx, cam.ny, 255); // .ppm header
  preparation_steps();
#pragma omp parallel shared(fp)
  {
    assert (!(cam.nx*cam.ny % npe()));
    unsigned char px[3];
    foreach_ray() {
      int jjj = cam.ny - jj;
      if ((jjj*cam.nx + ii) % npe() == tid()) {
	px[0]= BGR, px[1] = BGG, px[2] = BGB;
	get_color_ray (_r, px);
#if WATCH_ALONG
	XSketchPx (ii, jj, px);
#endif
      }
      if ((jjj*cam.nx + ii) % npe() == (npe() - 1)) 
	for (int i = 0; i < npe(); i++) {
	  if (tid() == i) {
	    fwrite (px, 1, 3, fp);
	    
	  }
#pragma omp barrier
	}
    }
  }  
  fflush (fp);
  aftercare();
  return true;
}
// Overload calls to store
#define store store_OPENMP
#endif //_OPENMP

/**
## Store adaptive

With [`raycaster`](raycaster.c) (or [this one](raycasterv.c))
installed, we can adaptively sample the scene (lossy compresion). It
can be effecient to generate high resolution movies.
 */
# if 0
#include <sys/wait.h>

#define ParentRead      read_pipe[0]
#define ParentWrite     write_pipe[1]
#define ChildRead       write_pipe[0]
#define ChildWrite      read_pipe[1]
#define BGRa (0)
#define BGGa (25.*(b[0] + b[1]) + 120)
#define BGBa (50.*b[1] + 120)

struct _sa {
  FILE * fp;  
  int ml;     // max level of output image (cam.nx/ny are ignored)
  double tol; // tolerance
};

trace
bool store_adaptive (struct _sa inp) {
  double tol = 25;
  int ml = 10;
  if (inp.tol)
    tol = inp.tol;
  if (inp.ml)
    ml = inp.ml;
  if (lights[0].type == 0) // No lights
    default_lights();
  preparation_steps();
  // Organize pipe and fork
  int read_pipe[2];
  int write_pipe[2];  
  pipe (read_pipe);
  pipe (write_pipe);
  pid_t p;
  p = fork();
  if (p == 0) { //child calls `raycaster`
    close (ParentRead);
    close (ParentWrite);
    if (system("which ./raycaster > /dev/null 2>&1")) {
      fprintf (stderr, "Error: Adaptive raycaster is not installed\n");
      assert(0);
    } else {
      dup2 (ChildRead, STDIN_FILENO);
      dup2 (ChildWrite, STDOUT_FILENO);
      dup2 (fileno(inp.fp), STDERR_FILENO);
      char m1[9], m2[15];
      sprintf (m1, "%d", ml);
      sprintf (m2, "%g", tol);
      execl ("raycaster", "raycaster", m1, m2, NULL);
    }
  } else { // Parent computes color codes
    close (ChildRead);
    close (ChildWrite);
    int a = sq(1 << (ml - 2)) + 1;
    int al = a;
    double * pos = malloc (2*al*sizeof(double));
    unsigned char * pxs = malloc (3*al*sizeof(unsigned char));
    while (a) {
      read (ParentRead, &a, sizeof(int));
      if (a > al) {
	al = 2*a;
	pos = realloc (pos, 2*al*sizeof(double));
	pxs = realloc (pxs, 3*al*sizeof(unsigned char));
      }
      for (int i = 0; i < a; i++)
      	read (ParentRead, &pos[2*i], 2*sizeof(double));
      coord cr, proj, vert, hori;
      foreach_dimension() 
	cr.x = cam.poi.x - cam.O.x;
      normalize (&cr);
      foreach_dimension() {
	proj.x = cam.up.x*cr.x;
	vert.x = cam.up.x - proj.x;
      }
      hori = cross (cr, vert);
      normalize (&vert);
      normalize (&hori);
#if _OPENMP
#pragma omp parallel
#endif
      for (int i = tid(); i < a; i += npe()) {
	coord apoint;
	double b[2];
	b[0] = pos[2*(i)];
	b[1] = pos[2*(i) + 1];
	foreach_dimension() 
	  apoint.x = (cam.poi.x
		      + cam.fov*(b[0] - 0.5)*hori.x
		      + cam.fov*(b[1] - 0.5)*vert.x);  
	// New ray
	ray _r;
	_r.depth = 0;
	_r.O = cam.O;
	foreach_dimension()
	  _r.dir.x = apoint.x - cam.O.x;
	normalize (&_r.dir);
	// Get and write color
	pxs[3*i] = BGRa;
	pxs[3*i + 1] = BGGa;
	pxs[3*i + 2] = BGBa;
	get_color_ray (_r, &pxs[3*i]);
      }
      write (ParentWrite, &pxs[0], 3*a*sizeof(unsigned char));
    }
    wait (NULL); //wait until the child is done writing the  image
    free (pos); pos = NULL;
    free (pxs); pxs = NULL;
  }
  aftercare();
  return true;
}
#endif
