
/**
# Isolines in 3D

In 3D, `isoline()` draws isolines on the plane,
 
$$\mathbf{n_p} \cdot \mathbf{x} = \alpha$$

Furthermore, `cross_section()` draw the cross section of a 3D
interface with a plane.
*/

struct _cross {
  char * f;   //fraction field
  char * fs;  // face fraction
  coord np;     //normal vector to plane for 3D
  double alpha; //define plane for 3D
  
  // all fields below must be identical to struct _draw_vof 
  char * c;
  char * s;
  bool edges;
  double larger;
  int filled;

  char * color;
  double min, max, spread;  // min and max icw n
  bool linear;
  Colormap map;
  float fc[3], lc[3], lw;
};

/**
## Helper functions:

A line with parameterization; 

$$\mathbf{x}(t) = \mathbf{P} + t \mathbf{n_i},$$

could cross through a cell. We find the corresponding $t$ values for
the possible points of entry and exit.
 */
void get_t (coord P, coord ni, double ts[2]) {
  coord t_in, t_out;
  foreach_dimension() {
    if (fabs(ni.x) > 1e-6) {
      t_in.x = (-P.x - sign(ni.x)/2.)/ni.x;
      t_out.x = t_in.x + fabs(1/ni.x);
    } else {
      t_in.x = -HUGE; 
      t_out.x = HUGE;
    }
  }
  ts[0] = max(t_in.x,  max(t_in.y,  t_in.z ));
  ts[1] = min(t_out.x, min(t_out.y, t_out.z));
}
/**
   A function that draws the parametric line between $t =$`t[0]` and
   $t=$`t[1]`.
*/
void plot_param_line (bview * view, coord P, double ts[2],
		      coord ni, Point point) {
  glBegin (GL_LINES);
  for (int i = 0; i < 2 ; i++)
    glvertex3d (view, (P.x + ts[i]*ni.x)*Delta + x,
		(P.y + ts[i]*ni.y)*Delta + y,
		(P.z + ts[i]*ni.z)*Delta + z);
  glEnd ();
}

bool cross_section (struct _cross p) {
  assert (dimension == 3); //3D only
  scalar f = lookup_field (p.f);
  face vector fs;
  if (p.fs) {
    struct { char x, y, z; } index = {'x', 'y', 'z'};
    foreach_dimension() {
      char name[80];
      sprintf (name, "%s.%c", p.fs, index.x);
      fs.x = lookup_field (name);
    }
  } else
    fs.x.i = 0;
  float alphap = 0;
  float lw = 2., lc[3] = {1., 1., 1.};
  coord np = {0., 0., 1.};
  if (p.alpha)
    alphap = p.alpha;
  if (p.lw)
    lw = p.lw;
  if (p.lc[0] || p.lc[1] || p.lc[2]) 
    lc[0] = p.lc[0], lc[1] = p.lc[1], lc[2] = p.lc[2]; 
  if (p.np.x || p.np.y ||  p.np.x)
    np.x = p.np.x, np.y = p.np.y, np.z = p.np.z;
  bview * view = draw();
  draw_lines (view, lc, lw) {
    foreach_visible_plane (view, np, alphap) {
      if (f[] > 1e-6 && f[] < (1. - 1e-6)) {
	/**
	 Next, the parameterizations of the planes, and the
	 intersecting line vector (`ni`) are computed. see:
	 [geom-algorithms'
	 website](http://geomalgorithms.com/a05-_intersect-1.html).
	*/
	coord nf = {0};
	if (fs.x.i)
	  nf  = facet_normal (point, f, fs);
	else
	  nf = interface_normal (point, f);
	double alphaf = -plane_alpha (f[], nf);
	double alphan = -(alphap - np.x*x - np.y*y - np.z*z)/Delta;
	coord ni;
	foreach_dimension() //cross product
	  ni.x = nf.y*np.z - nf.z*np.y;
	coord P;
	/**
	   We select the most robust Cartesian direction to describe
	   the line,
	*/
	double max_comp = 0.;
	foreach_dimension() 
	  if (fabs(ni.x) > max_comp)
	    max_comp = fabs(ni.x);
	foreach_dimension() {
	  if (fabs(ni.z) == max_comp) {
	    /**
	       and follow the recipy from the aforementioned
	       geometric algorithm to find a point on the line.
	    */
	    double det = np.x*nf.y - nf.x*np.y;
	    P.x = (np.y*alphaf - nf.y*alphan)/det;
	    P.y = (alphan*nf.x - alphaf*np.x)/det; 
	  P.z = 0.;
	  /**
	     Finally, the intersection is infinitely long, and we
	     only wish to plot the section inside the current
	     cell.
	  */
	  double ts[2];
	  get_t (P, ni, ts);
	  if (ts[0] < ts[1])  //Intersection must be within cell
	    plot_param_line (view, P, ts, ni, point);
	  } // Favorite direction
	} // Each direction
      } // Cells on isosurface
    } // Cells close to plane
  } //draw lines loop
  return true;
}

struct _isoline2 {
  char * phi;   //Scalar field
  double val;   //Value of scalar field on isoline
  int n;        //Number of lines (exclusive with val)
  coord np;     //normal vector to plane for 3D
  double alpha; //define plane for 3D
  
  // all fields below must be identical to struct _draw_vof above
  char * c;
  char * s;
  bool edges;
  double larger;
  int filled;

  char * color;
  double min, max, spread;  // min and max icw n
  bool linear;
  Colormap map;
  float fc[3], lc[3], lw;
};
/**
## The userfunction

It is based on the existing `isoline()` function.
*/
trace
bool isoline2 (struct _isoline2 p) {
  /**
In 2D, `draw_vof()` is useful.
   */
#if dimension == 2
  if (!p.color) p.color = p.phi;
  colorize_args (p);
  scalar phi = col, fiso[];
  face vector siso[];
  p.c = "fiso", p.s = "siso";
  struct _draw_vof a = *((struct _draw_vof *)&p.c);
  if (p.n < 2) {
    fractions (phi, fiso, siso, p.val);
    draw_vof (a);
  }
  else if (p.max > p.min) {
    double dv = (p.max - p.min)/(p.n - 1);
    for (p.val = p.min; p.val <= p.max; p.val += dv) {
      fractions (phi, fiso, siso, p.val);
      draw_vof (a);
    }
  }
#else // dimension == 3
  /**
     In 3D, we compute the intersection of an isosurface of $\phi$
     and the plane. Frist we set some default values and overload them
     if they are provided by the user.
   */
  scalar s = lookup_field (p.phi), f[];
  face vector fs[];
  double alphap = 0., val = 0., dv = 1.;
  float lw = 2., lc[3] = {1., 1., 1.};
  int n = 1;
  coord np = {0., 0., 1.};
  if (p.alpha)
    alphap = p.alpha;
  if (p.val)
    val = p.val;
  if (p.lw)
    lw = p.lw;
  if (p.lc[0] || p.lc[1] || p.lc[2]) 
    lc[0] = p.lc[0], lc[1] = p.lc[1], lc[2] = p.lc[2]; 
  if (p.np.x || p.np.y ||  p.np.x)
    np.x = p.np.x, np.y = p.np.y, np.z = p.np.z;
  if (p.n)
    n = p.n;
  if (n > 1 && !(p.min < p.max)) {//A guess
    p.min = statsf(s).min;
    p.max = statsf(s).max;
    dv = (p.max - p.min)/((double)n - 1.);
  }
  /**
     Second, the scalar field is transformed into a vertex field for
     the computation of the isosurface' facets.
   */
  vertex scalar phif[];
  foreach_vertex()
    phif[] = (s[] + s[-1] +  + s[0,-1,-1] + s[-1,-1,-1] +
	      s[0,-1] + s[0,0,-1] + s[-1,0,-1] + s[-1,-1])/8.;
  for (int j = 0 ; j < n ; j++) {
    if (n > 1)
      val = p.min + (double)j*dv;
    fractions (phif, f, fs, val); //Volume and face fractions
    cross_section ("f", "fs", np, alphap, lw = lw, lc = {lc[0],lc[1],lc[2]});
  } // Draw n cross sections
#endif // dimension == 3
  return true;
}
