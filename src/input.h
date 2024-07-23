/**
# *input_pgm()*: Importing Portable Gray Map (PGM) images

This function reads in the
[PGM](http://en.wikipedia.org/wiki/Netpbm_format) file in file *fp*
and imports the corresponding values into field *s*. The grayscale is
normalised and inverted so that the maximum value in field *s* is one
(black) and the minimum value is zero (white). 

By default the origin of the image (lower-left corner) is assumed to
be at (0,0) and the width of the image is set to `L0`. This can be
changed using the optional *ox*, *oy* and *width* parameters. */

#include "utils.h"

void input_pgm (scalar s, FILE * fp,
		double ox = 0., double oy = 0., double width = L0)
{
  char line[81];
  if (!fgets (line, 81, fp)) {
    fprintf (stderr, "input_pgm: could not read magic number\n");
    exit (1);
  }
  if (strcmp (line, "P2\n") && strcmp (line, "P5\n")) {
    fprintf (stderr, "input_pgm: magic number '%s' does not match PGM\n", 
	     line);
    exit (1);
  }
  int binary = !strcmp (line, "P5\n");
  if (!fgets (line, 81, fp)) {
    fprintf (stderr, "input_pgm: could not read width and height\n");
    exit (1);
  }
  int W, H;
  while (line[0] == '#' && fgets (line, 81, fp));
  if (line[0] == '#' || sscanf (line, "%d %d", &W, &H) != 2) {
    fprintf (stderr, "input_pgm: could not read width and height\n");
    exit (1);
  }
  if (!fgets (line, 81, fp)) {
    fprintf (stderr, "input_pgm: could not read maxval\n");
    exit (1);
  }
  int maxval;
  if (sscanf (line, "%d", &maxval) != 1) {
    fprintf (stderr, "input_pgm: could not read maxval\n");
    exit (1);
  }
  if (maxval < 256) {
    unsigned char * a = qmalloc (W*H, unsigned char);
    size_t n = 0;
    if (binary)
      n = fread (a, 1, W*H, fp);
    else {
      int v;
      while (n < W*H && fscanf (fp, "%d ", &v) == 1)
	a[n++] = v;
    }
    if (n != W*H) {
      fprintf (stderr, "input_pgm: read only %ld values\n", n);
      exit (1);
    }
    foreach() {
      int i = (x - ox)*W/width, j = (y - oy)*W/width;
      if (i >= 0 && i < W && j >= 0 && j < H)
	s[] = 1. - a[(H - 1 - j)*W + i]/(double)maxval;
      else
	s[] = 0.;
    }
    free (a);
  }
  else {
    unsigned short * a = qmalloc (W*H, unsigned short);
    size_t n = 0;
    if (binary)
      n = fread (a, 2, W*H, fp);
    else {
      int v;
      while (n < W*H && fscanf (fp, "%d ", &v) == 1)
	a[n++] = v;
    }
    if (n != W*H) {
      fprintf (stderr, "input_pgm: read only %ld values\n", n);
      exit (1);
    }
    foreach() {
      int i = (x - ox)*W/width, j = (y - oy)*W/width;
      if (i >= 0 && i < W && j >= 0 && j < H)
	s[] = 1. - a[(H - 1 - j)*W + i]/(double)maxval;
      else
	s[] = 0.;
    }
    free (a);
  }
}

static void next_char (FILE * fp, int target)
{
  int c = fgetc(fp), para = 0;
  while (c != EOF && (c != target || para > 0)) {
    if (c == '{') para++;
    if (c == '}') para--;
    c = fgetc(fp);
  }
  if (c != target) {
    fprintf (stderr, "input_gfs(): error: expecting '%c'\n", target);
    exit (1);
  }
}

static int next_string (FILE * fp, const char * target)
{
  int slen = strlen (target), para = 0;
  char s[slen + 1];
  s[slen] = '\0';
  int len = 0, c = fgetc (fp);
  while (c != EOF && len < slen) {
    if (c == '{') para++;
    if (c == '}') para--;
    s[len++] = c;
    c = fgetc (fp);
  }
  while (c != EOF && para >= 0) {
    if (!strcmp (s, target) && para == 0)
      break;
    if (c == '{') para++;
    if (c == '}') para--;
    for (int i = 0; i < slen - 1; i++)
      s[i] = s[i+1];
    s[slen - 1] = c;
    c = fgetc (fp);
  }
  if (strcmp (s, target))
    c = -1;
  return c;
}

/**
# *input_gfs()*: Gerris simulation format

The function reads simulation data in the format used in
[Gerris](http://gerris.dalembert.upmc.fr) simulation files. This is
the reciprocal function of [*output_gfs()*](output.h#output_gfs).

The arguments and their default values are:

*fp*
: a file pointer. Default is *file* or stdin.

*list*
: a list of scalar fields to read. Default is *all*.

*file*
: the name of the file to read from (mutually exclusive with *fp*). */

trace
void input_gfs (FILE * fp = stdin,
		scalar * list = NULL,
		char * file = NULL)
{
  not_mpi_compatible();
  
  if (file && !(fp = fopen (file, "r"))) {
    perror (file);
    exit (1);
  }

  bool input_all = (list == all);
  if (!list) list = all;

  #if TREE
  init_grid (1);
  #endif

  next_char (fp, '{');
  
  char * s = qmalloc (1, char);
  int len = 0;
  int c = fgetc(fp);
  while (c != EOF && c != '}') {
    s[len++] = c;
    qrealloc (s, len + 1, char);
    s[len] = '\0';
    c = fgetc(fp);
  }
  if (c != '}') {
    fprintf (stderr, "input_gfs(): error: expecting '}'\n");
    exit (1);
  }
  
  char * s1 = strstr (s, "variables");
  if (!s1) {
    fprintf (stderr, "input_gfs(): error: expecting 'variables'\n");
    exit (1);
  }

  s1 = strstr (s1, "=");
  if (!s1) {
    fprintf (stderr, "input_gfs(): error: expecting '='\n");
    exit (1);
  }
  s1++;

  while (strchr (" \t", *s1))
    s1++;

  scalar * input = NULL;
  s1 = strtok (s1, ", \t");
  while (s1) {
    char * name = replace (s1, '_', '.', false);
    bool found = false;
    for (scalar s in list)
      if (!is_constant(s) && s.name && !strcmp (s.name, name)) {
	input = list_append (input, s);
	found = true; break;
      }
    if (!found) {
      if (input_all) {
	scalar s = new scalar;
	free (s.name);
	s.name = strdup (name);
	input = list_append (input, s);
      }
      else
	input = list_append (input, (scalar){INT_MAX});
    }
    free (name);
    s1 = strtok (NULL, ", \t");
  }
  free (s);

  next_char (fp, '{');
  double t1 = 0.;
  if (next_string (fp, "Time") >= 0) {
    next_char (fp, '{');
    next_char (fp, 't');
    next_char (fp, '=');
    if (fscanf (fp, "%lf", &t1) != 1) {
      fprintf (stderr, "input_gfs(): error: expecting 't'\n");
      exit (1);
    }
    next_char (fp, '}');
    next_char (fp, '}');
  }

  if (next_string (fp, "Box") < 0) {
    fprintf (stderr, "input_gfs(): error: expecting 'GfsBox'\n");
    exit (1);
  }

  next_char (fp, '{');
  next_char (fp, '{');
  next_char (fp, '\n');

  scalar * listm = {cm,fm};
  scalar * listr = !is_constant(cm) ? listm : NULL;
  NOT_UNUSED (listr);

  foreach_cell() {
    unsigned flags;
    if (fread (&flags, sizeof (unsigned), 1, fp) != 1) {
      fprintf (stderr, "input_gfs(): error: expecting 'flags'\n");
      exit (1);
    }
    if (!(flags & (1 << 4)) && is_leaf(cell))
      refine_cell (point, listr, 0, NULL);
    double a;
    if (fread (&a, sizeof (double), 1, fp) != 1 || a != -1) {
      fprintf (stderr, "input_gfs(): error: expecting '-1'\n");
      exit (1);
    }
    for (scalar s in input) {
      if (fread (&a, sizeof (double), 1, fp) != 1) {
	fprintf (stderr, "input_gfs(): error: expecting a scalar\n");
	exit (1);
      }
      if (s.i != INT_MAX) {
	if (s.v.x.i >= 0) {
	  // this is a vector component, we need to rotate from
	  // Z-ordering (Gerris) to N-ordering (Basilisk)
#if dimension >= 2
	  if (s.v.x.i == s.i) {
	    s = s.v.y;
	    s[] = a;
	  }
	  else if (s.v.y.i == s.i) {
	    s = s.v.x;
	    s[] = - a;
	  }
#endif
#if dimension >= 3
	  else
	    s[] = a;
#endif
	}
	else	
	  s[] = a;
      }
    }
    if (is_leaf(cell))
      continue;
  }
  for (scalar s in listm)
    if (!is_constant(s))
      s.dirty = true;
  for (scalar s in input)
    if (!is_constant(s))
      s.dirty = true;

  free (input);
  if (file)
    fclose (fp);

  // the events are advanced to catch up with the time
  while (t < t1 && events (false))
    t = tnext;
  events (false);
}

/**
# *input_grd()*: Raster format (Esri grid)

This function reads a scalar field from a [Raster
file](https://en.wikipedia.org/wiki/Esri_grid). This is the reciprocal
function of [*output_grd()*](output.h#output_grd).

The arguments and their default values are:

*s* 
: the scalar where the data will be stored. No default value. You
must specify this parameter

*fp*
: a file pointer. Default is *file* or stdin.

*file*
: the name of the file to read from (mutually exclusive with *fp*).

*nodatavalue*
: the value of the NoDataValue. Default is the same as that defined in
the raster file. 

*linear*
: if true, the raster data is bilinearly interpolated. Default is false.

*periodic*
: if true, the x-axis is treated as periodic. Default is false.

*zero*
: if true NoDataValue are replaced by zero. Default is false.

*smooth*
: the number of Laplacian smoothing passes applied to the data before
interpolation. Default is zero.
*/

void input_grd (scalar s,
		FILE * fp = stdin, const char * file = NULL,
		double nodatavalue = 0.,
		bool linear = false, bool periodic = false, bool zero = false,
		int smooth = 0)
{
  scalar input = s;

  if (file && !(fp = fopen (file, "r"))) {
    perror (file);
    exit (1);
  }
  
  // Variables for the Raster data
  double DeltaGRD;
  int nx, ny;
  double XG0, YG0, ndv;

  // header
  char waste[100];
  if (fscanf (fp, "%s %d", waste, &nx) != 2) {
    fprintf (stderr, "input_grd(): error reading 'nx'\n");
    if (file) fclose (fp);
    return;
  }
  if (fscanf (fp, "%s %d", waste, &ny) != 2) {
    fprintf (stderr, "input_grd(): error reading 'ny'\n");
    if (file) fclose (fp);
    return;
  }
  if (fscanf (fp, "%s %lf", waste, &XG0) != 2) {
    fprintf (stderr, "input_grd(): error reading 'XG0'\n");
    if (file) fclose (fp);
    return;    
  }
  if (fscanf (fp, "%s %lf", waste, &YG0) != 2) {
    fprintf (stderr, "input_grd(): error reading 'YG0'\n");
    if (file) fclose (fp);
    return;    
  }
  if (fscanf (fp, "%s %lf", waste, &DeltaGRD) != 2) {
    fprintf (stderr, "input_grd(): error reading 'DeltaGRD'\n");
    if (file) fclose (fp);
    return;    
  }
  if (fscanf (fp, "%s %lf", waste, &ndv) != 2) {
    fprintf (stderr, "input_grd(): error reading 'ndv'\n");
    if (file) fclose (fp);
    return;    
  }

  //default value of NoData value
  if (!nodatavalue)
    nodatavalue = ndv;

  // read the data
  double * value = qmalloc (nx*ny, double);
  for (int i = ny - 1; i >= 0; i--)
    for (int j = 0 ; j < nx; j++) {
      if (fscanf (fp, "%lf ", &value[j + i*nx]) != 1) {
	fprintf (stderr, "input_grd(): error reading value %d,%d\n", i, j);
	if (file) fclose (fp);
	free (value);
	return;
      }
      if (zero && value[j + i*nx] == ndv)
	value[j + i*nx] = 0.;
    }

  // Laplacian smoothing
  if (smooth > 0) {
    double * smoothed = qmalloc (nx*ny, double);
    for (int s = 0; s < smooth; s++) {
      for (int i = 0; i < ny; i++)
	for (int j = 0 ; j < nx; j++) {
	  int n = 0;
	  smoothed[j + i*nx] = 0.;
	  for (int k = -1; k <= 1; k++)
	    for (int l = -1; l <= 1; l++)
	      if ((l != 0 || k != 0) &&
		  i + k >= 0 && i + k < ny &&
		  j + l >= 0 && j + l < nx &&
		  value[j + l + (i + k)*nx] != ndv)
		smoothed[j + i*nx] += value[j + l + (i + k)*nx], n++;
	  if (n == 0)
	    smoothed[j + i*nx] = zero ? 0. : ndv;
	  else
	    smoothed[j + i*nx] /= n;
	}
      swap (double *, value, smoothed);
    }
    free (smoothed);
  }
  
  bool warning = false;  
  foreach (serial) {
    if (periodic || input.boundary[right] == periodic_bc) {
      if (x > XG0 + nx*DeltaGRD)
	x -= nx*DeltaGRD;
      else if (x < XG0)
	x += nx*DeltaGRD;
    }
    // Test if the point in the Basilisk area is included in the raster area
    int j = (x - XG0 + DeltaGRD/2.)/DeltaGRD;
    int i = (y - YG0 + DeltaGRD/2.)/DeltaGRD;    
    if (i >= 0 && i < ny && j >= 0 && j < nx) {
      double val;
      // Test if we are on the ring of data around the raster grid
      int j1 = (x - XG0)/DeltaGRD;
      int i1 = (y - YG0)/DeltaGRD;
      if (linear && i1 >= 0 && j1 >= 0 && i1 < ny - 1 && j1 < nx - 1 &&
	  value[j1 + i1*nx] != ndv && value[j1 + 1 + i1*nx] != ndv &&
	  value[j1 + (i1 + 1)*nx] != ndv && value[j1 + 1 + (i1 + 1)*nx] != ndv) {
	// bi-linear interpolation
	double dx = x - (j1*DeltaGRD + XG0);
	double dy = y - (i1*DeltaGRD + YG0);
	val = (value[j1 + i1*nx] +
	       dx*(value[j1 + 1 + i1*nx] - value[j1 + i1*nx])/DeltaGRD +
	       dy*(value[j1 + (i1 + 1)*nx] - value[j1 + i1*nx])/DeltaGRD +
	       dx*dy*(value[j1 + i1*nx] + value[j1 + 1 + (i1 + 1)*nx] -
		      value[j1 + (i1 + 1)*nx] - value[j1 + 1 + i1*nx])
	       /sq(DeltaGRD));
      }
      else
	val = value[j + i*nx];
      if (val == ndv)
	input[] = nodata;
      else
	input[] = val;
    }
    else {
      input[] = nodata;
      warning = true;
    }
  }
  free (value);
  
  if (warning)
    fprintf (stderr,
	     "input_grd(): Warning: Raster data is not covering all"
	     " the simulation area\n");

  if (file)
    fclose (fp);
}
