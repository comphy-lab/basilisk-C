/**
# Output functions

## *output_field()*: Multiple fields interpolated on a regular grid (text format)

This function interpolates a *list* of fields on a *n x n* regular
grid. The resulting data are written in text format in the file
pointed to by *fp*. The correspondance between column numbers and
variables is summarised in the first line of the file. The data are
written row-by-row and each row is separated from the next by a
blank line. This format is compatible with the *splot* command of
*gnuplot* i.e. one could use something like

~~~bash
gnuplot> set pm3d map
gnuplot> splot 'fields' u 1:2:4
~~~

The arguments and their default values are:

*list*
: list of fields to output. Default is *all*.

*fp*
: file pointer. Default is *stdout*.

*n*
: number of points along each dimension. Default is *N*.

*linear*
: use first-order (default) or bilinear interpolation. */

struct OutputField {
  scalar * list;
  FILE * fp;
  int n;
  bool linear;
};

trace
void output_field (struct OutputField p)
{
/* @if _QCCACC */
/*   /\* Download data from device for output *\/ */
/*   ACC_UPDATE_HOST */
/* @endif */

  if (!p.list) p.list = all;
  if (p.n == 0) p.n = N;
  if (!p.fp) p.fp = stdout;

  int len = list_len(p.list);
  double ** field = matrix_new (p.n, p.n, len*sizeof(double));
  
  double Delta = L0/p.n;
@if _QCCACC
  ACC(acc parallel copy(p) copyin(p.list[0:(len+1)]) copyout(field[0:p.n][0:(p.n*len)]) present(acc_dataarr))
  ACC(acc loop independent)
@endif
  for (int i = 0; i < p.n; i++) {
    double xp = Delta*i + X0 + Delta/2.;
@if _QCCACC
  ACC(acc loop independent)
@endif
    for (int j = 0; j < p.n; j++) {
      double yp = Delta*j + Y0 + Delta/2.;
      if (p.linear) {
	int k = 0;
	for (scalar s in p.list)
@if _QCCACC
          field[i][len*j + k++] = interpolate2Dacc (s.i, xp, yp, acc_dataarr);
@else
	  field[i][len*j + k++] = interpolate (s, xp, yp);
@endif
      }
      else {
	Point point;
@if _QCCACC
        locate2Dacc (xp, yp, &point);
	int k = 0;
	for (scalar s in p.list)
	  field[i][len*j + k++] = point.level >= 0 ? s[] : nodata;
@else
	point = locate (xp, yp);
	int k = 0;
	for (scalar s in p.list)
	  field[i][len*j + k++] = point.level >= 0 ? s[] : nodata;
@endif
      }
    }
  }

  if (pid() == 0) { // master
@if _MPI
    MPI_Reduce (MPI_IN_PLACE, field[0], len*p.n*p.n, MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
@endif
    fprintf (p.fp, "# 1:x 2:y");
    int i = 3;
    for (scalar s in p.list)
      fprintf (p.fp, " %d:%s", i++, s.name);
    fputc('\n', p.fp);
    for (int i = 0; i < p.n; i++) {
      double xp = Delta*i + X0 + Delta/2.;
      for (int j = 0; j < p.n; j++) {
	double yp = Delta*j + Y0 + Delta/2.;
	fprintf (p.fp, "%g %g", xp, yp);
	int k = 0;
	for (scalar s in p.list)
	  fprintf (p.fp, " %g", field[i][len*j + k++]);
	fputc ('\n', p.fp);
      }
      fputc ('\n', p.fp);
    }
    fflush (p.fp);
  }
@if _MPI
  else // slave
    MPI_Reduce (field[0], NULL, len*p.n*p.n, MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
@endif

  matrix_free (field);
}

/**
## *output_matrix()*: Single field interpolated on a regular grid (binary format)

This function writes a binary representation of a single field
interpolated on a regular *n x n* grid. The format is compatible with
the binary matrix format of gnuplot i.e. one could use

~~~bash
gnuplot> set pm3d map
gnuplot> splot 'matrix' binary u 2:1:3
~~~

The arguments and their default values are:

*f*
: a scalar field (compulsory).

*fp*
: file pointer. Default is *stdout*.

*n*
: number of points along each dimension. Default is *N*.

*linear*
: use first-order (default) or bilinear interpolation. */

struct OutputMatrix {
  scalar f;
  FILE * fp;
  int n;
  bool linear;
};

trace
void output_matrix (struct OutputMatrix p)
{
@if _QCCACC
  /* Download data from device for output */
  ACC_UPDATE_HOST
@endif

  if (p.n == 0) p.n = N;
  if (!p.fp) p.fp = stdout;
  float fn = p.n;
  float Delta = L0/fn;
  fwrite (&fn, sizeof(float), 1, p.fp);
  for (int j = 0; j < p.n; j++) {
    float yp = Delta*j + X0 + Delta/2.;
    fwrite (&yp, sizeof(float), 1, p.fp);
  }
  for (int i = 0; i < p.n; i++) {
    float xp = Delta*i + X0 + Delta/2.;
    fwrite (&xp, sizeof(float), 1, p.fp);
    for (int j = 0; j < p.n; j++) {
      float yp = Delta*j + Y0 + Delta/2., v;
      if (p.linear)
	v = interpolate (p.f, xp, yp);
      else {
	Point point = locate (xp, yp);
	assert (point.level >= 0);
	v = val(p.f);
      }
      fwrite (&v, sizeof(float), 1, p.fp);
    }
  }
  fflush (p.fp);
}

/**
## Colormaps

Colormaps are arrays of (127) red, green, blue triplets. */

#define NCMAP 127

typedef void (* colormap) (double cmap[NCMAP][3]);

void jet (double cmap[NCMAP][3])
{
  for (int i = 0; i < NCMAP; i++) {
    cmap[i][0] = 
      i <= 46 ? 0. : 
      i >= 111 ? -0.03125*(i - 111) + 1. :
      i >= 78 ? 1. : 
      0.03125*(i - 46);
    cmap[i][1] = 
      i <= 14 || i >= 111 ? 0. : 
      i >= 79 ? -0.03125*(i - 111) : 
      i <= 46 ? 0.03125*(i - 14) : 
      1.;
    cmap[i][2] =
      i >= 79 ? 0. :
      i >= 47 ? -0.03125*(i - 79) :
      i <= 14 ? 0.03125*(i - 14) + 1.:
      1.;
  }
}

void cool_warm (double cmap[NCMAP][3])
{
  /* diverging cool-warm from:
   *  http://www.sandia.gov/~kmorel/documents/ColorMaps/CoolWarmFloat33.csv
   * see also:
   *  Diverging Color Maps for Scientific Visualization (Expanded)
   *  Kenneth Moreland
   */
  static double basemap[33][3] = {
    {0.2298057,   0.298717966, 0.753683153},
    {0.26623388,  0.353094838, 0.801466763},
    {0.30386891,  0.406535296, 0.84495867},
    {0.342804478, 0.458757618, 0.883725899},
    {0.38301334,  0.50941904,  0.917387822},
    {0.424369608, 0.558148092, 0.945619588},
    {0.46666708,  0.604562568, 0.968154911},
    {0.509635204, 0.648280772, 0.98478814},
    {0.552953156, 0.688929332, 0.995375608},
    {0.596262162, 0.726149107, 0.999836203},
    {0.639176211, 0.759599947, 0.998151185},
    {0.681291281, 0.788964712, 0.990363227},
    {0.722193294, 0.813952739, 0.976574709},
    {0.761464949, 0.834302879, 0.956945269},
    {0.798691636, 0.849786142, 0.931688648},
    {0.833466556, 0.860207984, 0.901068838},
    {0.865395197, 0.86541021,  0.865395561},
    {0.897787179, 0.848937047, 0.820880546},
    {0.924127593, 0.827384882, 0.774508472},
    {0.944468518, 0.800927443, 0.726736146},
    {0.958852946, 0.769767752, 0.678007945},
    {0.96732803,  0.734132809, 0.628751763},
    {0.969954137, 0.694266682, 0.579375448},
    {0.966811177, 0.650421156, 0.530263762},
    {0.958003065, 0.602842431, 0.481775914},
    {0.943660866, 0.551750968, 0.434243684},
    {0.923944917, 0.49730856,  0.387970225},
    {0.89904617,  0.439559467, 0.343229596},
    {0.869186849, 0.378313092, 0.300267182},
    {0.834620542, 0.312874446, 0.259301199},
    {0.795631745, 0.24128379,  0.220525627},
    {0.752534934, 0.157246067, 0.184115123},
    {0.705673158, 0.01555616,  0.150232812}	
  };
  
  for (int i = 0; i < NCMAP; i++) {
    double x = i*(32 - 1e-10)/(NCMAP - 1);
    int j = x; x -= j;
    for (int k = 0; k < 3; k++)
      cmap[i][k] = (1. - x)*basemap[j][k] + x*basemap[j+1][k];
  }
}

void gray (double cmap[NCMAP][3])
{
  for (int i = 0; i < NCMAP; i++)
    for (int k = 0; k < 3; k++)
      cmap[i][k] = i/(NCMAP - 1.);
}

void randomap (double cmap[NCMAP][3])
{
  srand(0);
  for (int i = 0; i < NCMAP; i++)
    for (int k = 0; k < 3; k++)
      cmap[i][k] = (noise() + 1.)/2.;
}

/**
Given a colormap and a minimum and maximum value, this function
returns the red/green/blue triplet corresponding to *val*. */

typedef struct {
  unsigned char r, g, b;
} color;

@if _QCCACC
ACC(acc routine seq)
void colormap_color_acc (double * cmap, double val, double min, double max, color * c)
{
  if (val == nodata) {
    c -> r = c -> g = c -> b = 0; // nodata is black
  }
  else {
    val = val <= min ? 0. : val >= max ? 0.9999 : (val - min)/(max - min);
    int i = val*(NCMAP - 1);
    double coef = val*(NCMAP - 1) - i;
    unsigned char * c1 = (unsigned char *) c;
    for (int j = 0; j < 3; j++)
      c1[j] = 255*(cmap[i*3+j]*(1. - coef) + cmap[(i + 1)*3+j]*coef);
  }
}
@endif

color colormap_color (double cmap[NCMAP][3], 
		      double val, double min, double max)
{
  color c;
  if (val == nodata) {
    c.r = c.g = c.b = 0; // nodata is black
    return c;
  }    
  val = val <= min ? 0. : val >= max ? 0.9999 : (val - min)/(max - min);
  int i = val*(NCMAP - 1);
  double coef = val*(NCMAP - 1) - i;
  assert (i < NCMAP - 1);
  unsigned char * c1 = (unsigned char *) &c;
  for (int j = 0; j < 3; j++)
    c1[j] = 255*(cmap[i][j]*(1. - coef) + cmap[i + 1][j]*coef);
  return c;
}

/**
## *output_ppm()*: Portable PixMap (PPM) image output

Given a field, this function outputs a colormaped representation as a
[Portable PixMap](http://en.wikipedia.org/wiki/Netpbm_format) image.

If [ImageMagick](http://www.imagemagick.org/) is installed on the
system, this image can optionally be converted to any image format
supported by ImageMagick.

The arguments and their default values are:

*f*
: a scalar field (compulsory).

*fp*
: a file pointer. Default is stdout.

*n*
: number of pixels. Default is *N*.

*file*
: sets the name of the file used as output for
ImageMagick. This allows outputs in all formats supported by
ImageMagick. For example, one could use

~~~c
output_ppm (f, file = "f.png");
~~~

to get a [PNG](http://en.wikipedia.org/wiki/Portable_Network_Graphics)
image.

*min, max*
: minimum and maximum values used to define the
colorscale. By default these are set automatically using the *spread*
parameter. 

*spread*
: if not specified explicitly, *min* and *max* are set to the average 
of the field minus (resp. plus) *spread* times the standard deviation. 
By default *spread* is five. 

*linear*
: whether to use bilinear or first-order interpolation. Default is 
first-order.

*box*
: the lower-left and upper-right coordinates of the domain to consider.
 Default is the entire domain.

*mask*
: if set, this field will be used to mask out (in black), the regions 
of the domain for which *mask* is negative. 

*map*
: the colormap: *jet*, *cool_warm* or *gray*. Default is *jet*.
*/

struct OutputPPM {
  scalar f;
  FILE * fp;
  int n;
  char * file;
  double min, max, spread, z;
  bool linear;
  double box[2][2];
  scalar mask;
  colormap map;
};

trace
void output_ppm (struct OutputPPM p)
{
  //@if _QCCACC
  /* Download data from device for output */
  //ACC_UPDATE_HOST
  //@endif

  // default values
  if (p.n == 0) p.n = N;
  if (p.min == 0 && p.max == 0) {
    stats s = statsf (p.f);
    double avg = s.sum/s.volume, spread = (p.spread ? p.spread : 5.)*s.stddev;
    p.min = avg - spread; p.max = avg + spread;
  }
  if (p.box[0][0] == 0. && p.box[0][1] == 0. && 
      p.box[1][0] == 0. && p.box[1][1] == 0.) {
    p.box[0][0] = X0;      p.box[0][1] = Y0;
    p.box[1][0] = X0 + L0; p.box[1][1] = Y0 + L0;
  }
  if (!p.map)
    p.map = jet;

  double fn = p.n;
  double Delta = (p.box[1][0] - p.box[0][0])/fn;
  int ny = (p.box[1][1] - p.box[0][1])/Delta;
  
  color ** ppm = matrix_new (ny, p.n, sizeof(color));
  double cmap[NCMAP][3];
  p.map (cmap);

@if _QCCACC
  double cmapacc[NCMAP*3];
  for(int i = 0; i < NCMAP; i++){
    for(int j = 0; j < 3; j++) {
      cmapacc[i*3+j] = cmap[i][j];
    }
  }
@endif

  OMP_PARALLEL()
  OMP(omp for schedule(static))
@if _QCCACC
  ACC(acc parallel copy(p) copyout(ppm[0:ny][0:p.n]) copyin(cmapacc[0:(NCMAP*3)]) present(acc_dataarr))
  ACC(acc loop independent)
@endif
  for (int j = 0; j < ny; j++) {
    double yp = Delta*j + p.box[0][1] + Delta/2.;
@if _QCCACC
  ACC(acc loop independent)
@endif
    for (int i = 0; i < p.n; i++) {
      double xp = Delta*i + p.box[0][0] + Delta/2., v;
      if (p.mask.i) { // masking
	if (p.linear) {
@if _QCCACC
          double m = interpolate2Dacc (p.mask.i, xp, yp, acc_dataarr);
@else
	  double m = interpolate (p.mask, xp, yp, p.z);
@endif
	  if (m < 0.)
	    v = nodata;
	  else {
@if _QCCACC
            v = interpolate2Dacc (p.f.i, xp, yp, acc_dataarr);
@else
            v = interpolate (p.f, xp, yp, p.z);
@endif
          }
	}
	else {
	  Point point;
@if _QCCACC
          locate2Dacc (xp, yp, &point);
@else
	  point = locate (xp, yp, p.z);
@endif
	  if (point.level < 0 || val(p.mask) < 0.)
	    v = nodata;
	  else
	    v = val(p.f);
	}
      }
      else if (p.linear) {
@if _QCCACC
        v = interpolate2Dacc (p.f.i, xp, yp, acc_dataarr);
@else
        v = interpolate (p.f, xp, yp, p.z);
@endif
      }
      else {
	Point point;
@if _QCCACC
        locate2Dacc (xp, yp, &point);
@else
	point = locate (xp, yp, p.z);
@endif
	v = point.level >= 0 ? val(p.f) : nodata;
      }
@if _QCCACC
      colormap_color_acc (cmapacc, v, p.min, p.max, &(ppm[ny - 1 - j][i]));
@else
      ppm[ny - 1 - j][i] = colormap_color (cmap, v, p.min, p.max);
@endif
    }
  }
  OMP_END_PARALLEL()
  
  if (pid() == 0) { // master
@if _MPI
    MPI_Reduce (MPI_IN_PLACE, ppm[0], 3*ny*p.n, MPI_UNSIGNED_CHAR, MPI_MAX, 0,
		MPI_COMM_WORLD);
@endif
    if (!p.fp) p.fp = stdout;
    if (p.file) {
      char * command = malloc (strlen ("convert ppm:- ") + strlen (p.file) + 1);
      strcpy (command, "convert ppm:- ");
      strcat (command, p.file);
      p.fp = popen (command, "w");
      free (command);
    }
    
    fprintf (p.fp, "P6\n%u %u 255\n", p.n, ny);
    fwrite (((void **) ppm)[0], sizeof(color), ny*p.n, p.fp);
    
    if (p.file)
      pclose (p.fp);
    else
      fflush (p.fp);
  }
@if _MPI
  else // slave
    MPI_Reduce (ppm[0], NULL, 3*ny*p.n, MPI_UNSIGNED_CHAR, MPI_MAX, 0,
		MPI_COMM_WORLD);
@endif
    
  matrix_free (ppm);
}

/**
We also define a convenience macro which calls *avconv* to generate a
MPEG movie from the PPM stream. */

#define mpeg(c,name,...) {					\
    static FILE * fp = popen ("ppm2mpeg > " name, "w");		\
    output_ppm (c, fp, __VA_ARGS__);				\
  }

/**
## *output_grd()*: ESRI ASCII Grid format

The [ESRI GRD format](http://en.wikipedia.org/wiki/Esri_grid) is a
standard format for importing raster data into [GIS
systems](http://en.wikipedia.org/wiki/Geographic_information_system).

The arguments and their default values are:

*f*
: a scalar field (compulsory).

*fp*
: a file pointer. Default is stdout.

$\Delta$
: size of a grid element. Default is 1/N.

*linear*
: whether to use bilinear or first-order interpolation. Default is 
first-order.

*box*
: the lower-left and upper-right coordinates of the domain to consider.
 Default is the entire domain.

*mask*
: if set, this field will be used to mask out, the regions 
of the domain for which *mask* is negative. */

struct OutputGRD {
  scalar f;
  FILE * fp;
  double Delta;
  bool linear;
  double box[2][2];
  scalar mask;
};

trace
void output_grd (struct OutputGRD p)
{
@if _QCCACC
  /* Download data from device for output */
  ACC_UPDATE_HOST
@endif

  // default values
  if (!p.fp) p.fp = stdout;
  if (p.box[0][0] == 0. && p.box[0][1] == 0. && 
      p.box[1][0] == 0. && p.box[1][1] == 0.) {
    p.box[0][0] = X0;      p.box[0][1] = Y0;
    p.box[1][0] = X0 + L0; p.box[1][1] = Y0 + L0;
    if (p.Delta == 0) p.Delta = L0/N;
  }

  double Delta = p.Delta;
  int nx = (p.box[1][0] - p.box[0][0])/Delta;
  int ny = (p.box[1][1] - p.box[0][1])/Delta;

  // header
  fprintf (p.fp, "ncols          %d\n", nx);
  fprintf (p.fp, "nrows          %d\n", ny);
  fprintf (p.fp, "xllcorner      %g\n", p.box[0][0]);
  fprintf (p.fp, "yllcorner      %g\n", p.box[0][1]);
  fprintf (p.fp, "cellsize       %g\n", Delta);
  fprintf (p.fp, "nodata_value   -9999\n");
  
  // data
  for (int j = ny-1; j >= 0; j--) {
    double yp = Delta*j + p.box[0][1] + Delta/2.;
    for (int i = 0; i < nx; i++) {
      double xp = Delta*i + p.box[0][0] + Delta/2., v;
      if (p.mask.i) { // masking
	if (p.linear) {
	  double m = interpolate (p.mask, xp, yp);
	  if (m < 0.)
	    v = nodata;
	  else
	    v = interpolate (p.f, xp, yp);
	}
	else {
	  Point point = locate (xp, yp);
	  if (point.level < 0 || val(p.mask) < 0.)
	    v = nodata;
	  else
	    v = val(p.f);
	}
      }
      else if (p.linear)
	v = interpolate (p.f, xp, yp);
      else {
	Point point = locate (xp, yp);
	v = point.level >= 0 ? val(p.f) : nodata;
      }
      if (v == nodata)
	fprintf (p.fp, "-9999 ");
      else
	fprintf (p.fp, "%f ", v);
    }
    fprintf (p.fp, "\n");
  }

  fflush (p.fp);
}

#if MULTIGRID

/**
## *output_gfs()*: Gerris simulation format

The function writes simulation data in the format used in
[Gerris](http://gerris.dalembert.upmc.fr) simulation files. These files can be read
with GfsView.

The arguments and their default values are:

*fp*
: a file pointer. Default is *name* or stdout.

*list*
: a list of scalar fields to write. Default is *all*. 

*t*
: the physical time. Default is zero. 

*file*
: the name of the file to write to (mutually exclusive with *fp*).

*translate*
: whether to replace "well-known" Basilisk variables with their Gerris
equivalents.
*/

struct OutputGfs {
  FILE * fp;
  scalar * list;
  double t;
  char * file;
  bool translate;
};

static char * replace (const char * input, int target, int with,
		       bool translate)
{
  if (translate) {
    if (!strcmp (input, "u.x"))
      return strdup ("U");
    if (!strcmp (input, "u.y"))
      return strdup ("V");
    if (!strcmp (input, "u.z"))
      return strdup ("W");
  }
  char * name = strdup (input), * i = name;
  while (*i != '\0') {
    if (*i == target)
      *i = with;
    i++;
  }
  return name;
}

trace
void output_gfs (struct OutputGfs p)
{
@if _QCCACC
  /* Download data from device for output */
  ACC_UPDATE_HOST
@endif

  bool opened = false;
  if (p.fp == NULL) {
    if (p.file == NULL)
      p.fp = stdout;
    else if (!(p.fp = fopen (p.file, "w"))) {
      perror (p.file);
      exit (1);
    }
    else
      opened = true;
  }
  scalar * list = p.list ? p.list : list_copy (all);

  fprintf (p.fp, 
	   "1 0 GfsSimulation GfsBox GfsGEdge { binary = 1"
	   " x = %g y = %g ",
	   0.5 + X0/L0, 0.5 + Y0/L0);
#if dimension == 3
  fprintf (p.fp, "z = %g ", 0.5 + Z0/L0);
#endif

  if (list != NULL && list[0].i != -1) {
    scalar s = list[0];
    char * name = replace (s.name, '.', '_', p.translate);
    fprintf (p.fp, "variables = %s", name);
    free (name);
    for (int i = 1; i < list_len(list); i++) {
      scalar s = list[i];
      if (s.name) {
	char * name = replace (s.name, '.', '_', p.translate);
	fprintf (p.fp, ",%s", name);
	free (name);
      }
    }
    fprintf (p.fp, " ");
  }
  fprintf (p.fp, "} {\n");
  if (p.t > 0.)
    fprintf (p.fp, "  Time { t = %g }\n", p.t);
  if (L0 != 1.)
    fprintf (p.fp, "  PhysicalParams { L = %g }\n", L0);
  fprintf (p.fp, "}\nGfsBox { x = 0 y = 0 z = 0 } {\n");

@if _MPI
  long header;
  if ((header = ftell (p.fp)) < 0) {
    perror ("output_gfs(): error in header");
    exit (1);
  }
  int cell_size = sizeof(unsigned) + sizeof(double);
  for (scalar s in list)
    if (s.name)
      cell_size += sizeof(double);
  scalar index = new scalar;
  size_t total_size = header + (z_indexing (index, false) + 1)*cell_size;
@endif
  
  // see gerris/ftt.c:ftt_cell_write()
  //     gerris/domain.c:gfs_cell_write()
  foreach_cell() {
@if _MPI // fixme: this won't work when combining MPI and mask()
    if (is_local(cell))
@endif
    {
@if _MPI
      if (fseek (p.fp, header + index[]*cell_size, SEEK_SET) < 0) {
	perror ("output_gfs(): error while seeking");
	exit (1);
      }
@endif
      unsigned flags = 
	level == 0 ? 0 :
#if dimension == 1
	child.x == 1;
#elif dimension == 2
      child.x == -1 && child.y == -1 ? 0 :
	child.x == -1 && child.y ==  1 ? 1 :
	child.x ==  1 && child.y == -1 ? 2 : 
	3;
#else // dimension == 3
      child.x == -1 && child.y == -1 && child.z == -1  ? 0 :
	child.x == -1 && child.y == -1 && child.z ==  1  ? 1 :
	child.x == -1 && child.y ==  1 && child.z == -1  ? 2 : 
	child.x == -1 && child.y ==  1 && child.z ==  1  ? 3 : 
	child.x ==  1 && child.y == -1 && child.z == -1 ? 4 :
	child.x ==  1 && child.y == -1 && child.z ==  1 ? 5 :
	child.x ==  1 && child.y ==  1 && child.z == -1 ? 6 : 
	7;
#endif
      if (is_leaf(cell))
	flags |= (1 << 4);
      fwrite (&flags, sizeof (unsigned), 1, p.fp);
      double a = -1;
      fwrite (&a, sizeof (double), 1, p.fp);
      for (scalar s in list)
	if (s.name) {
	  if (s.v.x.i >= 0) {
	    // this is a vector component, we need to rotate from
	    // N-ordering (Basilisk) to Z-ordering (Gerris)
#if dimension >= 2
	    if (s.v.x.i == s.i) {
	      s = s.v.y;
	      a = is_local(cell) && !isnan(s[]) && s[] != nodata ? s[] : DBL_MAX;
	    }
	    else if (s.v.y.i == s.i) {
	      s = s.v.x;
	      a = is_local(cell) && !isnan(s[]) && s[] != nodata ?
		- s[] : DBL_MAX;
	    }
#endif
#if dimension >= 3
	    else
	      a = is_local(cell) && !isnan(s[]) && s[] != nodata ? s[] : DBL_MAX;
#endif
	  }
	  else
	    a = is_local(cell) && !isnan(s[]) && s[] != nodata ? s[] : DBL_MAX;
	  fwrite (&a, sizeof (double), 1, p.fp);
	}
    }
    if (is_leaf(cell))
      continue;
  }
  
@if _MPI
  delete ({index});
  if (!pid() && fseek (p.fp, total_size, SEEK_SET) < 0) {
    perror ("output_gfs(): error while finishing");
    exit (1);
  }
  if (!pid())
@endif  
    fputs ("}\n", p.fp);
  fflush (p.fp);

  if (!p.list)
    free (list);
  if (opened)
    fclose (p.fp);
}

/**
## *dump()*: Basilisk snapshots

This function (together with *restore()*) can be used to dump/restore
entire simulations.

The arguments and their default values are:

*fp*
: a file pointer. Default is *name* or stdout.

*list*
: a list of scalar fields to write. Default is *all*. 

*t*
: the physical time. Default is zero. 

*file*
: the name of the file to write to (mutually exclusive with *fp*).
*/

struct Dump {
  FILE * fp;
  scalar * list;
  double t;
  char * file;
};

struct DumpHeader {
  double t;
  long len;
  int depth;
};

trace
void dump (struct Dump p)
{
@if _QCCACC
  /* Download data from device for output */
  ACC_UPDATE_HOST 
@endif

  FILE * fp = p.fp;
  scalar * lista = p.list ? p.list : all, * list = NULL;
  char * file = p.file;

  for (scalar s in lista)
    if (!s.face && s.i != cm.i)
      list = list_add (list, s);
  
  if (file && (fp = fopen (file, "w")) == NULL) {
    perror (file);
    exit (1);
  }
  assert (fp);

  struct DumpHeader header = { p.t, list_len(list), depth() };

  if (pid() == 0 && fwrite (&header, sizeof(header), 1, fp) < 1) {
    perror ("dump(): error while writing header");
    exit (1);
  }

  scalar index = {-1};
  
@if _MPI // Parallel
  index = new scalar;
  z_indexing (index, false);
  int cell_size = sizeof(unsigned) + header.len*sizeof(double);
@endif
  
  foreach_cell() {
@if _MPI // fixme: this won't work when combining MPI and mask()
    if (is_local(cell))
@endif
    {
@if _MPI
      if (fseek (fp, sizeof(header) + index[]*cell_size, SEEK_SET) < 0) {
	perror ("dump(): error while seeking");
	exit (1);
      }
@endif
      unsigned flags = is_leaf(cell) ? leaf : 0;
      if (fwrite (&flags, sizeof(unsigned), 1, fp) < 1) {
	perror ("dump(): error while writing flags");
	exit (1);
      }
      for (scalar s in list)
	if (fwrite (&s[], sizeof(double), 1, fp) < 1) {
	  perror ("dump(): error while writing scalars");
	  exit (1);
	}
    }
    if (is_leaf(cell))
      continue;
  }

  delete ({index});
  
  free (list);
  if (file)
    fclose (fp);
}

trace
bool restore (struct Dump p)
{
  FILE * fp = p.fp;
  scalar * lista = p.list ? p.list : all, * list = NULL;
  char * file = p.file;
  if (file && (fp = fopen (file, "r")) == NULL)
    return false;
  assert (fp);

  for (scalar s in lista)
    if (!s.face && s.i != cm.i)
      list = list_add (list, s);

  struct DumpHeader header;
  
  double t = 0.;
  if (fread (&header, sizeof(header), 1, fp) < 1) {
    fprintf (stderr, "restore(): error: expecting header\n");
    exit (1);
  }
  if (header.len != list_len (list)) {
    fprintf (stderr,
	     "restore(): error: the list lengths don't match: %ld != %d\n",
	     header.len, list_len (list));
    exit (1);
  }

#if TREE
  init_grid (1);

  foreach_cell() {
    cell.pid = pid();
    cell.flags |= active;
  }
  tree->dirty = true;
#else
  init_grid (1 << header.depth);
#endif

@if _QCCACC
  /* Download data from device for input */
  ACC_UPDATE_HOST
@endif

  foreach_cell() {
    unsigned flags;
    if (fread (&flags, sizeof(unsigned), 1, fp) != 1) {
      fprintf (stderr, "restore(): error: expecting 'flags'\n");
      exit (1);
    }
    for (scalar s in list)
      if (fread (&s[], sizeof(double), 1, fp) != 1) {
	fprintf (stderr, "restore(): error: expecting '%s'\n", s.name);
	exit (1);
      }/**
# Output functions

## *output_field()*: Multiple fields interpolated on a regular grid (text format)

This function interpolates a *list* of fields on a *n x n* regular
grid. The resulting data are written in text format in the file
pointed to by *fp*. The correspondance between column numbers and
variables is summarised in the first line of the file. The data are
written row-by-row and each row is separated from the next by a
blank line. This format is compatible with the *splot* command of
*gnuplot* i.e. one could use something like

~~~bash
gnuplot> set pm3d map
gnuplot> splot 'fields' u 1:2:4
~~~

The arguments and their default values are:

*list*
: list of fields to output. Default is *all*.

*fp*
: file pointer. Default is *stdout*.

*n*
: number of points along each dimension. Default is *N*.

*linear*
: use first-order (default) or bilinear interpolation. */

struct OutputField {
  scalar * list;
  FILE * fp;
  int n;
  bool linear;
};

trace
void output_field (struct OutputField p)
{
/* @if _QCCACC */
/*   /\* Download data from device for output *\/ */
/*   ACC_UPDATE_HOST */
/* @endif */

  if (!p.list) p.list = all;
  if (p.n == 0) p.n = N;
  if (!p.fp) p.fp = stdout;

  int len = list_len(p.list);
  double ** field = matrix_new (p.n, p.n, len*sizeof(double));
  
  double Delta = L0/p.n;
@if _QCCACC
  ACC(acc parallel copy(p) copyin(p.list[0:(len+1)]) copyout(field[0:p.n][0:(p.n*len)]) present(acc_dataarr))
  ACC(acc loop independent)
@endif
  for (int i = 0; i < p.n; i++) {
    double xp = Delta*i + X0 + Delta/2.;
@if _QCCACC
  ACC(acc loop independent)
@endif
    for (int j = 0; j < p.n; j++) {
      double yp = Delta*j + Y0 + Delta/2.;
      if (p.linear) {
	int k = 0;
	for (scalar s in p.list)
@if _QCCACC
          field[i][len*j + k++] = interpolate2Dacc (s.i, xp, yp, acc_dataarr);
@else
	  field[i][len*j + k++] = interpolate (s, xp, yp);
@endif
      }
      else {
	Point point;
@if _QCCACC
        locate2Dacc (xp, yp, &point);
	int k = 0;
	for (scalar s in p.list)
	  field[i][len*j + k++] = point.level >= 0 ? s[] : nodata;
@else
	point = locate (xp, yp);
	int k = 0;
	for (scalar s in p.list)
	  field[i][len*j + k++] = point.level >= 0 ? s[] : nodata;
@endif
      }
    }
  }

  if (pid() == 0) { // master
@if _MPI
    MPI_Reduce (MPI_IN_PLACE, field[0], len*p.n*p.n, MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
@endif
    fprintf (p.fp, "# 1:x 2:y");
    int i = 3;
    for (scalar s in p.list)
      fprintf (p.fp, " %d:%s", i++, s.name);
    fputc('\n', p.fp);
    for (int i = 0; i < p.n; i++) {
      double xp = Delta*i + X0 + Delta/2.;
      for (int j = 0; j < p.n; j++) {
	double yp = Delta*j + Y0 + Delta/2.;
	fprintf (p.fp, "%g %g", xp, yp);
	int k = 0;
	for (scalar s in p.list)
	  fprintf (p.fp, " %g", field[i][len*j + k++]);
	fputc ('\n', p.fp);
      }
      fputc ('\n', p.fp);
    }
    fflush (p.fp);
  }
@if _MPI
  else // slave
    MPI_Reduce (field[0], NULL, len*p.n*p.n, MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
@endif

  matrix_free (field);
}

/**
## *output_matrix()*: Single field interpolated on a regular grid (binary format)

This function writes a binary representation of a single field
interpolated on a regular *n x n* grid. The format is compatible with
the binary matrix format of gnuplot i.e. one could use

~~~bash
gnuplot> set pm3d map
gnuplot> splot 'matrix' binary u 2:1:3
~~~

The arguments and their default values are:

*f*
: a scalar field (compulsory).

*fp*
: file pointer. Default is *stdout*.

*n*
: number of points along each dimension. Default is *N*.

*linear*
: use first-order (default) or bilinear interpolation. */

struct OutputMatrix {
  scalar f;
  FILE * fp;
  int n;
  bool linear;
};

trace
void output_matrix (struct OutputMatrix p)
{
@if _QCCACC
  /* Download data from device for output */
  ACC_UPDATE_HOST
@endif

  if (p.n == 0) p.n = N;
  if (!p.fp) p.fp = stdout;
  float fn = p.n;
  float Delta = L0/fn;
  fwrite (&fn, sizeof(float), 1, p.fp);
  for (int j = 0; j < p.n; j++) {
    float yp = Delta*j + X0 + Delta/2.;
    fwrite (&yp, sizeof(float), 1, p.fp);
  }
  for (int i = 0; i < p.n; i++) {
    float xp = Delta*i + X0 + Delta/2.;
    fwrite (&xp, sizeof(float), 1, p.fp);
    for (int j = 0; j < p.n; j++) {
      float yp = Delta*j + Y0 + Delta/2., v;
      if (p.linear)
	v = interpolate (p.f, xp, yp);
      else {
	Point point = locate (xp, yp);
	assert (point.level >= 0);
	v = val(p.f);
      }
      fwrite (&v, sizeof(float), 1, p.fp);
    }
  }
  fflush (p.fp);
}

/**
## Colormaps

Colormaps are arrays of (127) red, green, blue triplets. */

#define NCMAP 127

typedef void (* colormap) (double cmap[NCMAP][3]);

void jet (double cmap[NCMAP][3])
{
  for (int i = 0; i < NCMAP; i++) {
    cmap[i][0] = 
      i <= 46 ? 0. : 
      i >= 111 ? -0.03125*(i - 111) + 1. :
      i >= 78 ? 1. : 
      0.03125*(i - 46);
    cmap[i][1] = 
      i <= 14 || i >= 111 ? 0. : 
      i >= 79 ? -0.03125*(i - 111) : 
      i <= 46 ? 0.03125*(i - 14) : 
      1.;
    cmap[i][2] =
      i >= 79 ? 0. :
      i >= 47 ? -0.03125*(i - 79) :
      i <= 14 ? 0.03125*(i - 14) + 1.:
      1.;
  }
}

void cool_warm (double cmap[NCMAP][3])
{
  /* diverging cool-warm from:
   *  http://www.sandia.gov/~kmorel/documents/ColorMaps/CoolWarmFloat33.csv
   * see also:
   *  Diverging Color Maps for Scientific Visualization (Expanded)
   *  Kenneth Moreland
   */
  static double basemap[33][3] = {
    {0.2298057,   0.298717966, 0.753683153},
    {0.26623388,  0.353094838, 0.801466763},
    {0.30386891,  0.406535296, 0.84495867},
    {0.342804478, 0.458757618, 0.883725899},
    {0.38301334,  0.50941904,  0.917387822},
    {0.424369608, 0.558148092, 0.945619588},
    {0.46666708,  0.604562568, 0.968154911},
    {0.509635204, 0.648280772, 0.98478814},
    {0.552953156, 0.688929332, 0.995375608},
    {0.596262162, 0.726149107, 0.999836203},
    {0.639176211, 0.759599947, 0.998151185},
    {0.681291281, 0.788964712, 0.990363227},
    {0.722193294, 0.813952739, 0.976574709},
    {0.761464949, 0.834302879, 0.956945269},
    {0.798691636, 0.849786142, 0.931688648},
    {0.833466556, 0.860207984, 0.901068838},
    {0.865395197, 0.86541021,  0.865395561},
    {0.897787179, 0.848937047, 0.820880546},
    {0.924127593, 0.827384882, 0.774508472},
    {0.944468518, 0.800927443, 0.726736146},
    {0.958852946, 0.769767752, 0.678007945},
    {0.96732803,  0.734132809, 0.628751763},
    {0.969954137, 0.694266682, 0.579375448},
    {0.966811177, 0.650421156, 0.530263762},
    {0.958003065, 0.602842431, 0.481775914},
    {0.943660866, 0.551750968, 0.434243684},
    {0.923944917, 0.49730856,  0.387970225},
    {0.89904617,  0.439559467, 0.343229596},
    {0.869186849, 0.378313092, 0.300267182},
    {0.834620542, 0.312874446, 0.259301199},
    {0.795631745, 0.24128379,  0.220525627},
    {0.752534934, 0.157246067, 0.184115123},
    {0.705673158, 0.01555616,  0.150232812}	
  };
  
  for (int i = 0; i < NCMAP; i++) {
    double x = i*(32 - 1e-10)/(NCMAP - 1);
    int j = x; x -= j;
    for (int k = 0; k < 3; k++)
      cmap[i][k] = (1. - x)*basemap[j][k] + x*basemap[j+1][k];
  }
}

void gray (double cmap[NCMAP][3])
{
  for (int i = 0; i < NCMAP; i++)
    for (int k = 0; k < 3; k++)
      cmap[i][k] = i/(NCMAP - 1.);
}

void randomap (double cmap[NCMAP][3])
{
  srand(0);
  for (int i = 0; i < NCMAP; i++)
    for (int k = 0; k < 3; k++)
      cmap[i][k] = (noise() + 1.)/2.;
}

/**
Given a colormap and a minimum and maximum value, this function
returns the red/green/blue triplet corresponding to *val*. */

typedef struct {
  unsigned char r, g, b;
} color;

@if _QCCACC
ACC(acc routine seq)
void colormap_color_acc (double * cmap, double val, double min, double max, color * c)
{
  if (val == nodata) {
    c -> r = c -> g = c -> b = 0; // nodata is black
  }
  else {
    val = val <= min ? 0. : val >= max ? 0.9999 : (val - min)/(max - min);
    int i = val*(NCMAP - 1);
    double coef = val*(NCMAP - 1) - i;
    unsigned char * c1 = (unsigned char *) c;
    for (int j = 0; j < 3; j++)
      c1[j] = 255*(cmap[i*3+j]*(1. - coef) + cmap[(i + 1)*3+j]*coef);
  }
}
@endif

color colormap_color (double cmap[NCMAP][3], 
		      double val, double min, double max)
{
  color c;
  if (val == nodata) {
    c.r = c.g = c.b = 0; // nodata is black
    return c;
  }    
  val = val <= min ? 0. : val >= max ? 0.9999 : (val - min)/(max - min);
  int i = val*(NCMAP - 1);
  double coef = val*(NCMAP - 1) - i;
  assert (i < NCMAP - 1);
  unsigned char * c1 = (unsigned char *) &c;
  for (int j = 0; j < 3; j++)
    c1[j] = 255*(cmap[i][j]*(1. - coef) + cmap[i + 1][j]*coef);
  return c;
}

/**
## *output_ppm()*: Portable PixMap (PPM) image output

Given a field, this function outputs a colormaped representation as a
[Portable PixMap](http://en.wikipedia.org/wiki/Netpbm_format) image.

If [ImageMagick](http://www.imagemagick.org/) is installed on the
system, this image can optionally be converted to any image format
supported by ImageMagick.

The arguments and their default values are:

*f*
: a scalar field (compulsory).

*fp*
: a file pointer. Default is stdout.

*n*
: number of pixels. Default is *N*.

*file*
: sets the name of the file used as output for
ImageMagick. This allows outputs in all formats supported by
ImageMagick. For example, one could use

~~~c
output_ppm (f, file = "f.png");
~~~

to get a [PNG](http://en.wikipedia.org/wiki/Portable_Network_Graphics)
image.

*min, max*
: minimum and maximum values used to define the
colorscale. By default these are set automatically using the *spread*
parameter. 

*spread*
: if not specified explicitly, *min* and *max* are set to the average 
of the field minus (resp. plus) *spread* times the standard deviation. 
By default *spread* is five. 

*linear*
: whether to use bilinear or first-order interpolation. Default is 
first-order.

*box*
: the lower-left and upper-right coordinates of the domain to consider.
 Default is the entire domain.

*mask*
: if set, this field will be used to mask out (in black), the regions 
of the domain for which *mask* is negative. 

*map*
: the colormap: *jet*, *cool_warm* or *gray*. Default is *jet*.
*/

struct OutputPPM {
  scalar f;
  FILE * fp;
  int n;
  char * file;
  double min, max, spread, z;
  bool linear;
  double box[2][2];
  scalar mask;
  colormap map;
};

trace
void output_ppm (struct OutputPPM p)
{
  //@if _QCCACC
  /* Download data from device for output */
  //ACC_UPDATE_HOST
  //@endif

  // default values
  if (p.n == 0) p.n = N;
  if (p.min == 0 && p.max == 0) {
    stats s = statsf (p.f);
    double avg = s.sum/s.volume, spread = (p.spread ? p.spread : 5.)*s.stddev;
    p.min = avg - spread; p.max = avg + spread;
  }
  if (p.box[0][0] == 0. && p.box[0][1] == 0. && 
      p.box[1][0] == 0. && p.box[1][1] == 0.) {
    p.box[0][0] = X0;      p.box[0][1] = Y0;
    p.box[1][0] = X0 + L0; p.box[1][1] = Y0 + L0;
  }
  if (!p.map)
    p.map = jet;

  double fn = p.n;
  double Delta = (p.box[1][0] - p.box[0][0])/fn;
  int ny = (p.box[1][1] - p.box[0][1])/Delta;
  
  color ** ppm = matrix_new (ny, p.n, sizeof(color));
  double cmap[NCMAP][3];
  p.map (cmap);

@if _QCCACC
  double cmapacc[NCMAP*3];
  for(int i = 0; i < NCMAP; i++){
    for(int j = 0; j < 3; j++) {
      cmapacc[i*3+j] = cmap[i][j];
    }
  }
@endif

  OMP_PARALLEL()
  OMP(omp for schedule(static))
@if _QCCACC
  ACC(acc parallel copy(p) copyout(ppm[0:ny][0:p.n]) copyin(cmapacc[0:(NCMAP*3)]) present(acc_dataarr))
  ACC(acc loop independent)
@endif
  for (int j = 0; j < ny; j++) {
    double yp = Delta*j + p.box[0][1] + Delta/2.;
@if _QCCACC
  ACC(acc loop independent)
@endif
    for (int i = 0; i < p.n; i++) {
      double xp = Delta*i + p.box[0][0] + Delta/2., v;
      if (p.mask.i) { // masking
	if (p.linear) {
@if _QCCACC
          double m = interpolate2Dacc (p.mask.i, xp, yp, acc_dataarr);
@else
	  double m = interpolate (p.mask, xp, yp, p.z);
@endif
	  if (m < 0.)
	    v = nodata;
	  else {
@if _QCCACC
            v = interpolate2Dacc (p.f.i, xp, yp, acc_dataarr);
@else
            v = interpolate (p.f, xp, yp, p.z);
@endif
          }
	}
	else {
	  Point point;
@if _QCCACC
          locate2Dacc (xp, yp, &point);
@else
	  point = locate (xp, yp, p.z);
@endif
	  if (point.level < 0 || val(p.mask) < 0.)
	    v = nodata;
	  else
	    v = val(p.f);
	}
      }
      else if (p.linear) {
@if _QCCACC
        v = interpolate2Dacc (p.f.i, xp, yp, acc_dataarr);
@else
        v = interpolate (p.f, xp, yp, p.z);
@endif
      }
      else {
	Point point;
@if _QCCACC
        locate2Dacc (xp, yp, &point);
@else
	point = locate (xp, yp, p.z);
@endif
	v = point.level >= 0 ? val(p.f) : nodata;
      }
@if _QCCACC
      colormap_color_acc (cmapacc, v, p.min, p.max, &(ppm[ny - 1 - j][i]));
@else
      ppm[ny - 1 - j][i] = colormap_color (cmap, v, p.min, p.max);
@endif
    }
  }
  OMP_END_PARALLEL()
  
  if (pid() == 0) { // master
@if _MPI
    MPI_Reduce (MPI_IN_PLACE, ppm[0], 3*ny*p.n, MPI_UNSIGNED_CHAR, MPI_MAX, 0,
		MPI_COMM_WORLD);
@endif
    if (!p.fp) p.fp = stdout;
    if (p.file) {
      char * command = malloc (strlen ("convert ppm:- ") + strlen (p.file) + 1);
      strcpy (command, "convert ppm:- ");
      strcat (command, p.file);
      p.fp = popen (command, "w");
      free (command);
    }
    
    fprintf (p.fp, "P6\n%u %u 255\n", p.n, ny);
    fwrite (((void **) ppm)[0], sizeof(color), ny*p.n, p.fp);
    
    if (p.file)
      pclose (p.fp);
    else
      fflush (p.fp);
  }
@if _MPI
  else // slave
    MPI_Reduce (ppm[0], NULL, 3*ny*p.n, MPI_UNSIGNED_CHAR, MPI_MAX, 0,
		MPI_COMM_WORLD);
@endif
    
  matrix_free (ppm);
}

/**
We also define a convenience macro which calls *avconv* to generate a
MPEG movie from the PPM stream. */

#define mpeg(c,name,...) {					\
    static FILE * fp = popen ("ppm2mpeg > " name, "w");		\
    output_ppm (c, fp, __VA_ARGS__);				\
  }

/**
## *output_grd()*: ESRI ASCII Grid format

The [ESRI GRD format](http://en.wikipedia.org/wiki/Esri_grid) is a
standard format for importing raster data into [GIS
systems](http://en.wikipedia.org/wiki/Geographic_information_system).

The arguments and their default values are:

*f*
: a scalar field (compulsory).

*fp*
: a file pointer. Default is stdout.

$\Delta$
: size of a grid element. Default is 1/N.

*linear*
: whether to use bilinear or first-order interpolation. Default is 
first-order.

*box*
: the lower-left and upper-right coordinates of the domain to consider.
 Default is the entire domain.

*mask*
: if set, this field will be used to mask out, the regions 
of the domain for which *mask* is negative. */

struct OutputGRD {
  scalar f;
  FILE * fp;
  double Delta;
  bool linear;
  double box[2][2];
  scalar mask;
};

trace
void output_grd (struct OutputGRD p)
{
@if _QCCACC
  /* Download data from device for output */
  ACC_UPDATE_HOST
@endif

  // default values
  if (!p.fp) p.fp = stdout;
  if (p.box[0][0] == 0. && p.box[0][1] == 0. && 
      p.box[1][0] == 0. && p.box[1][1] == 0.) {
    p.box[0][0] = X0;      p.box[0][1] = Y0;
    p.box[1][0] = X0 + L0; p.box[1][1] = Y0 + L0;
    if (p.Delta == 0) p.Delta = L0/N;
  }

  double Delta = p.Delta;
  int nx = (p.box[1][0] - p.box[0][0])/Delta;
  int ny = (p.box[1][1] - p.box[0][1])/Delta;

  // header
  fprintf (p.fp, "ncols          %d\n", nx);
  fprintf (p.fp, "nrows          %d\n", ny);
  fprintf (p.fp, "xllcorner      %g\n", p.box[0][0]);
  fprintf (p.fp, "yllcorner      %g\n", p.box[0][1]);
  fprintf (p.fp, "cellsize       %g\n", Delta);
  fprintf (p.fp, "nodata_value   -9999\n");
  
  // data
  for (int j = ny-1; j >= 0; j--) {
    double yp = Delta*j + p.box[0][1] + Delta/2.;
    for (int i = 0; i < nx; i++) {
      double xp = Delta*i + p.box[0][0] + Delta/2., v;
      if (p.mask.i) { // masking
	if (p.linear) {
	  double m = interpolate (p.mask, xp, yp);
	  if (m < 0.)
	    v = nodata;
	  else
	    v = interpolate (p.f, xp, yp);
	}
	else {
	  Point point = locate (xp, yp);
	  if (point.level < 0 || val(p.mask) < 0.)
	    v = nodata;
	  else
	    v = val(p.f);
	}
      }
      else if (p.linear)
	v = interpolate (p.f, xp, yp);
      else {
	Point point = locate (xp, yp);
	v = point.level >= 0 ? val(p.f) : nodata;
      }
      if (v == nodata)
	fprintf (p.fp, "-9999 ");
      else
	fprintf (p.fp, "%f ", v);
    }
    fprintf (p.fp, "\n");
  }

  fflush (p.fp);
}

#if MULTIGRID

/**
## *output_gfs()*: Gerris simulation format

The function writes simulation data in the format used in
[Gerris](http://gerris.dalembert.upmc.fr) simulation files. These files can be read
with GfsView.

The arguments and their default values are:

*fp*
: a file pointer. Default is *name* or stdout.

*list*
: a list of scalar fields to write. Default is *all*. 

*t*
: the physical time. Default is zero. 

*file*
: the name of the file to write to (mutually exclusive with *fp*).

*translate*
: whether to replace "well-known" Basilisk variables with their Gerris
equivalents.
*/

struct OutputGfs {
  FILE * fp;
  scalar * list;
  double t;
  char * file;
  bool translate;
};

static char * replace (const char * input, int target, int with,
		       bool translate)
{
  if (translate) {
    if (!strcmp (input, "u.x"))
      return strdup ("U");
    if (!strcmp (input, "u.y"))
      return strdup ("V");
    if (!strcmp (input, "u.z"))
      return strdup ("W");
  }
  char * name = strdup (input), * i = name;
  while (*i != '\0') {
    if (*i == target)
      *i = with;
    i++;
  }
  return name;
}

trace
void output_gfs (struct OutputGfs p)
{
@if _QCCACC
  /* Download data from device for output */
  ACC_UPDATE_HOST
@endif

  bool opened = false;
  if (p.fp == NULL) {
    if (p.file == NULL)
      p.fp = stdout;
    else if (!(p.fp = fopen (p.file, "w"))) {
      perror (p.file);
      exit (1);
    }
    else
      opened = true;
  }
  scalar * list = p.list ? p.list : list_copy (all);

  fprintf (p.fp, 
	   "1 0 GfsSimulation GfsBox GfsGEdge { binary = 1"
	   " x = %g y = %g ",
	   0.5 + X0/L0, 0.5 + Y0/L0);
#if dimension == 3
  fprintf (p.fp, "z = %g ", 0.5 + Z0/L0);
#endif

  if (list != NULL && list[0].i != -1) {
    scalar s = list[0];
    char * name = replace (s.name, '.', '_', p.translate);
    fprintf (p.fp, "variables = %s", name);
    free (name);
    for (int i = 1; i < list_len(list); i++) {
      scalar s = list[i];
      if (s.name) {
	char * name = replace (s.name, '.', '_', p.translate);
	fprintf (p.fp, ",%s", name);
	free (name);
      }
    }
    fprintf (p.fp, " ");
  }
  fprintf (p.fp, "} {\n");
  if (p.t > 0.)
    fprintf (p.fp, "  Time { t = %g }\n", p.t);
  if (L0 != 1.)
    fprintf (p.fp, "  PhysicalParams { L = %g }\n", L0);
  fprintf (p.fp, "}\nGfsBox { x = 0 y = 0 z = 0 } {\n");

@if _MPI
  long header;
  if ((header = ftell (p.fp)) < 0) {
    perror ("output_gfs(): error in header");
    exit (1);
  }
  int cell_size = sizeof(unsigned) + sizeof(double);
  for (scalar s in list)
    if (s.name)
      cell_size += sizeof(double);
  scalar index = new scalar;
  size_t total_size = header + (z_indexing (index, false) + 1)*cell_size;
@endif
  
  // see gerris/ftt.c:ftt_cell_write()
  //     gerris/domain.c:gfs_cell_write()
  foreach_cell() {
@if _MPI // fixme: this won't work when combining MPI and mask()
    if (is_local(cell))
@endif
    {
@if _MPI
      if (fseek (p.fp, header + index[]*cell_size, SEEK_SET) < 0) {
	perror ("output_gfs(): error while seeking");
	exit (1);
      }
@endif
      unsigned flags = 
	level == 0 ? 0 :
#if dimension == 1
	child.x == 1;
#elif dimension == 2
      child.x == -1 && child.y == -1 ? 0 :
	child.x == -1 && child.y ==  1 ? 1 :
	child.x ==  1 && child.y == -1 ? 2 : 
	3;
#else // dimension == 3
      child.x == -1 && child.y == -1 && child.z == -1  ? 0 :
	child.x == -1 && child.y == -1 && child.z ==  1  ? 1 :
	child.x == -1 && child.y ==  1 && child.z == -1  ? 2 : 
	child.x == -1 && child.y ==  1 && child.z ==  1  ? 3 : 
	child.x ==  1 && child.y == -1 && child.z == -1 ? 4 :
	child.x ==  1 && child.y == -1 && child.z ==  1 ? 5 :
	child.x ==  1 && child.y ==  1 && child.z == -1 ? 6 : 
	7;
#endif
      if (is_leaf(cell))
	flags |= (1 << 4);
      fwrite (&flags, sizeof (unsigned), 1, p.fp);
      double a = -1;
      fwrite (&a, sizeof (double), 1, p.fp);
      for (scalar s in list)
	if (s.name) {
	  if (s.v.x.i >= 0) {
	    // this is a vector component, we need to rotate from
	    // N-ordering (Basilisk) to Z-ordering (Gerris)
#if dimension >= 2
	    if (s.v.x.i == s.i) {
	      s = s.v.y;
	      a = is_local(cell) && !isnan(s[]) && s[] != nodata ? s[] : DBL_MAX;
	    }
	    else if (s.v.y.i == s.i) {
	      s = s.v.x;
	      a = is_local(cell) && !isnan(s[]) && s[] != nodata ?
		- s[] : DBL_MAX;
	    }
#endif
#if dimension >= 3
	    else
	      a = is_local(cell) && !isnan(s[]) && s[] != nodata ? s[] : DBL_MAX;
#endif
	  }
	  else
	    a = is_local(cell) && !isnan(s[]) && s[] != nodata ? s[] : DBL_MAX;
	  fwrite (&a, sizeof (double), 1, p.fp);
	}
    }
    if (is_leaf(cell))
      continue;
  }
  
@if _MPI
  delete ({index});
  if (!pid() && fseek (p.fp, total_size, SEEK_SET) < 0) {
    perror ("output_gfs(): error while finishing");
    exit (1);
  }
  if (!pid())
@endif  
    fputs ("}\n", p.fp);
  fflush (p.fp);

  if (!p.list)
    free (list);
  if (opened)
    fclose (p.fp);
}

/**
## *dump()*: Basilisk snapshots

This function (together with *restore()*) can be used to dump/restore
entire simulations.

The arguments and their default values are:

*fp*
: a file pointer. Default is *name* or stdout.

*list*
: a list of scalar fields to write. Default is *all*. 

*t*
: the physical time. Default is zero. 

*file*
: the name of the file to write to (mutually exclusive with *fp*).
*/

struct Dump {
  FILE * fp;
  scalar * list;
  double t;
  char * file;
};

struct DumpHeader {
  double t;
  long len;
  int depth;
};

trace
void dump (struct Dump p)
{
@if _QCCACC
  /* Download data from device for output */
  ACC_UPDATE_HOST 
@endif

  FILE * fp = p.fp;
  scalar * lista = p.list ? p.list : all, * list = NULL;
  char * file = p.file;

  for (scalar s in lista)
    if (!s.face && s.i != cm.i)
      list = list_add (list, s);
  
  if (file && (fp = fopen (file, "w")) == NULL) {
    perror (file);
    exit (1);
  }
  assert (fp);

  struct DumpHeader header = { p.t, list_len(list), depth() };

  if (pid() == 0 && fwrite (&header, sizeof(header), 1, fp) < 1) {
    perror ("dump(): error while writing header");
    exit (1);
  }

  scalar index = {-1};
  
@if _MPI // Parallel
  index = new scalar;
  z_indexing (index, false);
  int cell_size = sizeof(unsigned) + header.len*sizeof(double);
@endif
  
  foreach_cell() {
@if _MPI // fixme: this won't work when combining MPI and mask()
    if (is_local(cell))
@endif
    {
@if _MPI
      if (fseek (fp, sizeof(header) + index[]*cell_size, SEEK_SET) < 0) {
	perror ("dump(): error while seeking");
	exit (1);
      }
@endif
      unsigned flags = is_leaf(cell) ? leaf : 0;
      if (fwrite (&flags, sizeof(unsigned), 1, fp) < 1) {
	perror ("dump(): error while writing flags");
	exit (1);
      }
      for (scalar s in list)
	if (fwrite (&s[], sizeof(double), 1, fp) < 1) {
	  perror ("dump(): error while writing scalars");
	  exit (1);
	}
    }
    if (is_leaf(cell))
      continue;
  }

  delete ({index});
  
  free (list);
  if (file)
    fclose (fp);
}

trace
bool restore (struct Dump p)
{
  FILE * fp = p.fp;
  scalar * lista = p.list ? p.list : all, * list = NULL;
  char * file = p.file;
  if (file && (fp = fopen (file, "r")) == NULL)
    return false;
  assert (fp);

  for (scalar s in lista)
    if (!s.face && s.i != cm.i)
      list = list_add (list, s);

  struct DumpHeader header;
  
  double t = 0.;
  if (fread (&header, sizeof(header), 1, fp) < 1) {
    fprintf (stderr, "restore(): error: expecting header\n");
    exit (1);
  }
  if (header.len != list_len (list)) {
    fprintf (stderr,
	     "restore(): error: the list lengths don't match: %ld != %d\n",
	     header.len, list_len (list));
    exit (1);
  }

#if TREE
  init_grid (1);

  foreach_cell() {
    cell.pid = pid();
    cell.flags |= active;
  }
  tree->dirty = true;
#else
  init_grid (1 << header.depth);
#endif

@if _QCCACC
  /* Download data from device for input */
  ACC_UPDATE_HOST
@endif

  foreach_cell() {
    unsigned flags;
    if (fread (&flags, sizeof(unsigned), 1, fp) != 1) {
      fprintf (stderr, "restore(): error: expecting 'flags'\n");
      exit (1);
    }
    for (scalar s in list)
      if (fread (&s[], sizeof(double), 1, fp) != 1) {
	fprintf (stderr, "restore(): error: expecting '%s'\n", s.name);
	exit (1);
      }
    if (!(flags & leaf) && is_leaf(cell))
      refine_cell (point, NULL, 0, NULL);
    if (is_leaf(cell))
      continue;
  }

@if _QCCACC
  /* Upload data to device */
  ACC_UPDATE_DEVICE
@endif

  boundary (list);
  if (file)
    fclose (fp);

  // the events are advanced to catch up with the time
  double t1 = 0.;
  while (t1 < t && events (0, t1, false))
    t1 = tnext;

  free (list);
  return true;
}

#endif // MULTIGRID

    if (!(flags & leaf) && is_leaf(cell))
      refine_cell (point, NULL, 0, NULL);
    if (is_leaf(cell))
      continue;
  }

@if _QCCACC
  /* Upload data to device */
  ACC_UPDATE_DEVICE
@endif

  boundary (list);
  if (file)
    fclose (fp);

  // the events are advanced to catch up with the time
  double t1 = 0.;
  while (t1 < t && events (0, t1, false))
    t1 = tnext;

  free (list);
  return true;
}

#endif // MULTIGRID
