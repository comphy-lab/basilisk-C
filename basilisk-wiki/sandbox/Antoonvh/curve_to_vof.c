/**
# Find particle-based interface in cells.

An interface defined by a particle loop crosses some cells. We wish to
find which cells are cut by which section of the loop. We achieve
this using multigrid acceleration. 

~~~gnuplot Paricle index in blue and cell data
set term svg size 800, 800
set size ratio -1
set key off
dx = .03
plot 'out' w l t 'cells' ,				\
'log' u ($1-dx):($2+dx):3 with labels t 'cell data',	\
'log' u ($1):($2+dx):4 with labels t 'cell data',	\
'log' u ($1+dx):($2+dx):5 with labels t 'cell data',	\
'log' u ($1-dx):($2-dx):6 with labels t 'cell data',	\
'log' u ($1):($2-dx):7 with labels t 'cell data',	\
'log' u ($1+dx):($2-dx):8 with labels t 'cell data',	\
'particles' with labels textcolor 'blue' font ",15" ,\
'particles' w l lc 'blue'
~~~
*/

#include "particle.h"

typedef union {
  int * p;
  double d;
} datic;

int np = 11;

// Print data reference by data field (maybe multiple pairs)
void print_data (scalar data) {
  foreach() {
    if (data[]) {
      fprintf (stderr, "%g %g ", x, y);
      int * p = (datic){.d = data[]}.p;
      int i = 0;
      while (p[i] >= 0) {
	fprintf (stderr, "%d - %d ", p[i] , ((p[i] + 1)%np));
	i++;
      }
      fprintf (stderr, "\n");
    }
  }
}

// Compute intersection
bool segment_cell (coord a, coord b, Point point);
// Transfer indices to `data` reference field
void loop_to_grid (Particles interface, scalar data);

int main() {
  assert (sizeof (int *) == sizeof (double));
  L0 = 2.5;
  X0 = Y0 = -L0/2.0123;
  init_grid (16);
  output_cells();

  // Initialize particles
  Particles loop = new_particles (np);
  FILE * fp = fopen ("particles", "w");
  foreach_particle() {
    p().x = sin(2.*pi*j/np);
    p().y = cos(2.*pi*j/np);
    fprintf (fp, "%g %g %d\n", p().x, p().y, j);
  }
  fclose (fp);
  // declare and initialze data field
  scalar data[];
  foreach()
    data[] = 0.;
  multigrid_restriction ({data});
  loop_to_grid (loop, data);
  print_data (data);
  //cleanup
  free_p();
  foreach_cell() {
    free ((datic){.d = data[]}.p);
  }
}


// segment_cell intersection
bool segment_cell (coord a, coord b, Point point) { // From foreach_segment()
  coord t = {b.x - a.x, b.y - a.y};
  double norm = sqrt(sq(t.x) + sq(t.y));
  assert (norm > 0.);
  t.x = t.x/norm + 1e-6, t.y = t.y/norm - 1.5e-6;
  double alpha = (a.x*(b.y - a.y) -
		  a.y*(b.x - a.x))/norm;
  if (fabs(t.y*x - t.x*y - alpha) < 0.708*Delta) {
    coord o = {x,y}, p[2];
    int n = 0;
    foreach_dimension()
      if (t.x)
	for (int _i = -1; _i <= 1 && n < 2; _i += 2) {
	  p[n].x = o.x + _i*Delta/2.;
	  double ap = (p[n].x - a.x)/t.x;
	  p[n].y = a.y + ap*t.y;
	  if (fabs(p[n].y - o.y) <= Delta/2.) {
	    ap = clamp (ap, 0., norm);
	    p[n].x = a.x + ap*t.x, p[n].y = a.y + ap*t.y;
	    if ((2 - 1e-6)*fabs(p[n].x - o.x) <= Delta &&
		(2 - 1e-6)*fabs(p[n].y - o.y) <= Delta)
	      n++;
	  }
	}
    if (n == 2)
      return true;
  }
  return false;
}


void loop_to_grid (Particles interface, scalar data) {
  int np = pn[interface];
  // MG accelerate
  for (int _l = 0; _l <= depth(); _l++) {
    foreach_level(_l) {
      if (level == 0) {
	int n = 0;
	for (int i = 0; i < np; i++) {
	  coord a1 = {pl[interface][i].x, pl[interface][i].y};
	  coord a2 = {pl[interface][(i + 1)%np].x, pl[interface][(i + 1)%np].y};
	  if (segment_cell (a1, a2, point)) 
	    n++;
	}
	int * indi = malloc(sizeof(int)*(n + 1));
	int ind = 0;
	for (int i = 0; i < np; i++) {
	  coord a1 = {pl[interface][i].x, pl[interface][i].y};
	  coord a2 = {pl[interface][(i + 1)%np].x, pl[interface][(i + 1)%np].y};
	  if (segment_cell (a1, a2, point)) {
	    indi[ind++] = i; 
	  }
	}
	indi[ind] = -1;
	data[] = (datic){.p = indi}.d;
	
      } else {// (level > 0, check parent)
	if (coarse(data,0,0,0)) {
	  int * par = (datic){.d = coarse(data,0,0,0)}.p;
	  int i = 0, n = 0;
	  while (par[i] >= 0) {
	    coord a1 = {pl[interface][par[i]].x, pl[interface][par[i]].y};
	    coord a2 = {pl[interface][(par[i] + 1)%np].x, pl[interface][(par[i] + 1)%np].y};
	    if (segment_cell (a1, a2, point)) 
	      n++;
	    i++;
	  }
	  if (n) {
	    // printf ("%d %d\n", level, n);
	    int * indi = malloc(sizeof(int)*(n + 1));
	    int ind = 0;
	    i = 0;
	    while (par[i] >= 0) {
	      coord a1 = {pl[interface][par[i]].x, pl[interface][par[i]].y};
	      coord a2 = {pl[interface][(par[i] + 1)%np].x, pl[interface][(par[i] + 1)%np].y};
	      if (segment_cell (a1, a2, point)) 
		indi[ind++] = par[i]; 
	      i++;
	    }
	    indi[ind] = -1;
	    data[] = (datic){.p = indi}.d;
	  }
	}
      }
    }
  }
}

