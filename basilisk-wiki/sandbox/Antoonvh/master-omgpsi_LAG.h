/**
# Dual grid $\omega-\psi$ solver

This is used in combination with [slave-omgpsi_LAG.c]().

The master is concerned with time integration of the vorticity ($\omega$) field.
*/
#pragma autolink slave-omgpsi-LAG.o
#define ADD_PART_MEM double s; coord u; double ti;
#include "fpm.h"
#include "poisson.h"
#include "run.h"
#include "higher-order.h"
// Global fields
scalar omega[];

Particles p_omg;
#include "utils.h"
// A function for copying the vorticity field to the slave
Point locater (double xp = 0., double yp = 0., double zp = 0.)
{
  for (int l = depth(); l >= 0; l--) {
    Point point = {0};
    point.level = l;
    int n = 1 << point.level;
    point.i = (xp - X0)/L0*n + GHOSTS;
#if dimension >= 2
    point.j = (yp - Y0)/L0*n + GHOSTS;
#endif
#if dimension >= 3
    point.k = (zp - Z0)/L0*n + GHOSTS;
#endif
    if (point.i >= 0 && point.i < n + 2*GHOSTS
#if dimension >= 2
	&& point.j >= 0 && point.j < n + 2*GHOSTS
#endif
#if dimension >= 3
	&& point.k >= 0 && point.k < n + 2*GHOSTS
#endif
	) {
      if (allocated(0) && is_local(cell) && is_leaf(cell))
	return point;
    }
    else
      break;
  }
  Point point = {0};
  point.level = -1;
  return point;
}

double master_value_p_omg (double xp, double yp, int lev) {
  coord X = {xp, yp, 0};
  double coefs[10] = {0};
  double SCALE = L0/(1 << lev);
  least_squares_poly (X, coefs, p_omg);
  return coefs[0] + sq(SCALE)*(coefs[4] + coefs[5])/12.;
}

double master_value (const char * name, double xp, double yp, int i, int j, int lev)
{
  if (!grid) {
    fprintf (stderr, "slave_interpolate: error: no grid! this may be a "
	     "master/slave synchronization issue\n");
    exit (1);
  }
  scalar s = lookup_field (name);
  if (s.i < 0) {
    fprintf (stderr, "slave_interpolate: error: unknown field '%s'\n", name);
    exit (1);
  }
  Point point = {.i = i, .j = j, .level = lev}; // foreach_point() is not an alternative...
  Point p1 = locater (xp, yp);
  double val = 0;
  if (p1.level >= lev) {
    
    val = s[]; // more accurate and faster
  }
  else { // ?
    val = interpolate (s, xp, yp, 0, true);
  }
  return val;
}

// A function prototype for interpolating the stream function from the slave
extern double slave_interpolate(const char * name, double xp = 0, double yp = 0, double zp = 0, bool linear);

trace
void copy_field_slave (scalar s){
  foreach() {
    double zp = 0;
    if (dimension == 3)
      zp = z;
    s[] = slave_interpolate (s.name, x, y, zp, true);
  }
}

// Slave handling
extern void slave_solve_psi();
extern void slave_level();
extern int slave_init();
extern void slave_free();

#ifndef RKORDER
#define RKORDER (4)
#endif
#if (RKORDER == 3)
// Williamson, J. H.: Low-Storage Runge-Kutta schemes, J.
// Comput.Phys., 35, 48–56, 1980.
#define STAGES (3)
double An[STAGES] = {0., -5./9., -153./128.};
double Bn[STAGES] = {1./3., 15./16., 8./15.};
#else
// Carpenter, M.  H.  and Kennedy, C.  A.: Fourth-order
// 2N-storageRunge-Kutta schemes, Tech. Rep. TM-109112, NASA
// LangleyResearch Center, 1994
#define STAGES (5)
double An[STAGES] = {0.,
		     -567301805773. /1357537059087.,
		     -2404267990393./2016746695238.,
		     -3550918686646./2091501179385.,
		     -1275806237668./842570457699.};
double Bn[STAGES] = {1432997174477./9575080441755. ,
		     5161836677717./13612068292357.,
		     1720146321549./2090206949498. ,
		     3134564353537./4481467310338. ,
		     2277821191437./14882151754819.};
#endif

trace
void compute_omega (Particles p, scalar omg) {
  foreach() {
    coord X = {x,y,z};
    double coefs[10] = {0};
    least_squares_poly (X, coefs, p);
    omg[] = coefs[0] + sq(Delta)*(coefs[4] + coefs[5])/12.;
  }
  //"master_value" does not trigger automatic BC
  multigrid_restriction({omg});
  boundary({omg});
}

void advance_rk (Particles p, double dt) {
  for (int Stp = 0; Stp < STAGES; Stp++) {
    compute_omega (p_omg, omega);
    slave_solve_psi();
    foreach_particle_in(p) {
      p().u.x =  An[Stp]*p().u.x - slave_interpolate("uf.y", x, y, 0, true); 
      p().u.y =  An[Stp]*p().u.y + slave_interpolate("uf.x", x, y, 0 ,true); 
      foreach_dimension() 
	p().x += dt*Bn[Stp]*p().u.x;
    }
    ref_outdated = true;
  }
  particle_boundary(p);
}

event timestep (i++, last) {
  dt = dtnext (DT);
}


event advance (i++, last) {
  advance_rk(p_omg, dt);
}

event defaults (i = 0) {
  max_particles = 15;
  
  //2nd order suffices
  order_barrier[1] = 12;
  order_barrier[2] = 20;
}

event init (t = 0);

event call_timestep (t = 0) {
  event ("timestep"); 
}



/**
We let the Poisson-solver grid know when the simulation is finished.
*/

event stop (t = end) {
  slave_free();
}

#ifdef VIEW
void scatter_color (Particles p, float s = 3, Colormap map = jet, double minv = -1, double maxv = 1) {
  double cmap[NCMAP][3];
  (*map)(cmap);
  bview * view = draw();
#if dimension == 2
  glTranslatef (0., 0., view->lc*view->fov/24.); //from draw_lines()
#endif
  glPointSize(s*view->samples);
  glBegin (GL_POINTS);
  foreach_particle_in(p) {
    Color b = colormap_color (cmap, p().s, minv, maxv);		
    glColor3f (b.r/255., b.g/255., b.b/255.);			
    
#if dimension == 2    
    glvertex2d (view, x, y);
#else // dimension == 3
    glvertex3d (view, x, y, z);
#endif
  }
  glEnd();
  view->ni++; 
}
#endif


int no_more_particles_than (int maxp) {
  int rm = 0;

  if (ref_outdated) {
    free_scalar_data(reference);
    assign_particles(p_omg, reference);
    ref_outdated = false;
  }
  
  foreach(reduction(+:rm)) {
    int pn = 0;
    foreach_particle_point(reference, point)
      pn++;
    if (pn > maxp) {
      double err_max = 0;
      foreach_particle_point(reference, point) {
	coord X = {p().x, p().y};
	double coefs[10] = {1};
	least_squares_poly (X, coefs, p_omg, self = false);
	if (fabs(p().s - coefs[0]) > err_max)
	  err_max = fabs(p().s - coefs[0]);
      }
      {
	foreach_particle_point(reference, point) {
	  coord X = {p().x, p().y};
	  double coefs[10] = {1};
	  least_squares_poly (X, coefs, p_omg, self = false);
	  if (fabs(p().s - coefs[0]) == err_max) {
	    p().z = HUGE;
	    rm++;
	  }
	}
      }
    }
  }
  if (rm) {
    remove_particles (p_omg, p().z == HUGE);
    ref_outdated = true;
  }
  return rm;
}


int no_less_particles_than (int minp, int order_req = 2) {
  int ap = 0;
   
  if (ref_outdated) {
    free_scalar_data(reference);
    assign_particles(p_omg, reference);
    ref_outdated = false;
  }
  
  foreach(reduction(+:ap)) {
    int pnr = 0;
    foreach_particle_point(reference, point)
      pnr++;
    if (pnr < minp) {
      coord cc = {x, y, z};
      particle pn;
      // packing problem
      if (pnr <= 0) { // cell centre
	pn.x  = cc.x;
	pn.y = cc.y;
      } else { // center of oposing quadrant of "last" particle (? minp = 2)
	particle pc;
	foreach_particle_point(reference, point)
	  pc = p();
	foreach_dimension()
	  pn.x = cc.x + Delta/4.*(pc.x > cc.x ? -1 : 1);
      }
      coord X = {pn.x, pn.y};
      double coefs[10] = {1};
      int j = least_squares_poly (X, coefs, p_omg);
      if (j >= order_req)
	pn.s = coefs[0];
      else 
	pn.s = interpolate (omega, pn.x, pn.y);
      //pn.s = 2.*sin(2.*pi*pn.x)*sin(2.*pi*pn.y);
      pn.ti = t;
      add_particle (pn, p_omg);
      ap++;
    }
  }
  if (ap)
    ref_outdated = true;
  return ap;
}

int adapt_number(int minp = 2, int maxp = 5) {
  boundary({reference});
  int rem = no_more_particles_than(maxp);
  boundary({reference});
  int add = no_less_particles_than(minp);
  return (add - rem);
}
