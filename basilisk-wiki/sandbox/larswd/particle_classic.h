/**
Particles can be removed based on their position or another condition
 */

int remove_particle (particle p, Particles Plist) {
  long unsigned int * ind = malloc (sizeof(long int));
  int i = 0, a = 1; 
  foreach_particle_in(Plist) {
    foreach_dimension()
      if (p().x != p.x)
	      continue;
    if (i + 1 >= a) {
      ind = realloc (ind, a*2*sizeof(long int));
      a *= 2;
    }
    ind[i++] = _j_particle;
  }
  remove_particles_index (Plist, i, ind);
  free (ind); ind = NULL;
  return i;
}

/**
## Utility functions

### boundary conditions

peridoc box and `_MPI` boundaries.

Periodicity in vertical does not make sense
in the multilayer solver, and hence this function remains currently unchanged. 
 */

void particle_boundary (Particles P) {
  coord mind = {X0, Y0, Z0};
  foreach_particle_in(P) {
    foreach_dimension() {
      if (p().x < mind.x) 
	      p().x += L0;
      else if (p().x > (mind.x + L0))
	      p().x -= L0;
    }
  }
#if _MPI
  update_mpi (P);
#endif
}


/**
### Place `particle`s

For already allocated `particle`s referenced by `Particles P`:

The foreach() and foreach_dimension()-iterators do not see z-axis in multilayer,
and hence the support for this axis must be added manually where relevant. 
*/


void place_in_cells (Particles P) {
  particles pp = pl[P];                 
  long unsigned int np = 0;
  foreach(serial, reduction(+:np)) {
    coord loc = {x, y,z};
    foreach_dimension()
      pp[np].x = loc.x;
    np++;
  }
}

/**
### Simple particles statistics

The function computes, Average location, min, max and standard
deviation vectors. The correct statistics are only available for
`pid() == 0`.
 */


pstats statsp (Particles P) {
  coord avg = {0 ,0, 0}, min = {0,0,0},
    max = {0, 0, 0}, stddev = {0, 0, 0};
  foreach_dimension(3) {
    avg.x = stddev.x = 0;
    min.x = HUGE;
    max.x = -HUGE;
  }
  long unsigned int np = pn[P];
#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, &np, 1, MPI_UNSIGNED_LONG,
		 MPI_SUM, MPI_COMM_WORLD);
#endif
  if (np) {
    foreach_dimension() { //reduction of coord members does not work
      double avgx = 0;
      double minx = HUGE;
      double maxx = -HUGE;
      foreach_particle_in(P, reduction(+:avgx) reduction(max:maxx)
			  reduction (min:minx)) {
	  avgx += p().x;
	  if (p().x < minx)
	    minx = p().x;
	  if (p().x > maxx)
	    maxx = p().x;
      }
      avg.x = avgx/(double)np;
      min.x = minx;
      max.x = maxx;
    }
    foreach_dimension() {
      double stddevx = 0;
      foreach_particle_in(P, reduction(+:stddevx)) {
	stddevx += sq(avg.x - p().x);
      }
      stddev.x = sqrt(stddevx/(double)np);
    }
  }
  
  pstats s;
  s.max = max, s.min = min, s.avg = avg, s.stddev = stddev;
  return s;
}


/**
### Propability density function

Obtain a scalar PDF from particle locations 
 */

void particle_pdf (Particles P, scalar s) {
  foreach()
    s[] = 0;
  particle_boundary (P);
  foreach_particle_in(P) {
    Point point = locate (x,y,z);
    if (point.level > 0)
      s[]++;
  }
  foreach()
    s[] /= (pn[P]*dv());
  boundary ({s});
}


/**
### Random step

Particles displace a certain distance (`step`) in a random direction
*/

void random_step (Particles P, double step) {
  foreach_particle_in(P) {
    double theta = noise()*pi;
    #if (dimension == 1)
      coord d = {sign(theta)};
    #elif (dimension < 3)
      coord d = {sin(theta), cos(theta)};
    #else
      double phi = acos(noise());
      coord d = {sin(phi)*cos(theta), sin(phi)*sin(theta), cos(phi)};
    #endif
    foreach_dimension()
      p().x += d.x*step;
  }
}


/**
### Length of particle loop

A function that computes the line length of particles places in a loop
(mind the ordering). Special care is taken to obtain 4th order
accuracy.
 */
double plength (Particles P) {
  double lr = 0;
  int np = pn[P];
  particles pli = pl[P];
  foreach_particle_in(P) {
    int il = _j_particle - 1 < 0 ? np - 1: _j_particle - 1;
    int ir = _j_particle + 1 >= np ? 0 : _j_particle + 1;
    coord pl = (coord){pli[il].x, pli[il].y, pli[il].z};
    coord pm = (coord){pli[_j_particle].x, pli[_j_particle].y, pli[_j_particle].z};
    coord pr = (coord){pli[ir].x, pli[ir].y, pli[ir].z};
    coord bv = {0, 0, 0}, a1 = {0,0,0};
    coord a2 = {0,0,0};
    double bl = 0, ka = 0;
    foreach_dimension() {
      a1.x = pm.x - pl.x;
      bv.x = pr.x - pl.x;
      bl += sq(bv.x);
    }

    bl = sqrt(bl);
    foreach_dimension()
      ka += a1.x*bv.x/bl;

    foreach_dimension()
      a2.x = a1.x - bv.x*ka/bl;

    normalize (&a2);
    double al = 0, am = 0, ar = 0;
    foreach_dimension() {
      al -= a2.x*pl.x;
      am -= a2.x*pm.x;
      ar -= a2.x*pr.x;
    }
    double C = fabs(am - al);
    double xt = 0;
    double b = 0;
    foreach_dimension() {
      xt += sq(pl.x - pm.x);
      b += sq(pl.x - pr.x);
    }
    xt = sqrt (xt - sq(C));
    b = sqrt(b);
    double A =  C/((xt*(b - xt)));
    double xp1 = 0.5*b*(1 - 1/sqrt(3));
    double xp2 = 0.5*b*(1 + 1/sqrt(3));
    double l1 = b*sqrt(1. + sq(-2*A*xp1 + A*b));
    double l2 = b*sqrt(1. + sq(-2*A*xp2 + A*b));
    lr += (l1 + l2)/4;
  }
  return lr;
}

