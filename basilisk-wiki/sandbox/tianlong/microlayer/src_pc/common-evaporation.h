/**
# Common Evaporation

We define tolerances and functions useful
for evaporation models. */


#include "mycurvature.h"
#include "adapt_wavelet_leave_interface.h"

#ifndef F_ERR
# define F_ERR 1.e-10
#endif
#ifndef T_ERR
# define T_ERR 0.
#endif
#include "mass_refine_prolongation.h"

#ifndef I_TOL
# define I_TOL 1.e-8
#endif

/**
We define a custom for that simplifies
loops over the chemical species. */

#define foreach_elem(list, index) \
for (int index=0; index<list_len (list); index++)

/**
## Compute the normal in an interfacial cell
*/

coord normal (Point point, scalar c)
{
  coord n = interface_normal(point, c);
  double nn = 0.;
  foreach_dimension()
    nn += sq(n.x);
  nn = sqrt(nn);
  foreach_dimension()
    n.x /= nn;
  return n;
}

/**
## *mass2molefrac()*: Compute mole fractions from mass fractions

* *X*: vector filled with mole fractions
* *W*: vector with the mass fractions
* *MW*: vector with the molecular weights of each species
* *NS*: total number of species (vectors length)
*/

void mass2molefrac (double * X, const double * W, const double * MW, const int NS)
{
  double MWmix = 0.;
  for (int i=0; i<NS; i++) {
    MWmix += W[i]/MW[i];
  }
  for (int i=0; i<NS; i++) {
    X[i] = W[i]/MW[i]/(MWmix + 1.e-10);
  }
}

/**
## *mole2massfrac()*: Compute mass fractions from mole fractions

* *W*: vector filled with mole fractions
* *X*: vector with the mass fractions
* *MW*: vector with the molecular weights of each species
* *NS*: total number of species (vectors length)
*/

void mole2massfrac (double * W, const double * X, const double * MW, const int NS)
{
  double MWmix = 0.;
  for (int i=0; i<NS; i++) {
    MWmix += X[i]*MW[i];
  }
  for (int i=0; i<NS; i++) {
    W[i] = X[i]*MW[i]/(MWmix + 1.e-10);
  }
}

/**
## *mass2mw()*: Compute mixture molecular weight from mass fractions

* *W*: vector with the mass fractions
* *MW*: vector with the molecular weights of each species
* *NS*: total number of species (vectors length)
*/

double mass2mw (const double * W, const double * MW, const int NS)
{
  double MWmix = 0.;
  for (int i=0; i<NS; i++) {
    MWmix += W[i]/MW[i];
  }
  return 1./(MWmix + 1.e-10);
}

/**
## *mole2mw()*: Compute mixture molecular weight from mole fractions

* *X*: vector with the mass fractions
* *MW*: vector with the molecular weights of each species
* *NS*: total number of species (vectors length)
*/

double mole2mw (const double * X, const double * MW, const int NS)
{
  double MWmix = 0.;
  for (int i=0; i<NS; i++) {
    MWmix += X[i]*MW[i];
  }
  return MWmix;
}

/**
## *avg_neighbor()*: Compute the average value of a scalar field Y in a 3x3 stencil around the current cell
* *point*: current cell location
* *Y*: field to average
* *f*: vof volume fraction field
*/

double avg_neighbor (Point point, scalar Y, scalar f) {
  double fYnei = 0., fnei = 0.;
  foreach_neighbor(1) {
    fYnei += f[]*Y[];
    fnei  += f[];
  }
  return fYnei / fnei;
}

/**
## *avg_interface()*: Compute the average value of a scalar field Y in a region around the interface
* *Y*: fields to average
* *f*: vof volume fraction field
*/

double avg_interface (scalar Y, scalar f) {
  double Yavg = 0.;
  int counter = 0;
  foreach(reduction(+:Yavg) reduction(+:counter)) {
    if (f[] > F_ERR && f[] < 1.-F_ERR) {
      counter++;
      Yavg += Y[];
    }
  }
  return Yavg / counter;
}

/**
## *smooth_field()*: Smooth a discontinuous field, from sf calculation in two-phase.h

* *sf*: smoothed scalar field
* *f*: initial scalar field
* *ncycles*: number of smoothing cycles
*/

//void smooth_field (scalar sf, scalar f, int ncycles) {
//  for (int i=0; i<ncycles; i++) {
//#if dimension <= 2
//    foreach()
//      sf[] = (4.*f[] + 
//        2.*(f[0,1] + f[0,-1] + f[1,0] + f[-1,0]) +
//        f[-1,-1] + f[1,-1] + f[1,1] + f[-1,1])/16.;
//#else // dimension == 3
//    foreach()
//      sf[] = (8.*f[] +
//        4.*(f[-1] + f[1] + f[0,1] + f[0,-1] + f[0,0,1] + f[0,0,-1]) +
//        2.*(f[-1,1] + f[-1,0,1] + f[-1,0,-1] + f[-1,-1] + 
//      f[0,1,1] + f[0,1,-1] + f[0,-1,1] + f[0,-1,-1] +
//      f[1,1] + f[1,0,1] + f[1,-1] + f[1,0,-1]) +
//        f[1,-1,1] + f[-1,1,1] + f[-1,1,-1] + f[1,1,1] +
//        f[1,1,-1] + f[-1,-1,-1] + f[1,-1,-1] + f[-1,-1,1])/64.;
//#endif
//    foreach()
//      f[] = sf[];
//  }
//}

/**
## *shift_field()*: Shift a field localized at the interface toward the closest pure gas or liquid cells

* *fts*: field to shift
* *f*: vof volume fraction field
* *dir*: shifting direction: 1 liquid, gas otherwise
*/

void shift_field (scalar fts, scalar f, int dir) {

  scalar avg[];
#if TREE
  avg.refine = avg.prolongation = refinement_avg;
  avg.restriction = no_restriction;
  avg.dirty = true;
#endif

#if USE_MY_SOLID
  // Compute avg
  foreach() {
    avg[] = 0.;
    if (is_solid[] == 0 && f[] > F_ERR && f[] < 1. - F_ERR) {
      if (dir == 1) {
        int count = 0;
        foreach_neighbor (1) {
          if (is_solid[] == 0 && f[] > 1.-F_ERR) // Number of pure-liquid cells close to the interfacial cell
            count ++;
        }
        avg[] = count;
      }
      else {
        int count = 0;
        foreach_neighbor (1) {
          if (is_solid [] == 0 && f[] < F_ERR) // Number of pure-gas cells close to the interfacial cell
            count ++;
        }
        avg[] = count;
      }
    }
  }

  scalar sf0[];
  foreach() {
    sf0[] = fts[];
    fts[] = 0.;
  }
  boundary({sf0,fts});

  // Compute m
  foreach() {
    if (dir == 1) {
      if (is_solid[] == 0 && f[] > 1.-F_ERR) { // Move toward pure-liquid
        double val = 0.;
        foreach_neighbor (1) {
          if (is_solid[] == 0 && f[] > F_ERR && f[] < 1. - F_ERR && avg[] > 0) {
            val += sf0[]/avg[];
          }
        }
        fts[] = val;
      }
      else {
        fts[] = 0.;
      }
    }
    else {
      if (f[] < F_ERR) { // Move toward pure-gas
        double val = 0.;
        foreach_neighbor (1) {
          if (is_solid[] == 0 && f[] > F_ERR && f[] < 1. - F_ERR && avg[] > 0) {
            val += sf0[]/avg[];
          }
        }
        fts[] = val;
      }
      else {
        fts[] = 0.;
      }
    }
  }

#else
  // Compute avg
  foreach() {
    avg[] = 0.;
    if (f[] > F_ERR && f[] < 1. - F_ERR) {
      if (dir == 1) {
        int count = 0;
        foreach_neighbor (1) {
          if (f[] > 1.-F_ERR) // Number of pure-liquid cells close to the interfacial cell
            count ++;
        }
        avg[] = count;
      }
      else {
        int count = 0;
        foreach_neighbor (1) {
          if (f[] < F_ERR) // Number of pure-gas cells close to the interfacial cell
            count ++;
        }
        avg[] = count;
      }
    }
  }

  scalar sf0[];
  foreach() {
    sf0[] = fts[];
    fts[] = 0.;
  }
  boundary({sf0,fts});

  // Compute m
  foreach() {
    if (dir == 1) {
      if (f[] > 1.-F_ERR) { // Move toward pure-liquid
        double val = 0.;
        foreach_neighbor (1) {
          if (f[] > F_ERR && f[] < 1. - F_ERR && avg[] > 0) {
            val += sf0[]/avg[];
          }
        }
        fts[] = val;
      }
      else {
        fts[] = 0.;
      }
    }
    else {
      if (f[] < F_ERR) { // Move toward pure-gas
        double val = 0.;
        foreach_neighbor (1) {
          if (f[] > F_ERR && f[] < 1. - F_ERR && avg[] > 0) {
            val += sf0[]/avg[];
          }
        }
        fts[] = val;
      }
      else {
        fts[] = 0.;
      }
    }
  }
  #endif
}

/**
## *distribute_field()*: Distribute a field localized at the interface in the closest pure gas or liquid cells (5x5 stencil)
*/

void distribute_field (scalar ftd, scalar f, int dir) {
  scalar mdisttot[], df[];
  vector mf[];
  foreach() {
    mdisttot[] = 0.;
    if (f[] > F_ERR && f[] < 1.-F_ERR) {
      coord m = mycs (point, f);
      coord n = normal (point, f);
      double alpha = plane_alpha (f[], m);
      coord prel;
      plane_area_center (m, alpha, &prel);
      foreach_dimension()
        mf.x[] = m.x;

      double mdisttothere = 0.;
      foreach_neighbor (2) {
        if ( (f[] > 1.-F_ERR && dir == 1) || (f[] < F_ERR && dir != 1 ) ) {
          coord dist;
          foreach_dimension()
            dist.x = prel.x - x;
          double magdist = 0.;
          foreach_dimension()
            magdist += sq (dist.x);
          magdist = sqrt (magdist);
          double epsi = 0.;
          foreach_dimension()
            epsi += fabs (dist.x*n.x);
          mdisttothere += epsi/magdist;
        }
      }
      mdisttot[] = mdisttothere;
    }
  }

  foreach() {
    df[] = 0.;
    if ( (f[] > 1.-F_ERR && dir == 1) || (f[] < F_ERR && dir != 1 ) ) {
      double dfhere = 0.;
      foreach_neighbor (2) {
        if (f[] > F_ERR && f[] < 1.-F_ERR) {
          coord m, n;
          foreach_dimension() {
            m.x = mf.x[];
            n.x = mf.x[];
          }
          double nnorm = 0.;
          foreach_dimension()
            nnorm += sq (n.x);
          nnorm = sqrt (nnorm);
          foreach_dimension()
            n.x /= nnorm;
          double alpha = plane_alpha (f[], m);
          coord prel;
          plane_area_center (m, alpha, &prel);

          coord dist;
          foreach_dimension()
            dist.x = prel.x - x;
          double magdist = 0.;
          foreach_dimension()
            magdist += sq (dist.x);
          magdist = sqrt (magdist);
          double epsi = 0.;
          foreach_dimension()
            epsi += fabs (dist.x*n.x);
          double mdist = epsi/magdist;

          dfhere += ftd[]*mdist/mdisttot[];
        }
      }
      df[] += dfhere;
    }
  }

  foreach()
    ftd[] = df[];
}

