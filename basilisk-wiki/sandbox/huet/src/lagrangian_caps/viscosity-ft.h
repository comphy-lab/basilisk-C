/**
# Variable viscosity across a membrane
*/
#include "fractions.h"

#ifndef MUP
  #define MUP 1.;
#endif
#ifndef MUC
  #define MUC 1.;
#endif
#define CAPS_VISCOSITY 1.

/** We define the "grid gradient" $\bm{G}$, according to Tryggvason, JCP 2001.*/
vector G[];
void construct_divG(scalar divG, lagMesh* mesh) {
  comp_normals(mesh);
  #if dimension < 3
  compute_lengths(mesh);
  for(int i=0; i<mesh->nle; i++) {
    // compute the grid gradient on the midpoint of the edge
    coord gg; // grid gradient
    int en[2];
    en[0] = mesh->edges[i].node_ids[0];
    en[1] = mesh->edges[i].node_ids[1];
    foreach_dimension() gg.x = mesh->edges[i].normal.x*mesh->edges[i].length;
    for(int j=0; j<2; j++) {
      foreach_cache(mesh->nodes[en[j]].stencil) {
        //spread half the grid gradient of the edge from each of its nodes
        if (point.level >= 0) {
        coord dist;
          dist.x = GENERAL_1DIST(x, mesh->nodes[en[j]].pos.x);
          dist.y = GENERAL_1DIST(y, mesh->nodes[en[j]].pos.y);
          if (sq(dist.x) <= sq(2*Delta) && sq(dist.y) <= sq(2*Delta)) {
            double weight =
              (1 + cos(.5*pi*dist.x/Delta))*(1 + cos(.5*pi*dist.y/Delta))/
              (16.*sq(Delta));
            foreach_dimension() G.x[] -= weight*.5*gg.x;
          }
        }
      }
    }
  }
  #else
  for(int i=0; i<mesh->nlt; i++) {
    /** compute the grid gradient on the midpoint of the edge */
    coord gg; // gg for "grid gradient"
    int tn[3]; // tn for "triangle's nodes"
    tn[0] = mesh->triangles[i].node_ids[0];
    tn[1] = mesh->triangles[i].node_ids[1];
    tn[2] = mesh->triangles[i].node_ids[2];
    foreach_dimension()
      gg.x = mesh->triangles[i].normal.x*mesh->triangles[i].area;
    for(int j=0; j<3; j++) {
      foreach_cache(mesh->nodes[tn[j]].stencil) {
        /** Spread one third of the grid gradient of the triangle to each of its vertices */
        if (point.level >= 0) {
        coord dist;
          dist.x = GENERAL_1DIST(x, mesh->nodes[tn[j]].pos.x);
          dist.y = GENERAL_1DIST(y, mesh->nodes[tn[j]].pos.y);
          dist.z = GENERAL_1DIST(z, mesh->nodes[tn[j]].pos.z);
          if (sq(dist.x) <= sq(2*Delta) && sq(dist.y) <= sq(2*Delta)
            && sq(dist.z) <= sq(2*Delta)) {
            double weight = (1 + cos(.5*pi*dist.x/Delta))
              *(1 + cos(.5*pi*dist.y/Delta))*(1 + cos(.5*pi*dist.z/Delta))
              /(cube(4*Delta));
            foreach_dimension() G.x[] -= weight*gg.x/3.;
          }
        }
      }
    }
  }
  #endif
  foreach()
    if (cm[] > 1.e-20)
      foreach_dimension() divG[] += (G.x[1] - G.x[-1])/(2.*Delta);
}

double muc, mup;
scalar I[];
scalar prevI[];
scalar divG[];
event defaults (i = 0) {
  mu = new face vector;
  mup = MUP;
  muc = MUC;
  foreach() {
    prevI[] = 0.;
    I[] = 0.;
    foreach_dimension() divG[] = 0.;
  }
  /** We define below the homogeneous Dirichlet boundary conditions for the
  grid-gradient, and the indicator function on all walls. A consequence of this
  is that in the case of bi/tri-periodic boundary conditions in 2D/3D the
  Poisson solver has no Dirichlet BC to ensure the indicator function lies
  between 0 and 1. As a result, tri-periodic boxes (or bi-periodic squares)
  are not recommended with the current implementation. */
  foreach_dimension() {
    if (u.x.boundary[left] != periodic_bc) {
      I[left] = dirichlet(0.);
      I[right] = dirichlet(0.);
      prevI[left] = dirichlet(0.);
      prevI[right] = dirichlet(0.);
      divG[left] = dirichlet(0.);
      divG[right] = dirichlet(0.);
    }
  }
  #if EMBED
    I[embed] = dirichlet(0.);
    I.third = true;
  #endif
}

event properties (i++) {
  foreach() {
    if (cm[] > 1.e-20) {
      foreach_dimension() G.x[] = 0.;
      divG[] = 0.;
    }
  }

  for (int k=0; k<NCAPS; k++)
    if (CAPS(k).isactive) 
      construct_divG(divG, &CAPS(k));
  poisson(I, divG, tolerance = 1.e-6, minlevel = 4);

  // Simple clamping of I:
  foreach() {
    if (cm[] > 1.e-20) {
      if (fabs(divG[]) > 1.e-10) {
        I[] = clamp(I[], 0, 1);
        prevI[] = I[];
      }
      else {
        prevI[] = round(prevI[]);
        I[] = prevI[];
      }
      I[] = clamp(I[], 0, 1);
    }
  }
}

event properties (i++) {
  face vector muv = mu;
  foreach_face()
    if (fm.x[] > 1.e-20)
      muv.x[] = (mup + (muc - mup)*.5*(I[] + I[-1]))*fm.x[];
}
