/**
# Definition of the IBM regularized Dirac

The Lagrangian and Eulerian meshes communicate thanks to regularized Dirac-delta
functions: this is the core of the Immersed Boundary Method (IBM) as introduced
by [Peskin](Peskin1997).

This implementation of the IBM in Basilisk relies on the Cache structure to
efficiently loop through the Eulerian cells in the vicinity of the Lagrangian
nodes. Since the Cache structure is only defined for trees, the current
implementation is \textbf{only compatible with quadtree and octree grids, and
not with Cartesian nor multigrids}.
*/

#define BGHOSTS 2 // Having two layers of ghost cells should not be mandatory since the implementation relies on Caches

#define POS_PBC_X(X) ((u.x.boundary[left] != periodic_bc) ? (X) : (((X - (X0 +\
  L0/2)) > L0/2.) ? (X) - L0 : (X)))
#define POS_PBC_Y(Y) ((u.x.boundary[top] != periodic_bc) ? (Y) : (((Y - (Y0 +\
  L0/2)) > L0/2.) ? (Y) - L0 : (Y)))
#define POS_PBC_Z(Z) ((u.x.boundary[front] != periodic_bc) ? (Z) : (((Z - (Z0 +\
  L0/2)) > L0/2.) ? (Z) - L0 : (Z)))

struct _generate_lag_stencils_one_caps {
  lagMesh* mesh;
  bool no_warning;
};

/**
The function below loops through the Lagrangian nodes and "caches" the Eulerian
cells in a 5x5(x5) stencil around each node. In case of parallel simulations,
the cached cells are tagged with the process id.
*/
trace
void generate_lag_stencils_one_caps(struct _generate_lag_stencils_one_caps p) {
  lagMesh* mesh = p.mesh;
  bool no_warning = p.no_warning;
  for(int i=0; i<mesh->nln; i++) {
    mesh->nodes[i].stencil.n = 0;
    /**
    The current implementation assumes that the Eulerian cells around Lagrangian
    node are all at the maximum level.
    */
    double delta = (L0/(1 << grid->maxdepth));
    for(int ni=-2; ni<=2; ni++) {
      for(int nj=-2; nj<=2; nj++) {
        #if dimension < 3
        Point point = locate(POS_PBC_X(mesh->nodes[i].pos.x + ni*delta),
          POS_PBC_Y(mesh->nodes[i].pos.y + nj*delta));
        #else
        for(int nk=-2; nk<=2; nk++) {
          Point point = locate(POS_PBC_X(mesh->nodes[i].pos.x + ni*delta),
            POS_PBC_Y(mesh->nodes[i].pos.y + nj*delta),
            POS_PBC_Z(mesh->nodes[i].pos.z + nk*delta));
        #endif
        if (!(no_warning) && point.level >= 0 && point.level != grid->maxdepth)
          fprintf(stderr, "Warning: Lagrangian stencil not fully resolved.\n");
        cache_append(&(mesh->nodes[i].stencil), point, 0);
        #if _MPI
        #if dimension < 3
        if (ni == 0 && nj == 0) {
        #else
        if (ni == 0 && nj == 0 && nk == 0) {
        #endif
          if (point.level >= 0) mesh->nodes[i].pid = cell.pid;
          else mesh->nodes[i].pid = -1;
        }
        #endif
        #if dimension == 3
        }
        #endif
      }
    }
  }
}


struct _generate_lag_stencils {
  bool no_warning;
};

trace
void generate_lag_stencils(struct _generate_lag_stencils p) {
  for(int k=0; k<NCAPS; k++)
    if (CAPS(k).isactive)
      generate_lag_stencils_one_caps(mesh = &CAPS(k), no_warning = p.no_warning);
}


/**
The function below fills the scalar $forcing$ with the regularized source term
of magnitude $force$ located at $x0$. This assumes that $forcing$ is "clean",
i.e. that it is zero everywhere (or that if there is already some non-zero terms
the intention is to include them in the forcing).
*/
trace
void lag2eul(vector forcing, lagMesh* mesh) {
  for(int i=0; i<mesh->nln; i++) {
    foreach_cache(mesh->nodes[i].stencil) {
      if (point.level >= 0) {
        #if EMBED
          if (cs[] > 1.e-10) {
        #endif
        coord dist;
        dist.x = GENERAL_1DIST(x, mesh->nodes[i].pos.x);
        dist.y = GENERAL_1DIST(y, mesh->nodes[i].pos.y);
        #if dimension > 2
        dist.z = GENERAL_1DIST(z, mesh->nodes[i].pos.z);
        #endif
        #if dimension < 3
        if (fabs(dist.x) <= 2*Delta && fabs(dist.y) <= 2*Delta) {
          double weight =
            (1 + cos(.5*pi*dist.x/Delta))*(1 + cos(.5*pi*dist.y/Delta))
            /(sq(4*Delta));
        #else
        if (fabs(dist.x) <= 2*Delta && fabs(dist.y) <= 2*Delta &&
          fabs(dist.z) <= 2*Delta) {
          double weight =
            (1 + cos(.5*pi*dist.x/Delta))*(1 + cos(.5*pi*dist.y/Delta))
            *(1 + cos(.5*pi*dist.z/Delta))/(cube(4*Delta));
        #endif
          foreach_dimension() forcing.x[] += weight*mesh->nodes[i].lagForce.x;
        }
        #if EMBED
          }
        #endif
      }
    }
  }
}

/**
The function below interpolates the eulerian velocities onto the nodes of
the Lagrangian mesh.
*/
trace
void eul2lag(lagMesh* mesh) {
  for(int ii=0; ii<mesh->nln; ii++) {
    foreach_dimension() mesh->nodes[ii].lagVel.x = 0.;
    foreach_cache(mesh->nodes[ii].stencil) {
      if (point.level >= 0) {
        coord dist;
        dist.x = GENERAL_1DIST(x, mesh->nodes[ii].pos.x);
        dist.y = GENERAL_1DIST(y, mesh->nodes[ii].pos.y);
        #if dimension > 2
        dist.z = GENERAL_1DIST(z, mesh->nodes[ii].pos.z);
        #endif
        #if dimension < 3
        if (fabs(dist.x) <= 2*Delta && fabs(dist.y) <= 2*Delta) {
          double weight = (1 + cos(.5*pi*dist.x/Delta))*
            (1 + cos(.5*pi*dist.y/Delta))/16.;
        #else
        if (fabs(dist.x) <= 2*Delta && fabs(dist.y) <= 2*Delta
          && fabs(dist.z) <= 2*Delta) {
          double weight = (1 + cos(.5*pi*dist.x/Delta))*
            (1 + cos(.5*pi*dist.y/Delta))*(1 + cos(.5*pi*dist.z/Delta))/64.;
        #endif
        foreach_dimension() mesh->nodes[ii].lagVel.x += weight*u.x[];
        }
      }
    }
  }

  /**
  In case of parallel simulations, we communicate the Lagrangian velocity
  so that all processes have the same Lagrangian velocities.
  */
  #if _MPI
    if (mpi_npe > 1) reduce_lagVel(mesh);
  #endif
}

/**
The function below fills a scalar field "stencils" with noise in all "cached"
cells. Passing this scalar to the \textit{adapt_wavelet} function ensure all
the 5x5(x5) stencils around the Lagrangian nodes are at the same level.
*/
scalar stencils[];
trace
void tag_ibm_stencils_one_caps(lagMesh* mesh) {
  for(int i=0; i<mesh->nln; i++) {
    foreach_cache(mesh->nodes[i].stencil) {
      if (point.level >= 0) {
        coord dist;
        dist.x = GENERAL_1DIST(x, mesh->nodes[i].pos.x);
        dist.y = GENERAL_1DIST(y, mesh->nodes[i].pos.y);
        #if dimension > 2
        dist.z = GENERAL_1DIST(z, mesh->nodes[i].pos.z);
        #endif
        #if dimension < 3
        if (fabs(dist.x) <= 2*Delta && fabs(dist.y) <= 2*Delta) {
          stencils[] = sq(dist.x + dist.y)/sq(2.*Delta)*(2.+noise());
        #else
        if (fabs(dist.x) <= 2*Delta && fabs(dist.y) <= 2*Delta
          && fabs(dist.z) <= 2*Delta) {
          stencils[] = sq(dist.x + dist.y + dist.z)/cube(2.*Delta)*(2.+noise());
        #endif
        }
      }
    }
  }
}

trace
void tag_ibm_stencils() {
  foreach() stencils[] = 0.;
  for(int k=0; k<NCAPS; k++)
    if (CAPS(k).isactive) tag_ibm_stencils_one_caps(&CAPS(k));
}

/**
## References

~~~bib
@Article{Peskin1977,
  author    = {Peskin, C.S.},
  title     = {{Numerical analysis of blood flow in the heart}},
  journal   = {Journal of Computational Physics},
  year      = {1977},
  volume    = {25},
  number    = {3},
  pages     = {220--252},
  file      = {:files/Peskin1977 - Numerical Analysis of Blood Flow in the Heart.pdf:PDF},
  groups    = {FluidSolid flows},
  timestamp = {2013.07.18},
}
~~~
*/
