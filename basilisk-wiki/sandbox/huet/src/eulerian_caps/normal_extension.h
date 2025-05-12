/**
The function below assumes accurate normal vectors in a neighbourhood of the
interface and extends the scalar a defined in interfacial cells only to the
whole neighbourhood. Note that in the neighbourhood of the interface, initially
the scalar a needs to be initialized with a guess.

We use a Hamilton-Jacobi equation which characteristics are the normal vectors:
$$
\frac{\partial q}{\partial \tau} + S(\bm{x}) \bm{n} \cdot \nabla q = 0
$$
with $\tau$ a pseudo-time.
*/

#include "curvature.h"

#ifndef MAX_ITER_EXTENSION
  #define MAX_ITER_EXTENSION (10)
#endif
#ifndef MAX_ITER_NV_EXTENSION
  #define MAX_ITER_NV_EXTENSION (3*MAX_ITER_EXTENSION)
#endif
#ifndef RENORMALIZE_INTERFACIAL_NORMALS
  #define RENORMALIZE_INTERFACIAL_NORMALS 0
#endif

struct Extension {
  // Compulsory argument
  scalar* q;
  // Optional arguments
  scalar* q_prev; // FIXME: some tests have resulted in a
  // segfault when len(q_prev)<len(q)... As such, q_prev is not
  // really optional until this is fixed.
  int max_iter; // default is MAX_ITER_EXTENSION
};

trace
void normal_scalar_extension(struct Extension p) {
  scalar* q = p.q;
  int max_iter = (p.max_iter) ? p.max_iter : MAX_ITER_EXTENSION;
  double dtau = .25*L0/((1 << grid->maxdepth));
  int lenq = list_len(q);
  scalar* q_prev = NULL;
  if (p.q_prev && list_len(p.q_prev)>=lenq) q_prev = p.q_prev;
  else {
    for (scalar s in q) {
      scalar s_prev = new scalar;
      q_prev = list_add(q_prev, s_prev);
    }
  }
  int k = 0; // k counts the number of pseudo-time iterations
  boundary(q);
  while (k < max_iter) {
    foreach() {
      if (GAMMA) {
        for (int i=0; i<lenq; i++) {
          val(q_prev[i],0,0,0) = val(q[i],0,0,0);
        }
      }
    }
    boundary(q_prev);
    foreach() {
      if (GAMMA && !(IS_INTERFACE_CELL(point, f))) {
        for (int i=0; i<lenq; i++) {
          /** Unfortunately, the syntax "q[i][]" does not compile, so as far as
          I know the only workaround is to use the function val() and
          forget about using foreach_dimension().
          In some rare configuration, the upwind scheme still implies to
          access q_prev outside gamma. We choose to not consider these
          degenerated configurations at all, hence the precense of
          ngcaps>GRAD_THRESHOLD, testing the location of the cell. */
          val(q[i],0,0,0) -= dtau*(max(-sign(f[])*extended_n.x[],\
          0.)*(val(q_prev[i],0,0,0) - val(q_prev[i],-1,0,0)*(val(ngcaps,-1,0,0)\
          > GRAD_THRESHOLD)) + min(-sign(f[])*extended_n.x[],0.)\
          *(val(q_prev[i],1,0,0)*(val(ngcaps,1,0,0) > GRAD_THRESHOLD)\
          - val(q_prev[i],0,0,0)) + max(-sign(f[])*extended_n.y[],\
          0.)*(val(q_prev[i],0,0,0) - val(q_prev[i],0,-1,0)*\
          (val(ngcaps,0,-1,0) > GRAD_THRESHOLD)) + \
          min(-sign(f[])*extended_n.y[],0.)*(val(q_prev[i],0,1,0)*\
          (val(ngcaps,0,1,0) > GRAD_THRESHOLD) - val(q_prev[i],0,0,0))
          #if dimension==2
            )/Delta;
          #else
            + max(-sign(f[])*extended_n.z[], 0.)*(val(q_prev[i],0,0,0) -\
            val(q_prev[i],0,0,-1)*(val(ngcaps,0,0,-1) > GRAD_THRESHOLD)) + \
            min(-sign(f[])*extended_n.z[],0.)*(val(q_prev[i],0,0,1)*\
            (val(ngcaps,0,0,1) > GRAD_THRESHOLD) - val(q_prev[i],0,0,0)))/Delta;
          #endif
        }
      }
    }
    boundary(q);
    k++;
  }
  if (!(p.q_prev) || list_len(p.q_prev)<lenq) delete(q_prev);
}

trace
void initialize_normals_for_extension(vector normals) {
  coord n_coord;
  // Get accurate normals *on* the interface from the height functions
  // FIXME: we would like to have normal vectors at the center of the cells
  // instead. The current implementation for normal vectors is therefore
  // order 1 in space.

  vector fh = f.height, h = automatic (fh);
  if (!fh.x.i)
    heights (f, h);

  foreach() {
    if (GAMMA) {
      if (IS_INTERFACE_CELL(point, f)) {
        n_coord = height_normal(point, f, h);
        foreach_dimension() {
          normals.x[] = n_coord.x;
        }
        #if RENORMALIZE_INTERFACIAL_NORMALS
          double normn = sqrt(sq(normals.x[]) + sq(normals.y[]) + sq(normals.z[]));
          if (normn > 1.e-30) {
            foreach_dimension() {
              normals.x[] /= normn;
            }
          }
          else foreach_dimension() normals.x[] = 0.;
        #endif
      }
      else {
        foreach_dimension() {
          normals.x[] = - grad_caps.x[]/ngcaps[];
        }
      }
    }
    else {
      foreach_dimension() {
        normals.x[] = 0.;
      }
    }
  }
}

void upwind_extension(Point point, scalar my_field, vector my_normal, double dtau) {
  double a = 0.;
  foreach_dimension() {
    /** In some rare configurations, the upwind scheme needs to
    access q_prev outside gamma. We choose to not consider these
    degenerated configurations at all, hence the precense of
    NBG_GAMMA(*,*,*), testing the location of the cell. */
    a += dtau*(max(-sign(f[])*my_normal.x[], 0.)*(my_field[] -
        my_field[-1]*NBG_GAMMA(-1,0,0)) + min(-sign(f[])*my_normal.x[], 0.)*
        (my_field[1]*NBG_GAMMA(1,0,0) - my_field[]))/Delta;
  }
  my_field[] -= a;
}

trace
void normal_vector_extension(vector qv) {
  double dtau = .25*L0/((1 << grid->maxdepth));
  vector qv_prev[];
  int k = 0;
  while (k < MAX_ITER_NV_EXTENSION) {
    foreach() {
      if (GAMMA) foreach_dimension() qv_prev.x[] = qv.x[];
    }
    boundary((scalar *){qv_prev, qv});

    foreach() {
      if (GAMMA && !(IS_INTERFACE_CELL(point, f))) {
        foreach_dimension() {
          upwind_extension(point, qv.x, qv_prev, dtau);
        }
      }
    }
    #if RENORMALIZE_EXTENDED_NORMALS
    foreach() {
      if (GAMMA && !(IS_INTERFACE_CELL(point, f))) {
        double mnorm = sqrt(sq(qv.x[]) + sq(qv.y[]) + sq(qv.z[]));
        if (mnorm > 1.e-30) foreach_dimension() qv.x[] /= mnorm;
        else foreach_dimension() qv.x[] = 0.;
      }
    }
    #endif
    boundary((scalar *){qv});
    k++;
  }
}

coord youngs_normals(Point point, scalar c) {
  double m1, m2, mnorm;
  coord mxyz;
  m1 = c[-1,-1,-1] + c[-1,1,-1] + c[-1,-1,1] + c[-1,1,1] +
       2.*(c[-1,-1,0] + c[-1,1,0] + c[-1,0,-1] + c[-1,0,1]) +
       4.*c[-1,0,0];
  m2 = c[1,-1,-1] + c[1,1,-1] + c[1,-1,1] + c[1,1,1] +
       2.*(c[1,-1,0] + c[1,1,0] + c[1,0,-1] + c[1,0,1]) +
       4.*c[1,0,0];
  mxyz.x = m1 - m2;

  m1 = c[-1,-1,-1] + c[-1,-1,1] + c[1,-1,-1] + c[1,-1,1] +
       2.*( c[-1,-1,0] + c[1,-1,0] + c[0,-1,-1] + c[0,-1,1]) +
       4.*c[0,-1,0];
  m2 = c[-1,1,-1] + c[-1,1,1] + c[1,1,-1] + c[1,1,1] +
       2.*(c[-1,1,0] + c[1,1,0] + c[0,1,-1] + c[0,1,1]) +
       4.*c[0,1,0];
  mxyz.y = m1 - m2;

  m1 = c[-1,-1,-1] + c[-1,1,-1] + c[1,-1,-1] + c[1,1,-1] +
       2.*(c[-1,0,-1] + c[1,0,-1] + c[0,-1,-1] + c[0,1,-1]) +
       4.*c[0,0,-1];
  m2 = c[-1,-1,1] + c[-1,1,1] + c[1,-1,1] + c[1,1,1] +
       2.*(c[-1,0,1] + c[1,0,1] + c[0,-1,1] + c[0,1,1]) +
       4.*c[0,0,1];
  mxyz.z = m1 - m2;

  mnorm = fabs(mxyz.x) + fabs(mxyz.y) + fabs(mxyz.z);
  if (mnorm < 1e-30) {
    mxyz.x = 1.;
    mxyz.y = 0.;
    mxyz.z = 0.;
    return mxyz;
  }

  mxyz.x /= mnorm;
  mxyz.y /= mnorm;
  mxyz.z /= mnorm;
  return mxyz;
}
