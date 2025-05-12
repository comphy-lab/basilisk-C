struct vortex_filament{
  int     n_seg;
  double  a;
  double* t;
  coord*  c;
  coord*  tvec;
  coord*  nvec;
  coord*  bvec;
  coord   pcar;
};

void draw_space_curve(int n, coord *p){
  for (int i = 0; i < n-1; i++){
    glBegin(GL_LINES);
      glVertex3f(p[i  ].x, p[i  ].y, p[i  ].z);
      glVertex3f(p[i+1].x, p[i+1].y, p[i+1].z);
    glEnd();
  }
}

void fd_derivative( int n, double delta_t0, coord shift, coord *X, coord *dX){
  for (int i = 1; i < n-1; i++){
    foreach_dimension()
    dX[i].x = (X[i+1].x - X[i-1].x)/(2*delta_t0);
  }
  foreach_dimension(){
    dX[0].x   = (X[1].x - X[n-2].x + shift.x)/(2*delta_t0);
    dX[n-1].x = (X[1].x - X[n-2].x + shift.x)/(2*delta_t0);
  }
}

#include <gsl/gsl_spline.h>
#pragma autolink -lgsl -lgslcblas

coord gsl_interp1d( int n, double* t0, coord * P0, double tq){
  coord Pq;
  gsl_interp_accel *acc = gsl_interp_accel_alloc();

  foreach_dimension(){
    double *P_x;
    P_x = malloc(sizeof(double)*n);

    for (int i = 0; i < n; i++)
      P_x[i] = P0[i].x;
    gsl_spline *spline_x = gsl_spline_alloc(gsl_interp_cspline, n);
    gsl_spline_init(spline_x, t0, P_x, n);
    Pq.x = gsl_spline_eval (spline_x, tq, acc);
    gsl_spline_free (spline_x);

    free(P_x);
  }
  gsl_interp_accel_free (acc);
  return Pq;
}

double frenet_projection (double pos, void *params){
  struct vortex_filament *p = (struct vortex_filament *) params;

  int n_seg = p->n_seg;
  double* t   = p->t;
  coord*  c   = p->c;
  coord*  tvec = p->tvec;
  coord*  nvec = p->nvec;
  coord*  bvec = p->bvec;
  coord   pcar = p->pcar;

  coord ccar, frenet[3];
  ccar = gsl_interp1d( n_seg, t, c, pos);

  frenet[0] = gsl_interp1d( n_seg, t, tvec, pos);
  frenet[1] = gsl_interp1d( n_seg, t, nvec, pos);
  frenet[2] = gsl_interp1d( n_seg, t, bvec, pos);

  return vecdot(vecdiff(pcar, ccar), frenet[0]);
}

#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
double frenet_projection_min ( int n_seg, double as, double* t0, coord* c, coord* tvec, coord* nvec, coord* bvec, coord pcar, double r){
  struct vortex_filament params = {n_seg, as, t0, c, tvec, nvec, bvec, pcar};

  int status, verbose = 0;
  int iter = 0, max_iter = 100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;

  double x_lo = r - 0.25, x_hi = r + 0.25;
  gsl_function F;

  F.function = &frenet_projection;
  F.params = &params;

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_set_error_handler_off();
  gsl_root_fsolver_set (s, &F, x_lo, x_hi);

  if (verbose == 1) {
    printf ("using %s method\n", gsl_root_fsolver_name (s));
    printf ("%5s [%9s, %9s] %9s %10s %9s\n", "iter", "lower", "upper", "root", "err", "err(est)");
  }

  do {
    iter++;
    status = gsl_root_fsolver_iterate (s);
    r = gsl_root_fsolver_root (s);
    x_lo = gsl_root_fsolver_x_lower (s);
    x_hi = gsl_root_fsolver_x_upper (s);
    status = gsl_root_test_interval (x_lo, x_hi, 0, 1e-8);

    if ((status == GSL_SUCCESS) && (verbose == 1)){
      printf ("Converged:\n");
      printf ("%5d [%.7f, %.7f] %.7f %.7f\n", iter, x_lo, x_hi, r, x_hi - x_lo);
    }

  }
  while (status == GSL_CONTINUE && iter < max_iter);
  gsl_root_fsolver_free (s);

  return r;
}
