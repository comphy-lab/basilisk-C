#define BGHOSTS 2
#define VTOL 1e-10 //threshold for determining if cell is interfacial in flux calculation
#define THETA_LIMIT 0.1 //used in computing timestep restriction in diffusion solver
#define MIN_IT 2.0 //minimum number of iterations of diffusion solver per global timestep
#define MAX_IT 50.0 //maximum number of iterations of diffusion solver per global timestep
#define SORPTION_SCALING 3000.0 //used for limiting timestep restriction represents what is the maximum ratio between the subgrid and grid level gradient at interface a rule of thumb would be 100 times smaller than ratio between expected boundary layer thickness and grid size

#include "poisson.h"
#include "curvature.h"


attribute {
  scalar vol_frac;
  double phase_diffusivity;
  double alpha;
  double sorption_scaling;
}

double FTCS_diffusion_stability(scalar c){
  assert(c.phase_diffusivity > 0.0);
  double tmin = 0.5*THETA_LIMIT/dimension*sq(L0/(1 << grid->maxdepth))/c.phase_diffusivity;
  if (c.sorption_scaling > 0)
    tmin = tmin/clamp(c.sorption_scaling, 1.0, SORPTION_SCALING);
  return tmin;
}

foreach_dimension()
static double interface_fraction_x (coord m, double alpha, bool right){
#if dimension == 2
  alpha += (m.x + m.y)/2;
  coord n = m;
  double xo = (right ? 1. : 0.);
  if (n.y < 0.) {
    alpha -= n.y;
    n.y = - n.y;
  }
  if (n.y < 1e-4)
    return (n.x*(right ? 1 : -1) < 0. ? 1. : 0.);
  return clamp((alpha - n.x*xo)/n.y, 0., 1.);
#else // dimension == 3

  if (fabs(m.y) < 1e-4 && fabs(m.z) < 1e-4)
    return right ? (m.x < 0.) : (m.x > 0.);

  double n1, n2;
  double j;
  n1 = m.y/(fabs(m.y) + fabs(m.z));
  n2 = m.z/(fabs(m.y) + fabs(m.z));
  j = right ? 0.5 : -0.5;
  alpha -= j*m.x;
  alpha /= (fabs(m.y) + fabs(m.z));
  return clamp(line_area(n1, n2, alpha), 0., 1.);
#endif
}


foreach_dimension()
static void face_fraction_refine_2_x(Point point, scalar s){
  vector fs = s.v;
  scalar cs2 = fs.x.vol_frac;
  /**
  If the cell is empty or full, simple injection from the coarse cell
  value is used. */

  if (cs2[] <= 0. || cs2[] >= 1.) {

    /**
    We need to make sure that the fine cells face fractions match
    those of their neighbours. */

    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
  fine(fs.x,1,j,k) = cs2[];
    for (int i = 0; i <= 1; i++)
      if (!is_refined(neighbor(2*i-1)) && neighbor(2*i-1).neighbors &&
    (is_local(cell) || is_local(neighbor(2*i-1))))
  for (int j = 0; j <= 1; j++)
    for (int k = 0; k <= 1; k++)
      fine(fs.x,2*i,j,k) = fs.x[i];
  }
  else {

    /**
    If the cell contains the embedded boundary, we reconstruct the
    boundary using VOF linear reconstruction and a normal estimated
    from the surface fractions. */

    coord n = mycs(point, cs2);
    double alpha = plane_alpha (cs2[], n);

    /**
    We need to reconstruct the face fractions *fs* for the fine cells.

    For the fine face fractions contained within the coarse cell,
    we compute the intersections directly using the VOF
    reconstruction. */

#if dimension == 2

    /**
    In 2D, we obtain the face fractions by taking into
    account the orientation of the normal. */

    if (2.*fabs(alpha) < fabs(n.y)) {
      double yc = alpha/n.y;
      int i = yc > 0.;
      fine(fs.x,1,1 - i) = n.y < 0. ? 1. - i : i;
      fine(fs.x,1,i) = n.y < 0. ? i - 2.*yc : 1. - i + 2.*yc;
    }
    else
      fine(fs.x,1,0) = fine(fs.x,1,1) = alpha > 0.;

#else // dimension == 3

    /**
    in 3D, we use the 2D projection of the reconstruction. */

    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
  if (!fine(cs2,0,j,k) || !fine(cs2,1,j,k))
    fine(fs.x,1,j,k) = 0.;
  else {
    static const coord a = {0.,0.,0.}, b = {.5,.5,.5};
    coord nc;
    nc.x = 0., nc.y = (2.*j - 1.)*n.y, nc.z = (2.*k - 1.)*n.z;
    fine(fs.x,1,j,k) = rectangle_fraction (nc, alpha, a, b);
  }

#endif // dimension == 3

    /**
    For the fine face fractions coincident with the faces of the
    coarse cell, we obtain the intersection position from the
    coarse cell face fraction. */

    for (int i = 0; i <= 1; i++)
      if (neighbor(2*i-1).neighbors &&
    (is_local(cell) || is_local(neighbor(2*i-1)))) {
  if (!is_refined(neighbor(2*i-1))) {
    if (fs.x[i] <= 0. || fs.x[i] >= 1.)
      for (int j = 0; j <= 1; j++)
        for (int k = 0; k <= 1; k++)
    fine(fs.x,2*i,j,k) = fs.x[i];
    else {
#if dimension == 2

      /**
      In 2D the orientation is obtained by looking at the values
      of face fractions in the transverse direction. */

      double a = fs.y[0,1] <= 0. || fs.y[2*i-1,1] <= 0. ||
        fs.y[] >= 1. || fs.y[2*i-1] >= 1.;
      if ((2.*a - 1)*(fs.x[i] - 0.5) > 0.) {
        fine(fs.x,2*i,0) = a;
        fine(fs.x,2*i,1) = 2.*fs.x[i] - a;
      }
      else {
        fine(fs.x,2*i,0) = 2.*fs.x[i] + a - 1.;
        fine(fs.x,2*i,1) = 1. - a;
      }

#else  // dimension == 3

      /**
      In 3D we reconstruct the face fraction from the projection
      of the cell interface reconstruction, as above. */

      for (int j = 0; j <= 1; j++)
        for (int k = 0; k <= 1; k++) {
    static const coord a = {0.,0.,0.}, b = {.5,.5,.5};
    coord nc;
    nc.x = 0., nc.y = (2.*j - 1.)*n.y, nc.z = (2.*k - 1.)*n.z;
    fine(fs.x,2*i,j,k) =
      rectangle_fraction (nc, alpha - n.x*(2.*i - 1.)/2., a, b);
        }

#endif // dimension == 3
    }
  }

  /**
  The face fractions of empty children cells must be zero. */

  for (int j = 0; j <= 1; j++)
  #if dimension > 2
    for (int k = 0; k <= 1; k++)
  #endif
      if (fine(fs.x,2*i,j,k) && !fine(cs2,i,j,k))
        fine(fs.x,2*i,j,k) = 0.;
      }
  }
}

void face_fraction(scalar f, vector nn, scalar alpha, face vector s){
  foreach_face() {
    if ((f[-1] < VTOL) || (f[] < VTOL)){ // some cell is empty
      s.x[] = 0.;}
    else{ if ((f[-1] > (1. - VTOL)) && (f[] > (1. - VTOL))){ // both cells are full
      s.x[] = 1.;}
    else {
      double vleft = 1., vright = 1.;
      if (f[] < (1. - VTOL)) {
        coord m;
        m.x = nn.x[];
        m.y = nn.y[];
  
        #if dimension >= 3
        m.z = nn.z[];
        #endif

        vleft = interface_fraction_x (m, alpha[], false);

      }
      if (f[-1] < (1. - VTOL)) {
        coord m;
        m.x = nn.x[-1];
        m.y = nn.y[-1];

        #if dimension >= 3
        m.z = nn.z[-1];
        #endif

        vright = interface_fraction_x (m, alpha[-1], true);
  
      }
      s.x[] = sqrt(vleft*vright);
    }
   }
  }
  boundary((scalar*){s});
  restriction((scalar*){s});
}

/*
foreach_dimension()
static inline coord face_barycentre_x(Point point, scalar f, face vector fs2, vector nn, scalar alpha, int i){

  coord n_temp = {0.0,0.0,0.0};

  if (fs2.x[i,0,0] <= 0.0 || fs2.x[i,0,0] >= 1.0)
    return (coord){0.0,0.0,0.0};

  for (int side = -1; side <= 0; side++){
        if (((fabs(nn.y[side+i,0,0]) + fabs(nn.z[side+i,0,0]))) > 0.0){
          n_temp.y += nn.y[side+i,0,0]/((fabs(nn.y[side+i,0,0]) + fabs(nn.z[side+i,0,0])));
          n_temp.z += nn.z[side+i,0,0]/((fabs(nn.y[side+i,0,0]) + fabs(nn.z[side+i,0,0])));
      }
  }

  double nsum = 0.0;
  foreach_dimension()
    nsum += fabs(n_temp.x);

  if (!nsum)
    return (coord){0.0,0.0,0.0};

  foreach_dimension()
    n_temp.x /= nsum;

  coord p = {0.0,0.0,0.0}, p1;
  coord n;
  ((double *)&n)[0] = n_temp.y, ((double *)&n)[1] = n_temp.z;

  double alpha_temp= line_alpha (fs2.x[i,0,0], n);
  line_center (n, alpha_temp, fs2.x[i], &p1);
  p.x = 0.0, p.y = ((double *)&p1)[0], p.z = ((double *)&p1)[1];
  return p;
}

#define face_condition(fs, cs)            \
  (fs.x[i,j,k] > 0.5 && (fs.x[i,j,0] > 0.5 && fs.x[i,0,k] > 0.5) && \
   fs.y[i,j + (j < 0),0] && fs.y[i-1,j + (j < 0),0] &&      \
   fs.y[i,j + (j < 0),k] && fs.y[i-1,j + (j < 0),k] &&      \
   fs.z[i,0,k + (k < 0)] && fs.z[i-1,0,k + (k < 0)] &&      \
   fs.z[i,j,k + (k < 0)] && fs.z[i-1,j,k + (k < 0)] &&      \
   (cs[i-1,j,0] > VTOL) && (cs[i-1,0,k] > VTOL) && (cs[i-1,j,k] > VTOL) &&       \
   (cs[i,j,0] > VTOL) && (cs[i,0,k] > VTOL) && (cs[i,j,k] > VTOL))

foreach_dimension()
static inline double face_gradient_2_x(Point point, scalar c, scalar f, face vector fs2, vector n, scalar alpha, int i){

  if (f[i] > VTOL && f[i-1] > VTOL){
    if (f[i] < 1-VTOL || f[i-1] < 1.0-VTOL){

      coord p = face_barycentre_x(point, f, fs2, n, alpha, i);
      int j = sign(p.y), k = sign(p.z);

      if (face_condition(fs2, f)) {
        p.y = fabs(p.y), p.z = fabs(p.z);
        double gradient = (((c[i,0,0] - c[i-1,0,0])*(1. - p.y) +
           (c[i,j,0] - c[i-1,j,0])*p.y)*(1. - p.z) +
          ((c[i,0,k] - c[i-1,0,k])*(1. - p.y) +
           (c[i,j,k] - c[i-1,j,k])*p.y)*p.z)/Delta;

        return (c[i] - c[i-1])/Delta;
        //return (gradient*(c[i] - c[i-1]) > 0.0) ? gradient : (c[i] - c[i-1])/Delta;
      }
    }
    return (c[i] - c[i-1])/Delta;
  }
  return 0.0;
}*/

double embed_interpolate2 (Point point, scalar s, coord * p_in, scalar cs2) {

  coord p;
  foreach_dimension()
    p.x = (*p_in).x;

  int i = sign(p.x), j = sign(p.y);
#if dimension == 2
  if ((cs2[i] >= VTOL) && (cs2[0,j] >= VTOL) && (cs2[i,j] >= VTOL))
    // bilinear interpolation when all neighbors are defined
    return ((s[]*(1. - fabs(p.x)) + s[i]*fabs(p.x))*(1. - fabs(p.y)) + 
      (s[0,j]*(1. - fabs(p.x)) + s[i,j]*fabs(p.x))*fabs(p.y));
#else //dimension == 3, see cartesian-common.h
  int k = sign(p.z);
  x = fabs(p.x); y = fabs(p.y); z = fabs(p.z);
  /* trilinear interpolation */
  if ((cs2[i] >= VTOL) && (cs2[0,j] >= VTOL) && (cs2[i,j] >= VTOL) && (cs2[0,0,k] >= VTOL) &&
      (cs2[i,0,k] >= VTOL) && (cs2[0,j,k] >= VTOL) && (cs2[i,j,k] >= VTOL)) {
    return (((s[]*(1. - x) + s[i]*x)*(1. - y) + 
       (s[0,j]*(1. - x) + s[i,j]*x)*y)*(1. - z) +
      ((s[0,0,k]*(1. - x) + s[i,0,k]*x)*(1. - y) + 
       (s[0,j,k]*(1. - x) + s[i,j,k]*x)*y)*z);
  }
#endif
  else {
    // linear interpolation with gradients biased toward the
    // cells which are defined
    double val = s[];
    foreach_dimension() {
      int i = sign(p.x);
      if ((cs2[i] >= VTOL))
  val += fabs(p.x)*(s[i] - s[]);
      else if ((cs2[-i] >= VTOL))
  val += fabs(p.x)*(s[] - s[-i]);
    }
    return val;
  }
}


#define quadratic(x,a1,a2,a3) \
  clamp((((a1)*((x) - 1.) + (a3)*((x) + 1.))*(x)/2. - (a2)*((x) - 1.)*((x) + 1.)), min(min(a1,a2), a3), max(max(a1,a2), a3))


foreach_dimension()
static inline double dirichlet_gradient2_x(Point point, scalar s, scalar cs, face vector fs, coord * n_in, double alpha, coord * p_in, coord * p2_in, double * bc_in)
{
  double bc = *bc_in;

  coord n, p, p2;

  foreach_dimension(){
    n.x = (*n_in).x;
    p.x = (*p_in).x;
    p2.x = (*p2_in).x;
  }

  coord n2;
  foreach_dimension()
    n2.x = n.x;

  #if dimension == 2
    n2.z = 0.0;
    p2.z = 0.0;
  #endif

  normalize (&n);
  foreach_dimension()
    n.x = - n.x;
  double d[2], v[2] = {nodata,nodata};
  bool defined = true;
  foreach_dimension()
    if (defined && !fs.x[(n.x > 0.)])
      defined = false;
  if (defined)
    for (int l = 0; l <= 1; l++) {
      int i = (l + 1)*sign(n.x);
      d[l] = (i - p.x)/n.x;
      double y1 = p.y + d[l]*n.y;
      int j = y1 > 0.5 ? 1 : y1 < -0.5 ? -1 : 0;
      y1 -= j;
#if dimension == 2
      if (fs.x[i + (i < 0),j] && fs.y[i,j] && fs.y[i,j+1] &&
    (cs[i,j-1] >= VTOL) && (cs[i,j] >= VTOL) && (cs[i,j+1] >= VTOL))
  v[l] = quadratic (y1, (s[i,j-1]), (s[i,j]), (s[i,j+1]));
#else // dimension == 3
      double z = p.z + d[l]*n.z;
      int k = z > 0.5 ? 1 : z < -0.5 ? -1 : 0;
      z -= k;
      bool defined = fs.x[i + (i < 0),j,k];
      for (int m = -1; m <= 1 && defined; m++)
  if (!fs.y[i,j,k+m] || !fs.y[i,j+1,k+m] ||
      !fs.z[i,j+m,k] || !fs.z[i,j+m,k+1] ||
      !(cs[i,j+m,k-1] >= VTOL) || !(cs[i,j+m,k] >= VTOL) || !(cs[i,j+m,k+1] >= VTOL))
    defined = false;
      if (defined)
  // bi-quadratic interpolation
  v[l] =
    quadratic (z,
         quadratic (y1,
        (s[i,j-1,k-1]), (s[i,j,k-1]), (s[i,j+1,k-1])),
         quadratic (y1,
        (s[i,j-1,k]),   (s[i,j,k]),   (s[i,j+1,k])),
         quadratic (y1,
        (s[i,j-1,k+1]), (s[i,j,k+1]), (s[i,j+1,k+1])));
#endif // dimension == 3
      else
  break;
    }
  if (v[0] == nodata) {  
    d[0] = max(fabs(n2.x*p2.x + n2.y*p2.y + n2.z*p2.z - alpha)/sqrt(sq(n2.x)+sq(n2.y)+sq(n2.z)), VTOL/10);
      return (bc - s[])/(d[0]*Delta);
  }
  
  if (v[1] != nodata) // third order gradient
      //if third order gradient has same sign as 1st order gradient use it
      if (((d[1]*(bc - v[0])/d[0] - d[0]*(bc - v[1])/d[1])/((d[1] - d[0])*Delta))*((bc - s[])/(d[0]*Delta)) > 0.0){
        //fprintf(fout, "%g %g\n", Delta, (d[1]*(bc - v[0])/d[0] - d[0]*(bc - v[1])/d[1])/((d[1] - d[0])*Delta));
        return (d[1]*(bc - v[0])/d[0] - d[0]*(bc - v[1])/d[1])/((d[1] - d[0])*Delta);
      }
  //else if second order gradient has same sign as 1st order gradient use it
  if (((bc - v[0])/(d[0]*Delta))*((bc - s[])/(d[0]*Delta)) > 0)
    return (bc - v[0])/(d[0]*Delta); // second-order gradient
  // else use 1st order gradient
  d[0] = max(fabs(n2.x*p2.x + n2.y*p2.y + n2.z*p2.z - alpha)/sqrt(sq(n2.x)+sq(n2.y)+sq(n2.z)), VTOL/10);
      return (bc - s[])/(d[0]*Delta);
}


double dirichlet_gradient2 (Point point, scalar s, scalar cs, face vector fs, coord * n, double alpha, coord * p, coord * p2, double * bc){
#if dimension == 2
  foreach_dimension()
    if (fabs((*n).x) >= fabs((*n).y))
      return dirichlet_gradient2_x (point, s, cs, fs, n, alpha, p, p2, bc);
#else // dimension == 3
  if (fabs((*n).x) >= fabs((*n).y)) {
    if (fabs((*n).x) >= fabs((*n).z))
      return dirichlet_gradient2_x (point, s, cs, fs, n, alpha, p, p2, bc);
  }
  else if (fabs((*n).y) >= fabs((*n).z))
    return dirichlet_gradient2_y (point, s, cs, fs, n, alpha, p, p2, bc);
  return dirichlet_gradient2_z (point, s, cs, fs, n, alpha, p, p2, bc);
#endif // dimension == 3
  return nodata;
}


double dirichlet_gradient_first_order(Point point, scalar s, coord * n, double alpha, coord * p2, double bc){
  double dist = max(fabs((*n).x*(*p2).x + (*n).y*(*p2).y + (*n).z*(*p2).z - alpha)/sqrt(sq((*n).x)+sq((*n).y)+sq((*n).z)), VTOL/10);
  return (bc - s[])/(dist*Delta);
}



static inline double bilinear_weighted2(Point point, scalar s, scalar q)
{
  #if dimension == 1
    return (3.*coarse(s)*coarse(q) + coarse(s,child.x)*coarse(q,child.x))/(3*coarse(q) + 1*coarse(q,child.x));
  #elif dimension == 2
    return (9.*coarse(s)*coarse(q)
             + 3.*(coarse(s,child.x)*coarse(q, child.x) + coarse(s,0,child.y)*coarse(q, 0,child.y))
                         + coarse(s,child.x,child.y)*coarse(q, child.x,child.y))/
                         (9*coarse(q) + 3*(1*coarse(q, child.x) + 1*coarse(q, 0,child.y)) + 1*coarse(q, child.x,child.y));
  #else // dimension == 3
    if ((27.*coarse(q)
                                + 9.*(coarse(q,child.x) + coarse(q,0,child.y) + coarse(q,0,0,child.z))
                                + 3.*(coarse(q,child.x,child.y) + coarse(q,child.x,0,child.z) + coarse(q,0, child.y,child.z))
                                + coarse(q, child.x,child.y,child.z)) <= 0.0)
      return 0.0;
    return (27.*coarse(s)*coarse(q)
                  + 9.*(coarse(s,child.x)*coarse(q,child.x) + coarse(s,0,child.y)*coarse(q,0,child.y) + coarse(s,0,0,child.z)*coarse(q,0,0,child.z))
                        + 3.*(coarse(s,child.x,child.y)*coarse(q,child.x,child.y)
                        + coarse(s,child.x,0,child.z)*coarse(q,child.x,0,child.z) + coarse(s,0,child.y,child.z)*coarse(q,0, child.y,child.z))
                        + coarse(s,child.x,child.y,child.z)*coarse(q, child.x,child.y,child.z))/
                        (27.*coarse(q)
                                + 9.*(coarse(q,child.x) + coarse(q,0,child.y) + coarse(q,0,0,child.z))
                                + 3.*(coarse(q,child.x,child.y) + coarse(q,child.x,0,child.z) + coarse(q,0, child.y,child.z))
                                + coarse(q, child.x,child.y,child.z));
  #endif
}


void mg_cycle_weighted2 (scalar * a, scalar * res, scalar * da,
               void (* relax) (scalar * da, scalar * res, 
                               int depth, void * data),
               void * data,
               int nrelax, int minlevel, int maxlevel, scalar q)
{

  /**
  We first define the residual on all levels. */

  restriction (res);

  /**
  We then proceed from the coarsest grid (*minlevel*) down to the
  finest grid. */

  minlevel = min (minlevel, maxlevel);
  for (int l = minlevel; l <= maxlevel; l++) {

    /**
    On the coarsest grid, we take zero as initial guess. */

    if (l == minlevel)
      foreach_level_or_leaf (l)
        for (scalar s in da)
          foreach_blockf (s)
            s[] = 0.;

    /**
    On all other grids, we take as initial guess the approximate solution
    on the coarser grid bilinearly interpolated onto the current grid. */

    else
      foreach_level (l)
      for (scalar s in da)
          foreach_blockf (s)
            s[] = bilinear_weighted2 (point, s, q);
    
    /**
    We then apply homogeneous boundary conditions and do several
    iterations of the relaxation function to refine the initial guess. */

    boundary_level (da, l);
    for (int i = 0; i < nrelax; i++) {
      relax (da, res, l, data);
      boundary_level (da, l);
    }
  }

  /**
  And finally we apply the resulting correction to *a*. */

  foreach() {
    scalar s, ds;
    for (s, ds in a, da)
      foreach_blockf (s)
        s[] += ds[];
  }
}

struct MGSolve_weighted2 {
  scalar * a, * b;
  double (* residual) (scalar * a, scalar * b, scalar * res,
           void * data);
  void (* relax) (scalar * da, scalar * res, int depth, 
      void * data);
  void * data;
  
  int nrelax;
  scalar * res;
  int minlevel;
  double tolerance;
  scalar cs2;
};


int NITERMAX_DIFFUSION = 10;


mgstats mg_solve_weighted2(struct MGSolve_weighted2 p){

  /**
  We allocate a new correction and residual field for each of the scalars
  in *a*. */

  scalar * da = list_clone (p.a), * res = p.res;
  if (!res)
    res = list_clone (p.b);

  /**
  The boundary conditions for the correction fields are the
  *homogeneous* equivalent of the boundary conditions applied to
  *a*. */

  for (int b = 0; b < nboundary; b++)
    for (scalar s in da)
      s.boundary[b] = s.boundary_homogeneous[b];
  
  /**
  We initialise the structure storing convergence statistics. */

  mgstats s = {0};
  double sum = 0.;
  foreach (reduction(+:sum))
    for (scalar s in p.b)
      sum += s[];
  s.sum = sum;
  s.nrelax = p.nrelax > 0 ? p.nrelax : 4;
  
  /**
  Here we compute the initial residual field and its maximum. */

  double resb;
  resb = s.resb = s.resa = p.residual (p.a, p.b, res, p.data);

  /**
  We then iterate until convergence or until *NITERMAX_DIFFUSION* is reached. Note
  also that we force the solver to apply at least one cycle, even if the
  initial residual is lower than *TOLERANCE*. */

  if (p.tolerance == 0.)
    p.tolerance = TOLERANCE;
  for (s.i = 0;
       s.i < NITERMAX_DIFFUSION && (s.i < NITERMIN || s.resa > p.tolerance);
       s.i++) {
    mg_cycle_weighted2 (p.a, res, da, p.relax, p.data,
        s.nrelax,
        p.minlevel,
        grid->maxdepth, p.cs2);
    s.resa = p.residual (p.a, p.b, res, p.data);

    /**
    We tune the number of relaxations so that the residual is reduced
    by between 2 and 20 for each cycle. This is particularly useful
    for stiff systems which may require a larger number of relaxations
    on the finest grid. */

#if 1
    if (s.resa > p.tolerance) {
      if (resb/s.resa < 1.2 && s.nrelax < 100)
  s.nrelax++;
      else if (resb/s.resa > 10 && s.nrelax > 2)
  s.nrelax--;
    }
#else
    if (s.resa == resb) /* convergence has stopped!! */
      break;
    if (s.resa > resb/1.1 && p.minlevel < grid->maxdepth)
      p.minlevel++;
#endif

    resb = s.resa;
  }
  s.minlevel = p.minlevel;
  
  /**
  If we have not satisfied the tolerance, we warn the user. */

  if (s.resa > p.tolerance) {
    scalar v = p.a[0];
    fprintf (ferr, 
       "WARNING: convergence for %s not reached after %d iterations\n"
       "  res: %g sum: %g nrelax: %d\n", v.name,
       s.i, s.resa, s.sum, s.nrelax), fflush (ferr);
  }
    
  /**
  We deallocate the residual and correction fields and free the lists. */

  if (!p.res)
    delete (res), free (res);
  delete (da), free (da);

  return s;
}

/**
The equivalent residual function is obtained in a similar way in the
case of a Cartesian grid, however the case of the tree mesh
requires more careful consideration... */

/*
struct SGS_Diffusion {
  scalar a, b;
  (const) face vector alpha;
  (const) scalar lambda;
  (const) scalar cs2;
  (const) face vector fs2;
  (const) vector nn;
  (const) scalar gamma;
  double tolerance;
  int nrelax, minlevel;
  scalar * res;
#if EMBED
  double (* embed_flux) (Point, scalar, vector, double *);
#endif
};


static void SGS_relax (scalar * al, scalar * bl, int l, void * data){
  scalar a = al[0], b = bl[0];
  struct SGS_Diffusion * p = (struct SGS_Diffusion *) data;
  (const) face vector alpha = p->alpha;
  (const) scalar lambda = p->lambda;


#if JACOBI
  scalar c[];
#else
  scalar c = a;
#endif


  foreach_level_or_leaf (l) {
    double n = - sq(Delta)*b[], d = - lambda[]*sq(Delta);
    foreach_dimension() {
      n += alpha.x[1]*a[1] + alpha.x[]*a[-1];
      d += alpha.x[1] + alpha.x[];
    }
      c[] = n/d;
  }


#if JACOBI
  foreach_level_or_leaf (l)
    a[] = (a[] + 2.*c[])/3.;
#endif
  
#if TRASH
  scalar a1[];
  foreach_level_or_leaf (l)
    a1[] = a[];
  trash ({a});
  foreach_level_or_leaf (l)
    a[] = a1[];
#endif
}*/
/**
The equivalent residual function is obtained in a similar way in the
case of a Cartesian grid, however the case of the tree mesh
requires more careful consideration... */

/*static double SGS_residual(scalar * al, scalar * bl, scalar * resl, void * data)
{
  scalar a = al[0], b = bl[0], res = resl[0];
  struct SGS_Diffusion * p = (struct SGS_Diffusion *) data;
  (const) face vector alpha = p->alpha;
  (const) scalar lambda = p->lambda;
  (const) scalar cs2 = p->cs2;
  (const) face vector fs2 = p->fs2;
  (const) vector nn = p->nn;
  (const) scalar gamma = p->gamma;
  double maxres = 0.;

  face vector g[];
  foreach_face()
    g.x[] = alpha.x[]*face_gradient_2_x(point, a, cs2, fs2, nn, gamma, 0);
  foreach (reduction(max:maxres), nowarning) {
     res[] = b[] - lambda[]*a[];
    foreach_dimension()
      res[] -= (g.x[1] - g.x[])/Delta;

    //add line that if outside domain set residual to 0
    res[] = cs2[] > VTOL ? res[] : 0.0;
    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);
  }
#else // !TREE
   foreach (reduction(max:maxres), nowarning) {
    res[] = b[] - lambda[]*a[];
    foreach_dimension()
      res[] += (alpha.x[0]*face_gradient_2_x(point, a, cs2, fs2, nn, gamma, 0) -
    alpha.x[1]*face_gradient_2_x(point, a, cs2, fs2, nn, gamma, 1))/Delta;
    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);
  }
#endif // !TREE
  return maxres;
} */

 /** The gradient of is computed using a standard three-point scheme if we are far enough from the interface (as controlled by cmin), 
  * otherwise a two-point scheme biased away from the interface is used. 
  This is the same as vof concentration refine but it assumes that tracer is already normalized by the volume fraction**/

foreach_dimension()
static double concentration_gradient_x(Point point, scalar c, scalar t)
{
  static const double cmin = 0.0;
  double cl = c[-1], cc = c[], cr = c[1];
  if (t.inverse)
    cl = 1. - cl, cc = 1. - cc, cr = 1. - cr;
  if (cc >= cmin && t.gradient != zero) {
    if (cr >= cmin) {
      if (cl >= cmin) {
  if (t.gradient)
    return t.gradient (t[-1], t[], t[1])/Delta;
  else
    return (t[1] - t[-1])/(2.*Delta);
      }
      else
  return (t[1] - t[])/Delta;
    }
    else if (cl >= cmin)
      return (t[] - t[-1])/Delta;
  }
  return 0.;
}

static void concentration_refine (Point point, scalar s){
  scalar f = s.c;
  if (cm[] == 0. || (!s.inverse && f[] <= 0.) || (s.inverse && f[] >= 1.))
    foreach_child()
      s[] = 0.;
  else {
    coord g;
    foreach_dimension()
      g.x = Delta*concentration_gradient_x (point, f, s);
    double sc = s[], cmc = 4.*cm[];
    foreach_child() {
      s[] = sc;
      foreach_dimension()
        s[] += child.x*g.x*cm[-child.x]/cmc;
    }
  }
}

void diffusion_implicit(scalar c, scalar f, double dt_global, double bc, double tolerance_diff){

  vector n[];
  scalar alpha[];

  if (f.height.x.i){
    heights (f, f.height);
  }
  boundary({f});
  restriction({f});


  //modified vof reconstruction algorithm that computes normals with height function

  reconstruction(f, n, alpha);
  //boundary((scalar *){n});
  //boundary({alpha});


  face vector fs2[];
  face vector Diff[];
  scalar theta_idt[]; //volume correction (it takes into account the effective volume of species)

  foreach_dimension(){
    fs2.x.refine = fs2.x.prolongation = face_fraction_refine_2_x;
    fs2.x.vol_frac = f;}

  face_fraction(f, n, alpha, fs2);
  

  //MPI_Barrier(MPI_COMM_WORLD);
  //fprintf(ferr, "test_0\n");
  //MPI_Barrier(MPI_COMM_WORLD);


  foreach_face(){

    coord n_temp, p;
    double gamma_1 = 0.0, gamma_2 = 0.0;

    if ((f[] < (1.-VTOL)) && (f[] > VTOL)){
      n_temp.x = n.x[0,0,0];
      n_temp.y = n.y[0,0,0];
      n_temp.z = n.z[0,0,0];


      //normally not necessary but there were edge cases producing floating point errors
      if ((fabs(n_temp.x) + fabs(n_temp.y) + fabs(n_temp.z)) == 0.0){
        foreach_dimension()
          n_temp.x = 1.0;
        double nn = 0.;
        foreach_dimension()
          nn += fabs(n_temp.x);
        foreach_dimension()
          n_temp.x /= nn;
      }  

      plane_center (n_temp, plane_alpha (f[], n_temp), f[], &p);
      gamma_1 =  p.x;}

    if ((f[-1] < (1.-VTOL)) && (f[-1] > VTOL)){
      n_temp.x = n.x[-1,0,0];
      n_temp.y = n.y[-1,0,0];
      n_temp.z = n.z[-1,0,0];

      //normally not necessary but there were edge cases producing floating point errors
      if ((fabs(n_temp.x) + fabs(n_temp.y) + fabs(n_temp.z)) == 0.0){
        foreach_dimension()
          n_temp.x = 1.0;
        double nn = 0.;
        foreach_dimension()
          nn += fabs(n_temp.x);
        foreach_dimension()
          n_temp.x /= nn;
      }

      plane_center (n_temp, plane_alpha (f[-1], n_temp), f[-1], &p);
      gamma_2 = p.x;}

    //double gamma_centroid = 1.0;

    double gamma_centroid = 1./clamp((1. + gamma_1 - gamma_2), 0.1, 2.0); //minimum = 0.5, maximum = 10

    Diff.x[] = gamma_centroid*fm.x[]*fs2.x[]*c.phase_diffusivity;
    //fprintf(fout, "%g %g %g \n", Diff.x[], gamma_centroid, fs2.x[]);
  }

  //MPI_Barrier(MPI_COMM_WORLD);
  //fprintf(ferr, "test_1\n");
  //MPI_Barrier(MPI_COMM_WORLD);



 scalar transfer_beta[];
 scalar transfer_r[];


 scalar tr_c[];
 tr_c.gradient = minmod2;
 tr_c.c = f;
 tr_c.refine = tr_c.prolongation = concentration_refine;

 foreach(){
    tr_c[] = ((f[] >= VTOL) ? c[]/f[] : 0.);}

//MPI_Barrier(MPI_COMM_WORLD);

//  fprintf(ferr, "test_2\n");
//  MPI_Barrier(MPI_COMM_WORLD);

 foreach(){
   if ((f[] >= VTOL) && (f[] <= 1.0 - VTOL)){
     coord n_temp, p, p2;
     foreach_dimension()
       n_temp.x = n.x[];

     double area = plane_area_center(n_temp, alpha[], &p);

     #if dimension == 2
       line_center(n_temp, alpha[], f[], &p2);
     #else
       plane_center(n_temp, alpha[], f[], &p2);
     #endif

     assert(sq(n.x[])+sq(n.y[])+sq(n.z[]) > 0.0);

     double dist = max(fabs(n.x[]*p2.x + n.y[]*p2.y + n.z[]*p2.z - alpha[])/sqrt(sq(n.x[])+sq(n.y[])+sq(n.z[])), VTOL/10);
     //fprintf(fout, "scaling %g %g\n", t, scaling);

     assert(dist > 0.0);


     transfer_beta[] = -c.phase_diffusivity*area/dist/pow(Delta,2);
     transfer_r[] = bc*c.phase_diffusivity*area/dist/pow(Delta,2);
   }
   else{
     transfer_beta[] = 0.0;
     transfer_r[] = 0.0;
   }
 }

//fprintf(ferr, "test_3\n");
//MPI_Barrier(MPI_COMM_WORLD);

c.sorption_scaling = 1.0;

double dt_diffusion = FTCS_diffusion_stability(c);
double num_it = clamp(ceil(dt_global/dt_diffusion), MIN_IT, MAX_IT);
dt_diffusion = dt_global/num_it;

//fprintf(ferr, "test_4\n");
//MPI_Barrier(MPI_COMM_WORLD);

 theta_idt.gradient = zero;
 theta_idt.c = f;
 theta_idt.refine = theta_idt.prolongation = vof_concentration_refine;


 foreach() {
  theta_idt[] = - 1./dt_diffusion*cm[]*max(f[], VTOL);
 }
 restriction({theta_idt});

 scalar r[];
 scalar lambda[];




 for (int ii = 1; ii<=(int)num_it; ii++){

   foreach(){

     double scaling = 1.0;

     if ((f[] >= VTOL) && (f[] <= 1.0 - VTOL)){
      coord n_temp, p, p2;
      foreach_dimension()
        n_temp.x = n.x[];

      double area = plane_area_center(n_temp, alpha[], &p);

      #if dimension == 2
        line_center(n_temp, alpha[], f[], &p2);
      #else
        plane_center(n_temp, alpha[], f[], &p2);
      #endif

      double grad_first_order = dirichlet_gradient_first_order(point, tr_c, &n_temp, alpha[], &p2, bc);

      //create a scaling factor between higher order interface gradients and 1st order gradient. 
      scaling = (grad_first_order != 0.0) ? dirichlet_gradient2(point, tr_c, f, fs2, &n_temp, alpha[], &p, &p2, &bc)/grad_first_order : 1.0;
      scaling = clamp(scaling, 0, 10); //improves robustness of system by avoiding problematic cases
      //fprintf(fout, "scaling %g\n", scaling);
     }



     r[] = theta_idt[]*tr_c[] - scaling*transfer_r[];
     lambda[] = scaling*transfer_beta[] + theta_idt[];

  }

  restriction ((scalar *){Diff, r, lambda});

  struct Poisson q;
  q.a = tr_c, q.b = r; //variables in poisson equation
  q.alpha = Diff, q.lambda = lambda; //coefficients in poisson equation
  //q.cs2 = f, q.fs2 = fs2; //volume fraction and face fraction
  //q.nn = n, q.gamma = alpha; //interace recontstruction

  mgstats mg_diffusion = mg_solve_weighted2({tr_c}, {r}, residual, relax, &q, tolerance = tolerance_diff, cs2 = f);
  fprintf(ferr, "%g %g %d %d %d\n", t, num_it, ii, mg_diffusion.i, mg_diffusion.nrelax);

  foreach(){
    tr_c[] = f[] >= VTOL ? tr_c[] : 0.0;}

  //if there is a case where convergence is not reached, do all of the rest of the iterations in the subloop in one iteration to save time  
  //normally this should not be needed
  if (mg_diffusion.i == NITERMAX_DIFFUSION){
    fprintf(ferr, "save time by increasing tolerance\n");
    tolerance_diff = min(1e-4, mg_diffusion.resa*100);
  }
 }

 foreach(){
    c[] = tr_c[]*f[];}
}



void diffusion_implicit_axi(scalar c, scalar f, double dt_global, double bc, double tolerance_diff){

  vector n[];
  scalar alpha[];

  if (f.height.x.i){
    heights (f, f.height);
  }
  boundary({f});
  restriction({f});


  //modified vof reconstruction algorithm that computes normals with height function

  reconstruction(f, n, alpha);
  //boundary((scalar *){n});
  //boundary({alpha});


  face vector fs2[];
  face vector Diff[];
  scalar theta_idt[]; //volume correction (it takes into account the effective volume of species)

  foreach_dimension(){
    fs2.x.refine = fs2.x.prolongation = face_fraction_refine_2_x;
    fs2.x.vol_frac = f;}

  face_fraction(f, n, alpha, fs2);


  foreach_face(){

    coord n_temp, p;
    double gamma_1 = 0.0, gamma_2 = 0.0;

    if ((f[] < (1.-VTOL)) && (f[] > VTOL)){
      n_temp.x = n.x[0,0,0];
      n_temp.y = n.y[0,0,0];


      //normally not necessary but there were edge cases producing floating point errors
      if ((fabs(n_temp.x) + fabs(n_temp.y)) == 0.0){
        foreach_dimension()
          n_temp.x = 1.0;
        double nn = 0.;
        foreach_dimension()
          nn += fabs(n_temp.x);
        foreach_dimension()
          n_temp.x /= nn;
      }  

      plane_center (n_temp, plane_alpha (f[], n_temp), f[], &p);
      gamma_1 =  p.x;}

    if ((f[-1] < (1.-VTOL)) && (f[-1] > VTOL)){
      n_temp.x = n.x[-1,0,0];
      n_temp.y = n.y[-1,0,0];

      //normally not necessary but there were edge cases producing floating point errors
      if ((fabs(n_temp.x) + fabs(n_temp.y)) == 0.0){
        foreach_dimension()
          n_temp.x = 1.0;
        double nn = 0.;
        foreach_dimension()
          nn += fabs(n_temp.x);
        foreach_dimension()
          n_temp.x /= nn;
      }

      plane_center (n_temp, plane_alpha (f[-1], n_temp), f[-1], &p);
      gamma_2 = p.x;}

    //double gamma_centroid = 1.0;

    double gamma_centroid = 1./clamp((1. + gamma_1 - gamma_2), 0.1, 2.0); //minimum = 0.5, maximum = 10

    Diff.x[] = gamma_centroid*fm.x[]*fs2.x[]*c.phase_diffusivity;
    //fprintf(fout, "%g %g %g \n", Diff.x[], gamma_centroid, fs2.x[]);
  }



 scalar transfer_beta[];
 scalar transfer_r[];


 scalar tr_c[];
 tr_c.gradient = minmod2;
 tr_c.c = f;
 tr_c.refine = tr_c.prolongation = concentration_refine;

 foreach(){
    tr_c[] = ((f[] >= VTOL) ? c[]/f[] : 0.);}

//MPI_Barrier(MPI_COMM_WORLD);

 foreach(){
   if ((f[] >= VTOL) && (f[] <= 1.0 - VTOL)){
     coord n_temp, p, p2;
     foreach_dimension()
       n_temp.x = n.x[];

     double area = plane_area_center(n_temp, alpha[], &p);

     #if dimension == 2
       line_center(n_temp, alpha[], f[], &p2);
     #else
       plane_center(n_temp, alpha[], f[], &p2);
     #endif

     assert(sq(n.x[])+sq(n.y[]) > 0.0);

     double dist = max(fabs(n.x[]*p2.x + n.y[]*p2.y - alpha[])/sqrt(sq(n.x[])+sq(n.y[])), VTOL/10);
     //fprintf(fout, "scaling %g %g\n", t, scaling);

     assert(dist > 0.0);


     transfer_beta[] = -cm[]*c.phase_diffusivity*area/dist/pow(Delta,2);
     transfer_r[] = bc*cm[]*c.phase_diffusivity*area/dist/pow(Delta,2);
   }
   else{
     transfer_beta[] = 0.0;
     transfer_r[] = 0.0;
   }
 }


c.sorption_scaling = 1.0;

double dt_diffusion = FTCS_diffusion_stability(c);
double num_it = clamp(ceil(dt_global/dt_diffusion), MIN_IT, MAX_IT);
dt_diffusion = dt_global/num_it;

 theta_idt.gradient = zero;
 theta_idt.c = f;
 theta_idt.refine = theta_idt.prolongation = vof_concentration_refine;


 foreach() {
  theta_idt[] = - 1./dt_diffusion*cm[]*max(f[], VTOL);
 }
 restriction({theta_idt});

 scalar r[];
 scalar lambda[];


 for (int ii = 1; ii<=(int)num_it; ii++){

   foreach(){

     double scaling = 1.0;

     if ((f[] >= VTOL) && (f[] <= 1.0 - VTOL)){
      coord n_temp, p, p2;
      foreach_dimension()
        n_temp.x = n.x[];

      double area = plane_area_center(n_temp, alpha[], &p);

      #if dimension == 2
        line_center(n_temp, alpha[], f[], &p2);
      #else
        plane_center(n_temp, alpha[], f[], &p2);
      #endif

      double grad_first_order = dirichlet_gradient_first_order(point, tr_c, &n_temp, alpha[], &p2, bc);

      //create a scaling factor between higher order interface gradients and 1st order gradient. 
      scaling = (grad_first_order != 0.0) ? dirichlet_gradient2(point, tr_c, f, fs2, &n_temp, alpha[], &p, &p2, &bc)/grad_first_order : 1.0;
      scaling = clamp(scaling, 0, 10); //improves robustness of system by avoiding problematic cases
      //fprintf(fout, "scaling %g\n", scaling);
     }



     r[] = theta_idt[]*tr_c[] - scaling*transfer_r[];
     lambda[] = scaling*transfer_beta[] + theta_idt[];

  }

  restriction ((scalar *){Diff, r, lambda});

  struct Poisson q;
  q.a = tr_c, q.b = r; //variables in poisson equation
  q.alpha = Diff, q.lambda = lambda; //coefficients in poisson equation
  //q.cs2 = f, q.fs2 = fs2; //volume fraction and face fraction
  //q.nn = n, q.gamma = alpha; //interace recontstruction

  mgstats mg_diffusion = mg_solve({tr_c}, {r}, residual, relax, &q, tolerance = tolerance_diff);
  fprintf(ferr, "%g %g %d %d %d\n", t, num_it, ii, mg_diffusion.i, mg_diffusion.nrelax);

  foreach(){
    tr_c[] = f[] >= VTOL ? tr_c[] : 0.0;}

  //if there is a case where convergence is not reached, do all of the rest of the iterations in the subloop in one iteration to save time  
  //normally this should not be needed
  if (mg_diffusion.i == NITERMAX_DIFFUSION){
    fprintf(ferr, "save time by increasing tolerance\n");
    tolerance_diff = min(1e-4, mg_diffusion.resa*100);
  }
 }

 foreach(){
    c[] = tr_c[]*f[];}
}



void diffusion_bulk_implicit_transfer_explicit(scalar c, scalar f, double dt_global, double bc, double tolerance_diff) {

  vector n[];
  scalar alpha[];

  if (f.height.x.i){
    heights (f, f.height);
  }

  //modified vof reconstruction algorithm that computes normals with height function

  reconstruction(f, n, alpha);

  face vector fs2[];
  face vector Diff[];
  scalar theta_idt[]; //volume correction (it takes into account the effective volume of species)

  foreach_dimension(){
    fs2.x.refine = fs2.x.prolongation = face_fraction_refine_2_x;
    fs2.x.vol_frac = f;}

  face_fraction(f, n, alpha, fs2);
  restriction({f});

  foreach_face(){

    coord n_temp, p;
    double gamma_1 = 0.0, gamma_2 = 0.0;

    if (f[] < 1 && f[] > 0){
      n_temp.x = n.x[0,0,0];
      n_temp.y = n.y[0,0,0];
      n_temp.z = n.z[0,0,0];
      plane_center (n_temp, alpha[], f[], &p);
      gamma_1 =  p.x;}

    if (f[-1] < 1 && f[-1] > 0){
      n_temp.x = n.x[-1,0,0];
      n_temp.y = n.y[-1,0,0];
      n_temp.z = n.z[-1,0,0];
      plane_center (n_temp, alpha[-1], f[-1], &p);
      gamma_2 = p.x;}

    double gamma_centroid = clamp(1./(1. + gamma_1 - gamma_2), 0.0, 10.0);

    Diff.x[] = gamma_centroid*fm.x[]*fs2.x[]*c.phase_diffusivity;
    //fprintf(fout, "%g %g %g \n", Diff.x[], gamma_centroid, fs2.x[]);
  }


 scalar tr_c[];
 tr_c.gradient = minmod2;
 tr_c.c = f;
 tr_c.refine = tr_c.prolongation = concentration_refine;

 foreach(){
    tr_c[] = f[] >= VTOL ? c[]/f[] : 0.0;}

c.sorption_scaling = 1.0;

double dt_diffusion = FTCS_diffusion_stability(c);
double num_it = clamp(ceil(dt_global/dt_diffusion), MIN_IT, MAX_IT);
dt_diffusion = dt_global/num_it;
//fprintf(ferr, "num it %g %g %g %g\n", dt_global, FTCS_diffusion_stability(c, f, fs2, transfer_beta), dt_diffusion, num_it);

 theta_idt.gradient = zero;
 theta_idt.c = f;
 theta_idt.refine = theta_idt.prolongation = vof_concentration_refine;

 foreach() {
  theta_idt[] = - 1./dt_diffusion*cm[]*max(f[], VTOL);
 }
 restriction({theta_idt});

 scalar r[];
 scalar lambda[];

 for (int ii = 1; ii<=(int)num_it; ii++){

   foreach(){
     r[] = theta_idt[]*tr_c[];

     //compute transfer accross interface and add it as a source term
     if ((f[] >= VTOL) && (f[] <= 1.0 - VTOL)){
      coord n_temp, p, p2;
      foreach_dimension()
        n_temp.x = n.x[];

      double area = plane_area_center(n_temp, alpha[], &p);

      #if dimension == 2
        line_center(n_temp, alpha[], f[], &p2);
      #else
        plane_center(n_temp, alpha[], f[], &p2);
      #endif

      //forward euler integration of transfer across interface
      double grad = dirichlet_gradient2(point, tr_c, f, fs2, &n_temp, alpha[], &p, &p2, &bc);
      //use a flux limiter to maintain stability of explicit source term
      double flux = c.phase_diffusivity*area*grad*dt_diffusion/Delta;
      flux = min(fabs(flux), fabs(bc-tr_c[])*f[])*sign(flux); //flux_limiter to mantain stability in small cells
      //add source term
      r[] -= flux/dt_diffusion;
     }

     lambda[] = theta_idt[];
  }

  restriction ((scalar *){Diff, r, lambda});

  struct Poisson q;
  q.a = tr_c, q.b = r; //variables in poisson equation
  q.alpha = Diff, q.lambda = lambda; //coefficients in poisson equation
  //q.cs2 = f, q.fs2 = fs2; //volume fraction and face fraction
  //q.nn = n, q.gamma = alpha; //interace recontstruction

  mgstats mg_diffusion = mg_solve_weighted2({tr_c}, {r}, residual, relax, &q, tolerance = tolerance_diff, cs2 = f);
  fprintf(ferr, "%g %g %d %d %d\n", t, num_it, ii, mg_diffusion.i, mg_diffusion.nrelax);


  //if there is a case where convergence is not reached, do all of the rest of the iterations in the subloop in one iteration to save time  
  //normally this should not be needed
  if (mg_diffusion.i == NITERMAX_DIFFUSION){
    fprintf(ferr, "save time by increasing tolerance\n");
    tolerance_diff = min(1e-4, mg_diffusion.resa*100);
  }  

  foreach(){
    tr_c[] = f[] >= VTOL ? tr_c[] : 0.0;}
 }

  foreach(){
    c[] = tr_c[]*f[];}
}


//c1 = c2*henry_coef

void species_transfer_two_field_henry(scalar c1, scalar c2, scalar f, double dt_global, double henry_coef, double tolerance_diff){

  vector n1[], n2[];
  scalar alpha1[], alpha2[];

  if (f.height.x.i){
    heights (f, f.height);
  }

  //modified vof reconstruction algorithm that computes normals with height function
  scalar f2[];
  f2.refine = f2.prolongation = fraction_refine;

  reconstruction(f, n1, alpha1);

  //for opposite phase interface normal changes sign and alpha remains the same. f2 = 1-f[]
  foreach(){
    foreach_dimension()
      n2.x[] = -n1.x[];
      f2[] = 1. - f[];
  }
  foreach(){
    alpha2[] = -alpha1[];
  }

  alpha2.n = n2;
  alpha2.refine = alpha2.prolongation = alpha_refine;

  foreach_dimension()
    n2.x.refine = n2.x.prolongation = refine_injection;

  boundary((scalar *) {n2, f2});
  boundary({alpha2});


  face vector fs1[], fs2[];
  foreach_dimension(){
    fs1.x.refine = fs1.x.prolongation = face_fraction_refine_2_x;
    fs1.x.vol_frac = f;
    fs2.x.refine = fs2.x.prolongation = face_fraction_refine_2_x;
    fs2.x.vol_frac = f2;}

  face_fraction(f, n1, alpha1, fs1);
  face_fraction(f2, n2, alpha2, fs2);

  restriction({f, f2});


  face vector Diff1[], Diff2[];
  foreach_face(){

    coord n_temp, p;
    double gamma_1 = 0.0, gamma_2 = 0.0;

    if (f[] < 1 && f[] > 0){
      n_temp.x = n1.x[0,0,0];
      n_temp.y = n1.y[0,0,0];
      n_temp.z = n1.z[0,0,0];
      plane_center (n_temp, alpha1[], f[], &p);
      gamma_1 =  p.x;}

    if (f[-1] < 1 && f[-1] > 0){
      n_temp.x = n1.x[-1,0,0];
      n_temp.y = n1.y[-1,0,0];
      n_temp.z = n1.z[-1,0,0];
      plane_center (n_temp, alpha1[-1], f[-1], &p);
      gamma_2 = p.x;}

    double gamma_centroid = clamp(1./(1. + gamma_1 - gamma_2), 0.0, 2.0);


    Diff1.x[] = gamma_centroid*fm.x[]*fs1.x[]*c1.phase_diffusivity;

    //phase_2
    gamma_1 = 0.0, gamma_2 = 0.0;

    if (f2[] < 1 && f2[] > 0){
      n_temp.x = n2.x[0,0,0];
      n_temp.y = n2.y[0,0,0];
      n_temp.z = n2.z[0,0,0];
      plane_center (n_temp, alpha2[], f2[], &p);
      gamma_1 =  p.x;}

    if (f2[-1] < 1 && f2[-1] > 0){
      n_temp.x = n2.x[-1,0,0];
      n_temp.y = n2.y[-1,0,0];
      n_temp.z = n2.z[-1,0,0];
      plane_center (n_temp, alpha2[-1], f2[-1], &p);
      gamma_2 = p.x;}

    gamma_centroid = clamp(1./(1. + gamma_1 - gamma_2), 0.0, 2.0);

    Diff2.x[] = gamma_centroid*fm.x[]*fs2.x[]*c2.phase_diffusivity;
  }



 scalar tr_c1[];
 tr_c1.gradient = minmod2;
 tr_c1.c = f;
 tr_c1.refine = tr_c1.prolongation = concentration_refine;


 scalar tr_c2[];
 tr_c2.gradient = minmod2;
 tr_c2.c = f2;
 tr_c2.refine = tr_c2.prolongation = concentration_refine;

 foreach(){
    tr_c1[] = f[] >= VTOL ? c1[]/f[] : 0.0;
    tr_c2[] = f2[] >= VTOL ? c2[]/f2[] : 0.0;}

c1.sorption_scaling = 1.0;
c2.sorption_scaling = 1.0;

//timestep is limited by minimum timstep of both phases
double dt_diffusion = min(FTCS_diffusion_stability(c1), FTCS_diffusion_stability(c2));
double num_it = clamp(ceil(dt_global/dt_diffusion), MIN_IT, MAX_IT);
dt_diffusion = dt_global/num_it;
//fprintf(ferr, "num it %g %g %g %g\n", dt_global, FTCS_diffusion_stability(c, f, fs2, transfer_beta), dt_diffusion, num_it);

 scalar theta_idt1[], theta_idt2[];

 theta_idt1.gradient = zero;
 theta_idt1.c = f;
 theta_idt1.refine = theta_idt1.prolongation = vof_concentration_refine;

 theta_idt2.gradient = zero;
 theta_idt2.c = f2;
 theta_idt2.refine = theta_idt2.prolongation = vof_concentration_refine;

 foreach(){
  theta_idt1[] = - 1./dt_diffusion*cm[]*max(f[], VTOL);
  theta_idt2[] = - 1./dt_diffusion*cm[]*max(f2[], VTOL);
 }
 restriction((scalar *){theta_idt1, theta_idt2});

 scalar r1[], r2[];
 scalar lambda1[], lambda2[];


 for (int ii = 0; ii<(int)num_it; ii++){

   foreach(){

     r1[] = theta_idt1[]*tr_c1[];
     r2[] = theta_idt2[]*tr_c2[];

     //compute transfer across interface and add it as a source term algorithm 3 in Bothe 2013
     if ((f[] >= VTOL) && (f[] <= 1.0 - VTOL)){

        //Set the flux direction
        if (tr_c1[] - henry_coef*tr_c2[] > 0.0){

          //compute value at interface on donator side
           coord n_temp, p, p2;
            foreach_dimension()
              n_temp.x = n1.x[];

       plane_area_center(n_temp, alpha1[], &p);

       double donator_concentration_interface = clamp(embed_interpolate2(point, tr_c1, &p, f), tr_c2[]*henry_coef, tr_c1[]);

       double acceptor_concentration_interface = donator_concentration_interface/henry_coef;

       //compute interface normal derivative on acceptor side
       foreach_dimension()
          n_temp.x = n2.x[];

      double area = plane_area_center(n_temp, alpha2[], &p);

      #if dimension == 2
        line_center(n_temp, alpha2[], f2[], &p2);
      #else
        plane_center(n_temp, alpha2[], f2[], &p2);
      #endif

      //forward euler integration of transfer across interface
      double grad = dirichlet_gradient2(point, tr_c2, f2, fs2, &n_temp, alpha2[], &p, &p2, &acceptor_concentration_interface);

      double flux = -c2.phase_diffusivity*area*grad*dt_diffusion/Delta;
      double limit = -f2[]*f[]*(henry_coef*tr_c2[] - tr_c1[])/(f2[] + henry_coef*f[]);

      //enforce flux direction and use flux limiter with maximum flux equal to equilibriation in interfacial cell
      flux = clamp(flux, 0.0, limit);
      
      r1[] += flux/dt_diffusion;
      r2[] -= flux/dt_diffusion;

      //fprintf(ferr, "case 1 %12e %12e %12e %12e\n", grad, -c2.phase_diffusivity*area*grad*dt_diffusion/Delta, limit, flux/dt_diffusion);
        }
    else{

      //compute value at interface on donator side
        coord n_temp, p, p2;
        foreach_dimension()
          n_temp.x = n2.x[];

       plane_area_center(n_temp, alpha2[], &p);

       
       double donator_concentration_interface = clamp(embed_interpolate2(point, tr_c2, &p, f2), tr_c1[]/henry_coef, tr_c2[]);
       //fprintf(ferr, "donator_concentration_interface %g %g %g %g\n", embed_interpolate2(point, tr_c2, &p, f2), tr_c1[]/henry_coef, tr_c2[], donator_concentration_interface);
       double acceptor_concentration_interface = donator_concentration_interface*henry_coef;

       //compute interface normal derivative on acceptor side
       foreach_dimension()
          n_temp.x = n1.x[];

      double area = plane_area_center(n_temp, alpha1[], &p);

      #if dimension == 2
        line_center(n_temp, alpha1[], f[], &p2);
      #else
        plane_center(n_temp, alpha1[], f[], &p2);
      #endif

      //forward euler integration of transfer across interface
      double grad = dirichlet_gradient2(point, tr_c1, f, fs1, &n_temp, alpha1[], &p, &p2, &acceptor_concentration_interface);
      //enforce flux direction and use flux limiter with maximum flux equal to equilibriation in interfacial cell

      double flux = c1.phase_diffusivity*area*grad*dt_diffusion/Delta;
      double limit = f2[]*f[]*(henry_coef*tr_c2[] - tr_c1[])/(f2[] + henry_coef*f[]);
      flux = clamp(flux, 0.0, limit);

      r1[] -= flux/dt_diffusion;
      r2[] += flux/dt_diffusion;


      //fprintf(ferr, "case 2 %12e %12e %12e %12e\n", grad, c1.phase_diffusivity*area*grad*dt_diffusion/Delta,  (henry_coef*f[]*c2[] - f2[]*c1[])/(f2[] + henry_coef*f[]), -flux/dt_diffusion);

        }
     }

     lambda1[] = theta_idt1[];
     lambda2[] = theta_idt2[];
  }

 restriction ((scalar *){Diff1, Diff2, r1, r2, lambda1, lambda2});

  struct Poisson pp, pp2;
  pp.a = tr_c1, pp.b = r1; //variables in poisson equation
  pp.alpha = Diff1, pp.lambda = lambda1; //coefficients in poisson equation

  mgstats mg_diffusion = mg_solve_weighted2({tr_c1}, {r1}, residual, relax, &pp, tolerance = tolerance_diff, cs2 = f);



  pp2.a = tr_c2, pp2.b = r2; //variables in poisson equation
  pp2.alpha = Diff2, pp2.lambda = lambda2; //coefficients in poisson equation

  mgstats mg_diffusion_2 = mg_solve_weighted2({tr_c2}, {r2}, residual, relax, &pp2, tolerance = tolerance_diff, cs2 = f2);

  if (mg_diffusion_2.i == NITERMAX_DIFFUSION || mg_diffusion.i == NITERMAX_DIFFUSION){
    fprintf(ferr, "save time by increasing tolerance\n");
    tolerance_diff = min(1e-4, max(mg_diffusion.resa, mg_diffusion_2.resa)*100);
  }  

 fprintf(fout, "%g %g %d %d %d %d %d\n", dt_global, num_it, ii, mg_diffusion.i, mg_diffusion.nrelax, mg_diffusion_2.i, mg_diffusion_2.nrelax);

   foreach(){
    tr_c1[] = f[] >= VTOL ? tr_c1[] : 0.0;
    tr_c2[] = f2[] >= VTOL ? tr_c2[] : 0.0;}
 }

 foreach(){
    c1[] = tr_c1[]*f[];
    c2[] = tr_c2[]*f2[];}
}

/*
void compute_sgs_flux_scaling(face vector sgs_scaling, scalar sgs_sorption_scaling, scalar c, scalar f, vector n, scalar alpha, face vector fs2){
 foreach_face()
  sgs_scaling.x[] = 1.0;
 foreach()
  sgs_sorption_scaling[] = 1.0;
}*/


void compute_sgs_flux_scaling(face vector sgs_scaling, scalar sgs_sorption_scaling, scalar c, scalar f, vector n, scalar alpha, face vector fs2){

  scalar delta_b_tracer = c.delta_b;
  scalar c_boundary_tracer = c.c_boundary;

  //compute face bulk transfer scaling

 // fprintf(fout, "flux bulk\n");

  foreach_dimension(){
   foreach_face(x) {
    sgs_scaling.x[] = 1.0;
    if (f[] > VTOL && f[-1] > VTOL){
      double flux_sgs[2] = {0.0, 0.0};
      int num_interfacial = 0;
      double tag_sgs[2] = {0.0,0.0};

      for (int i = -1; i <= 0; i++){

          if (f[i] < 1-VTOL){

            num_interfacial++;

            if (delta_b_tracer[i] > 0.0){

              double eq_conc[2];

              coord n_temp2;
              n_temp2.x = sign(n.x[i,0,0])*max(fabs(n.x[i,0,0]), 1e-5);
              n_temp2.y = sign(n.y[i,0,0])*max(fabs(n.y[i,0,0]), 1e-5);
              n_temp2.z = sign(n.z[i,0,0])*max(fabs(n.z[i,0,0]), 1e-5);


              double summ = 0.0;
              foreach_dimension()
                summ += fabs(n_temp2.x);
              foreach_dimension()
                n_temp2.x /= summ;

              double alpha_temp2 = plane_alpha (f[i], n_temp2);


              coord n_out;
              double alpha_out;

              plane_transformation(n_temp2, alpha_temp2, &n_out, &alpha_out);


              double n_sort2[3] = {n_temp2.x, n_temp2.y, n_temp2.z};
              int idx2[3] = {0,1,2};

              for (int ii=0; ii<3; ii++){
                 for (int jj=ii+1; jj<3; jj++){
                   if (fabs(n_sort2[idx2[ii]]) < fabs(n_sort2[idx2[jj]])){
                     swap (double, idx2[ii], idx2[jj]);}
                 }
               }

              coord g;
              coord h;

              double g_temp1[3] = {((double)i)*(double)(sign(n_temp2.x)),0.0,0.0};
              //coord g_temp1 = (coord){0.0,0.0,0.0};

              g = (coord){g_temp1[idx2[0]], g_temp1[idx2[1]], g_temp1[idx2[2]]};
              h = (coord){g_temp1[idx2[0]]+1., g_temp1[idx2[1]]+1., g_temp1[idx2[2]]+1.};

              double cvof1 = plane_volume(n_temp2, alpha_temp2 + ((double)i)*n_temp2.x);

              double sgs1 = clamp(cvof1 > 0.0 ? compute_eta_sgs(n_out, alpha_out, delta_b_tracer[i]/Delta, cvof1, g,h) : 0.0, 0.0,1);

              eq_conc[0] = sgs1*(c.c_inf-c_boundary_tracer[i])+c_boundary_tracer[i];


              double g_temp2[3] = {((double)(i+1))*((double)(sign(n_temp2.x))),0.0,0.0};
              g = (coord){g_temp2[idx2[0]], g_temp2[idx2[1]], g_temp2[idx2[2]]};
              h = (coord){g_temp2[idx2[0]]+1., g_temp2[idx2[1]]+1., g_temp2[idx2[2]]+1.};


              double cvof2 = plane_volume(n_temp2, alpha_temp2 + ((double)(i+1))*n_temp2.x);

              double sgs2 = clamp(cvof2 > 0.0 ? compute_eta_sgs(n_out, alpha_out, delta_b_tracer[i]/Delta, cvof2, g,h) : 0.0,0.0,1);

              eq_conc[1] = sgs2*(c.c_inf-c_boundary_tracer[i])+c_boundary_tracer[i];

              double output_gradient = (eq_conc[0]-eq_conc[1])/Delta;

  
              double gamma_1 = 0.0;
              double gamma_2 = 0.0;
              coord p;

              //compute scaling factor
              if (cvof1 > 0.0 && cvof1 < 1.0){
                plane_center (n_temp2, alpha_temp2 + ((double)i)*n_temp2.x, cvof1, &p);
                gamma_1 =  p.x;
              }

              if (cvof2 > 0.0 && cvof2 < 1.0){
                plane_center (n_temp2, alpha_temp2 + ((double)(i+1))*n_temp2.x, cvof2, &p);
                gamma_2 =  p.x;
              }

              double gamma_centroid3 = 1./clamp((1. + gamma_1 - gamma_2), 0.1, 2.0);

              double a_temp[3] = {-((((double)(i)) + 0.5)*-sign(n_temp2.x) - 0.5),-1,-1};
              coord a = (coord){a_temp[idx2[0]], a_temp[idx2[1]], a_temp[idx2[2]]};

              //multiply by face area
              coord n_temp3;
              n_temp3.x = 0.0;
              n_temp3.y = n_temp2.y;
              n_temp3.z = n_temp2.z;

              double alpha_temp3 = alpha_temp2 + ((double)i + 0.5)*n_temp2.x;

              double face_area3 = plane_volume(n_temp3, alpha_temp3);

              //if (fabs(face_area3 - compute_zeta_sgs(n_out, alpha_out, 1e5, a)*1e5*sqrt(pi)/2) > 1e-4)
              //  fprintf (ferr, "%d %g %g %g\n", i, face_area3, compute_zeta_sgs(n_out, alpha_out, 1e5, a)*1e5*sqrt(pi)/2, fs2.x[]);
              //fprintf(ferr, "test %g %g\n",  delta_b_tracer[i]/Delta, compute_zeta_sgs(n_out, alpha_out, delta_b_tracer[i]/Delta, a));

              double sgs_flux = (c.c_inf - c_boundary_tracer[i])*clamp(compute_zeta_sgs(n_out, alpha_out, delta_b_tracer[i]/Delta, a), 0, face_area3*2/sqrt(M_PI)/(delta_b_tracer[i]/Delta))/Delta*-n.x[i,0,0]/sqrt(sq(n.x[i,0,0])+sq(n.y[i,0,0])+sq(n.z[i,0,0]));
              double pseudo_numerical_flux = output_gradient*face_area3*gamma_centroid3;
              flux_sgs[i+1] = fabs(pseudo_numerical_flux) > 0.0 ? max(sgs_flux/pseudo_numerical_flux, 0.0) : 0.0;
              //confirm that subgrid flux will have same sign as numerical flux if not do not use subgrid correction
              tag_sgs[i+1] = pseudo_numerical_flux*(c[]/f[]-c[-1]/f[-1]) >= 0.0 ? 1.0 : 0.0;
           }
        }
      }

      if ((tag_sgs[0] + tag_sgs[1]) > 0.0){
        sgs_scaling.x[] = (tag_sgs[0]*flux_sgs[0] + tag_sgs[1]*flux_sgs[1])/(tag_sgs[0] + tag_sgs[1]);
        //if face is in between interfacial and bulk cell limit scaling to be <= 1
        //sgs_scaling.x[] = ((num_interfacial == 1) ? min(sgs_scaling.x[], 1) : sgs_scaling.x[]);
        //if (num_interfacial == 1)
        //  fprintf(ferr, "test %g\n", sgs_scaling.x[]);
        assert(sgs_scaling.x[] >= 0.0);
        //fprintf(ferr, "%g\n", sgs_scaling.x[]);
      }
     }
   }
}

 // fprintf(fout, "flux interface\n");

  //compute interface transfer scaling
  foreach(){
    sgs_sorption_scaling[] = 1.0;
    if ((f[] > VTOL) && (f[] < (1.0 - VTOL)) && delta_b_tracer[] > 0.0){


      coord n_temp2, p2;
      foreach_dimension()
        n_temp2.x = sign(n.x[])*max(fabs(n.x[]), 1e-5);


      //L1 normalization
      double summ = 0.0;
      foreach_dimension()
        summ += fabs(n_temp2.x);
      foreach_dimension()
        n_temp2.x /= summ;

      double alpha_temp = plane_alpha (f[], n_temp2);

      #if dimension == 2
        line_center(n_temp2, alpha_temp, f[], &p2);
      #else
        plane_center(n_temp2, alpha_temp, f[], &p2);
      #endif

      double dist = max(fabs(n_temp2.x*p2.x + n_temp2.y*p2.y + n_temp2.z*p2.z - alpha_temp)/sqrt(sq(n_temp2.x)+sq(n_temp2.y)+sq(n_temp2.z)), VTOL/10);


      coord n_out;
      double alpha_out;
      plane_transformation(n_temp2, alpha_temp, &n_out, &alpha_out);

      double sgs_conc = clamp(compute_eta_sgs(n_out, alpha_out, delta_b_tracer[]/Delta, f[], (coord){0.,0.,0.},(coord){1.,1.,1.}),0.0,1)*(c.c_inf-c_boundary_tracer[])+c_boundary_tracer[];
    //fprintf(ferr, "%12e\n", ((1 - sgs_conc)/dist)/(2/sqrt(M_PI)/80));
    sgs_sorption_scaling[] = ((c_boundary_tracer[] - sgs_conc)/dist) != 0.0 ? max((-(c.c_inf - c_boundary_tracer[])*2/sqrt(M_PI)/delta_b_tracer[]*Delta)/((c_boundary_tracer[] - sgs_conc)/dist), 0.0) : 1.0;
    assert(sgs_sorption_scaling[] >= 0.0);
    }
  }
}

void diffusion_implicit_sgs(scalar c, scalar f, double dt_global, double tolerance_diff){

  f.refine = f.prolongation = fraction_refine;
  
  vector n[];
  scalar alpha[];

  scalar c_boundary_tracer = c.c_boundary;

 // fprintf(ferr, "heights\n");

  if (f.height.x.i){
    heights (f, f.height);
  }

  boundary({f});
  restriction({f});

  //modified vof reconstruction algorithm that computes normals with height function


 // fprintf(ferr, "reconstruction\n");

  reconstruction(f, n, alpha);
  boundary_layer(c);
  boundary({c});

  face vector fs2[];
  face vector Diff[];
  scalar theta_idt[]; //volume correction (it takes into account the effective volume of species)

  foreach_dimension(){
    fs2.x.refine = fs2.x.prolongation = face_fraction_refine_2_x;
    fs2.x.vol_frac = f;}

 // fprintf(ferr, "face_fraction\n");
 // MPI_Barrier(MPI_COMM_WORLD);

  face_fraction(f, n, alpha, fs2);
  

  //fprintf(ferr, "set scalars\n");
  //MPI_Barrier(MPI_COMM_WORLD);

foreach_face(){

    coord n_temp, p;
    double gamma_1 = 0.0, gamma_2 = 0.0;

    if ((f[] < (1.-VTOL)) && (f[] > VTOL)){
      n_temp.x = n.x[0,0,0];
      n_temp.y = n.y[0,0,0];
      n_temp.z = n.z[0,0,0];


      //normally not necessary but there were edge cases producing floating point errors
      if ((fabs(n_temp.x) + fabs(n_temp.y) + fabs(n_temp.z)) == 0.0){
        foreach_dimension()
          n_temp.x = 1.0;
        double nn = 0.;
        foreach_dimension()
          nn += fabs(n_temp.x);
        foreach_dimension()
          n_temp.x /= nn;
      }  

      plane_center (n_temp, plane_alpha (f[], n_temp), f[], &p);
      gamma_1 =  p.x;}

    if ((f[-1] < (1.-VTOL)) && (f[-1] > VTOL)){
      n_temp.x = n.x[-1,0,0];
      n_temp.y = n.y[-1,0,0];
      n_temp.z = n.z[-1,0,0];

      //normally not necessary but there were edge cases producing floating point errors
      if ((fabs(n_temp.x) + fabs(n_temp.y) + fabs(n_temp.z)) == 0.0){
        foreach_dimension()
          n_temp.x = 1.0;
        double nn = 0.;
        foreach_dimension()
          nn += fabs(n_temp.x);
        foreach_dimension()
          n_temp.x /= nn;
      }

      plane_center (n_temp, plane_alpha (f[-1], n_temp), f[-1], &p);
      gamma_2 = p.x;}

    //double gamma_centroid = 1.0;

    double gamma_centroid = 1./clamp((1. + gamma_1 - gamma_2), 0.1, 2.0); //minimum = 0.5, maximum = 10

    Diff.x[] = gamma_centroid*fm.x[]*fs2.x[]*c.phase_diffusivity;
    //fprintf(fout, "%g %g %g \n", Diff.x[], gamma_centroid, fs2.x[]);
  }


 scalar transfer_beta[];
 scalar transfer_r[];


 scalar tr_c[];
 tr_c.gradient = minmod2;
 tr_c.c = f;
 tr_c.refine = tr_c.prolongation = concentration_refine;

 foreach()
    tr_c[] = ((f[] >= VTOL) ? c[]/f[] : 0.);

 foreach(){
   if ((f[] >= VTOL) && (f[] <= 1.0 - VTOL)){
     coord n_temp, p, p2;
     foreach_dimension()
       n_temp.x = sign(n.x[])*max(fabs(n.x[]), 1e-5);

     double summ = 0.0;
     foreach_dimension()
        summ += fabs(n_temp.x);
     foreach_dimension()
        n_temp.x /= summ;

     double alpha_temp = plane_alpha (f[], n_temp);

     double area = plane_area_center(n_temp, alpha_temp, &p);

     #if dimension == 2
       line_center(n_temp, alpha_temp, f[], &p2);
     #else
       plane_center(n_temp, alpha_temp, f[], &p2);
     #endif

     double dist = max(fabs(n_temp.x*p2.x + n_temp.y*p2.y + n_temp.z*p2.z - alpha_temp)/sqrt(sq(n_temp.x)+sq(n_temp.y)+sq(n_temp.z)), VTOL/10);

     //fprintf(fout, "scaling %g %g\n", t, scaling);


     transfer_beta[] = -c.phase_diffusivity*area/dist/pow(Delta,2);
     transfer_r[] = c_boundary_tracer[]*c.phase_diffusivity*area/dist/pow(Delta,2);
   }
   else{
     transfer_beta[] = 0.0;
     transfer_r[] = 0.0;
   }
 }

face vector sgs_scaling[];
scalar sgs_sorption_scaling[]; 

//fprintf(ferr, "start flux scaling\n");
//MPI_Barrier(MPI_COMM_WORLD);

compute_sgs_flux_scaling(sgs_scaling, sgs_sorption_scaling, c, f, n, alpha, fs2);


//fprintf(ferr, "end flux scaling\n");
//MPI_Barrier(MPI_COMM_WORLD);

/*
 At the end of each diffusion step we compute what is the max scaling factor for the gradient at the interface. This is used to decide max diffusive timestep for the next iteration
 */
 double sorption_coef_max = 1.0;
 double count_interfacial_cells = 1e-6;
 double count_clipped = 0.0;
 c.sorption_scaling = 1.0;
 foreach(reduction(max:sorption_coef_max) reduction(+:count_interfacial_cells) reduction(+:count_clipped)){
    if (f[] >= 0.1 && f[] <= 1 - 0.1){
      sorption_coef_max = max(sorption_coef_max,sgs_sorption_scaling[]);
      count_interfacial_cells+= 1.0;
      count_clipped += sgs_sorption_scaling[] > SORPTION_SCALING ? 1.0 : 0.0;
    }
  }

//fprintf(ferr, "test 1\n");
//MPI_Barrier(MPI_COMM_WORLD);
c.sorption_scaling = clamp(sorption_coef_max, 1.0, SORPTION_SCALING);
if (count_clipped/count_interfacial_cells > 0.01) //if over 1 percent of interfacial cells are clipped send a warning
  fprintf(ferr, "Warning subgrid scale transfer clipped %g %g\n", count_clipped/count_interfacial_cells, sorption_coef_max);

//fprintf(ferr, "compute dt\n");
//MPI_Barrier(MPI_COMM_WORLD);
//assert(dt_global > 0);

double dt_diffusion = FTCS_diffusion_stability(c);
double num_it = clamp(ceil(dt_global/dt_diffusion), MIN_IT, MAX_IT);
dt_diffusion = dt_global/num_it;
//fprintf(ferr, "num it %g %g %g %g\n", dt_global, c.phase_diffusivity, dt_diffusion, num_it);
//MPI_Barrier(MPI_COMM_WORLD);

 theta_idt.gradient = zero;
 theta_idt.c = f;
 theta_idt.refine = theta_idt.prolongation = vof_concentration_refine;


 foreach() {
  theta_idt[] = - 1./dt_diffusion*cm[]*max(f[], VTOL);
 }
 restriction({theta_idt});

 scalar r[];
 scalar lambda[];
 face vector Diff2[];


// fprintf(ferr, "solve poisson\n");
//MPI_Barrier(MPI_COMM_WORLD);


 for (int ii = 1; ii<=(int)num_it; ii++){

  foreach(){
     r[] = theta_idt[]*tr_c[] - transfer_r[]*sgs_sorption_scaling[];
     lambda[] = transfer_beta[]*sgs_sorption_scaling[] + theta_idt[];}

  foreach_face()
    Diff2.x[] = Diff.x[]*sgs_scaling.x[];

  restriction ((scalar *){Diff, Diff2, r, lambda});

  struct Poisson q;
  q.a = tr_c, q.b = r; //variables in poisson equation
  q.alpha = Diff2, q.lambda = lambda; //coefficients in poisson equation
  //q.cs2 = f, q.fs2 = fs2; //volume fraction and face fraction
  //q.nn = n, q.gamma = alpha; //interace recontstruction

  mgstats mg_diffusion = mg_solve_weighted2({tr_c}, {r}, residual, relax, &q, tolerance = tolerance_diff, cs2 = f);
  fprintf(ferr, "%g %g %d %d %d\n", t, num_it, ii, mg_diffusion.i, mg_diffusion.nrelax);
  foreach(){
    tr_c[] = f[] >= VTOL ? tr_c[] : 0.0;
    c[] = tr_c[]*f[];}

  //if there is a case where convergence is not reached, do all of the rest of the iterations in the subloop in one iteration to save time  
  //normally this should not be needed
  if (mg_diffusion.i == NITERMAX_DIFFUSION){
    fprintf(ferr, "save time by increasing tolerance\n");
    tolerance_diff = min(1e-4, mg_diffusion.resa*100);
  }

  boundary_layer(c);
  boundary({c});
  compute_sgs_flux_scaling(sgs_scaling, sgs_sorption_scaling, c, f, n, alpha, fs2);
 }
}

