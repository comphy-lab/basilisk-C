/**
# Time-implicit discretisation of reaction--diffusion equations

This is an adaptation of [`diffusion.h`](/src/diffusion.h) for pairwise duffusion and flux partitioning

We want to discretise implicitly the reaction--diffusion equation
$$
\theta\partial_tf = \nabla\cdot(D\nabla f) + \beta f + r
$$ 
where $\beta f + r$ is a reactive term,  $D$ is the diffusion
coefficient and $\theta$ can be a density term.

Using a time-implicit backward Euler discretisation, this can be
written
$$
\theta\frac{f^{n+1} - f^{n}}{dt} = \nabla\cdot(D\nabla f^{n+1}) + \beta
f^{n+1} + r
$$
Rearranging the terms we get
$$
\nabla\cdot(D\nabla f^{n+1}) + (\beta - \frac{\theta}{dt}) f^{n+1} =
- \frac{\theta}{dt}f^{n} - r
$$
This is a Poisson--Helmholtz problem which can be solved with a
multigrid solver. */

#include "poisson-three_field.h"

/**
The parameters of the `diffusion()` function are a scalar field `f`,
scalar fields `r` and $\beta$ defining the reactive term, the timestep
`dt` and a face vector field containing the diffusion coefficient
`D`. If `D` or $\theta$ are omitted they are set to one. If $\beta$ is
omitted it is set to zero. Both `D` and $\beta$ may be constant
fields.

Note that the `r`, $\beta$ and $\theta$ fields will be modified by the solver.

The function returns the statistics of the Poisson solver. */

#define VTOL 1e-10 

static inline void restriction_embed_linear_i (Point point, scalar s)
{  
  // 0 children
  if (!(1-cs[])) {
    s[] = 0.;
    return;
  }

  /**
  We first try to interpolate "diagonally". If enough child cells are
  defined (i.e. have non-zero embedded fractions), we return the
  corresponding value. */

  double val = 0., nv = 0.;
  for (int i = 0; i <= 1; i++)
#if dimension > 2
    for (int j = 0; j <= 1; j++)
#endif
      if ((1-fine(cs,0,i,j)) && (1-fine(cs,1,!i,!j)))
  val += (fine(s,0,i,j) + fine(s,1,!i,!j))/2., nv++;
  if (nv > 0.) {
    s[] = val/nv;
    return;
  }

  /**
  Otherwise, we use the average of the child cells which are defined
  (there is at least one). */
  
  coord p = {0.,0.,0.};
  foreach_child()
    if (1-cs[])
      p.x += x, p.y += y, p.z += z, val += s[], nv++;
  assert (nv > 0.);
  s[] = val/nv;

}

/**
## Refinement/prolongation of cell-centered fields

For refinement, we use either bilinear interpolation, if the required
four coarse cell values are defined or trilinear interpolation if only
three coarse cell values are defined. If less that three coarse cell
values are defined ("pathological cases" below), we try to estimate
gradients in each direction and add the corresponding correction. */

static inline void refine_embed_linear_i (Point point, scalar s)
{
  foreach_child() {
    if (!(1-cs[]))
      s[] = 0.;
    else {
      assert ((1-coarse(cs)));
      int i = (child.x + 1)/2, j = (child.y + 1)/2;
#if dimension == 2
      if ((1-coarse(fs.x,i)) && (1-coarse(fs.y,0,j)) &&
    (coarse(cs) == 0. || coarse(cs,child.x) == 0. ||
     coarse(cs,0,child.y) == 0. || coarse(cs,child.x,child.y) == 0.)) {
  assert ((1-coarse(cs,child.x)) && (1-coarse(cs,0,child.y)));
  if ((1-coarse(fs.x,i,child.y)) && (1-coarse(fs.y,child.x,j))) {
    // bilinear interpolation
    assert ((1-coarse(cs,child.x,child.y)));
    s[] = (9.*coarse(s) + 
     3.*(coarse(s,child.x) + coarse(s,0,child.y)) + 
     coarse(s,child.x,child.y))/16.;
  }
  else
    // triangular interpolation   
    s[] = (2.*coarse(s) + coarse(s,child.x) + coarse(s,0,child.y))/4.;
      }
      else if ((1-coarse(cs,child.x,child.y)) &&
         (((1-coarse(fs.x,i)) && (1-coarse(fs.y,child.x,j))) ||
    ((1-coarse(fs.y,0,j)) && (1-coarse(fs.x,i,child.y))))) {
  // diagonal interpolation
  s[] = (3.*coarse(s) + coarse(s,child.x,child.y))/4.;
      }
#else // dimension == 3
      int k = (child.z + 1)/2;
      if (coarse(fs.x,i) < 0.25 && coarse(fs.y,0,j) < 0.25 &&
    coarse(fs.z,0,0,k) < 0.25 &&
    (coarse(cs) == 0. || coarse(cs,child.x) == 0. ||
     coarse(cs,0,child.y) == 0. || coarse(cs,child.x,child.y) == 0. ||
     coarse(cs,0,0,child.z) == 0. || coarse(cs,child.x,0,child.z) == 0. ||
     coarse(cs,0,child.y,child.z) == 0. ||
     coarse(cs,child.x,child.y,child.z) == 0.)) {
  assert ((1-coarse(cs,child.x)) && (1-coarse(cs,0,child.y)) &&
    (1-coarse(cs,0,0,child.z)));
  if ((1-coarse(fs.x,i,child.y)) && (1-coarse(fs.y,child.x,j)) &&
      (1-coarse(fs.z,child.x,child.y,k)) &&
      (1-coarse(fs.z,child.x,0,k)) && (1-coarse(fs.z,0,child.y,k))) {
    assert ((1-coarse(cs,child.x,child.y)) && (1-coarse(cs,child.x,0,child.z)) &&
      (1-coarse(cs,0,child.y,child.z)) &&
      (1-coarse(cs,child.x,child.y,child.z)));
    // bilinear interpolation
    s[] = (27.*coarse(s) + 
     9.*(coarse(s,child.x) + coarse(s,0,child.y) +
         coarse(s,0,0,child.z)) + 
     3.*(coarse(s,child.x,child.y) + coarse(s,child.x,0,child.z) +
         coarse(s,0,child.y,child.z)) + 
     coarse(s,child.x,child.y,child.z))/64.;
  }
  else
    // tetrahedral interpolation
    s[] = (coarse(s) + coarse(s,child.x) + coarse(s,0,child.y) +
     coarse(s,0,0,child.z))/4.;
      }
      else if ((1-coarse(cs,child.x,child.y,child.z)) &&
         (((1-coarse(fs.z,child.x,child.y,k)) &&
     (((1-coarse(fs.x,i)) && (1-coarse(fs.y,child.x,j))) ||
      ((1-coarse(fs.y,0,j)) && (1-coarse(fs.x,i,child.y)))))
    ||
    ((1-coarse(fs.z,0,0,k)) &&
     (((1-coarse(fs.x,i,0,child.z)) && (1-coarse(fs.y,child.x,j,child.z))) ||
      ((1-coarse(fs.y,0,j,child.z)) && (1-coarse(fs.x,i,child.y,child.z)))))
    ||
    ((1-coarse(fs.z,child.x,0,k)) &&
     (1-coarse(fs.x,i)) && (1-coarse(fs.y,child.x,j,child.z)))
    ||
    ((1-coarse(fs.z,0,child.y,k)) &&
     (1-coarse(fs.y,0,j)) && (1-coarse(fs.x,i,child.y,child.z)))
    ))
  // diagonal interpolation
  s[] = (3.*coarse(s) + coarse(s,child.x,child.y,child.z))/4.;
#endif // dimension == 3
      else {
  // Pathological cases, use 1D gradients.
  s[] = coarse(s);
  foreach_dimension() {
    if ((1-coarse(fs.x,(child.x + 1)/2)) && (1-coarse(cs,child.x)))
      s[] += (coarse(s,child.x) - coarse(s))/4.;
    else if ((1-coarse(fs.x,(- child.x + 1)/2)) && (1-coarse(cs,- child.x)))
      s[] -= (coarse(s,- child.x) - coarse(s))/4.;
  }
      }
    }
  }
}

foreach_dimension()
static double concentration_gradient_x(Point point, scalar c, scalar t)
{
  static const double cmin = 1.0;
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

static inline void restriction_concentration (Point point, scalar s)
{
  scalar f = s.c;
  double sum = 0.;
  double sumf = 0.0;
  foreach_child(){
    sum += f[]*s[];
    sumf += f[];
  }
  s[] = sumf > 0 ? sum/sumf : 0.0;
}

static void theta_refine (Point point, scalar s)
{
  scalar f = s.c;
  if (cm[] == 0. || (!s.inverse && f[] <= 0.) || (s.inverse && f[] >= 1.))
    foreach_child()
      s[] = coarse(s);
  else {
    coord g;
    foreach_dimension()
      g.x = Delta*vof_concentration_gradient_x (point, f, s);
    double sc = s.inverse ? s[]/(1. - f[]) : s[]/f[], cmc = 4.*cm[];
    foreach_child() {
      s[] = sc;
      foreach_dimension()
  s[] += child.x*g.x*cm[-child.x]/cmc;
      s[] *= s.inverse ? 1. - f[] : f[];
    }
  }
}






struct Diffusion_three_field {
  // mandatory
  scalar f1;
  scalar f2;
  scalar f3;
  scalar f;
  scalar cs2;
  scalar cs3;
  face vector fs2;
  double dt;
  // optional
  face vector D1;  // default 1
  face vector D2;  // default 1
  face vector D3;  // default 1

  scalar r1, beta1; // default 0
  scalar r2, beta2; // default 0
  scalar r3, beta3; // default 0

  scalar theta1;   // default 1
  scalar theta2;   // default 1
  scalar theta3;   // default 1
};


trace
mgstats diffusion_three_field (struct Diffusion_three_field p){

  
 // If *dt* is zero we don't do anything. 

  if (p.dt == 0.) {
    mgstats s = {0};
    return s;
  } 

 
 // We define $f$ and $r$ for convenience. 
  
  scalar f1 = p.f1, r1 = automatic (p.r1);
  scalar f2 = p.f2, r2 = automatic (p.r2);
  scalar f3 = p.f3, r3 = automatic (p.r3);

  scalar f = p.f;
  scalar cs2 = p.cs2;
  scalar cs3 = p.cs3;
  face vector fs2 = p.fs2;

  scalar theta1 = p.theta1;
  scalar theta2 = p.theta2;
  scalar theta3 = p.theta3;

 
//  We define a (possibly constant) field equal to $\theta/dt$

  const scalar idt[] = - 1./p.dt;
  (const) scalar theta_idt1 = theta1.i ? theta1 : idt;
  (const) scalar theta_idt2 = theta2.i ? theta2 : idt;
  (const) scalar theta_idt3 = theta3.i ? theta3 : idt;

  if (theta1.i) {
    scalar theta_idt1 = theta1;
    foreach()
      theta_idt1[] *= idt[];
  }
  if (theta2.i) {
    scalar theta_idt2 = theta2;
    foreach()
      theta_idt2[] *= idt[];
  }
  if (theta3.i) {
    scalar theta_idt3 = theta3;
    foreach()
      theta_idt3[] *= idt[];
  }


  
 // We use `r` to store the r.h.s. of the Poisson--Helmholtz solver.

  if (p.r1.i)
    foreach()
      r1[] = theta_idt1[]*f1[] - r1[];
  else // r was not passed by the user
    foreach()
      r1[] = theta_idt1[]*f1[];
  if (p.r2.i)
    foreach()
      r2[] = theta_idt2[]*f2[] - r2[];
  else // r was not passed by the user
    foreach()
      r2[] = theta_idt2[]*f2[];
  if (p.r3.i)
    foreach()
      r3[] = theta_idt3[]*f3[] - r3[];
  else // r was not passed by the user
    foreach()
      r3[] = theta_idt3[]*f3[];
  
  
 // If $\beta$ is provided, we use it to store the diagonal term $\lambda$.

  scalar lambda1 = theta_idt1;
  scalar lambda2 = theta_idt2;
  scalar lambda3 = theta_idt3;

  scalar beta1 = p.beta1;
  scalar beta2 = p.beta2;
  scalar beta3 = p.beta3;

  if (beta1.i) {
    foreach()
      beta1[] += theta_idt1[];
    lambda1 = beta1;
  }
  if (beta2.i){
    foreach()
      beta2[] += theta_idt2[];
    lambda2 = beta2;
  }
  if (beta3.i) {
    foreach()
      beta3[] += theta_idt3[];
    lambda3 = beta3;
  }
  
 // Finally we solve the system.

  return poisson_three_field ({f1, f2, f3}, {r1, r2, r3}, cs2, cs3, fs2, f, p.D1, p.D2, p.D3, lambda1 , lambda2, lambda3);
}



attribute {
  scalar vol_frac;
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


//field 1 corresponds to solid, field 2 corresponds to phase f = 1, field 3 corresponds to f = 0
void Diffusion_three_field_setup(scalar f, scalar cs2, scalar cs3, face vector fs2, scalar theta1, scalar theta2, scalar theta3, double rhoCp1, double rhoCp2, double rhoCp3){

  vector n[];
  scalar alpha[];

  vector n2[];
  scalar alpha2[];

  if (f.height.x.i){
    heights (f, f.height);
  }
  boundary({f});
  restriction({f});

  reconstruction(f, n, alpha);
  reconstruction(cs, n2, alpha2);

  cs2.refine = cs2.prolongation = fraction_refine;
  cs3.refine = cs3.prolongation = fraction_refine;

  theta1.c = cs;
  theta1.inverse = true;
  theta1.refine = theta1.prolongation = theta_refine;

  theta2.c = cs2;
  theta2.inverse = false;
  theta2.refine = theta2.prolongation = theta_refine;

  theta3.c = cs3;
  theta3.inverse = false;
  theta3.refine = theta3.prolongation = theta_refine;

  //set up cs2 and fs2
  foreach(){
    if (cs[] >= 1){
      cs2[] = f[];
    }
    else{
      if (cs[] <= 0){
        cs2[] = 0;
      }
      else{
        /*to avoid doing complicated geometry at the contact line we will perform and approximation and make the interface with the maximum absolute value normal component aligned with the grid. 
        This allows us to use rectangle fraction to compute the volume fraction of all three phase and the surface moments  */

        if (f[] > 0 && f[] < 1){

          double maxx_n1 = 0;
          double maxx_n2 = 0;

          foreach_dimension(){
            maxx_n1 = max(fabs(n.x[]), maxx_n1);
            maxx_n2 = max(fabs(n2.x[]), maxx_n2);
          }
          if (maxx_n2 >= maxx_n1){
            coord n_temp = (coord){0,0,0};
            coord nn = (coord){n2.x[],n2.y[],n2.z[]};
            #if dimension == 2
              if (fabs(nn.x) >= fabs(nn.y)){
                n_temp.x = -sign(nn.x);
              }
              else{
                n_temp.y = -sign(nn.y);
              }
            #else // dimension == 3
              if (fabs(nn.x) >= fabs(nn.y)) {
                if (fabs(nn.x) >= fabs(nn.z))
                  n_temp.x = -sign(nn.x);
              }
              else{
                if (fabs(nn.y) >= fabs(nn.z)){
                  n_temp.y = -sign(nn.y);
                }
                else{
                  n_temp.z = -sign(nn.z);
                }
              }
            #endif // dimension == 3

                coord aa = (coord){-0.5,-0.5,-0.5};
                coord bb = (coord){0.5,0.5,0.5};

                foreach_dimension(){
                  if (n_temp.x != 0){
                    if (n_temp.x > 0){
                      aa.x += (1-cs[]);
                    }
                    else{
                      bb.x -= (1-cs[]);
                    }
                  }
                }
                coord n_temp2;
                foreach_dimension()
                  n_temp2.x = n.x[];
                cs2[] = min(clamp(rectangle_fraction (n_temp2, alpha[], aa, bb), 0, 1)*cs[], cs[]);
          }
          else{
            coord n_temp = (coord){0,0,0};
            coord nn = (coord){n.x[],n.y[],n.z[]};
            #if dimension == 2
              if (fabs(nn.x) >= fabs(nn.y)){
                n_temp.x = -sign(nn.x);
              }
              else{
                n_temp.y = -sign(nn.y);
              }
            #else // dimension == 3
              if (fabs(nn.x) >= fabs(nn.y)) {
                if (fabs(nn.x) >= fabs(nn.z))
                  n_temp.x = -sign(nn.x);
              }
              else{
                if (fabs(nn.y) >= fabs(nn.z)){
                  n_temp.y = -sign(nn.y);
                }
                else{
                  n_temp.z = -sign(nn.z);
                }
              }
            #endif // dimension == 3

                coord aa = (coord){-0.5,-0.5,-0.5};
                coord bb = (coord){0.5,0.5,0.5};

                foreach_dimension(){
                  if (n_temp.x != 0){
                    if (n_temp.x > 0){
                      aa.x += (1-f[]);
                    }
                    else{
                      bb.x -= (1-f[]);
                    }
                  }
                }
                coord n_temp2;
                foreach_dimension()
                  n_temp2.x = n2.x[];
                cs2[] = min(clamp(rectangle_fraction (n_temp2, alpha2[], aa, bb), 0, 1)*f[], cs[]);

          }
        }
        else{
          cs2[] = f[]*cs[];
        }
      }
    }
  }

  foreach_dimension(){
    fs2.x.refine = fs2.x.prolongation = face_fraction_refine_2_x;
    fs2.x.vol_frac = f;}

  foreach_face() {
    if ((f[-1] < VTOL) || (f[] < VTOL)){ // some cell is empty
      fs2.x[] = 0.;}
    else{ 
      if ((f[-1] > (1. - VTOL)) && (f[] > (1. - VTOL))){ // both cells are full
        fs2.x[] = fs.x[];}
      else {
        if ((cs[-1] > VTOL) &&  (cs[-1] < (1-VTOL)) & (cs[0] > VTOL) &&  (cs[0] < (1-VTOL)) ){
            double vleft = 1., vright = 1.;
            if (f[] < (1. - VTOL)) {

              double maxx_n1 = 0;
              double maxx_n2 = 0;

              foreach_dimension(){
                maxx_n1 = max(fabs(n.x[]), maxx_n1);
                maxx_n2 = max(fabs(n2.x[]), maxx_n2);
              }

              if (maxx_n2 >= maxx_n1){

                coord n_temp = (coord){0,0,0};
                coord nn;
                nn.x = n2.x[];
                nn.y = n2.y[];
                #if dimension == 3
                nn.z = n2.z[];
                #endif
                #if dimension == 2
                  if (fabs(nn.x) >= fabs(nn.y)){
                    n_temp.x = -sign(nn.x);
                  }
                  else{
                    n_temp.y = -sign(nn.y);
                  }
                #else // dimension == 3
                  if (fabs(nn.x) >= fabs(nn.y)) {
                    if (fabs(nn.x) >= fabs(nn.z))
                      n_temp.x = -sign(nn.x);
                  }
                  else{
                    if (fabs(nn.y) >= fabs(nn.z)){
                      n_temp.y = -sign(nn.y);
                    }
                    else{
                      n_temp.z = -sign(nn.z);
                    }
                  }
                #endif // dimension == 3
                coord n_temp2;
                n_temp2.x = n.x[];
                n_temp2.y = n.y[];
                n_temp2.z = n.z[];
                double alpha_temp = alpha[];
                alpha_temp -= n_temp2.x*(-0.5);
                n_temp2.x = 0.0;
                double summ = 0;
                foreach_dimension()
                  summ += fabs(n_temp2.x);
                foreach_dimension()
                  n_temp2.x /= summ;
                alpha_temp /= summ;
                coord aa = (coord){-0.5,-0.5,-0.5};
                coord bb = (coord){0.5,0.5,0.5};
                foreach_dimension(){
                  if (n_temp.x != 0){
                    if (n_temp.x > 0){
                      aa.x += (1-cs[]);
                    }
                    else{
                      bb.x -= (1-cs[]);
                    }
                  }
                }
                vleft = rectangle_fraction (n_temp2, alpha_temp, aa, bb)*fs.x[];
              }
              else{

                coord n_temp = (coord){0,0,0};
                coord nn;
                nn.x = n.x[];
                nn.y = n.y[];
                #if dimension == 3
                nn.z = n.z[];
                #endif
                #if dimension == 2
                  if (fabs(nn.x) >= fabs(nn.y)){
                    n_temp.x = -sign(nn.x);
                  }
                  else{
                    n_temp.y = -sign(nn.y);
                  }
                #else // dimension == 3
                  if (fabs(nn.x) >= fabs(nn.y)) {
                    if (fabs(nn.x) >= fabs(nn.z))
                      n_temp.x = -sign(nn.x);
                  }
                  else{
                    if (fabs(nn.y) >= fabs(nn.z)){
                      n_temp.y = -sign(nn.y);
                    }
                    else{
                      n_temp.z = -sign(nn.z);
                    }
                  }
                #endif // dimension == 3
                coord n_temp2;
                n_temp2.x = n2.x[];
                n_temp2.y = n2.y[];
                n_temp2.z = n2.z[];
                double alpha_temp = alpha2[];
                alpha_temp -= n_temp2.x*(-0.5);
                n_temp2.x = 0.0;
                double summ = 0;
                foreach_dimension()
                  summ += fabs(n_temp2.x);
                foreach_dimension()
                  n_temp2.x /= summ;
                alpha_temp /= summ;
                coord aa = (coord){-0.5,-0.5,-0.5};
                coord bb = (coord){0.5,0.5,0.5};
                foreach_dimension(){
                  if (n_temp.x != 0){
                    if (n_temp.x > 0){
                      aa.x += (1-f[]);
                    }
                    else{
                      bb.x -= (1-f[]);
                    }
                  }
                }
                vleft = rectangle_fraction (n_temp2, alpha_temp, aa, bb)*f[];

              }
            }
            if (f[-1] < (1. - VTOL)) {

              double maxx_n1 = 0;
              double maxx_n2 = 0;

              maxx_n1 = max(fabs(n.x[-1]), maxx_n1);
              maxx_n1 = max(fabs(n.y[-1]), maxx_n1);
              #if dimension == 3
              maxx_n1 = max(fabs(n.z[-1]), maxx_n1);
              #endif
              maxx_n2 = max(fabs(n2.x[-1]), maxx_n2);
              maxx_n2 = max(fabs(n2.y[-1]), maxx_n2);
              #if dimension == 3
              maxx_n2 = max(fabs(n2.z[-1]), maxx_n2);
              #endif

              if (maxx_n2 >= maxx_n1){

                coord n_temp = (coord){0,0,0};
                coord nn;
                nn.x = n2.x[-1];
                nn.y = n2.y[-1];
                #if dimension == 3
                nn.z = n2.z[-1];
                #endif
                #if dimension == 2
                  if (fabs(nn.x) >= fabs(nn.y)){
                    n_temp.x = -sign(nn.x);
                  }
                  else{
                    n_temp.y = -sign(nn.y);
                  }
                #else // dimension == 3
                  if (fabs(nn.x) >= fabs(nn.y)) {
                    if (fabs(nn.x) >= fabs(nn.z))
                      n_temp.x = -sign(nn.x);
                  }
                  else{
                    if (fabs(nn.y) >= fabs(nn.z)){
                      n_temp.y = -sign(nn.y);
                    }
                    else{
                      n_temp.z = -sign(nn.z);
                    }
                  }
                #endif // dimension == 3
                coord n_temp2;
                n_temp2.x = n.x[-1];
                n_temp2.y = n.y[-1];
                #if dimension == 3
                n_temp2.z = n.z[-1];
                #endif
                double alpha_temp = alpha[-1];
                alpha_temp -= n_temp2.x*(0.5);
                n_temp2.x = 0.0;
                double summ = 0;
                foreach_dimension()
                  summ += fabs(n_temp2.x);
                foreach_dimension()
                  n_temp2.x /= summ;
                alpha_temp /= summ;
                coord aa = (coord){-0.5,-0.5,-0.5};
                coord bb = (coord){0.5,0.5,0.5};
                double cs_temp = cs[-1];
                foreach_dimension(){
                  if (n_temp.x != 0){
                    if (n_temp.x > 0){
                      aa.x += (1-cs_temp);
                    }
                    else{
                      bb.x -= (1-cs_temp);
                    }
                  }
                }
                vright = rectangle_fraction (n_temp2, alpha_temp, aa, bb)*fs.x[];
              }
              else{

                coord n_temp = (coord){0,0,0};
                coord nn;
                nn.x = n.x[-1];
                nn.y = n.y[-1];
                #if dimension == 3
                nn.z = n.z[-1];
                #endif
                #if dimension == 2
                  if (fabs(nn.x) >= fabs(nn.y)){
                    n_temp.x = -sign(nn.x);
                  }
                  else{
                    n_temp.y = -sign(nn.y);
                  }
                #else // dimension == 3
                  if (fabs(nn.x) >= fabs(nn.y)) {
                    if (fabs(nn.x) >= fabs(nn.z))
                      n_temp.x = -sign(nn.x);
                  }
                  else{
                    if (fabs(nn.y) >= fabs(nn.z)){
                      n_temp.y = -sign(nn.y);
                    }
                    else{
                      n_temp.z = -sign(nn.z);
                    }
                  }
                #endif // dimension == 3
                coord n_temp2;
                n_temp2.x = n2.x[-1];
                n_temp2.y = n2.y[-1];
                #if dimension == 3
                n_temp2.z = n2.z[-1];
                #endif
                double alpha_temp = alpha2[-1];
                alpha_temp -= n_temp2.x*(0.5);
                n_temp2.x = 0.0;
                double summ = 0;
                foreach_dimension()
                  summ += fabs(n_temp2.x);
                foreach_dimension()
                  n_temp2.x /= summ;
                alpha_temp /= summ;
                coord aa = (coord){-0.5,-0.5,-0.5};
                coord bb = (coord){0.5,0.5,0.5};
                double f_temp = f[-1];
                foreach_dimension(){
                  if (n_temp.x != 0){
                    if (n_temp.x > 0){
                      aa.x += (1-f_temp);
                    }
                    else{
                      bb.x -= (1-f_temp);
                    }
                  }
                }
                vright = rectangle_fraction (n_temp2, alpha_temp, aa, bb)*f[-1];
              }
              
            }

            fs2.x[] = min(sqrt(vleft*vright), fs.x[]); 

        }
        else{
           double vleft = 1., vright = 1.;
            if (f[] < (1. - VTOL)) {
              coord m;
              m.x = n.x[];
              m.y = n.y[];
        
              #if dimension >= 3
              m.z = n.z[];
              #endif

              vleft = interface_fraction_x (m, alpha[], false);

            }
            if (f[-1] < (1. - VTOL)) {
              coord m;
              m.x = n.x[-1];
              m.y = n.y[-1];

              #if dimension >= 3
              m.z = n.z[-1];
              #endif

              vright = interface_fraction_x (m, alpha[-1], true);
        
            }
            fs2.x[] = min(sqrt(vleft*vright), fs.x[]);
        } 
      }
   }
  }

  //solid theta
  foreach(){
    if ((cs[] > VTOL) && (cs[] < (1-VTOL))){
      #if AXI
        coord n_temp = facet_normal (point, cs, fs);
        foreach_dimension()
          n_temp.x = -n_temp.x;
        double alpha_temp = plane_alpha (1-cs[], n_temp);
        coord p;
        #if dimension == 2
          line_center(n_temp, alpha_temp, 1-cs[], &p);
        #else
          plane_center(n_temp, alpha_temp, 1-cs[], &p);
        #endif
          theta1[] = max(1-cs[], VTOL)*rhoCp1*(y + p.y*Delta);
      #else
        theta1[] = max(1-cs[], VTOL)*rhoCp1;
      #endif
    }
    else{
      #if AXI
        theta1[] = max(1-cs[], VTOL)*rhoCp1*(y);
      #else
        theta1[] = max(1-cs[], VTOL)*rhoCp1;
      #endif
    }
  }

    //f theta2
  foreach(){
    if ((cs2[] > VTOL) && (cs2[] < (1-VTOL))){
      #if AXI
        coord n_temp = interface_normal (cs2);
        double alpha_temp = plane_alpha (cs2[], n_temp);
        coord p;
        #if dimension == 2
          line_center(n_temp, alpha_temp, cs2[], &p);
        #else
          plane_center(n_temp, alpha_temp, cs2[], &p);
        #endif
          theta2[] = max(cs2[], VTOL)*rhoCp2*(y + p.y*Delta);
      #else
        theta2[] = max(cs2[], VTOL)*rhoCp2;
      #endif
    }
    else{
      #if AXI
        theta2[] = max(cs2[], VTOL)*rhoCp2*(y);
      #else
        theta2[] = max(cs2[], VTOL)*rhoCp2;
      #endif
    }
  }

  foreach()
    cs3[] = clamp(cs[] - cs2[], 0, 1);

  foreach(){
    if ((cs3[] > VTOL) && (cs3[] < (1-VTOL))){
      #if AXI
        coord n_temp = interface_normal (cs3);
        double alpha_temp = plane_alpha (cs3[], n_temp);
        coord p;
        #if dimension == 2
          line_center(n_temp, alpha_temp, cs3[], &p);
        #else
          plane_center(n_temp, alpha_temp, cs3[], &p);
        #endif
          theta3[] = max(cs3[], VTOL)*rhoCp3*(y + p.y*Delta);
      #else
        theta3[] = max(cs3[], VTOL)*rhoCp3;
      #endif
    }
    else{
      #if AXI
        theta3[] = max(cs3[], VTOL)*rhoCp3*(y);
      #else
        theta3[] = max(cs3[], VTOL)*rhoCp3;
      #endif
    }
  }

}