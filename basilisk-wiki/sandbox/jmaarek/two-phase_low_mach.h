/**
# Two-phase interfacial flows

This file helps setup simulations for flows of two fluids separated by
an interface (i.e. immiscible fluids). It is typically used in
combination with a [Navier--Stokes solver](navier-stokes/centered.h). 

The interface between the fluids is tracked with a Volume-Of-Fluid
method. The volume fraction in fluid 1 is $f=1$ and $f=0$ in fluid
2. The densities and dynamic viscosities for fluid 1 and 2 are *rho1*,
*mu1*, *rho2*, *mu2*, respectively. */

stats statsf2 (scalar f)
{
  double min = 1e100, max = -1e100, sum = 0., sum2 = 0., volume = 0.;
  foreach(reduction(+:sum) reduction(+:sum2) reduction(+:volume)
    reduction(max:max) reduction(min:min)) 
    if (dv() > 0. && cs[] > 0.5 && f[] != nodata) {
      volume += dv();
      sum    += dv()*f[];
      sum2   += dv()*sq(f[]);
      if (f[] > max) max = f[];
      if (f[] < min) min = f[];
    }
  stats s;
  s.min = min, s.max = max, s.sum = sum, s.volume = volume;
  if (volume > 0.)
    sum2 -= sum*sum/volume;
  s.stddev = sum2 > 0. ? sqrt(sum2/volume) : 0.;
  return s;
}

double damping_factor = 0.4;

#if !EMBED
const scalar cs[] = 1.;
#endif

static double interpolate_linear_coord (Point point, scalar v, coord p)
{
#if dimension == 1
  int i = sign(p.x);
  /* linear interpolation */
  return v[]*(1. - fabs(p.x)) + v[i]*fabs(p.x);
#elif dimension == 2
  int i =  sign(p.x), j = sign(p.y);
  /* bilinear interpolation */
  return ((v[]*(1. - fabs(p.x)) + v[i]*fabs(p.x))*(1. - fabs(p.y)) + 
    (v[0,j]*(1. - fabs(p.x)) + v[i,j]*fabs(p.x))*fabs(p.y));
#else // dimension == 3
  int i = sign(p.x), j = sign(p.y), k = sign(p.z);
  /* trilinear interpolation */
  return (((v[]*(1. - fabs(p.x)) + v[i]*fabs(p.x))*(1. - fabs(p.y)) + 
     (v[0,j]*(1. - fabs(p.x)) + v[i,j]*fabs(p.x))*fabs(p.y))*(1. - fabs(p.z)) +
    ((v[0,0,k]*(1. - fabs(p.x)) + v[i,0,k]*fabs(p.x))*(1. - fabs(p.y)) + 
     (v[0,j,k]*(1. - fabs(p.x)) + v[i,j,k]*fabs(p.x))*fabs(p.y))*fabs(p.z));
#endif  
}


#if USE_SUTHERLAND
double Sutherland_mu(double T){
  return mu_star*pow(T/273, 1.5)*(273+110.5)/(T + 110.5);
}
#else
double Sutherland_mu(double T){
  return mu_star;
}
#endif

double Sutherland_k(double T){
  return Sutherland_mu(T)*gamma_h*RR/((gamma_h-1)*Pr);
}

/*
#if USE_SUTHERLAND
double Sutherland_mu(double T){
  return 8.411e-5*pow(T/273, 1.5)*(273+97)/(T + 97);
}
#else
double Sutherland_mu(double T){
  return 8.411e-5;
}
#endif

double Sutherland_k(double T){
  return 0.168*pow(T/273, 1.5)*(273+120)/(T + 120);
}*/

#include "contact_embed/vof_contact.h"
#include "diffusion.h"


scalar f[], * interfaces = {f};

double rho1 = 1., mu1 = 0., rho2 = 1., mu2 = 0.;
scalar T2[], T1[], frho2[];

/**
Auxilliary fields are necessary to define the (variable) specific
volume $\alpha=1/\rho$ as well as the cell-centered density. */

face vector alphav[];
scalar rhov[];
face vector lambdavv[];
face vector muv[];


event defaults (i = 0) {
  alpha = alphav;
  rho = rhov;
  lambdav = lambdavv;
  mu = muv;

  foreach_dimension()
    u.x.gradient = minmod;

  f.tracers = list_copy ({T1, frho2, T2});


  T1.inverse = false;
  T1.gradient = minmod;
  T1.NO_1D_COMPRESSION = false;
  T1.depends = list_add (T1.depends, f);

  T1.refine = T1.prolongation = vof_concentration_refine;
  T1.restriction = restriction_volume_average;

  frho2.inverse = true;
  frho2.gradient = minmod;
  frho2.NO_1D_COMPRESSION = true;

  frho2.refine = frho2.prolongation = vof_concentration_refine;
  frho2.restriction = restriction_volume_average;
  frho2.depends = list_add (frho2.depends, f);

  T2.inverse = true;
  T2.gradient = minmod;
  T2.NO_1D_COMPRESSION = true;

  T2.refine = T2.prolongation = vof_concentration_refine;
  T2.restriction = restriction_volume_average;
  T2.depends = list_add (T2.depends, f);

  T2.depends = list_add (T2.depends, f);
  T2.depends = list_add (T2.depends, frho2);

  /**
  We add the interface to the default display. */

  display ("draw_vof (c = 'f');");
}

/**
The density and viscosity are defined using arithmetic averages by
default. The user can overload these definitions to use other types of
averages (i.e. harmonic). */

#ifndef rho
# define rho(f, rho2) (clamp(f,0.,1.)*(rho1 - rho2) + rho2)
#endif
#ifndef mu
# define mu(f, T)  (clamp(f,0.,1.)*(mu1 - Sutherland_mu(T))) + Sutherland_mu(T)
#endif

/**
We have the option of using some "smearing" of the density/viscosity
jump. */

#if FILTERED
scalar sf[];
#else
# define sf f
#endif

event tracer_advection (i++)
{
  
  /**
  When using smearing of the density jump, we initialise *sf* with the
  vertex-average of *f*. */

#ifndef sf
#if dimension <= 2
  foreach()
    sf[] = (4.*f[] + 
      2.*(f[0,1] + f[0,-1] + f[1,0] + f[-1,0]) +
      f[-1,-1] + f[1,-1] + f[1,1] + f[-1,1])/16.;
#else // dimension == 3
  foreach()
    sf[] = (8.*f[] +
      4.*(f[-1] + f[1] + f[0,1] + f[0,-1] + f[0,0,1] + f[0,0,-1]) +
      2.*(f[-1,1] + f[-1,0,1] + f[-1,0,-1] + f[-1,-1] + 
    f[0,1,1] + f[0,1,-1] + f[0,-1,1] + f[0,-1,-1] +
    f[1,1] + f[1,0,1] + f[1,-1] + f[1,0,-1]) +
      f[1,-1,1] + f[-1,1,1] + f[-1,1,-1] + f[1,1,1] +
      f[1,1,-1] + f[-1,-1,-1] + f[1,-1,-1] + f[-1,-1,1])/64.;
#endif
#endif // !sf

#if TREE
  sf.prolongation = refine_bilinear;
  sf.dirty = true; // boundary conditions need to be updated
#endif
}

#include "fractions.h"

event properties (i++)
{
  scalar T_temp[], rho_temp[];
  T_temp.third = false;
  rho_temp.third = false;
  foreach(){
    rho_temp[] = max(frho2[], 0)/clamp(1-f[],1e-16,1.);
    T_temp[] = max(T2[], 0)/clamp(1-f[],1e-16,1.);
  }

  foreach_face() {
    if (fm.x[] > 0){
      double ff = (sf[] + sf[-1])/2.;
      double rho_ff = 0, T_ff = 0;
      if (f[-1] < 1 && f[] < 1){
        T_ff = face_value(T_temp,0);
        rho_ff = face_value(rho_temp,0);
      }
      else{
        if (f[-1] < 1){
          rho_ff = rho_temp[-1];
          T_ff = T_temp[-1];
        }
        if (f[] < 1){
            rho_ff = rho_temp[];
            T_ff = T_temp[];
        }
      }

      //fprintf(ferr, "test prop %g %g\n", rho(ff, rho_ff), mu(ff, T_ff) );

      alphav.x[] = fm.x[] > 0 ? fm.x[]/rho(ff, rho_ff) : 0.0;
      muv.x[] = fm.x[] > 0 ? fm.x[]*mu(ff, T_ff) : 0.0;
      lambdav.x[] = -2./3.*fm.x[]*(clamp(ff,0.,1.)*(0.0 - Sutherland_mu(T_ff)) + Sutherland_mu(T_ff));
    }
    else{
      alphav.x[] = 0;
      muv.x[] = 0.0;
      lambdav.x[] = 0.0;
    }
  }
  
  foreach(){
    rhov[] = cm[]*rho(sf[], max(frho2[], 0)/clamp(1-sf[],1e-16,1.));
  }

#if TREE  
  sf.prolongation = fraction_refine;
  sf.dirty = true; // boundary conditions need to be updated
  
#endif
}

static void momentum_refine (Point point, scalar u) {
  refine_bilinear (point, u);
  double rhou = 0.;
  foreach_child()
    rhou += cm[]*rho(f[], max(frho2[], 0)/clamp(1-f[],1e-16,1.))*u[];
  double du = u[] - rhou/((1 << dimension)*(cm[] + SEPS)*rho(f[], max(frho2[], 0)/clamp(1-f[],1e-16,1.)));
  foreach_child()
    u[] += du;
}

static void momentum_restriction (Point point, scalar u)
{
  double rhou = 0.;
  foreach_child()
    rhou += cm[]*rho(f[], max(frho2[], 0)/clamp(1-f[],1e-16,1.))*u[];
  u[] = rhou/((1 << dimension)*(cm[] + SEPS)*rho(f[], max(frho2[], 0)/clamp(1-f[],1e-16,1.)));
}


event defaults (i = 0)
{
  stokes = true;

  foreach_dimension() {
    u.x.depends = list_add (u.x.depends, f);
    u.x.depends = list_add (u.x.depends, frho2);
    u.x.gradient = minmod;
  }

#if TREE

  /**
  On trees, the refinement and restriction functions above rely on the
  volume fraction field *f* being refined/restricted before the
  components of velocity. To ensure this, we move *f* to the front of
  the field list (*all*). */

  int i = 0;
  while (all[i].i != f.i) i++;
  while (i > 0 && all[i].i)
    all[i] = all[i-1], i--;
  all[i] = f;
    
  /**
  We then set the refinement and restriction functions for the
  components of the velocity field. The boundary conditions on
  $\mathbf{u}$ now depend on those on $f$. */

  frho2.restriction = restriction_volume_average;
  frho2.refine = frho2.prolongation = vof_concentration_refine;

  T2.restriction = restriction_volume_average;
  T2.refine = T2.prolongation = vof_concentration_refine;
  
  foreach_dimension() {
    u.x.refine = u.x.prolongation = momentum_refine;
    u.x.restriction = momentum_restriction;
    u.x.depends = list_add (u.x.depends, f);
    u.x.depends = list_add (u.x.depends, T2);
  }
#endif
}

foreach_dimension()
static double boundary_q1_x (Point neighbor, Point point, scalar q1, bool * data)
{
  return clamp(f[],0.,1.)*rho1*u.x[];
}

foreach_dimension()
static double boundary_q2_x (Point neighbor, Point point, scalar q2, bool * data)
{
  return frho2[]*u.x[];
}


#if TREE
foreach_dimension()
static void prolongation_q1_x (Point point, scalar q1) {
  foreach_child()
    q1[] = clamp(f[],0.,1.)*rho1*u.x[];
}

foreach_dimension()
static void prolongation_q2_x (Point point, scalar q2) {
  foreach_child()
    q2[] = frho2[]*u.x[];
}
#endif


static scalar * interfaces1 = NULL;
const vector u_zero[] = {0,0,0};
(const) vector u_solid = u_zero;

scalar P_loc[];

event stability(i++){
  double dt_temp = HUGE;
  foreach(reduction(min:dt_temp)) {
    double div = 0.;
    foreach_dimension()
      div += uf.x[1] - uf.x[];
    div /= Delta;
    dt_temp = min(dt_temp, 0.1/max(fabs(div), 1e-16));
  }
  dtmax = min(dtmax, dt_temp);
}

event vof (i++) {

  //stats s1 = statsf2 (T2);
  //stats s2 = statsf2 (frho2);

  //fprintf(ferr, "T2 stats %g %g %g %g\n", s1.min, s1.max, s1.stddev, s1.volume);
  //fprintf(ferr, "frho2 stats %g %g %g %g\n", s2.min, s2.max, s2.stddev, s2.volume);

  // no slip boundary condition
  foreach(){
    T2[] = ( (cs[] > 0.0 && (1-f[]) > 0) ? frho2[]*T2[]*cs[]/(1-f[]) : 0.0); //conserves total enthalpy integrated in conservative form right now we assume c_p is constant
    frho2[] = ( cs[] > 0.0 ? frho2[]*cs[] : 0.0); //conserves total mass
    foreach_dimension()
      u.x[] = ( cs[] > 0.0 ? u.x[] : 0.0) ;
  }

  /**
  We allocate two temporary vector fields to store the two components
  of the momentum and set the boundary conditions and prolongation
  functions. The boundary conditions on $q_1$ and $q_2$ depend on the
  boundary conditions on $f$. */
  
  vector q1[], q2[];
  scalar f_cs[];

  for (scalar s in {q1,q2}) {
    s.depends = list_add (s.depends, f);
    s.depends = list_add (s.depends, frho2);
    foreach_dimension()
      s.v.x.i = -1; // not a vector
  }
  for (int i = 0; i < nboundary; i++)
    foreach_dimension() {
      q1.x.boundary[i] = boundary_q1_x;
      q2.x.boundary[i] = boundary_q2_x;
    }
#if TREE
  foreach_dimension() {
    q1.x.prolongation = prolongation_q1_x;
    q2.x.prolongation = prolongation_q2_x;
  }
  frho2.prolongation = vof_concentration_refine;
  T2.prolongation = vof_concentration_refine;
#endif

  /**
  We split the total momentum $q$ into its two components $q1$ and
  $q2$ associated with $f$ and $1 - f$ respectively. */

  foreach(){
    coord u_p;
    foreach_dimension()
      u_p.x = u.x[];
    if (cs[] > 0 && cs[] < 1){
      coord nn = facet_normal (point, cs, fs);
      double alpha = plane_alpha (cs[], nn);
      coord pp;
      plane_center (nn, alpha, cs[], &pp);
      u_p.x =  interpolate_linear_coord(point, u.x, pp);
      u_p.y =  interpolate_linear_coord(point, u.y, pp);
      #if dimension == 3
      u_p.z =  interpolate_linear_coord(point, u.z, pp);
      #endif
    }
    f_cs[] = clamp(f[],0,1)*cs[];
    foreach_dimension() {
      double fc = clamp(f[],0,1);
      q1.x[] = fc*rho1*u_p.x*cs[];
      q2.x[] = frho2[]*u_p.x*cs[];
    }
  }

  /**
  Momentum $q2$ is associated with $1 - f$, so we set the *inverse*
  attribute to *true*. We use the same slope-limiting as for the
  velocity field. */

  foreach_dimension() {
    q2.x.inverse = true;
    q2.x.NO_1D_COMPRESSION = true;
    q1.x.gradient = q2.x.gradient = u.x.gradient;
  }

  f_cs.gradient = zero;
  f_cs.NO_1D_COMPRESSION = false;

  /**
  We associate the transport of $q1$ and $q2$ with $f$ and transport
  all fields consistently using the VOF scheme. */
  
  //T2.gradient = zero;

  scalar * tracers = f.tracers;
  f.tracers = list_concat (tracers, (scalar *){q1, q2, f_cs});

  vof_advection ({f}, i);
  free (f.tracers);
  f.tracers = tracers;
  
  /**
  We recover the advected velocity field using the total momentum and
  the density */

  foreach(){
    foreach_dimension(){
       u.x[] = 0;
       if (cs[] > 0){
        u.x[] = ((rho1*f_cs[] + frho2[]) > 0 ? (q1.x[] + q2.x[])/(rho1*f_cs[] + frho2[])*cs[] + u_solid.x[]*(1-cs[]) : u_solid.x[]);  
       }
    }
      P_loc[] = ((cs[] > 0) ? T2[]*RR/cs[] : nodata);
      frho2[] = ((cs[] > 0) ? frho2[]/cs[] : 0);
      T2[] = ((frho2[]*cs[] > 0) ? T2[]/frho2[]/cs[]*(1-f[]) : 0.0);
    }
    
    //T2.gradient = minmod;
       

  /**
  We set the list of interfaces to NULL so that the default *vof()*
  event does nothing (otherwise we would transport $f$ twice). */
  
  interfaces1 = interfaces, interfaces = NULL;
}

/**
We set the list of interfaces back to its default value. */

event tracer_advection (i++) {
  interfaces = interfaces1;
}


double dP_dt = 0.0;
double dt1 = 0.0;
double dP_dt1 = 0.0;

event stability(i++){
  //fprintf(ferr, "stability %g %g\n ", dP_dt, 0.1*P_ref*gamma_h/max(fabs((dP_dt)), 1e-16)/(gamma_h-1));
  dtmax = min(0.1*P_ref*gamma_h/max(fabs((dP_dt)), 1e-16)/(gamma_h-1), dtmax);
}


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

event tracer_diffusion (i++) {

  //stats s1 = statsf2 (T2);
  //stats s2 = statsf2 (frho2);

  //fprintf(ferr, "T2 stats %g %g %g %g\n", s1.min, s1.max, s1.stddev, s1.volume);
  //fprintf(ferr, "frho2 stats %g %g %g %g\n", s2.min, s2.max, s2.stddev, s2.volume);

  P_ref += 0; //if there is a mass flux into the cavity then you would increase the thermodynamic pressure associated with the advection step

  face vector Diff[];

  scalar rho_temp2[];

   T2.refine = T2.prolongation = concentration_refine;

  foreach(){
    f[] = clamp(f[], 0, 1);
    frho2[] = (1-f[]) > 0 ? max(frho2[], 0) : 0.0;
    rho_temp2[] = cm[]*c_p*frho2[];
    T2[] = (1-f[]) > 0 ? max(T2[], 0)/(1-f[]) : 0.0;
  }


  foreach_face(){
      Diff.x[] = Sutherland_k(face_value(T2,0))*fm.x[];
  }

  scalar T_old[];
  scalar rr[];

  foreach(){
    T_old[] = T2[];
    if ((dt1 > 0.0) && (dP_dt1 != 0.0)){
      rr[] = cm[]*(dP_dt + (dP_dt-dP_dt1)*dt/dt1);}
    else{
      rr[] = dP_dt*cm[];
    }
  }

  double temp_r = 0;
  if (dt1 > 0.0 && dP_dt1 != 0.0){
      temp_r = (dP_dt + (dP_dt-dP_dt1)*dt/dt1);}
  else{
      temp_r = (dP_dt != nodata ? dP_dt : 0.0);
    }

  fprintf(ferr, "test 2 dP_dt %g %g\n", dP_dt, temp_r);

  mgstats mgT = diffusion(T2,dt,Diff, theta = rho_temp2, r = rr);

  mgpf = mgT;

  fprintf(ferr, "diffusion stats %d %d\n", mgT.i, mgT.nrelax);


  //divergence of k grad T = rho cp dT/dt - rr. strange equation comes from the fact that theta and rr are modified in diffusion()
  foreach()
    T_old[] = cs[] > 0 ? (T2[]-T_old[])*(rho_temp2[]/(-1./dt))/dt - (rho_temp2[]*T_old[]-rr[]) : 0.0;

  double temp1 = 0;
  double temp2 = 0;

  foreach(reduction(+:temp1) reduction(+:temp2)){
    temp1 += cs[] > 0 ? T_old[]*dv()/cm[] : 0;
    temp2 += cs[] > 0 ? dv() : 0;
  }
  //double dP_dt_temp = (gamma_h-1)*temp1/temp2;

  dP_dt1 = dP_dt;
  dP_dt = (gamma_h-1)*temp1/temp2;
  dt1 = dt;

  P_ref += dP_dt*dt;
  
  double P_V = 0; 
  double v_sum = 0; 

  double P_max = -HUGE; 
  double P_min = HUGE; 

  foreach(reduction(+:P_V) reduction(+:v_sum) reduction(max:P_max) reduction(min:P_min)){
    P_V += ((cs[] > 0) ? frho2[]*T2[]*RR*dv() : 0);
    v_sum += (1-f[])*dv();
    P_max = (((1-f[]) > 0 && cs[] >= 1) ? max(P_max, frho2[]*T2[]/(1-f[])) : P_max);
    P_min = (((1-f[]) > 0 && cs[] >= 1) ? min(P_min, frho2[]*T2[]/(1-f[])) : P_min);
  }

  fprintf (ferr, "pressure error %g %g %g\n", (P_ref-P_V/v_sum)/P_ref, P_max, P_min);

  P_ref = P_V/v_sum;

  T2.refine = T2.prolongation = vof_concentration_refine;


  foreach(){
    drhodt[] = cs[] > 0 ? RR/P_ref/c_p*(T_old[] - cm[]*dP_dt/(gamma_h-1)) + damping_factor/gamma_h/P_ref*(frho2[]*T2[]*RR*cm[] - P_ref*(1-f[])*cm[])/dt  : 0.0;
    //drhodt[] = 0;
    T2[] *= (1-f[]);
  }

  //s1 = statsf2 (T2);
  //s2 = statsf2 (frho2);

  //fprintf(ferr, "T2 stats %g %g %g %g\n", s1.min, s1.max, s1.stddev, s1.volume);
  //fprintf(ferr, "frho2 stats %g %g %g %g\n", s2.min, s2.max, s2.stddev, s2.volume);



  }
