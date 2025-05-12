/**
# Planetary Geostrophic model
   
   This code follows the implementation of
   [J. Callies](https://github.com/joernc/pgcm) of the planetary
   geostrophic equations [Samelson and Vallis 1997]. The equations are
   written in terrain following coordinates.
   
## Variables declaration
   
   Although the equations are 3d, we use a multilayer implentation
   because we need to compute the hydorstatic pressure (not easy to do
   with octree).
*/

#include "predictor-corrector.h"
#include "poisson.h"
#include "timestep.h"

// layered variables
scalar * wl  = NULL;
scalar * bl  = NULL;
scalar * pl  = NULL;
scalar * kfsl  = NULL;
//scalar * hdfl  = NULL;

/* TODO: should be face vector?: how to declare face vector pointer */
vector * ul = NULL;
vector * kfl  = NULL;
vector * hdiffl  = NULL;

scalar * evolving = NULL;

// vertical grid
int nl = 1;
double * sc;
double * sf;
double ds;

// physical parameters
double r = 0.1;
double a = 0.2;
double ys = 0; // southern latitude
double tend = 1; // end time


// forcing
scalar b_surf[];
double tau_surf = 1e-2;

// topography
scalar hc[];
vector hpc[];
face vector hf[];
tensor hpf[];

// 2d variables
scalar gg[];
scalar jebar[];
scalar psibt[];
scalar wind_effect[];

face vector u[];
face vector hdiff[];

// diffusion coefficient
face vector kf[];

// elliptic solver
mgstats mgD;
face vector ronh[];
vector fonh[];
double omega = 0.3;

char * outdir = "./";

// user defined topography and topography gradient
double h (double x, double y);
double hx (double x, double y);
double hy (double x, double y);

// user defined diffusivity coeff
double k (double x, double y, double s);

// user defined wind and wind derivative
double taux   (double x, double y);
double tauy   (double x, double y);
double taux_y (double x, double y);
double tauy_x (double x, double y);


/**
## Elliptic Solver 

   The routines for the elliptic solver are almost identical to the one
   in poisson.h. The main difference is the pseudo SOR method with the
   relaxation factor $\omega$. We added this because if the frictional
   term is weak the matrix wee want to invert is no longer diagonally
   dominant. 
*/

struct Btsolver {
  scalar a, b;
  (const) face vector alpha;
  (const)  vector beta;
  (const) double omega;
  double tolerance;
  int nrelax;
  scalar * res;
//  double minlevel;
};

static double residual_bt (scalar * al, scalar * bl, scalar * resl, void * data)
{
  scalar a = al[0], b = bl[0], res = resl[0];
  struct Btsolver * p = data;
  (const) face vector alpha = p->alpha;
  (const) vector beta = p->beta;
  double maxres = 0.;
  struct { double x, y; } f = {-1.,1.};
#if TREE
  /* conservative coarse/fine discretisation (2nd order) */
  face vector g[];
  foreach_face()
    g.x[] = alpha.x[]*(a[] - a[-1])/Delta;
  boundary_flux ({g});
  foreach (reduction(max:maxres)) {
    res[] = b[];
    foreach_dimension() {
      res[] += (g.x[] - g.x[1])/Delta;
      res[] += f.x*beta.y[]*0.5*(a[1] - a[-1])/Delta;
    }
    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);
  }
#else
  /* "naive" discretisation (only 1st order on trees) */
  foreach (reduction(max:maxres)) {
    res[] = b[];
    foreach_dimension(){
      res[] += ((alpha.x[1] + alpha.x[])*a[]
  		- alpha.x[1]*a[1] - alpha.x[]*a[-1])/sq(Delta);
      res[] += f.x*beta.y[]*0.5*(a[1] - a[-1])/Delta;
    }

    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);
  }
#endif
  boundary (resl);
  return maxres;
}


static void relax_bt (scalar * al, scalar * bl, int l, void * data)
{
  scalar a = al[0], b = bl[0];
  struct Btsolver * p = data;
  (const) face vector alpha = p->alpha;
  (const) vector beta = p->beta;
  (const) double omega = p->omega;
  struct { double x, y; } f = {-1.,1.};

  /**
     We use either Jacobi (under)relaxation or we directly reuse values
     as soon as they are updated. For Jacobi, we need to allocate space
     for the new field *c*. Jacobi is useful mostly as it gives results
     which are independent of the order in which the cells are
     traversed. This is not the case for the simple traversal, which
     means for example that results will depend on whether a tree or
     a multigrid is used (because cells will be traversed in a different
     order). The same comment applies to OpenMP or MPI parallelism. In
     practice however Jacobi convergence tends to be slower than simple
     reuse. */
  
#if JACOBI
  scalar c[];
#else
  scalar c = a;
#endif
  
  /**
     We use the face values of $\alpha$ to weight the gradients of the
     5-points Laplacian operator. We get the relaxation function. */

  foreach_level_or_leaf (l) {

    double n = - sq(Delta)*b[], d = 0.;
    foreach_dimension() {
      n += alpha.x[1]*a[1] + alpha.x[]*a[-1];
      d += alpha.x[1] + alpha.x[];
      n -= f.x*beta.y[]*0.5*(a[1] - a[-1])*Delta;
    }
    c[] = (1-omega)*c[] + omega*n/d;
  }

  /**
     For weighted Jacobi we under-relax by using a weight of 2/3. */
  
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
}

trace
mgstats btsolver (struct Btsolver p)
{

  face vector alpha = p.alpha;
  vector beta = p.beta;
  restriction ((scalar *){alpha,beta});

  /**
     If *tolerance* is set it supersedes the default of the multigrid
     solver. */

  double defaultol = TOLERANCE;
  if (p.tolerance)
    TOLERANCE = p.tolerance;

//  if (p.minlevel)
//    MINLEVEL = p.minlevel;

  scalar a = p.a, b = p.b;
  mgstats s = mg_solve ({a}, {b}, residual_bt, relax_bt, &p, p.nrelax, p.res);

  /**
     We restore the default. */

  if (p.tolerance)
    TOLERANCE = defaultol;

  return s;
}


/**
##  3d advection    

   Vertical boundary conditions (ghost points) needed for advection
   and horizontal diffusion because of the metric term in sigma
   coordinates
*/

void vertbc  (scalar * bl) 
{
  scalar b0 = bl[0];
  scalar b1 = bl[1];
  scalar b2 = bl[2];
  scalar b  = bl[3];
  foreach()
    b0[] = b[] - 3*b2[] + 3*b1[];
  boundary({b0});

  b0 = bl[nl+1];
  b1 = bl[nl];
  foreach()
    b0[] = 2*b_surf[] - b1[];
  boundary({b0});
}

/**
   We compute the horizontal velocities with the frictional thermal wind
   equation and the verticalvelocity with the incompressibility condition
*/

@def foreach_layer(l)
  foreach()
    for (int l = 1; l < nl ; l++) {
@
@define end_foreach_layer() } end_foreach();
  
trace
double velocity  (scalar * bl, vector * ul, scalar * wl, 
		  scalar psibt, double dtmax)
{

  vertbc(bl);
  /**
     p is the hydrostatic pressure at cell center. 
     $$
     \frac{\partial p}{\partial \sigma} = bh \, .
     $$
     We also compute $\gamma$
     $$
     \gamma = -\int_{-1}^0 \sigma b d\sigma
     $$
  */

  scalar p = pl[0];
  scalar b = bl[1];
  foreach(){
    p[] = b[]*hc[]*0.5*ds;
    gg[] = -b[]*ds*sc[0];
  }
  
  foreach_layer (l) {
    b = bl[l+1];
    scalar b0 = bl[l];
    p = pl[l];
    scalar p0 = pl[l-1];
    p[] = p0[] + 0.5*(b0[] + b[])*hc[]*ds;
    gg[] -= b[]*ds*sc[l];
  }

  for (int l = 0; l < nl ; l++) {
    p = pl[l];
    boundary ({p});
  }

  boundary ({gg,psibt});

  /**
     rhs for the barotropic stream function solver */
  foreach(){
    jebar[] = - hpc.y[]*(gg[1] - gg[-1])/(2*Delta)
      + hpc.x[]*(gg[0,1] - gg[0,-1])/(2*Delta)
      + wind_effect[];
  }
  
  /**
     The barotropic stream function is the solution of the elliptic equation
     $$
     \nabla\cdot (\alpha\nabla \psi) + \beta \cdot (\nabla \psi) = JEBAR
     $$ 
     with 
     $$
     \alpha = \frac{r}{h} \, , \quad \textrm{and} \quad 
     \beta = \left(\frac{\partial}{\partial y} \frac{y}{h}, 
     - \frac{\partial}{\partial x} \frac{y}{h}\right)
     $$
     The right hand side is the joint effect of baroclinicity and relief (JEBAR)
     $$
     JEBAR = J (h,\gamma) \, , \quad \textrm{with} \quad 
     J(a,b) = \partial_x a \partial_y b - \partial_y a \partial_x b
     $$


  */
//  mgD = btsolver(psibt,jebar,ronh,fonh,omega, tolerance = 1e-5, minlevel = 4);
  mgD = btsolver(psibt,jebar,ronh,fonh,omega, tolerance = 1e-5);
  
  /**
     compute baroclinic velocity at faces times metric 

     The baroclinic velocity is given by the frictional geostrophic balance
     $$
     -f v = -\frac{\partial p}{\partial x} + \sigma b h_x - r u
     $$
     and
     $$
     f u = -\frac{\partial p}{\partial y} + \sigma b h_y - r v
     $$

  */
  struct { double x, y; } f = {-1.,1.};

  face vector u0[];
  foreach_face()
    u0.x[] = 0.;

  foreach_face(){
    for (int l = 0; l < nl ; l++) {
      b = bl[l+1];
      p = pl[l];
      u = ul[l];
      
      u.x[] = (-r*((p[] - p[-1])/Delta - sc[l]*hpf.x.x[]*0.5*(b[] + b[-1]))
	       + f.x*y*(0.25*(p[0,1] - p[0,-1] + p[-1,1] - p[-1,-1])/Delta
			- sc[l]*hpf.x.y[]*0.5*(b[] + b[-1])))
	/(sq(r) + sq(y))
	*hf.x[];
      u0.x[] += u.x[];
    }
  }
  
  foreach_face()
    u0.x[] = u0.x[]/nl;
  
  /**
     The barotropic velocity is incorect: we subtract it and add the
     correct barotropic velocity from the elliptic equation  (metric
     term included) */
  foreach_face(){
    for (int l = 0; l < nl ; l++) {
      u = ul[l];
      u.x[] += -u0.x[]
        + f.x*(psibt[0,1] + psibt[-1,1] - psibt[0,-1] - psibt[-1,-1])/(4.*Delta);
    }
  }

  for (int l = 0; l < nl ; l++) {
    u = ul[l];
    boundary ((scalar *){u});
    dtmax = timestep (u, dtmax);
  }
  
  /**
     We know the horizontal velocity; we now use the continuity
     equation (incompressibility) to get the vertical velocity. The
     incompressibility condition is
     $$
     \frac{\partial }{\partial x} (h u) + \frac{\partial }{\partial y} (h v) + 
     \frac{\partial }{\partial \sigma} (h w) = 0
     $$

     compute vertical velocity (metric included) */
  foreach(){
    for (int l = 0; l < nl ; l++) {
      scalar w0 = wl[l];
      scalar w = wl[l+1];
      u = ul[l];
      w[] = w0[] -(u.x[1] - u.x[] + u.y[0,1] - u.y[])*ds/Delta;
    }
  }

  return dtmax;
}


/**
   the advection is a 3d flux divergence: 
   $$
   \nabla (u b) = \frac{1}{h}\left( \frac{\partial }{\partial x} (h u b) 
   + \frac{\partial }{\partial y} (h v b) + 
   \frac{\partial }{\partial \sigma} (h w b)  \right)
   $$
*/

  // TODO: add upper Ekman layer

trace
void advection  (scalar * bl, vector * ul, scalar * wl, scalar * dbl)
{  

  vertbc(bl);
  foreach() {
    for (int l = 0; l < nl ; l++) {
            
      scalar b = bl[l+1];
      scalar b0 = bl[l];
      scalar b1 = bl[l+2];
      scalar db = dbl[l+1];
      scalar w0 = wl[l];
      scalar w = wl[l+1];
      u = ul[l];
      
      db[] += (((b[] + b[-1,0])*u.x[] -
      	       (b[] + b[1,0] )*u.x[1,0] +
      	       (b[] + b[0,-1])*u.y[] -
      	       (b[] + b[0,1] )*u.y[0,1]
      	       )/(2.*Delta) +
      	      ((b[] + b0[])*w0[] -
      	       (b[] + b1[])*w[])/(2.*ds)
      	      )/hc[];

    }
  }
}

/**
## diffusion

   The horizontal diffusion is computed explicitely and the vertical
   diffusion is computed implicitly. The treatment of vertical boundaries
   (air-sea interface and sea floor) is handled with ghots points which
   are adjusted in vertbc 

   The vertical diffustion flux is
   $$
   -\kappa \frac{1 + \sigma^2(h_x^2 + h_y^2)}{h^2} \frac{\partial b}{\partial \sigma}
   $$
   is treated here implicitely 

*/

trace
void vdiff_implicit(scalar * bl, double dt)
{
  double ad[nl], bd[nl], cd[nl], rhs[nl];
 
  /**
     surface BC */
  scalar b0 = bl[nl];
  scalar kfs = kfsl[nl];
  foreach()
    b0[] += dt*2*kfs[]/(sq(ds*hc[]))*b_surf[];
  
  foreach() {
    for (int l = 0; l < nl; l++) {
      b0 = bl[l+1];
      rhs[l] = b0[];
    }
    
    kfs = kfsl[1];
    double K0 = 0;
    double K1 = kfs[]*(1+sq(a)*sq(sf[1])*(sq(hpc.x[])+sq(hpc.y[])))/(sq(hc[]));
    
    for (int l = 0; l < nl-1; l++) {
      
      ad[l] = - dt*K0/sq(ds);
      cd[l] = - dt*K1/sq(ds);
      bd[l] = 1. - ad[l] - cd[l];
      
      kfs = kfsl[l+1];
      K0 = K1;
      K1 = kfs[]*(1+sq(a)*sq(sf[l+1])*(sq(hpc.x[])+sq(hpc.y[])))/(sq(hc[]));
      
    }
    // upper layer
    ad[nl-1] = - dt*K0/sq(ds);
    bd[nl-1] = 1 + dt*K0/sq(ds) + 2*dt*K1/sq(ds);
    
    /**
       We can now solve the tridiagonal system using the [Thomas
       algorithm](https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm). */
    
    for (int l = 1; l < nl; l++) {
      bd[l] -= ad[l]*cd[l-1]/bd[l-1];
      rhs[l] -= ad[l]*rhs[l-1]/bd[l-1];
    }
    
    b0 = bl[nl];
    b0[] = ad[nl-1] = rhs[nl-1]/bd[nl-1];
    for (int l = nl - 2; l >= 0; l--) {
      b0 = bl[l+1];
      b0[] = ad[l] = (rhs[l] - cd[l]*ad[l+1])/bd[l];
    }
  }
}


/**
   There are three terms left for the diffusion operator 
   on the $x$ coordinate:
   $$
   -\kappa \left( \frac{\partial b}{\partial x} 
   - \frac{\sigma h_x}{h} \frac{\partial b}{\partial \sigma}\right)
   $$
   on $y$:
   $$
   -\kappa \left( \frac{\partial b}{\partial y} 
   - \frac{\sigma h_y}{h} \frac{\partial b}{\partial \sigma}\right)
   $$
   and on $\sigma$:
   $$
   -\kappa \left( - \frac{\sigma h_x}{h} \frac{\partial b}{\partial x}
   - \frac{\sigma h_y}{h} \frac{\partial b}{\partial y} \right)
   $$
*/

trace
void hdiffusion  (scalar * bl, scalar * dbl) 
{
  vertbc(bl);

  foreach_face() {
    for (int l = 0; l < nl ; l++) {
      scalar b = bl[l+1];
      scalar b0 = bl[l];
      scalar b1 = bl[l+2];
      hdiff = hdiffl[l];
      kf = kfl[l];

      hdiff.x[] = sq(a)*kf.x[]*(hf.x[]*(b[] - b[-1])/Delta - 
				sc[l]*hpf.x.x[]*0.25*(b1[]-b0[] + b1[-1]-b0[-1])/ds);
    }
  }
  /**
     no flux side boundary conditions */

  for (int l = 0; l < nl ; l++) {
    hdiff = hdiffl[l];
    boundary ((scalar *){hdiff});
  }

  /**
     and divergence of the flux (horizontal part)*/
  foreach(){
    for (int l = 0; l < nl ; l++) {
      scalar db = dbl[l+1];
      hdiff = hdiffl[l];
      db[] += (hdiff.x[1,0]- hdiff.x[] + hdiff.y[0,1] - hdiff.y[])/(hc[]*Delta);
    }
  }

  /**
     cross terms in the vertical gradient */
  double vflux[nl+1];
  vflux[0] = 0.; // no flux bottom
  foreach() {
    for (int l = 1; l < nl+1 ; l++) {
      scalar b = bl[l];
      scalar b1 = bl[l+1];
      scalar kfs = kfsl[l];
//      scalar hdf = hdfl[l];

      vflux[l] = -sq(a)*kfs[]*sf[l]*(hpc.x[]*0.25*(b[1]-b[-1] + b1[1]-b1[-1])/Delta
				    + hpc.y[]*0.25*(b[0,1]-b[0,-1] + b1[0,1]-b1[0,-1])/Delta) ;
      scalar db = dbl[l];
      db[] += (vflux[l] - vflux[l-1])/(ds*hc[]);
    }
  }

  /* /\** */
  /*    and divergence of the flux (vertical part)*\/ */
  /* // keep hdfl global variable for lower and upper flux BC */
  /* foreach(){ */
  /*   for (int l = 0; l < nl ; l++) {     */
  /*     scalar db = dbl[l+1]; */
  /*     scalar hdf = hdfl[l]; */
  /*     scalar hdf0 = hdfl[l+1]; */
  /*     db[] += (hdf0[] - hdf[])/(ds*hc[]); */
  /*   } */
  /* } */
}

/**
## Convection
Restore unstable stratification to neutral profile */

void convection(scalar * bl)
{
  foreach() {
    /**
       Upper layer instability criterion convective instability occurs
       if b[l] < b[l-1]. If the upper layer is stable, we assume the
       entire water column is stable. */
    int l = nl;
    scalar b1 = bl[l];
    scalar b0 = bl[l-1];

    while( b1[] - b0[] < 0 && l > 2) {
      /**
  	 layer equaly spaced: standard averaging */
      b1[] = 0.5*(b0[] + b1[]);
      b0[] = b1[];

      l--;
      b1 = bl[l];
      b0 = bl[l-1];
    }
    /**
     Lower layer case*/
    if( b1[] - b0[] < 0 && l == 2) {
      b1[] = 0.5*(b0[] + b1[]);
      b0[] = b1[];
    }
  }
}

/**
## surface forcing
 */

void forcing_implicit(scalar * bl, double dt)
{

  scalar b = bl[nl];
  foreach() 
    b[] = (b_surf[]*dt + b[]*tau_surf)/(dt + tau_surf);
}

/**
## time stepping routines 

   We use the predictor corrector implementation */

static void advance_pg (scalar * output, scalar * input, 
			scalar * updates, double dt)
{
  foreach() {
    for (int l = 0; l < nl ; l++) {
      scalar bi = input[l+1];
      scalar bo = output[l+1];
      scalar db = updates[l+1];
      bo[] = bi[] + db[]*dt;
    }
  }
  forcing_implicit(output,dt);
  vdiff_implicit(output,dt);
  convection(output);

  for (int l = 0; l < nl ; l++) {
    scalar b  = output[l+1];
    boundary ({b});
  }
}

double update_pg (scalar * evolving, scalar * updates, double dtmax)
{
  foreach()
    for (scalar s in updates)
      s[] = 0.;

  dtmax = velocity(evolving, ul, wl, psibt, dtmax);
  advection(evolving, ul, wl, updates);
//  adjust_kv(evolving);
  hdiffusion(evolving,updates);

  return dtmax;
}

/**
## Declaration of layered variables and boundary conditions 

we put it in a function in order to be able to call it as an external
module */

void set_vars()
{
  assert (wl     == NULL);
  assert (bl     == NULL);
  assert (pl     == NULL);
  assert (ul     == NULL);
  assert (hdiffl == NULL);
//  assert (hdfl  == NULL);
  assert (nl > 0);

  //TODO: use nl+2 for everyone for simplicity
  for (int l = 0; l < nl; l++) {
    scalar p = new scalar;
    pl = list_append (pl, p);

    face vector u = new face vector;
    ul = vectors_append (ul, u);
    face vector kf = new face vector;
    kfl = vectors_append (kfl, kf);
    face vector hdiff = new face vector;
    hdiffl = vectors_append (hdiffl, hdiff);

    scalar w = new scalar;
    wl = list_append (wl, w);
    scalar kfs = new scalar;
    kfsl = list_append (kfsl, kfs);
//    scalar hdf = new scalar;
//    hdfl = list_append (hdfl, hdf);

    sprintf(u.x.name, "%s%d", "u",l);
    sprintf(u.y.name, "%s%d", "v",l);
    sprintf(w.name, "%s%d", "w",l);
  }
  // w: nl+ 1
  scalar w = new scalar;
  wl = list_append (wl, w);
  scalar kfs = new scalar;
  kfsl = list_append (kfsl, kfs);
//  scalar hdf = new scalar;
//  hdfl = list_append (hdfl, hdf);

  // b,dbl: nl+ 2
  for (int l = 0; l < nl+2; l++) {
    scalar b = new scalar;
    bl  = list_append (bl, b);
    sprintf(b.name, "%s%d", "b",l);
  }

  evolving = bl;

  // init vertical grid
  ds = 1./nl;
  sc = malloc (nl*sizeof(double));
  sf = malloc ((nl+1)*sizeof(double));
  
  sf[0] = -1.0;
  for (int l = 1; l < nl+1; l++)
    sf[l] = sf[l-1] + ds;

  sc[0] = -1.0 + 0.5*ds;
  for (int l = 1; l < nl; l++)
    sc[l] = sc[l-1] + ds;


  /**
     Initialize variables */
  foreach() {
    for (scalar w in wl)
      w[] = 0.0;
    for (scalar b in bl)
      b[] = 0.0;
//    for (scalar hdf in hdfl)
//      hdf[] = 0.0;
  }
  foreach_face() {
    for (vector u in ul)
      u.x[] = 0.0;
  }  


  /**
     topography at cell center and cell face and coefs fonh (beta) and
     ronh (alpha) for the elliptic equation */

  /* TODO: define refine and prolongations (cf axi.h exemple) of all
     these metric*/
  // cell center
  foreach(){
    hc[] = h(x,y);
    hpc.x[] = hx(x,y);
    hpc.y[] = hy(x,y);
  }
  // cell face
  foreach_face(){
    hf.x[] = h(x,y);
  }
  // derivative at cell face
  foreach_face(x) {
    hpf.x.x[] = hx(x,y) ;
    hpf.x.y[] = hy(x,y) ;
  }
  foreach_face(y) {
    hpf.y.x[] = hx(x,y) ;
    hpf.y.y[] = hy(x,y) ;
  }

  // elliptic solver
  foreach(){
    fonh.x[] = -hpc.x[]*y/sq(hc[]);
    fonh.y[] = (hc[] - hpc.y[]*y )/sq(hc[]);
  }
  foreach_face(){
    ronh.x[] = r/hf.x[];
  }

  foreach(){
    wind_effect[] = taux_y(x,y)/hc[] - taux(x,y)*hpc.y[]/sq(hc[])
      - tauy_x(x,y)/hc[] + tauy(x,y)*hpc.x[]/sq(hc[]);
  }

  // diffusivity coeff at faces
  foreach_face() {
    for (int l = 0; l < nl ; l++) {
      kf = kfl[l];
      kf.x[] = k(x,y,sc[l]);
    }
  }

  foreach(){
    for (int l = 0; l < nl+1 ; l++) {
      scalar kfs = kfsl[l];
      kfs[] = k(x,y,sf[l]);
    }
  }

  /**
     boundary conditions: we use the default BC (no flux) even though
     the code is written in terrain following coordinates. So strickly
     speaking the value of the scalar at the ghost point involves
     vertical derivatives. We can get away with it because the
     topography vanishes next to the boundary. For the barotropic stream
     function we use the BC $\psi = 0$. */

  psibt[right]  = dirichlet(0);
  psibt[left]   = dirichlet(0);
  psibt[top]    = dirichlet(0);
  psibt[bottom] = dirichlet(0);

  //TODO: these BC seem to be the default but it seems that I still
  //need to specify it

  for (vector u in ul) {
   u.n[right] = 0.0;
   u.n[left] = 0.0;
   u.n[top] = 0.0;
   u.n[bottom] = 0.0;
  }

  for (vector hdiff in hdiffl) {
   hdiff.n[right] = 0.0;
   hdiff.n[left] = 0.0;
   hdiff.n[top] = 0.0;
   hdiff.n[bottom] = 0.0;
  }

  /**
     We overload the default 'advance' and 'update' functions of the
     predictor-corrector scheme and (TODO) setup the prolongation and
     restriction methods on trees. */

  advance = advance_pg;
  update = update_pg;

}

event defaults (i = 0){
  set_vars();
}


/**
The event below will happen after all the other initial events to take
into account user-defined field initialisations. */

event init (i = 0)
{
  boundary (all);
}

/**
## Cleanup */

void trash_vars(){
  free (wl), wl = NULL;
  free (bl), bl = NULL;
  free (pl), pl = NULL;
  free (ul), ul = NULL;
  free (kfsl), kfsl = NULL;
  free (kfl), kfl = NULL;
//  free (hdfl), hdfl = NULL;
  free (hdiffl), hdiffl = NULL;
  free (sc);
  free (sf);
}

event cleanup (i = end, last) {
  trash_vars();
}

/* #ifdef H5FILE_NAME */
/* #include "output_xmf.h" */
/* void backup_fields (scalar * sfields, vector * vfields, int nf) */
/* { */
/*   FILE * fp ; */
/*   char name[80], stamp[1024];  */
/*   sprintf(stamp, "%6.6i", nf); */

/*   nf > 0 ? sprintf(name, "fields_%6.6d.xmf", nf) : sprintf(name, "fields.xmf"); */
/*   fp = fopen(name, "w"); output_xmf_h5_foreach (sfields, vfields, N, fp, stamp); fclose (fp); */
/* } */
/* #endif */
