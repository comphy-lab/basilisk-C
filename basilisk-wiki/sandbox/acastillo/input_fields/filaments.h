/** 
\renewcommand*{\vec}[1]{\boldsymbol{#1}}
\newcommand*{\uvec}[1]{\boldsymbol{\hat{#1}}}
# Curved vortex filaments

This section borrows from the development for a general curved vortex presented
by [Callegari and Ting (1978)](#callegari1978) and notes from the book by
Maurice Rossi, Stéphane Le Dizès, and others.

<center>
![Geometrical elements for a curved line (a) in the osculating plane (b) 
in the plane orthogonal to the vortex axis](rossi1.png)
</center>

Let us define a curvilinear orthogonal basis associated with a reference
**space-curve** $\mathcal{C}(s,t)$. Here $s$ is the **arc-length** along the
initial curve measured from a given reference point. 

Each point along $\mathcal{C}(s,t)$ is characterized by its **position vector**
$\vec{x}_c$, local **curvature** $\kappa(s,t)$, local **torsion** $\tau(s,t)$,
and the corresponding **Frenet–Serret** frame composed of three vectors
$(\uvec{t}, \uvec{n}, \uvec{b})$
$$
\begin{aligned}
  \uvec{t}(s) \equiv \frac{d\vec{x}_c}{ds}, \quad
  \uvec{n}(s) \equiv \frac{1}{\kappa}\frac{d\uvec{T}}{ds}, \quad
  \uvec{b}(s) \equiv \uvec{T}\times\uvec{N}
\end{aligned}
$$
which are related to one another as:
$$
\begin{aligned}
\frac{d}{ds}
\begin{bmatrix}
\uvec{t} \\ \uvec{n} \\ \uvec{b}
\end{bmatrix}
=
\begin{bmatrix}
0 & \kappa(s) & 0 \\
-\kappa(s) & 0 & \tau(s) \\
0 & -\tau(s) & 0 
\end{bmatrix}
\begin{bmatrix}
\uvec{t} \\ \uvec{n} \\ \uvec{b}
\end{bmatrix}
\end{aligned}
$$

The position of a given point $M$ can be expressed as two coordininates in the
plane $\mathcal{A}(s,t)$ locally perpendicular to $\mathcal{C}(s,t)$ and one
coordinate $s$ along $\mathcal{C}(s,t)$. This way, we may introduce a local
**radial** and **angular** coordinates $(\rho,\varphi)$ and write the position
vector of a point $M$ as:
$$
\begin{aligned}
  \vec{x} = \vec{x}_c (s) + \rho\vec{e}_\rho(s,\varphi)
\end{aligned}
$$
where the unit vectors $\vec{e}_\rho$ and $\vec{e}_\varphi$ are defined as
$$
\begin{aligned}
  \uvec{e}_\rho &=& \cos\varphi\uvec{n} + \sin\varphi\uvec{b}
  \\
  \uvec{e}_\varphi &=& -\sin\varphi\uvec{n} + \cos\varphi\uvec{b}
\end{aligned}
$$

Note that $(\rho,\varphi,s)$ are curvilinear coordinates but they are not
orthogonal in general. However, if we define an angle $\theta$
$$
\begin{aligned}
  \theta = \varphi - \theta_0(s), \quad \text{ where } \quad
  \frac{d\theta_0}{ds} = -\tau(s)
\end{aligned}
$$
then $(\rho,\theta,s)$ form a set of curvilinear orthogonal coordinates. In
this system, the infinitesimal increment reads as
$$
\begin{aligned}
  d\vec{x} &=&  d\rho ~\uvec{e}_\rho(s,\varphi)
      + \rho d\varphi ~\uvec{e}_\varphi(s,\varphi) + h_s ds ~\uvec{T}(s)
\end{aligned}
$$
where 
$$
\begin{aligned}
  h_s = (1-\kappa(s)\rho\cos(\theta+\theta_0(s)))
\end{aligned}
$$

\br

# Representation of vortex filaments in Basilisk
### *vortex_filament*: defines the properties of a vortex filament

This structure represents a vortex filament. It includes properties such as the
number of segments, a scalar quantity, and various vectors representing the
geometry and orientation of the filament.

*nseg*
: number of segments

*t*
: arbitrary of parameter, usually $\xi_n$

*C*
: Cartesian coordinates describing the curve, $\vec{x}_c \in \mathcal{C}(s(\xi_n))$

*Tvec*
: tangent unit vector, $\uvec{T}(s(\xi_n))$

*Nvec*
: normal unit vector, $\uvec{N}(s(\xi_n))$

*Bvec*
: binormal unit vector, $\uvec{B}(s(\xi_n))$

*s*
: arclenght coordinate, $s(\xi_n)$ or $\ell(\xi_n)$

*pcar*
: placeholder for a Cartesian coordinate 

*sigma*
: local stretching factor, $\sigma(s(\xi_n))$  (optional)

*kappa*
: local curvature, $\kappa(s(\xi_n))$  (optional)

*tau*
: local torsion, $\tau(s(\xi_n))$  (optional)

*theta0*
: cumulative torsion, $\varphi_0(s(\xi_n))$  (optional)

*/

#include "PointTriangle.h"
#pragma autolink -lgsl -lgslcblas
#include <gsl/gsl_spline.h>

struct filament_splines {
  gsl_interp_accel *acc;
  gsl_spline *C[3];
  gsl_spline *Tvec[3];
  gsl_spline *Nvec[3];
  gsl_spline *Bvec[3];
  gsl_spline *s;
  gsl_spline *sigma;
  gsl_spline *kappa;
  gsl_spline *tau;
  gsl_spline *theta0;
  gsl_spline *a;
};

struct vortex_filament{
  int     nseg;     // number of segments
  double* phi;    // arbitrary parametrisation of C  
  coord*  C;        // cartesian coords of filament
  coord*  Tvec;     // unit tangent vector
  coord*  Nvec;     // unit normal vector
  coord*  Bvec;     // unit binormal vector
  double* s;        // arc-length coordinate
  coord   pcar;     // current point in Cartesian coordinates  
  double* sigma;    // Stretching factor
  double* kappa;    // Curvature
  double* tau;      // Torsion
  double* theta0;  // Cumulative torsion
  double* a;        // Core size
  coord*  dC;       // 
  coord*  d2C;      // 
  coord*  d3C;      // 
  coord*  Ulocal;   //
  coord*  Uauto;    // 
  coord*  Umutual;  //
  coord*  U;        // 
  coord*  Uprev;    // 
  double* vol;      // Initial volume 
  void*   splines;  // Cached splines (filament_splines struct pointer)
};

// Function to free memory for the cached splines
void free_filament_splines(void* splines_ptr) {
  if (splines_ptr == NULL) return;
  struct filament_splines* s = (struct filament_splines*)splines_ptr;
  for (int d = 0; d < 3; d++) {
    if (s->C[d]) gsl_spline_free(s->C[d]);
    if (s->Tvec[d]) gsl_spline_free(s->Tvec[d]);
    if (s->Nvec[d]) gsl_spline_free(s->Nvec[d]);
    if (s->Bvec[d]) gsl_spline_free(s->Bvec[d]);
  }
  if (s->s) gsl_spline_free(s->s);
  if (s->sigma) gsl_spline_free(s->sigma);
  if (s->kappa) gsl_spline_free(s->kappa);
  if (s->tau) gsl_spline_free(s->tau);
  if (s->theta0) gsl_spline_free(s->theta0);
  if (s->a) gsl_spline_free(s->a);
  if (s->acc) gsl_interp_accel_free(s->acc);
  free(s);
}

// Function to allocate and initialize splines for the filament
void setup_filament_splines(struct vortex_filament* filament) {
  if (filament == NULL) return;
  if (filament->splines != NULL) {
    free_filament_splines(filament->splines);
    filament->splines = NULL;
  }

  struct filament_splines* s = malloc(sizeof(struct filament_splines));
  s->acc = gsl_interp_accel_alloc();

  double* temp = malloc(filament->nseg * sizeof(double));

  for (int d = 0; d < 3; d++) {
    // C
    for (int i = 0; i < filament->nseg; i++) {
      temp[i] = (d == 0 ? filament->C[i].x : (d == 1 ? filament->C[i].y : filament->C[i].z));
    }
    s->C[d] = gsl_spline_alloc(gsl_interp_cspline, filament->nseg);
    gsl_spline_init(s->C[d], filament->phi, temp, filament->nseg);

    // Tvec
    for (int i = 0; i < filament->nseg; i++) {
      temp[i] = (d == 0 ? filament->Tvec[i].x : (d == 1 ? filament->Tvec[i].y : filament->Tvec[i].z));
    }
    s->Tvec[d] = gsl_spline_alloc(gsl_interp_cspline, filament->nseg);
    gsl_spline_init(s->Tvec[d], filament->phi, temp, filament->nseg);

    // Nvec
    for (int i = 0; i < filament->nseg; i++) {
      temp[i] = (d == 0 ? filament->Nvec[i].x : (d == 1 ? filament->Nvec[i].y : filament->Nvec[i].z));
    }
    s->Nvec[d] = gsl_spline_alloc(gsl_interp_cspline, filament->nseg);
    gsl_spline_init(s->Nvec[d], filament->phi, temp, filament->nseg);

    // Bvec
    for (int i = 0; i < filament->nseg; i++) {
      temp[i] = (d == 0 ? filament->Bvec[i].x : (d == 1 ? filament->Bvec[i].y : filament->Bvec[i].z));
    }
    s->Bvec[d] = gsl_spline_alloc(gsl_interp_cspline, filament->nseg);
    gsl_spline_init(s->Bvec[d], filament->phi, temp, filament->nseg);
  }

  s->s = gsl_spline_alloc(gsl_interp_cspline, filament->nseg);
  gsl_spline_init(s->s, filament->phi, filament->s, filament->nseg);

  s->sigma = gsl_spline_alloc(gsl_interp_cspline, filament->nseg);
  gsl_spline_init(s->sigma, filament->phi, filament->sigma, filament->nseg);

  s->kappa = gsl_spline_alloc(gsl_interp_cspline, filament->nseg);
  gsl_spline_init(s->kappa, filament->phi, filament->kappa, filament->nseg);

  s->tau = gsl_spline_alloc(gsl_interp_cspline, filament->nseg);
  gsl_spline_init(s->tau, filament->phi, filament->tau, filament->nseg);

  s->theta0 = gsl_spline_alloc(gsl_interp_cspline, filament->nseg);
  gsl_spline_init(s->theta0, filament->phi, filament->theta0, filament->nseg);

  s->a = gsl_spline_alloc(gsl_interp_cspline, filament->nseg);
  gsl_spline_init(s->a, filament->phi, filament->a, filament->nseg);

  free(temp);
  filament->splines = (void*)s;
}

// Function to allocate memory for the *members* of a vortex_filament struct
// Assumes the struct itself has already been created (e.g., on the stack or
// heap)
void allocate_vortex_filament_members(struct vortex_filament* filament, int nseg) {
  
  if (filament == NULL) {
    perror("Failed to allocate memory for filament");
    return ;
  }

  filament->nseg  = nseg;
  filament->pcar  = (coord) {0.0, 0.0, 0.0};
  filament->splines = NULL;

  // Allocate memory for the double arrays
  filament->phi = malloc(nseg * sizeof(double));  
  filament->s     = malloc(nseg * sizeof(double));
  filament->sigma = malloc(nseg * sizeof(double));
  filament->kappa = malloc(nseg * sizeof(double));
  filament->tau   = malloc(nseg * sizeof(double));  
  filament->a     = malloc(nseg * sizeof(double));
  filament->vol   = malloc(nseg * sizeof(double));
  filament->theta0 = malloc(nseg * sizeof(double));
  
  // Allocate memory for the coord arrays
  filament->C     = malloc(nseg * sizeof(coord));
  filament->dC    = malloc(nseg * sizeof(coord));
  filament->d2C   = malloc(nseg * sizeof(coord));
  filament->d3C   = malloc(nseg * sizeof(coord));  
  filament->Tvec  = malloc(nseg * sizeof(coord));
  filament->Nvec  = malloc(nseg * sizeof(coord));
  filament->Bvec  = malloc(nseg * sizeof(coord));

  filament->Ulocal  = malloc(nseg * sizeof(coord));
  filament->Uauto   = malloc(nseg * sizeof(coord));
  filament->Umutual = malloc(nseg * sizeof(coord));
  filament->U       = malloc(nseg * sizeof(coord));
  filament->Uprev   = malloc(nseg * sizeof(coord));
}

// Function to free memory for the *members* of a vortex_filament struct
// Assumes the struct itself will NOT be freed by this function
void free_vortex_filament_members(struct vortex_filament* filament) {
  if (filament == NULL) {
    // Nothing to free
    return;
  }

  // Free the double arrays
  free(filament->phi);      filament->phi = NULL; // Set to NULL after freeing
  free(filament->s);        filament->s = NULL;
  free(filament->sigma);    filament->sigma = NULL;
  free(filament->kappa);    filament->kappa = NULL;
  free(filament->tau);      filament->tau = NULL;
  free(filament->a);        filament->a = NULL;
  free(filament->vol);      filament->vol = NULL;
  free(filament->theta0);   filament->theta0 = NULL;

  // Free the coord arrays
  free(filament->C);        filament->C = NULL;
  free(filament->dC);       filament->dC = NULL;
  free(filament->d2C);      filament->d2C = NULL;
  free(filament->d3C);      filament->d3C = NULL;
  free(filament->Tvec);     filament->Tvec = NULL;
  free(filament->Nvec);     filament->Nvec = NULL;
  free(filament->Bvec);     filament->Bvec = NULL;

  free(filament->Ulocal);   filament->Ulocal = NULL;
  free(filament->Uauto);    filament->Uauto = NULL;
  free(filament->Umutual);  filament->Umutual = NULL;
  free(filament->U);        filament->U = NULL;
  free(filament->Uprev);    filament->Uprev = NULL;

  if (filament->splines != NULL) {
    free_filament_splines(filament->splines);
    filament->splines = NULL;
  }
}

struct local_filament{  
  bool near;       // flag based on distance from C
  double phi;      // arbitrary parametrisation of C  
  coord  C;        // cartesian coords of filament
  coord  Tvec;     // unit tangent vector
  coord  Nvec;     // unit normal vector
  coord  Bvec;     // unit binormal vector
  double s;        // arc-length coordinate
  coord  pcar;     // current point in Cartesian coordinates  
  double sigma;    // Stretching factor
  double kappa;    // Curvature
  double tau;      // Torsion
  double theta0;   // Cumulative torsion
  double a;        // Core size
  coord Mcar;      //
  coord Mrad;      //
  double rho;      //
  double theta;    //
};


/**
## Interpolate Values Using Cubic Splines
### *gsl_interp1d_vector()*: performs 1D interpolation on a vector using cubic splines

This function performs 1D interpolation on a vector using cubic splines from the
GNU Scientific Library (GSL). It takes discretized coord values `V0` at given
spatial coordinate `phi0` and interpolates to find the corresponding value at 
a coordinate `phiq`.

The arguments and their descriptions are:

*nseg*
: number of segments

*phi0*
: array with the coordinates $\phi_0$

*V0*
: array with the coord values $V_0$

*phiq*
: target value at which to interpolate $V_q$

*/

coord gsl_interp1d_vector( int nseg, double* phi0, coord * V0, double phiq){
  coord Vq;
  gsl_interp_accel *acc = gsl_interp_accel_alloc();

  foreach_dimension(){
    double *V_x;
    V_x = malloc(sizeof(double)*nseg);
    for (int i = 0; i < nseg; i++)
      V_x[i] = V0[i].x;

    gsl_spline *spline_x = gsl_spline_alloc(gsl_interp_cspline, nseg);
    gsl_spline_init(spline_x, phi0, V_x, nseg);
    Vq.x = gsl_spline_eval (spline_x, phiq, acc);
    gsl_spline_free (spline_x);
    free(V_x);
  }
  gsl_interp_accel_free (acc);
  return Vq;
}

/**
### *gsl_interp1d_scalar()*: performs 1D interpolation on a scalar using cubic splines

This function performs 1D interpolation on a scalar using cubic splines from the
GNU Scientific Library (GSL). It takes discretized scalar values `P0` at given
spatial coordinate `phi0` and interpolates to find the corresponding value at 
a coordinate `phiq`.

The arguments and their descriptions are:

*nseg*
: number of segments

*phi0*
: array with the coordinates $\phi_0$

*P0*
: array with the scalar values $P_0$

*phiq*
: target value at which to interpolate $P_q$
*/

double gsl_interp1d_scalar( int nseg, double* phi0, double * P0, double phiq){
  double Pq;
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *spline_x = gsl_spline_alloc(gsl_interp_cspline, nseg);
  gsl_spline_init(spline_x, phi0, P0, nseg);
  Pq = gsl_spline_eval (spline_x, phiq, acc);
  gsl_spline_free (spline_x);
  gsl_interp_accel_free (acc);
  return Pq;
}

/**
## Project a Point onto the Frenet Frame of a Vortex Filament
### *frenet_projection()*: projects a point onto the Frenet frame

This function calculates the projection of a point $M$ onto the Frenet-Serret
frame of a vortex filament at a specified position $s(\xi_q)$. The Frenet frame
$(\uvec{T}, \uvec{N}, \uvec{B})$ is interpolated at this position.

The function returns the dot product of the vector from the origin
$\vec{O}(s(\xi_q))$ to the point $\vec{M}$ with the tangent vector
$\uvec{T}(s(\xi_q))$:

$$
(\vec{M} - \vec{O}(s(\xi_q))) \cdot \uvec{T}(s(\xi_q)).
$$

If the point $\vec{M}$ lies within the plane $\mathcal{A}(s(\xi_q))$, the output
of this function will be zero. This property will be utilized as the objective
function to minimize when projecting points onto the Frenet-Serret frame.

The arguments and their descriptions are:

*tq*
: double representing the position along the vortex filament where the
projection is to be computed.

*params*
: pointer to a `struct vortex_filament`, which contains the properties of the
vortex filament, including the number of segments, coordinates, and Frenet frame
vectors.

*/

double frenet_projection (double phiq, void *params){
  struct vortex_filament *p = (struct vortex_filament *) params;

  double phi_min = p->phi[0];
  double phi_max = p->phi[p->nseg - 1];

  // Detect periodicity by checking if the endpoints of C are close
  bool is_periodic = (vecdist2(p->C[0], p->C[p->nseg - 1]) < 1e-6);
  if (is_periodic) {
    double len = phi_max - phi_min;
    if (len > 0.0) {
      double wrapped = phiq - phi_min;
      wrapped = fmod(wrapped, len);
      if (wrapped < 0.0) wrapped += len;
      phiq = phi_min + wrapped;
    }
  } else {
    if (phiq < phi_min) phiq = phi_min;
    if (phiq > phi_max) phiq = phi_max;
  }

  coord ccar, Tvec;
  if (p->splines != NULL) {
    struct filament_splines* s = (struct filament_splines*)p->splines;
    ccar.x = gsl_spline_eval(s->C[0], phiq, s->acc);
    ccar.y = gsl_spline_eval(s->C[1], phiq, s->acc);
    ccar.z = gsl_spline_eval(s->C[2], phiq, s->acc);

    Tvec.x = gsl_spline_eval(s->Tvec[0], phiq, s->acc);
    Tvec.y = gsl_spline_eval(s->Tvec[1], phiq, s->acc);
    Tvec.z = gsl_spline_eval(s->Tvec[2], phiq, s->acc);
  } else {
    ccar = gsl_interp1d_vector( p->nseg, p->phi, p->C, phiq);
    Tvec = gsl_interp1d_vector( p->nseg, p->phi, p->Tvec, phiq);
  }

  return vecdot(vecdiff(p->pcar, ccar), Tvec); 
}

/**
### *frenet_projection_min()*: find the position with the minimum projection on the Frenet frame

This function solves a minimization problem to find the position along a vortex
filament where the projection of a point `pcar` onto the tangent vector of the
Frenet-Serret frame is minimized. In other words, we search for $\xi_q$ such
that :
$$
(\vec{M} - \vec{O}(s(\xi_q))) \cdot \uvec{T}(s(\xi_q)) = 0
$$

It uses the Brent method from the GNU Scientific
Library (GSL) to find the root of the projection function within a specified
interval.

The arguments and their descriptions are:

*params*
: pointer to a `struct vortex_filament`, which contains the properties of the
vortex filament, including the number of segments, coordinates, and Frenet frame
vectors.

*r*
: an initial guess for $\xi_q$


*/

#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
double frenet_projection_min(struct vortex_filament params, double r) {

  // Use a dynamic search window proportional to the parameter spacing
  double dphi = (params.nseg > 1) ? (params.phi[1] - params.phi[0]) : 0.25;
  double window = 4.0 * dphi;
  if (window < 0.01) window = 0.01;
  if (window > 0.5) window = 0.5;

  double x_lo = r - window, x_hi = r + window;
  gsl_function F;

  F.function = &frenet_projection;
  F.params = &params;

  // Evaluate endpoints to verify if they bracket a root
  double f_lo = frenet_projection(x_lo, &params);
  double f_hi = frenet_projection(x_hi, &params);

  // If not bracketed, return r (the closest segment point) directly to avoid solver failure
  if (f_lo * f_hi >= 0.0) {
    return r;
  }

  int status, verbose = 0;
  int iter = 0, max_iter = 100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;

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


/**
### *get_local_coordinates()*: computes the coordinates in the local frame

This function computes the local coordinates required for the vorticity field.
Each point $M$ is projected into the local Frenet-Serret frame to
obtain a set of local coordinates, such that:

$$
\begin{aligned}
\textit{i)} \quad\quad & 
(\vec{M} - \vec{O}(s(\xi_q))) \cdot \uvec{T}(s(\xi_q)) = 0
\\
\textit{ii)} \quad\quad & 
(\vec{M} - \vec{O}(s(\xi_q))) \cdot \uvec{N}(s(\xi_q)) = x_N
\\
\textit{iii)} \quad\quad & 
(\vec{M} - \vec{O}(s(\xi_q))) \cdot \uvec{B}(s(\xi_q)) = x_B
\end{aligned}
$$

This requires finding the value of $s(\xi_q)$ along each space curve that
verifies *i)* through a minization process.


Then, we use the local coordinates $(x_N, x_B)$ to define a local radial and 
angular coordinates:
$$
\begin{aligned}
x_\rho = \sqrt{x_N^2 + x_B^2}, &&
x_\phi = \arctan(x_B,x_N)
\end{aligned}
$$
*/

struct local_filament get_local_coordinates(int spatial_period, double max_distance=L0, struct vortex_filament *vortex) {
  struct local_filament local_coordinates = {0}; 

  coord cart_coord, radial_coord; 
  double min_distance = 1e30, min_phi = 0;
  double local_radial_coord, local_angular_coord;
  
  // Iterate through each segment to find the segment closest to the point
  // this should give us a good starting point
  for (int i = 0; i < vortex->nseg; i++) {
    double current_distance = vecdist2(vortex->pcar, vortex->C[i]);
    if (current_distance < min_distance) {
      min_distance = current_distance;
      min_phi = vortex->phi[i];
    }
  }

  // Adjust the initial phi guess if the curve has periodicity
  if (spatial_period != 0)
    min_phi = fmod(min_phi + spatial_period * 2 * pi, spatial_period * 2 * pi);

  
  // If the point is close enough to the vortex, refine the initial guess
  if (min_distance < max_distance) {
    
    // Find the value of xi
    double phi = frenet_projection_min(*vortex, min_phi);

    coord Ocar, Tvec, Nvec, Bvec;
    double torsion_angle;
    double s, sigma, kappa, tau, a;
    
    if (vortex->splines != NULL) {
      struct filament_splines* s_cache = (struct filament_splines*)vortex->splines;
      Ocar.x = gsl_spline_eval(s_cache->C[0], phi, s_cache->acc);
      Ocar.y = gsl_spline_eval(s_cache->C[1], phi, s_cache->acc);
      Ocar.z = gsl_spline_eval(s_cache->C[2], phi, s_cache->acc);

      Tvec.x = gsl_spline_eval(s_cache->Tvec[0], phi, s_cache->acc);
      Tvec.y = gsl_spline_eval(s_cache->Tvec[1], phi, s_cache->acc);
      Tvec.z = gsl_spline_eval(s_cache->Tvec[2], phi, s_cache->acc);

      Nvec.x = gsl_spline_eval(s_cache->Nvec[0], phi, s_cache->acc);
      Nvec.y = gsl_spline_eval(s_cache->Nvec[1], phi, s_cache->acc);
      Nvec.z = gsl_spline_eval(s_cache->Nvec[2], phi, s_cache->acc);

      Bvec.x = gsl_spline_eval(s_cache->Bvec[0], phi, s_cache->acc);
      Bvec.y = gsl_spline_eval(s_cache->Bvec[1], phi, s_cache->acc);
      Bvec.z = gsl_spline_eval(s_cache->Bvec[2], phi, s_cache->acc);

      torsion_angle = gsl_spline_eval(s_cache->theta0, phi, s_cache->acc);

      s     = gsl_spline_eval(s_cache->s, phi, s_cache->acc);
      sigma = gsl_spline_eval(s_cache->sigma, phi, s_cache->acc);
      kappa = gsl_spline_eval(s_cache->kappa, phi, s_cache->acc);
      tau   = gsl_spline_eval(s_cache->tau, phi, s_cache->acc);
      a     = gsl_spline_eval(s_cache->a, phi, s_cache->acc);
    } else {
      // Find the Cartesian coordinates of point O(xi)
      Ocar = gsl_interp1d_vector(vortex->nseg, vortex->phi, vortex->C, phi);

      // Compute the Frenet-Serret frame vectors at the projected phi
      Tvec = gsl_interp1d_vector(vortex->nseg, vortex->phi, vortex->Tvec, phi);
      Nvec = gsl_interp1d_vector(vortex->nseg, vortex->phi, vortex->Nvec, phi);
      Bvec = gsl_interp1d_vector(vortex->nseg, vortex->phi, vortex->Bvec, phi);

      // Compute the torsion angle
      torsion_angle = gsl_interp1d_scalar(vortex->nseg, vortex->phi, vortex->theta0, phi);

      // Then, we compute the other properties
      s       = gsl_interp1d_scalar( vortex->nseg, vortex->phi, vortex->s,       phi);
      sigma   = gsl_interp1d_scalar( vortex->nseg, vortex->phi, vortex->sigma,   phi);
      kappa   = gsl_interp1d_scalar( vortex->nseg, vortex->phi, vortex->kappa,   phi);
      tau     = gsl_interp1d_scalar( vortex->nseg, vortex->phi, vortex->tau,     phi);    
      a       = gsl_interp1d_scalar( vortex->nseg, vortex->phi, vortex->a,       phi);    
    }
    
    // Compute local coordinates in the Frenet-Serret frame
    cart_coord.x = vecdot(vecdiff(vortex->pcar, Ocar), Tvec); // x_T (must be zero) 
    cart_coord.y = vecdot(vecdiff(vortex->pcar, Ocar), Nvec); // x_N 
    cart_coord.z = vecdot(vecdiff(vortex->pcar, Ocar), Bvec); // x_B

    // Convert local coordinates to radial coordinates
    local_radial_coord = sqrt(vecdot(cart_coord, cart_coord));
    local_angular_coord = atan2(cart_coord.z, cart_coord.y);
    
    // Set radial coordinates
    radial_coord.x = cart_coord.x;
    radial_coord.y = local_radial_coord;
    radial_coord.z = local_angular_coord - torsion_angle;

    local_coordinates = (struct local_filament){1, phi, Ocar, Tvec, Nvec, Bvec, 
      s, vortex->pcar, sigma, kappa, tau, torsion_angle, a, cart_coord, radial_coord, 
      local_radial_coord, local_angular_coord - torsion_angle}; 
  }
  return local_coordinates;
}

double get_min_distance(int spatial_period, double max_distance=L0, struct vortex_filament *vortex) {
    
  double min_distance = 1e30;
  
  // Iterate through each segment to find the segment closest to the point
  // this should give us a good starting point
  for (int i = 0; i < vortex->nseg; i++) {
    double current_distance = vecdist2(vortex->pcar, vortex->C[i]);
    if (current_distance < min_distance) {
      min_distance = current_distance;      
    }
  }
  return min_distance;
}

/**
## Compute the Finite Difference Derivative of a Coordinate Array
### *fd_derivative()*: calculates the first derivative of a coordinate array using finite differences

This function computes the first derivative of a coordinate array `X` using a
central finite difference scheme. It handles the interior points using a
standard central difference formula and applies a periodic boundary condition
with an optional shift to the endpoints.

The arguments and their descriptions are:

*n*
: integer representing the number of points in the coordinate array.

*dphi*
: double representing the spacing between consecutive points in the array, used
as the step size in the finite difference calculation.

*shift*
: a `coord` structure representing an optional shift applied to the periodic
boundary condition.

*X*
: pointer to a `coord` array representing the input coordinate values.

*dX*
: pointer to a `coord` array where the computed derivatives will be stored.
*/

void fd_derivative( int n, double dphi, coord shift, coord *X, coord *dX){
  for (int i = 1; i < n-1; i++){
    foreach_dimension(){
      dX[i].x = (X[i+1].x - X[i-1].x)/(2*dphi);
    }
  }
  foreach_dimension(){
    dX[0].x   = (X[1].x - X[n-2].x + shift.x)/(2*dphi);
    dX[n-1].x = (X[1].x - X[n-2].x + shift.x)/(2*dphi);
  }
}

/** 
 Finally, we create a macro so we can initialize the vortex filaments more
 easily
*/
macro initialize_filaments (struct vortex_filament filament, int nseg, double dphi, double* phi, double* a, coord* C, coord xshift, coord dxshift)
{
  
  for (int i = 0; i < nseg; i++){   
    filament.a[i] = a[i]; 
    filament.phi[i] = phi[i]; 
    foreach_dimension(){            
      filament.C[i].x =  C[i].x;  
    }
  }  

  // Find the 1st, 2nd, and 3rd derivatives of C  
  fd_derivative(nseg, dphi,  xshift,   filament.C,  filament.dC);
  fd_derivative(nseg, dphi, dxshift,  filament.dC, filament.d2C);
  fd_derivative(nseg, dphi, dxshift, filament.d2C, filament.d3C);

  // Compute the Frenet-Serret frame, curvature, and torsion
  for (int i = 0; i < nseg; i++){   
    foreach_dimension(){                      
      filament.Tvec[i].x =  filament.dC[i].x/sqrt(vecdot( filament.dC[i],  filament.dC[i]));      
      filament.Nvec[i].x = filament.d2C[i].x/sqrt(vecdot(filament.d2C[i], filament.d2C[i]));      
    }
    filament.sigma[i] = sqrt(vecdot( filament.dC[i],  filament.dC[i]));
    filament.Bvec[i] = vecdotproduct(filament.Tvec[i], filament.Nvec[i]);   
      
    filament.kappa[i] = sqrt(vecdot(filament.d2C[i], filament.d2C[i]))/sq(filament.sigma[i]);
    coord var1 = vecdotproduct(filament.dC[i], filament.d2C[i]);    
    filament.tau[i] = vecdot(var1, filament.d3C[i])/vecdot(var1,var1);        
  }

  // Compute the arc-lenght coordinate and cumulative torsion  
  memset (filament.s,      0, nseg*sizeof (double));
  memset (filament.theta0, 0, nseg*sizeof (double));  
  for (int i = 0; i < nseg-1; i++){
    filament.s[i+1] = filament.s[i] + filament.sigma[i+1]*dphi;
    filament.theta0[i+1] = filament.theta0[i] - filament.sigma[i+1]*filament.tau[i+1]*dphi;
  }
}

/**
## Compute the curvature terms at next to leading order for a Batchelor vortex

At the next to leading order, corrections due to curvature are obtained by
solving a second- order differential equation
$$
\mathbf{L} \psi = R 
$$



In the functions below, we solve numerically $\psi$ and compute the corrections
at this order.

*/

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

/** Struct to pass parameters ($a$, $U_c$) to the GSL integration functions
*/ 
typedef struct {
  double a;
  double U_c;
} integration_params;

/** Struct to hold the results: $\psi(\rho)$ and its first two derivatives
*/ 
typedef struct {
  double A;    // Psi(rho)
  double A_p;  // First derivative Psi'(rho)
  double A_pp; // Second derivative Psi''(rho)
} integration_results;

/** 

$$
f(\rho,a) = \exp(-\rho^2/a^2)
$$
 
*/ 
double gauss(double s, const integration_params *params) {
  return exp(-s * s / (params->a * params->a));
}

/**
$$
R = u^{(0)}_\theta + 2\rho\omega^{(0)}_\xi - 2\rho\omega^{(0)}_\theta u^{(0)}_\xi/u^{(0)}_\theta
$$
Here, $\vec{u}^{(0)}$ and $\vec{\omega}^{(0)}$ are the leading order solutions.
Subscripts $(\rho,\theta,\xi)$ indicate the local **radial**, **azimuthal**, and
**axial** components, respectively. 
*/

double right_hand_side(double s, void *p) {
  integration_params *params = (integration_params *)p;
  const double a = params->a;
  const double U_c = params->U_c;
  const double tol = 1e-18;

  if (fabs(s) < tol) return 0.0;

  double gs = gauss(s, params);
  double u_th = (1.0 / s) * (1.0 - gs);
  double u_xi = U_c * gs;
  double w_th = U_c * (2.0 * s / sq(a)) * gs;
  double w_xi = (2.0 / sq(a)) * gs;
  double R = 2.0 * s * w_xi + u_th - 2.0 * s * w_th * u_xi / (u_th + tol);
  return R;
}

/**
The solution of $\psi$ is obtained using the method of variation of constants
$$
\psi = y_1 \int_0^\rho G(s)/(s y_1^2(s)) ds
$$
where the integrand $G(s)$ reads
$$
G(\rho) = \int_0^\rho s y_1(s) R(s) ds
$$
and $y_1$ is the fundamental solution 
$$
y_1 = u^{(0)}_\theta = \frac{1}{\rho} (1 - f(\rho,a))
$$
*/

double integrand_G(double s, void *p) {
  integration_params *params = (integration_params *)p;
  return (1.0 - gauss(s, params)) * right_hand_side(s, p);
}

double integrand_Psi(double s, void *p) {
  integration_params *params = (integration_params *)p;

  if (fabs(s) < 1e-9) return 0.0;

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
  double G_s, error;
  
  gsl_function F_G;
  F_G.function = &integrand_G;
  F_G.params = p;

  gsl_integration_qag(&F_G, 0, s, 0, 1e-7, 1000, GSL_INTEG_GAUSS61, w, &G_s, &error);
  gsl_integration_workspace_free(w);

  double gs = gauss(s, params);
  double denominator = sq(1.0 - gs);
  
  if (fabs(denominator) < 1e-18) return 0.0;
  
  return (s * G_s) / denominator;
}

/**

We perform the numerical integration using GLS and evaluate the first and 
second derivatives with respect to $\rho$

*/

void compute_A_with_derivatives(double rho, double a, double U_c, integration_results *results) {
  const double tol = 1e-18;
  if (rho < tol){
    results->A = 0.;
    results->A_p = 0.;
    results->A_pp = 0.; 
    return;
  }

  integration_params params = {a, U_c};

  // --- Step 1: Compute the necessary integrals ---
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
  double G_rho, I_psi_rho, error;
  
  gsl_function F_G, F_psi;
  F_G.function = &integrand_G;
  F_G.params = &params;
  F_psi.function = &integrand_Psi;
  F_psi.params = &params;
  
  // Temporarily turn off default GSL error handler
  gsl_error_handler_t *old_handler = gsl_set_error_handler_off();

  // Integral G(rho)
  gsl_integration_qag(&F_G, 0, rho, 0, 1e-7, 1000, GSL_INTEG_GAUSS61, w, &G_rho, &error);
  
  // Integral for A(rho), which we call I_psi(rho)
  gsl_integration_qag(&F_psi, 0, rho, 0, 1e-7, 1000, GSL_INTEG_GAUSS61, w, &I_psi_rho, &error);
  
  gsl_set_error_handler(old_handler); // Restore handler
  gsl_integration_workspace_free(w);

  // --- Step 2: Compute intermediate values at rho ---
  double g_rho = gauss(rho, &params);
  double R_rho = right_hand_side(rho, &params);

  // --- Step 3: Calculate A, A', and A'' using analytical formulas ---
  results->A = (1.0 / rho) * (1.0 - g_rho) * I_psi_rho;

  // First derivative A'(rho)
  double C1 = -(1.0 - g_rho) / sq(rho) + (2.0 * g_rho) / sq(a);
  results->A_p = C1 * I_psi_rho + G_rho / (1.0 - g_rho);

  // Second derivative A''(rho)
  double one_minus_g_sq = sq(1.0 - g_rho);
  double Ipsip_rho = (rho * G_rho) / one_minus_g_sq; // This is I_psi'(rho)
  
  double C1p = 2.0 * (1.0 - g_rho) / cube(rho) 
         - (2.0 * g_rho) / (sq(a) * rho) 
         - (4.0 * rho * g_rho) / (sq(a) * sq(a));
         
  double C2p = R_rho - (2.0 * rho * g_rho * G_rho) / (sq(a) * one_minus_g_sq);
  
  results->A_pp = C1p * I_psi_rho + C1 * Ipsip_rho + C2p;
}

/**
Having solved $\psi$ we may compute the velocity and vorticity terms
*/

// Struct to hold the final calculated output values.
typedef struct {
  double u0_r;  // Velocity at leading order
  double u0_th; 
  double u0_xi;
  double w0_r;  // Vorticity at leading order
  double w0_th;
  double w0_xi;
  double u1_r;  // Velocity at the following order
  double u1_th;
  double u1_xi;
  double w1_r;  // Vorticity at the following order
  double w1_th;
  double w1_xi;
} batchelor_vortex;

/**
 
At **leading-order**, the velocity reads 
$$
u^{(0)}_\rho = 0,
\quad
u^{(0)}_\theta = \frac{1}{\rho}(1 - f(\rho,a)),
\quad 
u^{(0)}_\xi = U_c f(\rho,a)
$$
while the vorticity reads 
$$
\omega^{(0)}_\rho = 0,
\quad
\omega^{(0)}_\theta = (2 U_c \rho / a^2) f(\rho,a),
\quad 
\omega^{(0)}_\xi = (2/a^2) f(\rho,a)
$$

while the **next-to-leading-order** terms are functions of $\psi$, $\kappa$, 
$\varphi$ and the leading order terms.

*/

batchelor_vortex calculate_vortex_flow(struct local_filament* vortex, double U_c, 
  const integration_results* corr) {
  
  batchelor_vortex results;

  // Extracting values for readability
  double rho = vortex->rho; 
  double a = vortex->a; 
  double kappa = vortex->kappa; 
  double phi = vortex->theta + vortex->theta0;

  double A = corr->A;
  double dA = corr->A_p;
  double d2A = corr->A_pp;
  

  // Intermediate calculations
  double g_rho = exp(-sq(rho / a));
  double rho2 = sq(rho);
  double a2 = sq(a);

  // --- Zeroth-Order Base Flow (u0, w0) ---
  results.u0_r = 0.0;
  results.u0_th = (1.0 / rho) * (1.0 - g_rho);
  results.u0_xi = U_c * g_rho;

  results.w0_r = 0.0;
  results.w0_th = U_c * (2.0 * rho / a2) * g_rho;
  results.w0_xi = (2.0 / a2) * g_rho;

  // --- Intermediate terms needed for first-order terms ---
  // --- Swirl number at leading order
  double S0 = rho * results.w0_th / results.u0_th; 
  double denom = 1.0 - g_rho;
  double dS0_dr = (3.0 * rho2 - 2.0 * sq(rho2) / a2) * g_rho / denom
                - (2.0 * sq(rho2) / a2) * sq(g_rho) / (denom * denom);
  dS0_dr *= 2.0 * U_c / a2;

  // --- First-Order Perturbations (u1, w1) ---
  double cos_phi = cos(phi);
  double sin_phi = sin(phi);

  results.u1_r  = -A/rho2;
  results.u1_th = results.u0_th - dA/rho;
  results.u1_xi = results.u0_xi + S0*A/rho2;

  results.w1_r  = -S0*A/cube(rho);
  results.w1_th = results.w0_th + S0*A/cube(rho) - dS0_dr*A/rho2 - S0*dA/rho2 ;
  //double dw0_xi_dr = -2.0 * rho * results.w0_xi / a2;
  //results.w1_xi = -results.w0_xi - dw0_xi_dr * A / ((results.u0_th + 1e-18) * rho) + 2.0 * S0 * results.u0_xi / rho;
  results.w1_xi = results.u0_th/rho + results.w0_xi - d2A/rho - dA/rho2 + A/cube(rho);

  
  results.u1_r  *= kappa * rho * sin_phi;
  results.u1_th *= kappa * rho * cos_phi; 
  results.u1_xi *= kappa * rho * cos_phi;
  results.w1_r  *= kappa * rho * sin_phi; 
  results.w1_th *= kappa * rho * cos_phi;
  results.w1_xi *= kappa * rho * cos_phi;

  return results;
}

/**
# References

~~~bib

@article{callegari1978,
  title={Motion of a curved vortex filament with decaying vortical core and axial velocity},
  author={Callegari, AJ and Ting, Lu},
  journal={SIAM Journal on Applied Mathematics},
  volume={35},
  number={1},
  pages={148--175},
  year={1978},
  publisher={SIAM}
}

@article{blanco2015internal,
  title={Internal structure of vortex rings and helical vortices},
  author={Blanco-Rodr{\'\i}guez, Francisco J and Le Diz{\`e}s, St{\'e}phane and Sel{\c{c}}uk, Can and Delbende, Ivan and Rossi, Maurice},
  journal={Journal of Fluid Mechanics},
  volume={785},
  pages={219--247},
  year={2015},
  publisher={Cambridge University Press}
}

~~~
*/