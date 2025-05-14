/** 
\renewcommand*{\vec}[1]{\boldsymbol{#1}}
\newcommand*{\uvec}[1]{\boldsymbol{\hat{#1}}}
# Curved vortex filaments

This section borrows from the development for a general curved vortex presented
by [Callegari and Ting (1978)](#callegari1978) and notes from the book by
Maurice Rossi, Ivan Delbende, Stéphane Le Dizès, and others.

## Curvilinear coordinate system

A curvilinear coordinate system can be defined starting from the orthogonal
Cartesian one. If $x$, $y$, $z$ are the Cartesian coordinates, the curvilinear
ones, $u$, $v$, $w$, can be expressed as smooth functions of $x$, $y$, $z$,
according to $u = u(x,y,z)$, $v = v(x,y,z)$, and $w = w(x,y,z)$, which can be
inverted as $x = x(u,v,w)$, $y = y(u,v,w)$, $z = z(u,v,w)$. The infinitesimal
increment $d\vec{r}$ is the same regardless the coordinate system. For instance,
starting from Cartesian coordinates, the infinitesimal increment
$d\vec{r}=(dx,dy,dz)$ can be written as:
$$
\begin{aligned}
  d\vec{r} =\frac{\partial\vec{r}}{\partial u} du 
  + \frac{\partial\vec{r}}{\partial v} dv 
  + \frac{\partial\vec{r}}{\partial w} dw
\end{aligned}
$$
or alternatively as
$$
\begin{aligned}
  d\vec{r} =
  h_u du \uvec{e}_u  + h_v dv \uvec{e}_v + h_w dw \uvec{e}_w
\end{aligned}
$$
where $h_u$, $h_v$, and $h_w$, are the scale (also known as Lamé) coefficients
$$
\begin{aligned}
  h_u \equiv \sqrt{ \frac{\partial\vec{r}}{\partial u} \cdot \frac{\partial\vec{r}}{\partial u} }, \quad
  h_v \equiv \sqrt{ \frac{\partial\vec{r}}{\partial v} \cdot \frac{\partial\vec{r}}{\partial v} }, \quad
  h_w \equiv \sqrt{ \frac{\partial\vec{r}}{\partial w} \cdot \frac{\partial\vec{r}}{\partial w} }
\end{aligned}
$$
and $\vec{e}_u$, $\vec{e}_v$, and $\vec{e}_w$, are the curvilinear
orthonormal basis vectors
$$
\begin{aligned}
  \uvec{e}_u = \frac{1}{h_u} \frac{\partial\vec{r}}{\partial u}, \quad
  \uvec{e}_v = \frac{1}{h_v} \frac{\partial\vec{r}}{\partial v}, \quad
  \uvec{e}_w = \frac{1}{h_w} \frac{\partial\vec{r}}{\partial w}
\end{aligned}
$$

For instance, in this system, the elementary arc length reads
$$
\begin{aligned}  
  ds \equiv
  \sqrt{d\vec{r}\cdot d\vec{r}} =
  \sqrt{ (h_u du)^2 + (h_v dv)^2 + (h_w dw)^2 }
\end{aligned}
$$
while the divergence of a vector field $\vec{v}$ reads as
$$
\begin{aligned}
  \nabla\cdot\vec{v} =
  \frac{1}{h_u h_v h_w}
  \left[
  \frac{ \partial(v_u h_v h_w) }{\partial u}  +
  \frac{ \partial(v_v h_u h_w) }{\partial v}  +
  \frac{ \partial(v_w h_v h_u) }{\partial w}
  \right]
\end{aligned}
$$

## Frenet–Serret basis 

Now, let us define a curvilinear orthogonal basis associated with a reference
space curve $\mathcal{C}(s)$, which lies in the filament. Here $s$ is the arc
length along the initial curve measured from a given reference point. Each point
along $\mathcal{C}(s)$ is characterized by its position $\vec{x}_c$, local
curvature $\kappa(s)$, local torsion $\tau(s)$, and the corresponding
Frenet–Serret frame composed of three vectors $(\uvec{T}, \uvec{N}, \uvec{B})$
$$
\begin{aligned}
  \uvec{T}(s) \equiv \frac{d\vec{x}_c}{ds}, \quad
  \uvec{N}(s) \equiv \frac{1}{\kappa}\frac{d\uvec{T}}{ds}, \quad
  \uvec{B}(s) \equiv \uvec{T}\times\uvec{N}
\end{aligned}
$$
which are related to one another as:
$$
\begin{aligned}
\frac{d}{ds}
\begin{bmatrix}
\uvec{T} \\ \uvec{N} \\ \uvec{B}
\end{bmatrix}
=
\begin{bmatrix}
0 & \kappa(s) & 0 \\
-\kappa(s) & 0 & \tau(s) \\
0 & -\tau(s) & 0 
\end{bmatrix}
\begin{bmatrix}
\uvec{T} \\ \uvec{N} \\ \uvec{B}
\end{bmatrix}
\end{aligned}
$$

Here $\kappa(s, t)$ and $\tau(s, t)$ are called the curvature and torsion at
point $O(s, t)$
$$
\begin{aligned}
\kappa(s, t) &\equiv 
\frac{\lVert \vec{r}'(s,t) \times \vec{r}''(s,t) \rVert}{\lVert \vec{r}'(s,t) \rVert^3}
\\
\tau(s, t) &\equiv \frac{\lVert (\vec{r}'(s,t) \times \vec{r}''(s,t)) \cdot \vec{r}'''(s,t) \rVert}{\lVert \vec{r}'(s,t) \times \vec{r}''(s,t) \rVert^2}
\end{aligned}
$$

## Position vector in the Frenet-Serret frame 

<center>![Geometrical elements for a curved line (a) in the osculating plane (b) 
in the plane orthogonal to the vortex axis](rossi1.png)</center>

### Introducing the two local Cartesian coordinates

The position of a given point $M$ can be expressed as
$$
\begin{aligned}
  \vec{x} = \vec{x}_c (s) + \vec{a}
\end{aligned}
$$
where $\vec{a}$ lies in the plane $A_n(s)$ perpendicular to $\mathcal{C}(s)$
$$
\begin{aligned}
  \vec{a} = a_2 \uvec{N}(s) + a_3 \uvec{B}(s)
\end{aligned}
$$

### Introducing the local radial and angular coordinates
We introduce a local radial and angular coordinates $(\rho,\phi)$ in the plane
$A(s)$
$$
\begin{aligned}
  \uvec{e}_\rho &=& \cos\phi\uvec{N} + \sin\phi\uvec{B}
  \\
  \uvec{e}_\phi &=& -\sin\phi\uvec{N} + \cos\phi\uvec{B}
\end{aligned}
$$
such that
$$
\begin{aligned}
  \vec{x} = \vec{x}_c (s) + \rho \vec{e}_\rho(s,\phi)
\end{aligned}
$$
Differentiating with respect to $s$ gives
$$
\begin{aligned}
  \frac{d\vec{x}}{ds} &=& \frac{d\rho}{ds} \uvec{e}_\rho(s,\phi)
    + \rho\left(\frac{d\phi}{ds} + \tau(s)\right)\uvec{e}_\phi(s,\phi)
    + (1-\kappa(s)\rho\cos\phi)\uvec{T}(s)
\end{aligned}
$$
from which we deduce the following expression:
$$
\begin{aligned}
  d\vec{x} = d\rho \uvec{e}_\rho(s,\phi)
    + \rho\left( d\phi + \tau(s)ds\right)\uvec{e}_\phi(s,\phi) + (1-\kappa(s)\rho\cos\phi)\uvec{T}(s) ds
\end{aligned}
$$

### A Locally Orthogonal basis

Note that $(\rho,\phi,s)$ are curvilinear coordinates but they are not
orthogonal in general. However, if we define an angle $\varphi$
$$
\begin{aligned}
  \varphi = \phi - \varphi_0(s), \quad \text{ where } \quad
  \frac{d\varphi_0}{ds} = \tau(s)
\end{aligned}
$$
then $(\rho,\varphi,s)$ form a set of curvilinear orthogonal coordinates. In
this system, the infinitesimal increment reads as
$$
\begin{aligned}
  d\vec{x} &=&  h_\rho d\rho ~\uvec{e}_\rho(s,\varphi)
      + h_\varphi d\varphi ~\uvec{e}_\varphi(s,\varphi) + h_s ds ~\uvec{T}(s)
\end{aligned}
$$
where the Lamé coefficients read:
$$
\begin{aligned}
  h_\rho = 1, \quad
  h_\varphi = \rho, \quad
  h_s = (1-\kappa(s)\rho\cos(\varphi+\varphi_0(s)))
\end{aligned}
$$
and the new orthogonal basis reads:
$$
\begin{aligned}
  \begin{bmatrix}
    \uvec{e}_\rho(\varphi,s) \\
    \uvec{e}_\varphi(\varphi,s) \\
    \uvec{T}(s)
  \end{bmatrix}
  =
  \begin{bmatrix}
  \cos(\varphi+\varphi_0(s)) & \sin(\varphi+\varphi_0(s)) & 0 \\
  -\sin(\varphi+\varphi_0(s)) & \cos(\varphi+\varphi_0(s)) & 0 \\
  0 & 0 & 1
  \end{bmatrix}
  \begin{bmatrix}
    \uvec{N}(s) \\
    \uvec{B}(s) \\
    \uvec{T}(s)
  \end{bmatrix}
\end{aligned}
$$

For instance, in this system, the divergence of a vector field
$\vec{v}=(v_\rho\uvec{e}_\rho + v_\varphi\uvec{e}_\varphi + v_t\uvec{T})$ would
read as: 
\begin{aligned}
  \nabla\cdot\vec{v} &=&
  \frac{1}{\rho h_s(\varphi,s)}\left(
  \frac{d(\rho h_s(\varphi,s) v_\rho)}{d\rho} + \frac{d(h_s(\varphi,s) v_\varphi)}{d\varphi} + \frac{d(\rho v_t)}{ds}
  \right)
\end{aligned}

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

*varphi0*
: cumulative torsion, $\varphi_0(s(\xi_n))$  (optional)

*/

#include "PointTriangle.h"
struct vortex_filament{
  int     nseg;     // number of segments
  double* theta;    // arbitrary parametrisation of C  
  coord*  C;        // cartesian coords of filament
  coord*  Tvec;     // unit tangent vector
  coord*  Nvec;     // unit normal vector
  coord*  Bvec;     // unit binormal vector
  double* s;        // arc-length coordinate
  coord   pcar;     // current point in Cartesian coordinates  
  double* sigma;    // Stretching factor
  double* kappa;    // Curvature
  double* tau;      // Torsion
  double* varphi0;  // Cumulative torsion
  double* a;        // Core size
  coord*  dC;       // 
  coord*  d2C;      // 
  coord*  d3C;      // 
  coord*  Ulocal;   //
  coord*  Uauto;    // 
  coord*  Umutual;  //
  coord*  Utotal;   // 
  double* vol;      // Initial volume 
  
};

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

  // Allocate memory for the double arrays
  filament->theta = malloc(nseg * sizeof(double));  
  filament->s     = malloc(nseg * sizeof(double));
  filament->sigma = malloc(nseg * sizeof(double));
  filament->kappa = malloc(nseg * sizeof(double));
  filament->tau   = malloc(nseg * sizeof(double));  
  filament->a     = malloc(nseg * sizeof(double));
  filament->vol   = malloc(nseg * sizeof(double));
  filament->varphi0 = malloc(nseg * sizeof(double));
  
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
  filament->Utotal  = malloc(nseg * sizeof(coord));
}

// Function to free memory for the *members* of a vortex_filament struct
// Assumes the struct itself will NOT be freed by this function
void free_vortex_filament_members(struct vortex_filament* filament) {
  if (filament == NULL) {
    // Nothing to free
    return;
  }

  // Free the double arrays
  free(filament->theta);    filament->theta = NULL; // Set to NULL after freeing
  free(filament->s);        filament->s = NULL;
  free(filament->sigma);    filament->sigma = NULL;
  free(filament->kappa);    filament->kappa = NULL;
  free(filament->tau);      filament->tau = NULL;
  free(filament->a);        filament->a = NULL;
  free(filament->vol);      filament->vol = NULL;
  free(filament->varphi0);  filament->varphi0 = NULL;

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
  free(filament->Utotal);   filament->Utotal = NULL;
}

struct local_filament{  
  bool near;       // flag based on distance from C
  double theta;    // arbitrary parametrisation of C  
  coord  C;        // cartesian coords of filament
  coord  Tvec;     // unit tangent vector
  coord  Nvec;     // unit normal vector
  coord  Bvec;     // unit binormal vector
  double s;        // arc-length coordinate
  coord  pcar;     // current point in Cartesian coordinates  
  double sigma;    // Stretching factor
  double kappa;    // Curvature
  double tau;      // Torsion
  double varphi0;  // Cumulative torsion
  double a;        // Core size
  coord Mcar;      //
  coord Mrad;      //
  double rho;      //
  double phi;      //
};


/**
## Interpolate Values Using Cubic Splines
### *gsl_interp1d_vector()*: performs 1D interpolation on a vector using cubic splines

This function performs 1D interpolation on a vector using cubic splines from the
GNU Scientific Library (GSL). It takes discretized coord values `V0` at given
spatial coordinate `theta0` and interpolates to find the corresponding value at 
a coordinate `thetaq`.

The arguments and their descriptions are:

*nseg*
: number of segments

*theta0*
: array with the coordinates $\theta_0$

*V0*
: array with the coord values $V_0$

*thetaq*
: target value at which to interpolate $V_q$

*/

#pragma autolink -lgsl -lgslcblas
#include <gsl/gsl_spline.h>
coord gsl_interp1d_vector( int nseg, double* theta0, coord * V0, double thetaq){
  coord Vq;
  gsl_interp_accel *acc = gsl_interp_accel_alloc();

  foreach_dimension(){
    double *V_x;
    V_x = malloc(sizeof(double)*nseg);
    for (int i = 0; i < nseg; i++)
      V_x[i] = V0[i].x;

    gsl_spline *spline_x = gsl_spline_alloc(gsl_interp_cspline, nseg);
    gsl_spline_init(spline_x, theta0, V_x, nseg);
    Vq.x = gsl_spline_eval (spline_x, thetaq, acc);
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
spatial coordinate `theta0` and interpolates to find the corresponding value at 
a coordinate `thetaq`.

The arguments and their descriptions are:

*nseg*
: number of segments

*theta0*
: array with the coordinates $\theta_0$

*P0*
: array with the scalar values $P_0$

*thetaq*
: target value at which to interpolate $P_q$
*/

double gsl_interp1d_scalar( int nseg, double* theta0, double * P0, double thetaq){
  double Pq;
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *spline_x = gsl_spline_alloc(gsl_interp_cspline, nseg);
  gsl_spline_init(spline_x, theta0, P0, nseg);
  Pq = gsl_spline_eval (spline_x, thetaq, acc);
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

double frenet_projection (double thetaq, void *params){
  struct vortex_filament *p = (struct vortex_filament *) params;

  coord ccar, frenet[3];
  ccar = gsl_interp1d_vector( p->nseg, p->theta, p->C, thetaq);

  frenet[0] = gsl_interp1d_vector( p->nseg, p->theta, p->Tvec, thetaq);
  frenet[1] = gsl_interp1d_vector( p->nseg, p->theta, p->Nvec, thetaq);
  frenet[2] = gsl_interp1d_vector( p->nseg, p->theta, p->Bvec, thetaq);

  return vecdot(vecdiff(p->pcar, ccar), frenet[0]); 
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
  double min_distance = 1e30, min_theta = 0;
  double local_radial_coord, local_angular_coord;
  
  // Iterate through each segment to find the segment closest to the point
  // this should give us a good starting point
  for (int i = 0; i < vortex->nseg; i++) {
    double current_distance = vecdist2(vortex->pcar, vortex->C[i]);
    if (current_distance < min_distance) {
      min_distance = current_distance;
      min_theta = vortex->theta[i];
    }
  }

  // Adjust the initial theta guess if the curve has periodicity
  if (spatial_period != 0)
    min_theta = fmod(min_theta + spatial_period * 2 * pi, spatial_period * 2 * pi);

  
  // If the point is close enough to the vortex, refine the initial guess
  if (min_distance < max_distance) {
    
    // Find the value of xi
    double theta = frenet_projection_min(*vortex, min_theta);

    // Find the Cartesian coordinates of point O(xi)
    coord Ocar = gsl_interp1d_vector(vortex->nseg, vortex->theta, vortex->C, theta);

    // Compute the Frenet-Serret frame vectors at the projected theta
    coord Tvec, Nvec, Bvec;
    Tvec = gsl_interp1d_vector(vortex->nseg, vortex->theta, vortex->Tvec, theta);
    Nvec = gsl_interp1d_vector(vortex->nseg, vortex->theta, vortex->Nvec, theta);
    Bvec = gsl_interp1d_vector(vortex->nseg, vortex->theta, vortex->Bvec, theta);
    
    // Compute local coordinates in the Frenet-Serret frame
    cart_coord.x = vecdot(vecdiff(vortex->pcar, Ocar), Tvec); // x_T (must be zero) 
    cart_coord.y = vecdot(vecdiff(vortex->pcar, Ocar), Nvec); // x_N 
    cart_coord.z = vecdot(vecdiff(vortex->pcar, Ocar), Bvec); // x_B

    // Convert local coordinates to radial coordinates
    local_radial_coord = sqrt(vecdot(cart_coord, cart_coord));
    local_angular_coord = atan2(cart_coord.z, cart_coord.y);
    
    // Compute the torsion angle
    double torsion_angle = gsl_interp1d_scalar(vortex->nseg, vortex->theta, vortex->varphi0, theta);

    // Set radial coordinates
    radial_coord.x = cart_coord.x;
    radial_coord.y = local_radial_coord;
    radial_coord.z = local_angular_coord - torsion_angle;

    // Then, we compute the other properties
    double s       = gsl_interp1d_scalar( vortex->nseg, vortex->theta, vortex->s,       theta);
    double sigma   = gsl_interp1d_scalar( vortex->nseg, vortex->theta, vortex->sigma,   theta);
    double kappa   = gsl_interp1d_scalar( vortex->nseg, vortex->theta, vortex->kappa,   theta);
    double tau     = gsl_interp1d_scalar( vortex->nseg, vortex->theta, vortex->tau,     theta);    
    double a       = gsl_interp1d_scalar( vortex->nseg, vortex->theta, vortex->a,       theta);    

    local_coordinates = (struct local_filament){1, theta, Ocar, Tvec, Nvec, Bvec, 
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

*dtheta*
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

void fd_derivative( int n, double dtheta, coord shift, coord *X, coord *dX){
  for (int i = 1; i < n-1; i++){
    foreach_dimension(){
      dX[i].x = (X[i+1].x - X[i-1].x)/(2*dtheta);
    }
  }
  foreach_dimension(){
    dX[0].x   = (X[1].x - X[n-2].x + shift.x)/(2*dtheta);
    dX[n-1].x = (X[1].x - X[n-2].x + shift.x)/(2*dtheta);
  }
}

/** 
 Finally, we create a macro so we can initialize the vortex filaments more
 easily
*/
macro initialize_filaments (struct vortex_filament filament, int nseg, double dtheta, double* theta, double* a, coord* C, coord xshift, coord dxshift)
{
  
  for (int i = 0; i < nseg; i++){   
    filament.a[i] = a[i]; 
    filament.theta[i] = theta[i]; 
    foreach_dimension(){            
      filament.C[i].x =  C[i].x;  
    }
  }  

  // Find the 1st, 2nd, and 3rd derivatives of C  
  fd_derivative(nseg, dtheta,  xshift,   filament.C,  filament.dC);
  fd_derivative(nseg, dtheta, dxshift,  filament.dC, filament.d2C);
  fd_derivative(nseg, dtheta, dxshift, filament.d2C, filament.d3C);

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
  memset (filament.s,       0, nseg*sizeof (double));
  memset (filament.varphi0, 0, nseg*sizeof (double));  
  for (int i = 0; i < nseg-1; i++){
    filament.s[i+1] = filament.s[i] + filament.sigma[i+1]*dtheta;
    filament.varphi0[i+1] = filament.varphi0[i] + filament.sigma[i+1]*filament.tau[i+1]*dtheta;
  }
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


~~~
*/