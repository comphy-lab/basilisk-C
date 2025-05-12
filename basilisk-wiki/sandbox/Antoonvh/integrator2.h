/**
## Solving for a line-integral field 

The line integral of $b$ from $\mathbf{x_o}$ in the (constant)
direction $\mathbf{n}$ is defined as $a$,

$$a(\mathbf{x_o}) = \int_{\mathbf{x_o}} b\ \mathrm{d}\mathbf{n},$$ 

using some intuitive vector calculus, $a$ statisfies

$$\frac{\partial a}{\partial \mathbf{n}} = -b,$$

Because solving stragies for both first and second-order accurate
discretizations of this problem have specific issues, we derive the
equation again,

$$\frac{\partial^2 a}{\partial \mathbf{n}^2} = -\frac{\partial b}{\partial \mathbf{n}},$$

or, the so-called "implicit integral equation"

$$\mathbf{n}\cdot\left(\nabla\left(\mathbf{n}\cdot \nabla a\right) \right) = -\mathbf{n} \cdot \nabla b.$$

i.e. in 2 dimensions,

$$n_x^2a_{xx} + n_y^2a_{yy} + 2n_xn_ya_{xy} = -n_xb_x - n_yb_y,$$

or in 3, 

$$n_x^2a_{xx} + n_y^2a_{yy} + n_z^2a_{zz} + 2n_xn_ya_{xy} + 2n_xn_za_{xz}+ 2n_yn_za_{yz} = -n_xb_x - n_yb_y - n_zb_z,$$

where the subscripts indicate components or derivatives where
sensible. Such an equation maybe solved iteratively, using multigrid
acceleration.
 */
#include "poisson.h"

struct Integrate_dn {
  scalar a, b;
  coord n;
  double tolerance;
  int nrelax, minlevel;
  scalar * res;
};

static double residual_int (scalar * al, scalar * rhsl, scalar * resl, void * data) {
  scalar a = al[0], rhs = rhsl[0], res = resl[0];
  struct Integrate_dn * p = (struct Integrate_dn *) data;
  coord n = p->n;
  double maxres = 0;
  boundary(al);
  foreach () {
    res[] = rhs[];
    foreach_dimension() {
      res[] -= sq(n.x)*(a[1] + a[-1] - 2.*a[])/(sq(Delta));  
    }
#if (dimension == 2)
    res[] -= 2.*n.x*n.y*((a[1, 1]    - a[1, -1])    - (a[-1, 1] - a[-1, -1]))/(4*sq(Delta));
#elif (dimension == 3)
    res[] -= 2.*(n.x*n.y*((a[1, 1, 0] - a[1, -1, 0]) - (a[-1, 1, 0] - a[-1, -1, 0])) +
		 n.x*n.z*((a[1, 0, 1] - a[1, 0, -1]) - (a[-1, 0, 1] - a[-1, 0, -1])) +
		 n.y*n.z*((a[0, 1, 1] - a[0, 1, -1]) - (a[0, -1, 1] - a[0, -1, -1]))
		 )/(4*sq(Delta));
#endif
    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);
  }
  boundary (resl);
  return maxres;
}

static void relax_int (scalar * al, scalar * rhsl, int l, void * data) {
  scalar a = al[0], rhs = rhsl[0];;
  struct Integrate_dn * p = (struct Integrate_dn *) data;
  coord n = p->n;
#if JACOBI
  scalar c[];
#else
  scalar c = a;
#endif
  foreach_level (l) {
    double nu = -rhs[]*sq(Delta), d = 0;
    foreach_dimension() {
      d += 2.*sq(n.x);
      nu += sq(n.x)*(a[1] + a[-1]);
    }
#if (dimension == 2)
    nu += 2.*n.x*n.y*((a[1, 1] - a[1, -1]) - (a[-1, 1] - a[-1, -1]))/4.;
#elif (dimension == 3)
    nu += 2.*(n.x*n.y*((a[1, 1, 0] - a[1, -1, 0]) - (a[-1, 1, 0] - a[-1, -1, 0])) +
	      n.x*n.z*((a[1, 0, 1] - a[1, 0, -1]) - (a[-1, 0, 1] - a[-1, 0, -1])) +
	      n.y*n.z*((a[0, 1, 1] - a[0, 1, -1]) - (a[0, -1, 1] - a[0, -1, -1]))
	      )/(4.);
#endif
    c[] = nu/d;
  }
#if JACOBI
  foreach_level (l) {
    a[] = (a[] + 2.*c[])/3.;
  }
#endif
}

/**
## Boundary conditions

The solution to the explicit equation is only relevant when subjected
to the proper boundary conditions. 
 */

#include "bwatch.h"
trace
double BC_int (coord O, coord n, scalar s) {
  double integral = 0; 
  ray r;
  r.O = O;
  r.dir = n;
  foreach_ray_cell_intersection(r reduction(+:integral)) {
    double l = 0;
    foreach_dimension()
      l += sq(_a[1].x - _a[0].x);
    l = l > 0 ? sqrt (l) : 0;
    integral += l*interpolate (s, (_a[0].x + _a[1].x)/2.,
			       (_a[0].y + _a[1].y)/2., (_a[0].z + _a[1].z)/2.); 
  }
  return integral;
}

/**
## User interface

`integrate_dn()` forms the user interface. For the moment, boundary conditions are not computed properly when using MPI.
*/

face vector _BC;

trace
mgstats integrate_dn (struct Integrate_dn p) {
  scalar a = p.a, b = p.b, temp[];
  _BC = new_face_vector ("_BC");
  foreach()
    temp[] = a[];
  boundary ({temp});
  normalize (&p.n);
  double tol = L0*1e-4;

  //We compute the boundary values
  temp[left]  = dirichlet (_BC.x[]);
  temp[right] = dirichlet (_BC.x[]);
  foreach_boundary (left)  {
    _BC.x[] = BC_int ((coord){x + tol, y, z}, p.n, b);
  }
  foreach_boundary (right) {
    _BC.x[] = BC_int ((coord){x - tol, y, z}, p.n, b);
  }
#if (dimension >= 2)
  temp[bottom] = dirichlet (_BC.y[]);
  temp[top]    = dirichlet (_BC.y[]);
  foreach_boundary (bottom) {
    _BC.y[] = BC_int ((coord){x, y + tol, z}, p.n, b);
  }
  foreach_boundary (top) {
    _BC.y[] = BC_int ((coord){x, y - tol, z}, p.n, b);
  }
#endif
#if (dimension == 3)
  temp[back]  = dirichlet (_BC.z[]);
  temp[front] = dirichlet (_BC.z[]);
  foreach_boundary (back) { 
    _BC.z[] = BC_int ((coord){x, y, z + tol}, p.n, b);
  }
  foreach_boundary (front) {
    _BC.z[] = BC_int ((coord){x, y, z - tol}, p.n, b);
  }
#endif
  boundary_flux ({_BC}); //Coarse level boundary values.
  scalar rhs[];
  foreach() {
    rhs[] = 0;
    foreach_dimension() 
      rhs[] += -p.n.x*(b[1] - b[-1])/(2.*Delta);
  }
  mgstats mgint = mg_solve ({temp}, {rhs}, residual_int, relax_int, &p,
			    p.nrelax, p.res, minlevel = max(1, p.minlevel));
  delete ((scalar*){_BC});
  foreach()
    a[] = temp[];
  boundary ({a});
  return mgint;
}
