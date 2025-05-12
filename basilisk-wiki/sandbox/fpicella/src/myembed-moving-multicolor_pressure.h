/**
# Advection of an arbitrary number of discrete rigid bodies

We detail here the first-order in time coupling algorithm we use to
advect an arbitrary number of discrete rigid bodies $\Gamma_{i,
\Delta}$, of boundary $\delta \Gamma_{i, \Delta}$, compatible only
with the [centered Navier-Stokes solver](mycentered.h). We consider
two advection scenarios:

* imposed embedded boundary motion;

* explicit weak fluid-solid coupling.

## Setup

In each scenario, each discrete rigid body $\Gamma_{i,\Delta}$ is
characterized by: the position of the center of mass
$\mathbf{x}_{\Gamma,i}$ stored in *c*, the translation velocity of the
center of mass $\mathbf{u}_{\Gamma,i}$ stored in *u*, the angular
velocity about the center of mass $\mathbf{\omega}_{\Gamma,i}$ stored
in *w* and their corresponding accelerations stored in *au, aw*.

The remaining unknown is the location and shape of the discrete rigid
boundary $\delta \Gamma_{i, \Delta}$, which we describe using a
user-defined distance function $\Phi_{i}$ given by the function
*phi()* along with a corresponding color field *col*.

Note that when using the explicit weak fluid-solid coupling, the
physical parameters describing the particle (density, volume, moment
of inertia) can be #defined via the hook *ADD_PARTICLE_PHY_PARAMS*. */

typedef struct {
  coord c;      // Center of mass
  coord u, w;   // Velocity
  coord au, aw; // Acceleration

#if dimension == 2
  double (* phi) (double x, double y);
#else // dimension == 3
  double (* phi) (double x, double y, double z);
#endif // dimension  

  scalar col;

#ifdef ADD_PARTICLE_PHY_PARAMS
  ADD_PARTICLE_PHY_PARAMS
#endif // ADD_PARTICLE_PHY_PARAMS
} particle;

/**
We define the type *particles* to avoid double pointer types and
mistakes. */

typedef particle * particles;

/**
By default we define only 1 particle. */

#ifndef p_n // number of particles
#define p_n (1)
#endif

particle pl[p_n] = {0};

/**
The following function *p_phi()* provides a default distance function
describing a circle or a sphere that the user can overload. */

#if dimension == 2
double p_phi (double x, double y)
{
  return sq (x) + sq (y) - sq (0.5);
}
#else // dimension == 3
double p_phi (double x, double y, double z)
{
  return sq (x) + sq (y) + sq (z) - sq (0.5);
}
#endif // dimension

/**
Then, by computing the intersection of all distance functions
$\Phi_{i}$ in the function *p_shape()*, we obtain the embedded
fractions representing all the embedded boundaries. Note that this
function must contain a call to the *fractions_cleanup()* function. */

#if PARTICLE_PERIOD_X
#define ADD_PARTICLE_PERIODICITY_X for (double xp = -(L0); xp <= (L0); xp += (L0))
#else
#define ADD_PARTICLE_PERIODICITY_X double xp = 0.;
#endif // PARTICLE_PERIOD_X

#if PARTICLE_PERIOD_Y
#define ADD_PARTICLE_PERIODICITY_Y for (double yp = -(L0); yp <= (L0); yp += (L0))
#else
#define ADD_PARTICLE_PERIODICITY_Y double yp = 0.;
#endif // PARTICLE_PERIOD_Y

#if dimension == 3 && PARTICLE_PERIOD_Z
#define ADD_PARTICLE_PERIODICITY_Z for (double zp = -(L0); zp <= (L0); zp += (L0))
#else
#define ADD_PARTICLE_PERIODICITY_Z double zp = 0.;
#endif // dimension == 3 && PARTICLE_PERIOD_Z

void p_shape (scalar c, face vector f, particles pl)
{
  vertex scalar phi[];
  foreach_vertex() {
    phi[] = HUGE;
    for (int i = 0; i < (p_n); i++) {
      ADD_PARTICLE_PERIODICITY_X
        ADD_PARTICLE_PERIODICITY_Y  
#if dimension == 2
	  phi[] = intersection (phi[], pl[i].phi ((x + xp - pl[i].c.x),
						  (y + yp - pl[i].c.y)));
#else // dimension == 3
          ADD_PARTICLE_PERIODICITY_Z
	    phi[] = intersection (phi[], pl[i].phi ((x + xp - pl[i].c.x),
						    (y + yp - pl[i].c.y),
						    (z + zp - pl[i].c.z)));
#endif // dimension
    }
  }
  boundary ({phi});
  fractions (phi, c, f);
  fractions_cleanup (c, f,
		     smin = 1.e-14, cmin = 1.e-14, opposite = true);
}

/**
Finally, the following function defines the color field *col*
associated with the ith-particle. */

void p_shape_col (particles pl)
{
  for (int i = 0; i < p_n; i++) {
    
    scalar col = automatic (pl[i].col);
    
    vertex scalar phi[];
    foreach_vertex() {
      phi[] = HUGE;
      ADD_PARTICLE_PERIODICITY_X
	ADD_PARTICLE_PERIODICITY_Y
#if dimension == 2
	  phi[] = intersection (phi[], pl[i].phi ((x + xp - pl[i].c.x),
						  (y + yp - pl[i].c.y)));
#else // dimension == 3
          ADD_PARTICLE_PERIODICITY_Z
	    phi[] = intersection (phi[], pl[i].phi ((x + xp - pl[i].c.x),
						    (y + yp - pl[i].c.y),
						    (z + zp - pl[i].c.z)));
#endif // dimension
    }
    boundary ({phi});
    fractions (phi, col);
    restriction ({col});
  }
}

/**
## Fluid-solid coupling

The following no-slip Dirichlet boundary condition for velocity and
hommogeneous Neumann boundary condition for pressure on the discrete
rigid boundary $\delta \Gamma_{i,\Delta}$ allow us to couple the
motion of the fluid and the discrete rigid body $\Gamma_{i,\Delta}$:
$$
\left\{
\begin{aligned}
&
\mathbf{u} = \mathbf{u}_{\Gamma,i} = \mathbf{u_{\Gamma,i}} +
\omega_{\Gamma,i} \times \left(\mathbf{x} - \mathbf{x}_{\Gamma,i}\right)
\\
&
{\nabla}_{\Gamma,i} p = 0.
\end{aligned}
\right.
$$
The homogeneous Neumann boundary condition is suitable for a fixed
rigid body and we have found that using it with a moving rigid body
does not significantly affect the computed solution.

#### No-slip Dirichlet boundary condition for velocity

The function *velocity_noslip_x()* computes the previously defined
no-slip boundary condition for the $x$-component of the velocity
$\mathbf{u}$. */

foreach_dimension()
static inline double velocity_noslip_x (Point point, particles pl,
					double xc, double yc, double zc)
{
  assert (cs[] > 0. && cs[] < 1.);

  for (int i = 0; i < (p_n); i++) {
    
    scalar col = automatic (pl[i].col);

    if (col[] > 0. && col[] < 1.) {

      /**
      We first compute the relative position $\mathbf{r_{i}} =
      \mathbf{x} - \mathbf{x}_{\Gamma,i}$. */ 
  
      // The coordinate x,y,z are not permuted with foreach_dimension()
      coord r = {xc, yc, zc};
      foreach_dimension() {
	r.x -= pl[i].c.x;
	if (Period.x) {
	  if (fabs (r.x) > fabs (r.x + (L0)))
	    r.x += (L0);
	  if (fabs (r.x) > fabs (r.x - (L0)))
	    r.x -= (L0);
	}
      }
      
      /**
      We then compute the veolcity (translation + rotation). */
    
#if dimension == 2
      coord sgn = {-1, 1};
      return (pl[i].u.x) + sgn.x*pl[i].w.x*(r.y);
#else // dimension == 3
      return (pl[i].u.x) + pl[i].w.y*(r.z) - pl[i].w.z*(r.y);
#endif // dimension
    }
  }
  // Other fixed embedded boundaries (u=0)
  return 0.;
}

/**
## Dump and restore a particle */

struct p_Dump {
  char * file;    // File name
  particles list; // List of particles
  FILE * fp;      // File pointer
  bool unbuffered;
};

void p_dump (struct p_Dump p)
{
  FILE * fp = p.fp;
  char def[] = "p_dump", * file = p.file ? p.file : p.fp ? NULL : def;

  char * name = NULL;
  if (file) {
    name = (char *) malloc (strlen(file) + 2);
    strcpy (name, file);
    if (!p.unbuffered)
      strcat (name, "~");
    if ((fp = fopen (name, "w")) == NULL) {
      perror (name);
      exit (1);
    }
  }
  assert (fp);

  // Get particle data
  particles p_list = p.list;

  // Dump particle data
  fwrite (p_list, sizeof(*p_list), (p_n), fp);
  
  /* free (p_list); */
  if (file) {
    fclose (fp);
    if (!p.unbuffered)
      rename (name, file);
    free (name);
  }
}

bool p_restore (struct p_Dump p)
{
  FILE * fp = p.fp;
  char * file = p.file;
  if (file && (fp = fopen (file, "r")) == NULL)
    return 0;
  assert (fp);

  // Read particle data
  particles p_list = p.list;

  if (fread (p_list, sizeof(*p_list), (p_n), fp) < 1) {
    fprintf (ferr, "#p_restore(): error reading particle data\n");
    exit (1);
  }

  if (file)
    fclose (fp);
  
  return true;
}

/**
## Initialization */

event defaults (i = 0)
{
  for (int i = 0; i < (p_n); i++) {

    scalar col = new scalar;
    foreach()
      col[] = 1.;
    boundary ({col});
    
#if TREE
    col.restriction = restriction_average;
    col.refine = col.prolongation = fraction_refine;
#endif

#ifdef ADD_PARTICLE_PHY_PARAMS
    pl[i] = (particle) {.c = {0., 0., 0.},
                        .u = {0., 0., 0.},
                        .w = {0., 0., 0.},
			.au = {0., 0., 0.},
			.aw = {0., 0., 0.},
			.phi = p_phi,
			.col = col,
			.r = 0.,
			.v = 0.,
			.i = {0., 0., 0.}};
#else
    pl[i] = (particle) {.c = {0., 0., 0.},
                        .u = {0., 0., 0.},
                        .w = {0., 0., 0.},
			.au = {0., 0., 0.},
			.aw = {0., 0., 0.},
			.phi = p_phi,
			.col = col};
#endif // ADD_PARTICLE_PHY_PARAMS
  }
}

event init (i = 0)
{
  /**
  We decrease the value of the *CFL* (same value as the one used with
  the VOF algorithm). */
  
  CFL = 0.5;

  if (!restore (file = "restart")) { // No restart
  
    /**
    We initialize the embedded boundary in the test case file. We also
    initialize the velocity in the test case file. */
  }
  else { // Restart

    /**
    We first restore the particle properties. */
    
    bool p_restart = p_restore ("p_restart", pl);
    assert (p_restart == true);

    /**
    We then need to initialize the face fraction *fs* since it is not
    dumped. */
  
    p_shape (cs, fs, pl);
  }

  /**
  We initialize, even when restarting, the color fields. */
  
  p_shape_col (pl);

  /**
  We then initialize the volume fraction at the previous timestep
  *csm1*. This needs to be done even when restarting the simulation as
  *csm1* is not dumped. */

  trash ({csm1});
  foreach()
    csm1[] = cs[];
  boundary    ({csm1});
  restriction ({csm1}); // Since restriction/prolongation depend on csm1

  /**
  Finally, we define the boundary conditions for the velocity, the
  pressure gradient *g* and the presssure *p* on the embedded
  boundaries. */
  
  p[embed]  = neumann (0);
  pf[embed] = neumann (0);

#if dimension == 2
  u.n[embed]   = dirichlet (velocity_noslip_x (point, pl, x, y, z));
  u.t[embed]   = dirichlet (velocity_noslip_y (point, pl, x, y, z));
  uf.n[embed]  = dirichlet (velocity_noslip_x (point, pl, x, y, z));
  uf.t[embed]  = dirichlet (velocity_noslip_y (point, pl, x, y, z));

  g.n[embed] = dirichlet (0);
  g.t[embed] = dirichlet (0);
#else // dimension == 3
  u.n[embed]   = dirichlet (velocity_noslip_x (point, pl, x, y, z));
  u.t[embed]   = dirichlet (velocity_noslip_y (point, pl, x, y, z));
  u.r[embed]   = dirichlet (velocity_noslip_z (point, pl, x, y, z));
  uf.n[embed]  = dirichlet (velocity_noslip_x (point, pl, x, y, z));
  uf.t[embed]  = dirichlet (velocity_noslip_y (point, pl, x, y, z));
  uf.r[embed]  = dirichlet (velocity_noslip_z (point, pl, x, y, z));
  
  g.n[embed] = dirichlet (0);
  g.t[embed] = dirichlet (0);
  g.r[embed] = dirichlet (0);
#endif // dimension

  /**
  As *rho* is used in the Neumann boundary condition for pressure, we
  need *rho* on all levels of the grid. */

  restriction ({rho});

  boundary ({p, u, g, pf, uf});
}

/**
## Timestep

We modify the maximun timestep *dtmax* to account for the velocity of
the discrete rigid body $\Gamma_{\Delta}$. Note that this event occurs
before moving the embedded boundaries to their $t^{n+1}$ position. We
therefore use the position and boundary conditions at time $t^n$. */

event stability (i++)
{
  foreach(reduction(min:dtmax)) {
    if (cs[] > 0. && cs[] < 1.) {

      // Barycenter and normal of the embedded fragment
      coord b, n;
      embed_geometry (point, &b, &n);

      // Local maximum velocity, in the direction of the normal
      double umax = 0.;
      foreach_dimension() {

	/**
	We use here the boundary condition on the embedded
	boundary. */
	
	bool dirichlet = true;
	double ub = (u.x.boundary[embed] (point, point,
					  u.x, &dirichlet));
	assert (dirichlet);	
	umax += (ub*n.x);
      }
      
      // Non-restrictive timestep (independent of *cs* and *fs*)
      double dte = Delta/(fabs (umax) + SEPS);
      if (dte < dtmax)
	dtmax = dte;
    }
  }
}

/**
## Prediction */

event advection_term (i++)
{
  /**
  In case of a periodic domain, we shift the coordinates of the center
  of mass *pi_p*. */

  coord p_o = {(X0), (Y0), (Z0)};
  for (int i = 0; i < (p_n); i++)
    foreach_dimension() {
      if (Period.x) {
	if (pl[i].c.x < p_o.x)
	  pl[i].c.x += (L0);
	if (pl[i].c.x > p_o.x + (L0))
	  pl[i].c.x -= (L0);
      }
    }
  
  /**
  #### Step 1

  We store the volume fraction defined at time *t*. */

  trash ({csm1});
  foreach()
    csm1[] = cs[];
  boundary    ({csm1});
  restriction ({csm1}); // Since restriction/prolongation depend on csm1

  /**
  #### Step 2 (prediction)

  We advance the embedded boundary to time *t+dt*. This step requires
  the user to define the quantities *p_p*, *p_u*, *p_w*, *p_au* and
  *p_aw* at time *t+dt*. */
  
  p_shape (cs, fs, pl);

  /**
  We then make sure not to use any values in newly emerged cells by
  setting the flag *emerged* to false and update all boundaries as the
  boundary conditions depend on the embedded boundaries. */

  emerged = false;
  boundary (all);

  /**
  We update the fluid properties to account for changes in the
  metric. */
  
  event ("properties");

  /**
  #### Step 3 (emerged cells)
  
  Since the advection event is explict, we define here the values of
  the centered velocity *u* and centered pressure gradient *g* in
  emerged cells. We also define the pressures *p* and *pf* in emerged
  cells to provide an improved initial guess to the multigrid
  projection solver.

  In the solid cells, we set all variables to 0. This is necessary to
  avoid mesh adaptation inside the solid boundaries. This might
  however lead to inaccuracies when using the default
  *restriction_average* operator. */

  for (scalar s in {p, u, g, pf})
    foreach()
      if (cs[] <= 0.)
	s[] = 0.;

  /**
  In the emerged cells, we use an extrapolation along the normal,
  using the function *embed_extrapolate*, to update the velocity *u*
  and the pressure gradient *g* (both use Dirichlet boundary
  conditions). */

  for (scalar s in {u, g})
    foreach() {
      if (csm1[] <= 0. && cs[] > 0.) {

  	// Normal emerged cell
  	if (cs[] < 1.) {

  	  // Cell centroid, barycenter and normal of the embedded fragment
	  coord c, b, n;
	  embed_geometry (point, &b, &n);
	  double alpha = plane_alpha (cs[], n);
	  plane_center (n, alpha, cs[], &c); // Different from line_area_center
	
  	  // Dirichlet boundary condition on the embedded boundary
  	  bool dirichlet = true;
	  double sb = (s.boundary[embed] (point, point, s, &dirichlet));
	  assert (dirichlet);
	  
  	  // Emerged cell value
  	  s[] = embed_extrapolate (point, s, cs, n, c, sb);
  	}
	
  	// Pathological emerged cell (cs = 1)
  	else {
	  int navg = 0;
	  double savg = 0.;
	  foreach_neighbor(1)
	    if (cs[] > 0. && (emerged || csm1[] > 0.)) {
	      navg += 1;
	      savg += s[];
	    }
	  s[] = savg/(navg + SEPS);
  	}
      }
    }

  /**
  As the pressure uses a Neumann boundary condition, we update the
  value of the pressures *p* and *pf* in emerged cells using their
  average over neighboring cells. Note that these values of pressure
  are used only as the initial condition for the Poisson solver. */
  
  for (scalar s in {p, pf})
    foreach() {
      if (csm1[] <= 0. && cs[] > 0.) {
	int navg = 0;
	double savg = 0.;
	foreach_neighbor(1)
	  if (cs[] > 0. && (emerged || csm1[] > 0.)) {
	    navg += 1;
	    savg += s[];
	  }
	s[] = savg/(navg + SEPS);
      }
    }
  
  /**
  Before using the *boundary* function, we set the *emerged* flag to
  true to indicate that all emerged cells have been updated and can
  now be used. */
  
  emerged = true;
  boundary (all);
  
  /**
  ## Step 4 (color function) */

  p_shape_col (pl);
  
  /**
  ## Step 5 (boundary conditions)

  As *rho* is used in the Neumann boundary condition for pressure, we
  need *rho* on all levels of the grid. */

  restriction ({rho});
}
