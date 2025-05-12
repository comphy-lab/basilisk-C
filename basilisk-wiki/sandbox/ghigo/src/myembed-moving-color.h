/**
# Advection of a discrete rigid body

We detail here the first-order in time coupling algorithm we use to
advect the discrete rigid body $\Gamma_{\Delta}$, of boundary $\delta
\Gamma_{\Delta}$, compatible only with the [centered Navier-Stokes
solver](mycentered.h). We assume that there is only one rigid body and
consider two advection scenarios:

* imposed embedded boundary motion;

* explicit weak fluid-solid coupling.

## Setup

In each scenario, the discrete rigid body $\Gamma_{\Delta}$ is
characterized by: the position of the center of mass
$\mathbf{x}_{\Gamma}$ stored in *p_p*, the translation velocity of the
center of mass $\mathbf{u}_{\Gamma}$ stored in *p_u*, the angular
velocity about the center of mass $\mathbf{\omega}_{\Gamma}$ stored in
*p_w* and their corresponding accelerations stored in *p_au, p_aw*. */

coord p_p;        // Position
coord p_u, p_w;   // Velocity
coord p_au, p_aw; // Acceleration

/**
The remaining unknown is the location and shape of the discrete rigid
boundary $\delta \Gamma_{\Delta}$, which we describe using a
user-defined distance function $\Phi$ given by the function
*p_shape()*. Note that this function must contain a call to the
*fractions_cleanup()* function. */

extern void p_shape (scalar c, face vector f, coord p);

/**
## Color field

We define here a color field that helps us track the particle only,
allowing us to define other fixed embedded boundaries such as
walls. */

scalar p_col[];

/**
The following function is case specific and defines the location and
shape of the embedded boundaries associated to the color field
*p_col*. */

extern void p_shape_col (scalar c, coord p); // Particle

/**
## Fluid-solid coupling

The following no-slip Dirichlet boundary condition for velocity and
Neumann boundary condition for pressure on the discrete rigid boundary
$\delta \Gamma_{\Delta}$ allow us to couple the motion of the fluid
and the discrete rigid body $\Gamma_{\Delta}$:

$$
\left\{
\begin{aligned}
&
\mathbf{u} = \mathbf{u}_{\Gamma} = \mathbf{u_{\Gamma}} +
\omega_{\Gamma} \times \left(\mathbf{x} - \mathbf{x}_{\Gamma}\right)
\\
&
{\nabla}_{\Gamma} p = -\rho \frac{\mathrm{D} \mathbf{u}}{\mathrm{D}t}
\cdot \mathbf{n}_{\Gamma} + \left(\nabla \cdot \left(\mu \nabla \mathbf{u}
\right)\right) \cdot \mathbf{n}_{\Gamma} = -\rho \left( \frac{\mathrm{d}
\mathbf{u_{\Gamma}}}{\mathrm{d}t} + \frac{\mathrm{d}
\mathbf{\omega_{\Gamma}}}{\mathrm{d}t} \times \left(\mathbf{x} - \mathbf{x}_{\Gamma}\right) +
\mathbf{\omega_{\Gamma}} \times \mathbf{\omega_{\Gamma}} \times
\left(\mathbf{x} - \mathbf{x}_{\Gamma}\right) \right) \cdot \mathbf{n}_{\Gamma} + \left(\nabla \cdot \left(\mu
\nabla \mathbf{u} \right)\right) \cdot \mathbf{n}_{\Gamma}.
\end{aligned}
\right.
$$

In practice, we however use the following homogeneous Neumann boundary condition for pressure:
$$
{\nabla}_{\Gamma} p = 0
$$
This condition is suitable for a fixed rigid body and we have found
that using it with a moving rigid body does not significantly affect
the computed solution.

#### No-slip Dirichlet boundary condition for velocity

The function *velocity_noslip_x()* computes the previously defined
no-slip boundary condition for the $x$-component of the velocity
$\mathbf{u}$. */

foreach_dimension()
static inline double velocity_noslip_x (Point point,
					coord up, coord wp, coord pp,
					double xc, double yc, double zc)
{
  assert (cs[] > 0. && cs[] < 1.);

  // Colored embedded boundaries
  if (p_col[] > 0. && p_col[] < 1.) {

    /**
    We first compute the relative position $\mathbf{r} = \mathbf{x} -
    \mathbf{x}_{\Gamma}$. */ 
  
    // The coordinate x,y,z are not permuted with foreach_dimension()
    coord r = {xc, yc, zc};
    foreach_dimension() {
      r.x -= pp.x;
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
    return (up.x) + sgn.x*wp.x*(r.y);
#else // dimension == 3
    return (up.x) + wp.y*(r.z) - wp.z*(r.y);
#endif // dimension
  }
  // Other fixed embedded boundaries (u=0)
  return 0.;
}

/**
#### Neumann boundary condition for pressure

The function *pressure_gradient()* returns in *da* the contribution of
all the components of the particle acceleration and viscous stresses
to the pressure gradient. */

static inline void pressure_acceleration (Point point,
					  coord dup, coord dwp,
					  coord wp, coord pp,
					  double xc, double yc, double zc,
					  coord * da)
{
  assert (cs[] > 0. && cs[] < 1.);
  
  // Colored embedded boundaries
  if (p_col[] > 0. && p_col[] < 1.) {

    /**
    We first compute the relative position $\mathbf{r} = \mathbf{x} -
    \mathbf{x}_{\Gamma}$. */ 
  
    // The coordinate x,y,z are not permuted with foreach_dimension()
    coord r = {xc,yc,zc};
    foreach_dimension() {
      r.x -= pp.x;
      if (Period.x) {
	if (fabs (r.x) > fabs (r.x + (L0)))
	  r.x += (L0);
	if (fabs (r.x) > fabs (r.x - (L0)))
	  r.x -= (L0);
      }
    }
    
    /**
    We first compute the acceleration of the particle *gu*. */
    
    coord gu;
    foreach_dimension()
      gu.x = 0.;
    
    // Rotational acceleration
    coord dw1, dw2;
#if dimension == 2
    coord sgn = {-1, 1};
    foreach_dimension() {
      dw1.x = sgn.x*dwp.x*(r.y);
      dw2.x = - sq(wp.x)*(r.x);
    } 
#else // dimension = 3
    double w1 = 0., w2 = 0.;
    foreach_dimension() {
      w1 += wp.x*(r.x);
      w2 += sq(wp.x);
    }
    foreach_dimension() {
      dw1.x = dwp.y*(r.z) - dwp.z*(r.y);
      dw2.x = w1*wp.x - w2*(r.x);
    }
#endif // dimension
    
    /**
    We include the translational and rotational accelerations. */

    foreach_dimension()
      gu.x = dup.x + dw1.x + dw2.x;
    
    /**
    We then compute the viscous contribution. */
    
    coord gmu;
    foreach_dimension()
      gmu.x = 0.;
      
    /* foreach_dimension() { */
    /*   scalar s = u.x; */
    /*   double a = 0.; */
    /*   foreach_dimension() */
    /* 	a += (mu.x[1]*face_gradient_x (s, 1) - mu.x[0]*face_gradient_x (s, 0))/Delta; */
    /*   double b, c = embed_flux (point, u.x, mu, &b); */
    /*   a += (b + c*u.x[]); */
    /*   gmu.x = -a; */
    /* } */
    
    /**
    Finally, we combine both contributions. */

    foreach_dimension()
      da->x = (gu.x + gmu.x);
  }
  // Other fixed embedded boundaries (gradp=0)
  else {
    foreach_dimension()
      da->x = 0.;
  }
}

foreach_dimension()
static inline double pressure_acceleration_x (Point point,
					      coord dup, coord dwp,
					      coord wp, coord pp,
					      double xc, double yc, double zc)
{
  coord da;
  pressure_acceleration (point, dup, dwp, wp, pp, xc, yc, zc, &da);
  return da.x;
}

/**
The following function finally computes the Neumann boundary condition
for pressure, where the inward unit normal $\mathbf{n}_{\Gamma}$ is
pointing from fluid to solid. To switch to a homogenous Neumann
boundary condition, the user can set the scalar attribute
*neumann_zero* to true. The default is false. */

attribute {
  bool neumann_zero;
}

static inline double pressure_neumann (Point point,
				       coord dup, coord dwp,
				       coord wp, coord pp,
				       double xc, double yc, double zc)
{
  coord da;
  pressure_acceleration (point, dup, dwp, wp, pp, xc, yc, zc, &da);
  
  coord b, n;
  embed_geometry (point, &b, &n);
  
  double dpdn = 0.;
  foreach_dimension()
    dpdn += (-da.x)*n.x;
  dpdn *= rho[]/(cs[] + SEPS);
  return dpdn;
}

/**
## Dump and restore a particle */

typedef struct {
  coord c;      // Center of mass
  coord u, w;   // Velocity
  coord au, aw; // Acceleration
} particle;

struct p_Dump {
  char * file; // File name
  particle * list; // List of particles
  FILE * fp; // File pointer
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

  // Get particle data (only 1 particle for now)
  int p_n = 1;
  particle * p_list = p.list;

  // Dump particle data
  fwrite(p_list, sizeof(*p_list), p_n, fp);
  
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

  // Read particle data (only 1 particle for now)
  int p_n = 1;
  particle * p_list = p.list;

  if (fread (p_list, sizeof(*p_list), (p_n), fp) < 1) {
    fprintf (ferr, "#p_restore(): error reading particle data\n");
    exit (1);
  }

  foreach_dimension() {
    p_p.x  = p_list->c.x;
    p_u.x  = p_list->u.x;
    p_w.x  = p_list->w.x;
    p_au.x = p_list->au.x;
    p_aw.x = p_list->aw.x;
  }
  
  /* free (p_list); */
  if (file)
    fclose (fp);
  
  return true;
}

/**
## Initialization */

event defaults (i = 0)
{
#if TREE
  p_col.restriction = restriction_average;
  p_col.refine = p_col.prolongation = fraction_refine;
#endif
  
  foreach_dimension() {
    p_p.x  = 0.;
    p_u.x  = 0.;
    p_w.x  = 0.;
    p_au.x = 0.;
    p_aw.x = 0.;
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
    
    particle pp_restore;
    bool p_restart = p_restore ("p_restart", &pp_restore);
    assert (p_restart == true);

    /**
    We then need to initialize the face fraction *fs* since it is not
    dumped. */
  
    p_shape (cs, fs, p_p);
  }

  /**
  We initialize, even when restarting, the color field. */
  
  p_shape_col (p_col, (p_p));
  boundary    ({p_col});
  restriction ({p_col}); // Since some BC might depend on p_col

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
  
  p[embed]  = neumann (p.neumann_zero ?  0. : pressure_neumann (point, (p_au), (p_aw), (p_w), (p_p), x, y, z));
  pf[embed] = neumann (pf.neumann_zero ? 0. : pressure_neumann (point, (p_au), (p_aw), (p_w), (p_p), x, y, z));

#if dimension == 2
  u.n[embed]   = dirichlet (velocity_noslip_x (point, (p_u), (p_w), (p_p), x, y, z));
  u.t[embed]   = dirichlet (velocity_noslip_y (point, (p_u), (p_w), (p_p), x, y, z));
  uf.n[embed]  = dirichlet (velocity_noslip_x (point, (p_u), (p_w), (p_p), x, y, z));
  uf.t[embed]  = dirichlet (velocity_noslip_y (point, (p_u), (p_w), (p_p), x, y, z));

  g.n[embed] = dirichlet (p.neumann_zero ? 0. : pressure_acceleration_x (point, (p_au), (p_aw), (p_w), (p_p), x, y, z)); 
  g.t[embed] = dirichlet (p.neumann_zero ? 0. : pressure_acceleration_y (point, (p_au), (p_aw), (p_w), (p_p), x, y, z));
#else // dimension == 3
  u.n[embed]   = dirichlet (velocity_noslip_x (point, (p_u), (p_w), (p_p), x, y, z));
  u.t[embed]   = dirichlet (velocity_noslip_y (point, (p_u), (p_w), (p_p), x, y, z));
  u.r[embed]   = dirichlet (velocity_noslip_z (point, (p_u), (p_w), (p_p), x, y, z));
  uf.n[embed]  = dirichlet (velocity_noslip_x (point, (p_u), (p_w), (p_p), x, y, z));
  uf.t[embed]  = dirichlet (velocity_noslip_y (point, (p_u), (p_w), (p_p), x, y, z));
  uf.r[embed]  = dirichlet (velocity_noslip_z (point, (p_u), (p_w), (p_p), x, y, z));
  
  g.n[embed] = dirichlet (p.neumann_zero ? 0. : pressure_acceleration_x (point, (p_au), (p_aw), (p_w), (p_p), x, y, z)); 
  g.t[embed] = dirichlet (p.neumann_zero ? 0. : pressure_acceleration_y (point, (p_au), (p_aw), (p_w), (p_p), x, y, z));
  g.r[embed] = dirichlet (p.neumann_zero ? 0. : pressure_acceleration_z (point, (p_au), (p_aw), (p_w), (p_p), x, y, z));
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
  of mass *p_p*. */

  coord p_o = {(X0), (Y0), (Z0)};
  foreach_dimension() {
    if (Period.x) {
      if (p_p.x < p_o.x)
	p_p.x += (L0);
      if (p_p.x > p_o.x + (L0))
	p_p.x -= (L0);
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
  
  p_shape (cs, fs, p_p);

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

  p_shape_col (p_col, (p_p));
  boundary    ({p_col});
  restriction ({p_col}); // Since some BC might depend on p_col
  
  /**
  ## Step 5 (boundary conditions)

  As *rho* is used in the Neumann boundary condition for pressure, we
  need *rho* on all levels of the grid. */

  restriction ({rho});
}
