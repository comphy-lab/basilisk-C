/**
# Fictitious-domain method with distributed Lagrangian multipliers 

We are looking for a numerical approximation for the solution of the
following equations 
$$
\rho\left(\partial_t{\mathbf{u}} +\mathbf{u}\cdot \mathbf{\nabla}\mathbf{u}\right) = -\mathbf{\nabla} p + \mathbf{\nabla} \cdot\left(2\mu\mathbf{D}\right) + \rho\mathbf{a} + \mathbf{\lambda}~\text{over}~\Omega
$$
$$
\left(1-\frac{\rho}{\rho_s}\right)\left(M\left(\partial_t{\mathbf{U}}-\mathbf{g}\right)\right)=-\int_{P(t)} {\mathbf{\lambda}} dv~\text{over}~P(t)
$$
$$
\left(1-\frac{\rho}{\rho_s}\right)\left(\mathbf{I}\partial_t{\mathbf{\omega}} + \mathbf{\omega} \times\mathbf{I}\mathbf{\omega}\right) = -\int_{P(t)}\mathbf{r}\times{\mathbf{\lambda}} dv~\text{over}~P(t)
$$
$$
\mathbf{u}-\left(\mathbf{U}+\mathbf{\omega}\times \mathbf{r}\right)=0~\text{over}~P(t)
$$
$$
\mathbf{\nabla}\cdot\mathbf{u} = \mathbf{0} ~\text{over}~\Omega
$$

From a numerical point of view, this is a complicated set of equations
to solve direclty. One has to deal with three main difficulties:

  * an advection-diffusion equation
  * the imcompressibility constraint with the pressure $p$ as unknown 
  * the particle's solid body's constraint with the Lagrange
multipliers as unknow.

The classic strategy is to split the problem into subproblems and
solve them successively. We chose here a two-steps time-spliting.
*/

# define BGHOSTS 2
# define BSIZE 128

# ifndef DLM_Moving_particle 
# define DLM_Moving_particle 0
# define TRANSLATION 0
# define ROTATION 0
# endif

# if DLM_Moving_particle
# ifndef TRANSLATION 
# define TRANSLATION 1
# endif
# ifndef ROTATION 
# define ROTATION 1
# endif
# endif

# ifndef NPARTICLES 
# define NPARTICLES 1
# endif

# ifndef DLM_alpha_coupling
# define DLM_alpha_coupling 0
# endif

#ifndef debugBD
# define debugBD 0     // set 1 to desactivate boundary points
#endif
#ifndef debugInterior
# define debugInterior 0  // set 1 to desactivate interior points
#endif

/** Functions and structures needed for the implementation of
    fictitious domain method. */
# include "save_data.h"
# include "dlmfd_functions.h"

/** Basilisk scalars and vectors needed for the implementation of the
    fictitious domain method. */
vector DLM_lambda[];
scalar flagfield[];
scalar flagfield_mailleur[];
vector index_lambda[];
vector DLM_periodic_shift[];

static FILE * pdata[NPARTICLES];
static FILE * fdata[NPARTICLES];
particle particles[NPARTICLES] = {0};

#if DLM_alpha_coupling
vector DLM_explicit[];
#endif
/**
## First sub-problem: Navier Stokes problem

We are using the centred solver for this problem.
$$
\rho\left(\partial_t{\mathbf{u}}+\mathbf{u}\cdot {\mathbf{\nabla}}\mathbf{u}\right) = -{\mathbf{\nabla}}p + {\mathbf{\nabla}}\cdot\left(2\mu\mathbf{D}\right) + \rho\mathbf{a}~\text{over}~\Omega
$$
$$
{\mathbf{\nabla}}\cdot\mathbf{u} = \mathbf{0} ~\text{over}~\Omega. 
$$
*/

# include "navier-stokes/centered.h"

/* Adding this small macro because dv() breaks the compability with embed.h */
#if dimension == 1
# define dlmfd_dv() (Delta)
#elif dimension == 2
# define dlmfd_dv() (sq(Delta))
#else // dimension == 3
# define dlmfd_dv() (cube(Delta))
#endif


/**
## Second sub-problem: Fictitious-domain problem
We are solving a fictitious problem given by:
$$
\rho \partial_t\mathbf{u} = {\mathbf{\lambda}}~\text{over}~\Omega
$$
$$
\left(1-\frac{\rho}{\rho_s}\right)\left(M\left(\partial_t\mathbf{U}-\mathbf{g}\right)\right)=-\int_{P(t)}{\mathbf{\lambda}} dv~\text{over}~P(t)
$$
$$
\left(1-\frac{\rho}{\rho_s}\right)\left(\mathbf{I}\partial_t {\mathbf{\omega}} + {\mathbf{\omega}} \times\mathbf{I}{\mathbf{\omega}}\right) = -\int_{P(t)}\mathbf{r}\times{\mathbf{\lambda}} dv~\text{over}~P(t)
$$
$$
\mathbf{u}-\left(\mathbf{U}+{\mathbf{\omega}}\times \mathbf{r}\right)=0~\text{over}~P(t),
$$
with unknowns $\mathbf{u}, \mathbf{U}, {\mathbf{\omega}}$ and ${\mathbf{\lambda}}$ being respectively the fluid velocity field, the particle's translational velocity, the particle's rotational velocity and the Lagrange multipliers. The particle occupies the domain $P(t)$ with density $\rho_s$, Inertia tensor $\mathbf{I}$ and mass $M$. The vector $\mathbf{g}$ is the gravity acceleration here.

This leads to a saddle-point problem which is solved with an iterative solver (Uzawa, conjugate gradient algorithm). It is implemented in the function below. 
*/


/* Timers and Timings */
timing dlmfd_globaltiming = {0.};

# include "dlmfd_perf.h"

void DLMFD_subproblem (particle * p, const int i, const double rho_f) {

  /* Timers and Timings */
  timer dlmfd_timer = timer_start ();
  double mpitimings[npe()];
  
  double DLM_alpha = 0., DLM_beta = 0.;
  double DLM_tol = 1.e-5, DLM_nr2 = 0., DLM_nr2_km1 = 0., DLM_wv = 0.;
  int ki = 0;
  int DLM_maxiter = 100;
  vector qu[], tu[];
  coord ppshift = {0, 0, 0};
  
  static FILE * converge = fopen ("converge-uzawa", "w");
  static FILE * cellvstime = fopen ("dlmfd-cells", "w");
 
  /* Allocate and initialize particles */
  allocate_and_init_particles (p, NPARTICLES, index_lambda, flagfield, flagfield_mailleur, DLM_periodic_shift);

#if debugInterior == 0
  Cache * Interior[NPARTICLES];
#endif

#if debugBD == 0
  Cache * Boundary[NPARTICLES];
  SolidBodyBoundary * sbm[NPARTICLES];
#endif
  GeomParameter * gci[NPARTICLES];

  
#if DLM_Moving_particle
#if TRANSLATION
  coord * qU[NPARTICLES];
  coord * U[NPARTICLES];
  coord * tU[NPARTICLES];
#endif
#if ROTATION
  coord * qw[NPARTICLES];
  coord * w[NPARTICLES];
  coord * tw[NPARTICLES];
#endif
#endif
  
  for (int k = 0; k < NPARTICLES; k++) {
    gci[k] = &(p[k].g);
    
#if debugInterior == 0
    Interior[k] = &(p[k].Interior);
#endif
#if debugBD == 0 
    Boundary[k] = &(p[k].reduced_domain);
    sbm[k] = &(p[k].s);
#endif
    
#if DLM_Moving_particle
#if TRANSLATION
    qU[k] = &(p[k].qU);
    U[k] = &(p[k].U);
    tU[k] = &(p[k].tU);
#endif
#if ROTATION
    qw[k] = &(p[k].qw);
    w[k] = &(p[k].w);
    tw[k] = &(p[k].tw);
#endif
#endif
  }

  if (NPARTICLES > 1)
    remove_too_close_multipliers (p, index_lambda);

#if debugBD == 0   
#if DLM_Moving_particle
  reverse_fill_flagfield (p, flagfield_mailleur, index_lambda, 0);
#endif
  
  /* Consider fictitious domain's boundary points only if they are far
     enough of the domain's solid boundaries (i.e a distance greater
     than 2*Delta, with Delta being the smallest cell size) otherwise
     the fictitious domain is pixalated. This does not apply when the
     box is periodic */
  double twodelta;

# if dimension == 2
  if (!Period.x && !Period.y) {
#elif dimension == 3 
  if (!Period.x && !Period.y && !Period.z) {
#endif
    foreach_level (depth()) {
      twodelta = 2.*Delta;

#if dimension == 2
      if (( x > (L0 - twodelta + X0)) || (x  < (twodelta + X0)) || ( y > (L0 - twodelta + Y0)) || (y < (twodelta + Y0)))
#elif dimension == 3
	if ((x  > (L0 - twodelta + X0)) ||  (x  < (twodelta + X0)) || (y > (L0 - twodelta + Y0)) || (y < (twodelta + Y0)) || (z > (L0 - twodelta +Z0)) || (z  < (twodelta + Z0)))
#endif
	  {
	    if (index_lambda.x[] > -1.) {
	      index_lambda.x[] = -1.;
	    }
	  }
    }
  }
  boundary((scalar*) {index_lambda});


  /* Tagg the cells belonging to the stencil of the boundary-multipliers points  */
  reverse_fill_flagfield (p, flagfield, index_lambda, 1);
#endif
  
  /* Alllocate arrays for CG/Uzawa. */
  vector DLM_r[];
  vector DLM_w[];
  vector DLM_v[];

  /* Create fictitious-domain's cache for the interior domain. */
#if debugInterior == 0 
  for (int k = 0; k < NPARTICLES; k++) {

    if (!(p[k].iswall) && !(p[k].iscube)) { 
      create_FD_Interior_Particle (&p[k], index_lambda, DLM_periodic_shift, flagfield);
    }

    if (p[k].iswall) {
      create_FD_Interior_Wall(&p[k], index_lambda, flagfield);
    }

    if (p[k].iscube) {
      create_FD_Interior_Cube_v2 (&p[k], index_lambda, DLM_periodic_shift);
    }
  }
#endif

  boundary ((scalar*) {DLM_periodic_shift});

  int allpts = 0; int lm = 0; int tcells = 0;
  foreach() {
    tcells ++; DLM_lambda.x[] = 0.; DLM_r.x[] = 0.; DLM_w.x[] = 0.; DLM_v.x[] = 0.; qu.x[] = 0.; tu.x[] = 0.;
    DLM_lambda.y[] = 0.; DLM_r.y[] = 0.; DLM_w.y[] = 0.; DLM_v.y[] = 0.; qu.y[] = 0.; tu.y[] = 0.;
    DLM_lambda.z[] = 0.; DLM_r.z[] = 0.; DLM_w.z[] = 0.; DLM_v.z[] = 0.; qu.z[] = 0.; tu.z[] = 0.;
  }
  
#if _MPI
  mpi_all_reduce (tcells, MPI_INTEGER, MPI_SUM);
#endif
  
  lm = total_dlmfd_multipliers (p, NPARTICLES);
  /* Get track of the number of multipliers for statistics */
  for (int k = 0; k < NPARTICLES; k++)
    p[k].tmultipliers += lm;

  allpts = total_dlmfd_cells (p, NPARTICLES);
  /* Get track of the number of cells involved in the dlmfd solver for statistics */
  for (int k = 0; k < NPARTICLES; k++)
    p[k].tcells += allpts;

  fprintf (stderr,"# Total Lagrange multipliers: %d, constraint cells: %d, Total cells: %d \n", lm, allpts, tcells);

  if (i == 0)
    fprintf (cellvstime,"# Time Iteration \t Lagrange multipliers \t Constraint cells \t Total cells\n");
  fprintf (cellvstime,"%d \t %d \t %d \t %d \n", i, lm, allpts, tcells);
  /* fflush (cellvstime); */
  

#if DLM_Moving_particle
#if TRANSLATION
  for (int k = 0; k < NPARTICLES; k++) {
    (*qU[k]).x = 0.; (*tU[k]).x = 0.;
    (*qU[k]).y = 0.; (*tU[k]).y = 0.;
    (*qU[k]).z = 0.; (*tU[k]).z = 0.;
  }
#endif
#if ROTATION
  for (int k = 0; k < NPARTICLES; k++) {
    (*qw[k]).x = 0.; (*tw[k]).x = 0.;
    (*qw[k]).y = 0.; (*tw[k]).y = 0.;
    (*qw[k]).z = 0.; (*tw[k]).z = 0.;
  }
#endif
#endif

  /** # Iterative Uzawa/conjugate-gradient solver 
  Implementation of the iterative Uzawa/conjugate gradient solver for
  the fictitious domain problem. The code below follows closely the
  steps given in the annex of the paper "Wachs et al J.Eng. Math,
  2011". (minus all the sign errors)*/
 
  /* Iterative CG/Uzawa */
  /* Initialization */
  /* ------------- */
  /* Compute qu = fu - (M_u^T)*lambda^0 with fu = rho*dV/dt*u^n+1/2 +
     alpha*lambda^n. */ 
  /* Note that (M_u^T)*lambda <=> <lambda,v>_P(t) (v is the test function
     for u)*/
  
  /* For moving particles add:  */
  /* qU = fU -(M_U^T)*lambda^0 with fU = (1 - rho_f /
     rho_s)*M*(U^{n+1/2}/dt + g) with M the particle's mass and g the
     gravity. */
  /* Note that (M_U^T)*lambda <=> - <lambda, V>_P(t) (V is the test
     function for U)*/
  /* qw = fw -(M_w^T)*lambda^0 with fw=(1-rho_f/rho_s)/dt*Ip*w^n+1/2
   * with Ip the Inertia tensor.*/
  /* Note that (M_w^T)*lambda <=> - <lambda, xi^r_GM>_P(t) (xi is the test
   * function for w)*/


  
  foreach() {
    foreach_dimension() {
      qu.x[] = rho_f*dlmfd_dv()*u.x[]/dt;

#if DLM_alpha_coupling
      qu.x[] += DLM_explicit.x[];
      DLM_explicit.x[] = 0.;
#endif
    }
  }

  boundary((scalar*){qu});
  
  /* Interior points qu = fu -(M_u^T)*lambda^0  */
  /* qu = fu -(M_u^T)*lambda^0 = fu -<lambda,v>_P(t) */

  /* Interior points qU = fU -(M_U^T)*lambda^0  */
  /* qU = fU -(M_U^T)*lambda^0 = fU  + <lambda,V>_P(t) */
  
  /* Interior points qw = fw -(M_w^T)*lambda^0  */
  /* qw = fw -(M_w^T)*lambda^0 = fw  + <lambda,xi^r_GM>_P(t) */
  
  /* For moving particles the fU and fw parts are added after the
     scalar product (for mpi purpose) */

#if debugInterior == 0
  for (int k = 0; k < NPARTICLES; k++) {
    foreach_cache((*Interior[k])) {
      if ((flagfield[]  < 1) && ((int)index_lambda.y[] == k)) {
	foreach_dimension() {
	  qu.x[] -= DLM_lambda.x[];

#if DLM_Moving_particle
#if TRANSLATION
	  (*qU[k]).x += DLM_lambda.x[];
#endif
#endif
	}
	
#if DLM_Moving_particle
#if ROTATION
	/* Modify temporerly the particle center position for periodic boundary condition */
	foreach_dimension()
	  (*gci[k]).center.x += DLM_periodic_shift.x[];
	  
	/* as a.(b^c) = c.(a^b) = b.(c^a) */
	/* <lambda,xi^r_GM>_P(t) = <r_GM, lambda^xi>_P(t) = <xi,r_GM^lambda>_P(t) */
	(*qw[k]).x += (DLM_lambda.z[]*(y - (*gci[k]).center.y) - DLM_lambda.y[]*(z - (*gci[k]).center.z)); // lambda_z*r_y - lambda_y*r_z 
	(*qw[k]).y += (DLM_lambda.x[]*(z - (*gci[k]).center.z) - DLM_lambda.z[]*(x - (*gci[k]).center.x)); // lambda_x*r_z - lambda_z*r_x 
	(*qw[k]).z += (DLM_lambda.y[]*(x - (*gci[k]).center.x) - DLM_lambda.x[]*(y - (*gci[k]).center.y)); // lambda_y*r_x - lambda_x*r_y

	foreach_dimension()
	  (*gci[k]).center.x -= DLM_periodic_shift.x[];
#endif
#endif
      }
    }
  }
#endif

  /* Boundary points qu = fu -(M_u^T)*lambda^0 */
  /* Boundary points qU = fU -(M_U^T)*lambda^0 */
  /* Boundary points qw = fw -(M_w^T)*lambda^0 */
  
# if debugBD == 0
  double weight = 0.;
  coord weightcellpos = {0, 0, 0};
  coord lambdacellpos = {0, 0, 0};
  coord lambdapos = {0, 0, 0};
  coord sum = {0, 0, 0};
  boundary ((scalar*) {DLM_lambda});

  for (int k = 0; k < NPARTICLES; k++) {
    particle * pp = &p[k];
    
    foreach_cache ((*Boundary[k])) {
   
      sum.x = 0; sum.y = 0; sum.z = 0; 
      weightcellpos.x = x; weightcellpos.y = y; weightcellpos.z = 0.; 
#if dimension == 3
      weightcellpos.z = z;
#endif
     
      foreach_neighbor() {
	if (((int)index_lambda.x[] > -1) && (level == depth()) && is_leaf(cell) && ((int)index_lambda.y[] == pp->pnum)) {
	  lambdacellpos.x = x; lambdacellpos.y = y; lambdacellpos.z = 0.;
	  lambdapos.x = (*sbm[k]).x[(int)index_lambda.x[]];
	  lambdapos.y = (*sbm[k]).y[(int)index_lambda.x[]];
	  lambdapos.z = 0.;
	  
#if dimension == 3
	  lambdacellpos.z = z;
	  lambdapos.z = (*sbm[k]).z[(int)index_lambda.x[]];
#endif
	  
	 
	  foreach_dimension() {
	    ppshift.x = 0.;
	    if (lambdacellpos.x < 0.) ppshift.x = -L0;
	    if (lambdacellpos.x > L0) ppshift.x = +L0;
	  }
	  
	  weight = reversed_weight (pp, weightcellpos, lambdacellpos, lambdapos, Delta, ppshift,1,1);
	      
	  foreach_dimension()
	    sum.x += weight*DLM_lambda.x[];
	}
      }
      
      qu.x[] -= sum.x; qu.y[] -= sum.y; qu.z[] -= sum.z;
    
#if DLM_Moving_particle
      if (index_lambda.x[] > -1 && ((int)index_lambda.y[] == pp->pnum)) {
	lambdapos.x = (*sbm[k]).x[(int)index_lambda.x[]];
	lambdapos.y = (*sbm[k]).y[(int)index_lambda.x[]];
	lambdapos.z = 0.;
#if dimension == 3
	lambdapos.z = (*sbm[k]).z[(int)index_lambda.x[]];
#endif

#if TRANSLATION
	foreach_dimension()
	  (*qU[k]).x += DLM_lambda.x[];
#endif
#if ROTATION
	
	/* Modify temporerly the particle center position for periodic boundary condition */
	foreach_dimension()
	  (*gci[k]).center.x += DLM_periodic_shift.x[];

	/* as a.(b^c) = c.(a^b) = b.(c^a) */
	/* <lambda,xi^r_GM>_P(t) = <r_GM, lambda^xi>_P(t) = <xi,r_GM^lambda>_P(t) */
	(*qw[k]).x += (DLM_lambda.z[]*(lambdapos.y - (*gci[k]).center.y) - DLM_lambda.y[]*(lambdapos.z - (*gci[k]).center.z)); //lambda_z*r_y - lambda_y*r_z
	(*qw[k]).y += (DLM_lambda.x[]*(lambdapos.z - (*gci[k]).center.z) - DLM_lambda.z[]*(lambdapos.x - (*gci[k]).center.x)); //lambda_x*r_z - lambda_z*r_x
	(*qw[k]).z += (DLM_lambda.y[]*(lambdapos.x - (*gci[k]).center.x) - DLM_lambda.x[]*(lambdapos.y - (*gci[k]).center.y)); //lambda_y*r_x - lambda_x*r_y 

	foreach_dimension()
	  (*gci[k]).center.x -= DLM_periodic_shift.x[];
#endif
      }
#endif
      
	
    }
  }
#endif

  /* Broadcast through mpi */
#if DLM_Moving_particle
#if _MPI
  for (int k = 0; k < NPARTICLES; k++) {
#if TRANSLATION
    mpi_all_reduce ((*qU[k]).x, MPI_DOUBLE, MPI_SUM);
    mpi_all_reduce ((*qU[k]).y, MPI_DOUBLE, MPI_SUM);
    mpi_all_reduce ((*qU[k]).z, MPI_DOUBLE, MPI_SUM);
#endif
#if ROTATION
    mpi_all_reduce ((*qw[k]).x, MPI_DOUBLE, MPI_SUM);
    mpi_all_reduce ((*qw[k]).y, MPI_DOUBLE, MPI_SUM);
    mpi_all_reduce ((*qw[k]).z, MPI_DOUBLE, MPI_SUM);
#endif
  }
#endif

  /* Adding here the fU and fw parts to qU and qw */
  for (int k = 0; k < NPARTICLES; k++) {
#if TRANSLATION
    /* (*qU[k]).x += (1. - (rho_f/p[k].rho_s))*(p[k].M)*(p[k].gravity.x + p[k].adforce.x + (*U[k]).x/dt) ; */
    /* (*qU[k]).y += (1. - (rho_f/p[k].rho_s))*(p[k].M)*(p[k].gravity.y + p[k].adforce.y + (*U[k]).y/dt) ; */
    /* (*qU[k]).z += (1. - (rho_f/p[k].rho_s))*(p[k].M)*(p[k].gravity.z + p[k].adforce.z + (*U[k]).z/dt) ; */

    /* Modification for the explicit added term  */
    (*qU[k]).x += (1. - (rho_f/p[k].rho_s))*(p[k].M)*(p[k].gravity.x + p[k].adforce.x ) + p[k].DLMFD_couplingfactor * p[k].M * (*U[k]).x / dt ;
    (*qU[k]).y += (1. - (rho_f/p[k].rho_s))*(p[k].M)*(p[k].gravity.y + p[k].adforce.y ) + p[k].DLMFD_couplingfactor * p[k].M * (*U[k]).y / dt ;
    (*qU[k]).z += (1. - (rho_f/p[k].rho_s))*(p[k].M)*(p[k].gravity.z + p[k].adforce.z ) + p[k].DLMFD_couplingfactor * p[k].M * (*U[k]).z / dt ;
    
#endif
#if ROTATION
    /* The inertia tensor is */
    /*  Ixx -Ixy -Ixz */
    /* -Iyx  Iyy -Iyz */
    /* -Izx -Izy  Izz */ 
    /* with */
    /* Ip[0] = Ixx */
    /* Ip[1] = Iyy */
    /* Ip[2] = Izz */
    /* Ip[3] = Ixy */
    /* Ip[4] = Ixz */
    /* Ip[5] = Iyz */
    /* (*qw[k]).x += (1. - (rho_f/p[k].rho_s))*( (p[k].Ip[0])*(*w[k]).x - (p[k].Ip[3])*(*w[k]).y - (p[k].Ip[4])*(*w[k]).z)/dt; */
    /* (*qw[k]).y += (1. - (rho_f/p[k].rho_s))*(-(p[k].Ip[3])*(*w[k]).x + (p[k].Ip[1])*(*w[k]).y - (p[k].Ip[5])*(*w[k]).z)/dt; */
    /* (*qw[k]).z += (1. - (rho_f/p[k].rho_s))*(-(p[k].Ip[4])*(*w[k]).x - (p[k].Ip[5])*(*w[k]).y + (p[k].Ip[2])*(*w[k]).z)/dt; */

     /* Modification for the explicit added term  */
    (*qw[k]).x += p[k].DLMFD_couplingfactor * ( (p[k].Ip[0])*(*w[k]).x - (p[k].Ip[3])*(*w[k]).y - (p[k].Ip[4])*(*w[k]).z) / dt;
    (*qw[k]).y += p[k].DLMFD_couplingfactor * (-(p[k].Ip[3])*(*w[k]).x + (p[k].Ip[1])*(*w[k]).y - (p[k].Ip[5])*(*w[k]).z) / dt;
    (*qw[k]).z += p[k].DLMFD_couplingfactor * (-(p[k].Ip[4])*(*w[k]).x - (p[k].Ip[5])*(*w[k]).y + (p[k].Ip[2])*(*w[k]).z) / dt;
   
#endif
  }
#endif


  /* Invert L*u^0 = qu with L=dv*rho/dt*Identity_Matrix */
  foreach() {
    foreach_dimension() {
      u.x[] = qu.x[]*dt/(rho_f*dlmfd_dv());
    }
  }
  
  /* Invert A_U U^0 = qU with A_U = ((1 - rho_f/rho_s)*M)/dt *
     Identity_Matrix */
  
  /* Invert A_w w^0 = qw with A_w = ((1 - rho_f/rho_s)*Ip)/dt *
     Identity_Matrix */
  
#if DLM_Moving_particle
  for (int k = 0; k < NPARTICLES; k++) {
#if TRANSLATION
    /* (*U[k]).x = (dt*(*qU[k]).x)/((1. - rho_f/p[k].rho_s)*p[k].M); */
    /* (*U[k]).y = (dt*(*qU[k]).y)/((1. - rho_f/p[k].rho_s)*p[k].M); */
    /* (*U[k]).z = (dt*(*qU[k]).z)/((1. - rho_f/p[k].rho_s)*p[k].M); */

    /* Modification for the explicit added term  */
    (*U[k]).x = (dt*(*qU[k]).x)/(p[k].DLMFD_couplingfactor*p[k].M);
    (*U[k]).y = (dt*(*qU[k]).y)/(p[k].DLMFD_couplingfactor*p[k].M);
    (*U[k]).z = (dt*(*qU[k]).z)/(p[k].DLMFD_couplingfactor*p[k].M);
#endif
#if ROTATION
    /* The inertia tensor is */
    /*  Ixx -Ixy -Ixz */
    /* -Iyx  Iyy -Iyz */
    /* -Izx -Izy  Izz */ 
    /* with */
    /* Ip[0] = Ixx */
    /* Ip[1] = Iyy */
    /* Ip[2] = Izz */
    /* Ip[3] = Ixy */
    /* Ip[4] = Ixz */
    /* Ip[5] = Iyz */

    /* Allocate temporly a matrix for Ip */
    double ** Imat = malloc(3 * sizeof(double*));
    for (int ii = 0; ii < 3; ii++) {
      Imat[ii] = malloc(3 * sizeof(double));
    }

    Imat[0][0] = p[k].Ip[0];
    Imat[1][1] = p[k].Ip[1];
    Imat[2][2] = p[k].Ip[2];

    Imat[0][1] = -p[k].Ip[3];
    Imat[0][2] = -p[k].Ip[4];
    Imat[1][2] = -p[k].Ip[5];
    
    Imat[1][0] = Imat[0][1]; 
    Imat[2][0] = Imat[0][2];
    Imat[2][1] = Imat[1][2];

    /* A_w = ((1 - rho_f/rho_s)*Ip)/dt */
    
    double Aw[3][3];
    for (int ll = 0; ll < 3; ll++) {
      for (int nn = 0; nn < 3; nn++) {
	/* Aw[ll][nn] = (1 - rho_f/p[k].rho_s)*Imat[ll][nn]/dt; */
	/* Modification for the explicit added term  */
	Aw[ll][nn] = p[k].DLMFD_couplingfactor * Imat[ll][nn] / dt;
      }
    }

    /* Takes in Aw and return Aw^-1 as Imat */
    inverse3by3matrix(Aw, Imat);
  
    (*w[k]).x = Imat[0][0]*(*qw[k]).x + Imat[0][1]*(*qw[k]).y + Imat[0][2]*(*qw[k]).z;
    (*w[k]).y = Imat[1][0]*(*qw[k]).x + Imat[1][1]*(*qw[k]).y + Imat[1][2]*(*qw[k]).z;
    (*w[k]).z = Imat[2][0]*(*qw[k]).x + Imat[2][1]*(*qw[k]).y + Imat[2][2]*(*qw[k]).z;

    for (int ii = 0; ii < 3; ii++) 
      free(Imat[ii]);
    free(Imat);
#endif
  }
#endif
  
  /* Compute residual r^0 = G - M_u*u^0 - M_U*U^0 - M_w*w^0 */
  /* Note that M_u*u^0 <=> <alpha, u>_P(t) */
  /* Note that M_U*U^0 <=> - <alpha, U>_P(t) */
  /* Note that M_w*w^0 <=> - <alpha, w^r_GM>_P(t) */
  /* alpha is he test function for lambda the Lagrange multipliers */

  /* So r^0 = G -<alpha, u>_P(t) + <alpha, U>_P(t) + <alpha, w^r_GM>_P(t) */

  /* In the 2D case when rotation is enabled, the z-component of the
   * residual becomes involved in the algorithm, this is why the
   * foreach_dimension() iterators are being removed from this
   * point. Otherwise the z components would be discarded */


  /* Interior points: r^0 = G - M_u*u^0 - M_U*U^0 - M_w*w^0 */
  /* So r^0 = G -<alpha, u>_P(t) + <alpha, U>_P(t) + <alpha, w^r_GM>_P(t) */
#if debugInterior == 0
  for (int k = 0; k < NPARTICLES; k++) {
    
#if !DLM_Moving_particle
    coord imposedU = (p[k]).imposedU;
    coord imposedw = (p[k]).imposedw;
#endif
    
    foreach_cache((*Interior[k])) {
      if ((flagfield[]  < 1) && ((int)index_lambda.y[] == k)) {
	foreach_dimension() {
	  DLM_r.x[] = - u.x[];
 
#if DLM_Moving_particle
#if TRANSLATION
	  DLM_r.x[] +=  (*U[k]).x;
#endif
#else
	  DLM_r.x[] += imposedU.x;
#endif
	}
	
#if DLM_Moving_particle
#if ROTATION
	/* Modify temporerly the particle center position for periodic boundary condition */
	foreach_dimension()
	  (*gci[k]).center.x += DLM_periodic_shift.x[];

	/* if (DLM_periodic_shift.x[] < 0) { */
	/*   printf ("position (%f,%f,%f), center shifted of %f\n", x, y, z, DLM_periodic_shift.x[]); */
	/*   printf ("new sphere position (%f,%f,%f)\n", (*gci[k]).center.x, (*gci[k]).center.y, (*gci[k]).center.z); */
	/* } */
	DLM_r.x[] += ((*w[k]).y*(z - (*gci[k]).center.z) - (*w[k]).z*(y - (*gci[k]).center.y)); // w_y*r_z - w_z*r_y
	DLM_r.y[] += ((*w[k]).z*(x - (*gci[k]).center.x) - (*w[k]).x*(z - (*gci[k]).center.z)); // w_z*r_x - w_x*r_z
	DLM_r.z[] += ((*w[k]).x*(y - (*gci[k]).center.y) - (*w[k]).y*(x - (*gci[k]).center.x)); // w_x*r_y - w_y*r_x

	foreach_dimension()
	  (*gci[k]).center.x -= DLM_periodic_shift.x[];
#endif
	
#else
	/* Modify temporerly the particle center position for periodic boundary condition */
	foreach_dimension()
	  (*gci[k]).center.x += DLM_periodic_shift.x[];

	DLM_r.x[] += (imposedw.y*(z - (*gci[k]).center.z) - imposedw.z*(y - (*gci[k]).center.y));
	DLM_r.y[] += (imposedw.z*(x - (*gci[k]).center.x) - imposedw.x*(z - (*gci[k]).center.z));
	DLM_r.z[] += (imposedw.x*(y - (*gci[k]).center.y) - imposedw.y*(x - (*gci[k]).center.x));

	foreach_dimension()
	  (*gci[k]).center.x -= DLM_periodic_shift.x[];
	
#endif
      }
    }
  }
#endif
  
  /* Boundary points: r^0 = G - M_u*u^0 - M_U*U^0 - M_w*w^0 */
  /* So r^0 = G -<alpha, u>_P(t) + <alpha, U>_P(t) + <alpha, w^r_GM>_P(t) */
#if debugBD == 0
  double testweight = 0.;
  boundary((scalar*){u});

  for (int k = 0; k < NPARTICLES; k++) {
    particle * pp = &p[k];
#if !DLM_Moving_particle
    coord imposedU = (p[k]).imposedU;
    coord imposedw = (p[k]).imposedw;
#endif
    foreach_cache((*Boundary[k])) {
      if (index_lambda.x[] > -1 && ((int)index_lambda.y[] == pp->pnum)) {
	lambdacellpos.x = x; lambdacellpos.y = y; lambdacellpos.z = 0.;
	lambdapos.x = (*sbm[k]).x[(int)index_lambda.x[]];
	lambdapos.y = (*sbm[k]).y[(int)index_lambda.x[]];
	lambdapos.z = 0.;

#if dimension == 3
	lambdacellpos.z = z;
	lambdapos.z = (*sbm[k]).z[(int)index_lambda.x[]];
#endif
	sum.x = 0; sum.y = 0; sum.z = 0;
	testweight = 0.;

	foreach_dimension()
	    ppshift.x = DLM_periodic_shift.x[];
	
	foreach_neighbor() {
	  if ((level == depth())) {
	    weightcellpos.x = x; weightcellpos.y = y; weightcellpos.z = 0.;
#if dimension == 3
	    weightcellpos.z = z;
#endif
    
	    weight = reversed_weight (pp, weightcellpos, lambdacellpos, lambdapos, Delta, ppshift, 1., 0.);

	    testweight += weight;
	    foreach_dimension()
	      sum.x += weight*u.x[];
	  }
	}

	if ((testweight < 1. - 0.000001) || (testweight > 1. + 0.000001))
	  printf("testweight = %f\n",testweight);
	
	  DLM_r.x[] =  -sum.x; DLM_r.y[] =  -sum.y; DLM_r.z[] =  -sum.z;
	
		
#if DLM_Moving_particle
#if TRANSLATION
	foreach_dimension()
	  DLM_r.x[] +=  (*U[k]).x;
#endif
#if ROTATION
	
	/* Modify temporerly the particle center position for periodic boundary condition */
	foreach_dimension()
	  (*gci[k]).center.x += DLM_periodic_shift.x[];
		
	DLM_r.x[] +=  ((*w[k]).y*(lambdapos.z - (*gci[k]).center.z) - (*w[k]).z*(lambdapos.y - (*gci[k]).center.y)); // w_y*r_z - w_z*r_y
	DLM_r.y[] +=  ((*w[k]).z*(lambdapos.x - (*gci[k]).center.x) - (*w[k]).x*(lambdapos.z - (*gci[k]).center.z)); // w_z*r_x - w_x*r_z
	DLM_r.z[] +=  ((*w[k]).x*(lambdapos.y - (*gci[k]).center.y) - (*w[k]).y*(lambdapos.x - (*gci[k]).center.x)); // w_x*r_y - w_y*r_x

	foreach_dimension()
	  (*gci[k]).center.x -= DLM_periodic_shift.x[];
	
#endif
#else
	foreach_dimension() {
	  DLM_r.x[] += imposedU.x;
	}

	/* Modify temporerly the particle center position for periodic boundary condition */
	foreach_dimension()
	  (*gci[k]).center.x += DLM_periodic_shift.x[];
	
	DLM_r.x[] +=  (imposedw.y*(lambdapos.z - (*gci[k]).center.z) - imposedw.z*(lambdapos.y - (*gci[k]).center.y));
	DLM_r.y[] +=  (imposedw.z*(lambdapos.x - (*gci[k]).center.x) - imposedw.x*(lambdapos.z - (*gci[k]).center.z));
	DLM_r.z[] +=  (imposedw.x*(lambdapos.y - (*gci[k]).center.y) - imposedw.y*(lambdapos.x - (*gci[k]).center.x));

	foreach_dimension()
	  (*gci[k]).center.x -= DLM_periodic_shift.x[];
	
#endif
      }
    }
  }
#endif
  boundary((scalar*){DLM_r});
  
  DLM_nr2 = 0.;
  /* set DLM_w = DLM_r */
  foreach() {
    DLM_w.x[] = DLM_r.x[];
    DLM_w.y[] = DLM_r.y[];
    DLM_w.z[] = DLM_r.z[];    
    DLM_nr2 += sq(DLM_r.x[]) + sq(DLM_r.y[]) + sq(DLM_r.z[]);
  }


#if _MPI
  mpi_all_reduce(DLM_nr2, MPI_DOUBLE, MPI_SUM);
#endif
  

  /* Iterative loop */
  /* -------------- */

  for (ki = 1; ki < DLM_maxiter && (sqrt(DLM_nr2) > DLM_tol); ki++) {

    /* (1) Initialize and compute qu = (M^T)*DLM_w */
    foreach() {
      qu.x[] = 0.; qu.y[] = 0.; qu.z[] = 0.;
      tu.x[] = 0.; tu.y[] = 0.; tu.z[] = 0.;
    }
    
#if DLM_Moving_particle
    for (int k = 0; k < NPARTICLES; k++) {
#if TRANSLATION
      (*qU[k]).x = 0.; (*tU[k]).x = 0.;
      (*qU[k]).y = 0.; (*tU[k]).y = 0.;
      (*qU[k]).z = 0.; (*tU[k]).z = 0.;
#endif
#if ROTATION
      (*qw[k]).x = 0.; (*tw[k]).x = 0.;
      (*qw[k]).y = 0.; (*tw[k]).y = 0.;
      (*qw[k]).z = 0.; (*tw[k]).z = 0.;
#endif
    }
#endif
     
    /* Interior points qu = (M_u^T)DLM_w =   <DLM_w, v>_P(t) */
    /* Interior points qU = (M_U^T)DLM_w = - <DLM_w, V>_P(t) */
    /* Interior points qw = (M_w^T)DLM_w = - <DLM_w, xi^r_GM>_P(t) */
    /* - <DLM_w, xi^r_GM>_P(t) = - <r_GM, DLM_w^xi>_P(t) = - <xi, r_GM^DLM_w>_P(t) */
    
#if debugInterior == 0
  for (int k = 0; k < NPARTICLES; k++) {
    foreach_cache((*Interior[k])) {
      if ((flagfield[]  < 1) && ((int)index_lambda.y[] == k)) {
	qu.x[] = DLM_w.x[];
	qu.y[] = DLM_w.y[];
	qu.z[] = DLM_w.z[];
	  
#if DLM_Moving_particle
#if TRANSLATION
	(*qU[k]).x += -DLM_w.x[];
	(*qU[k]).y += -DLM_w.y[];
	(*qU[k]).z += -DLM_w.z[];
#endif
#endif
			
#if DLM_Moving_particle
#if ROTATION

	/* Modify temporerly the particle center position for periodic boundary condition */
	foreach_dimension()
	  (*gci[k]).center.x += DLM_periodic_shift.x[];
	
	(*qw[k]).x +=  -(DLM_w.z[]*(y - (*gci[k]).center.y) - DLM_w.y[]*(z - (*gci[k]).center.z)); // -(w_z*r_y - w_y*r_z)
	(*qw[k]).y +=  -(DLM_w.x[]*(z - (*gci[k]).center.z) - DLM_w.z[]*(x - (*gci[k]).center.x)); // -(w_x*r_z - w_z*r_x)
	(*qw[k]).z +=  -(DLM_w.y[]*(x - (*gci[k]).center.x) - DLM_w.x[]*(y - (*gci[k]).center.y)); // -(w_y*r_x - w_x*r_y)

	foreach_dimension()
	  (*gci[k]).center.x -= DLM_periodic_shift.x[];
	
#endif
#endif
      }
    }
  }
#endif
  
  /* Boundary points qu = (M_u^T)DLM_w =   <DLM_w, v>_P(t) */
  /* Boundary points qU = (M_U^T)DLM_w = - <DLM_w, V>_P(t) */
  /* Boundary points qw = (M_w^T)DLM_w = - <DLM_w, xi^r_GM>_P(t) */
   /* - <DLM_w, xi^r_GM>_P(t) = - <r_GM, DLM_w^xi>_P(t) = - <xi, r_GM^DLM_w>_P(t) */
  
# if debugBD == 0
  boundary ((scalar*){DLM_w, qu});
  for (int k = 0; k < NPARTICLES; k++) {
    particle * pp = &p[k];
    foreach_cache((*Boundary[k])) {
      sum.x = 0; sum.y = 0; sum.z = 0;
      weightcellpos.x = x; weightcellpos.y = y; weightcellpos.z = 0.;
#if dimension == 3
      weightcellpos.z = z;
#endif
      foreach_neighbor() {
	if(((int)index_lambda.x[] > -1) && (level == depth()) && (is_leaf(cell)) && ((int)index_lambda.y[] == pp->pnum)) {
	  lambdacellpos.x = x;
	  lambdacellpos.y = y;
	  lambdacellpos.z = 0.;
	  lambdapos.x = (*sbm[k]).x[(int)index_lambda.x[]];
	  lambdapos.y = (*sbm[k]).y[(int)index_lambda.x[]];
	  lambdapos.z = 0.;

#if dimension == 3
	  lambdacellpos.z = z;
	  lambdapos.z = (*sbm[k]).z[(int)index_lambda.x[]];
#endif

	  foreach_dimension() {
	    ppshift.x = 0.;
	    if (lambdacellpos.x < 0.) ppshift.x = -L0;
	    if (lambdacellpos.x > L0) ppshift.x = +L0;
	  }
	  
	  weight = reversed_weight (pp, weightcellpos, lambdacellpos, lambdapos, Delta, ppshift, 1., 1.);

	  sum.x += weight*DLM_w.x[]; sum.y += weight*DLM_w.y[]; sum.z += weight*DLM_w.z[];
	}
      }

      qu.x[] += sum.x;
      qu.y[] += sum.y;
      qu.z[] += sum.z; //+= here as one fluid cell can be affected by multiples particle's boundary multipliers
	
#if DLM_Moving_particle
      if (index_lambda.x[] > -1 && ((int)index_lambda.y[] == pp->pnum)) {
	lambdapos.x = (*sbm[k]).x[(int)index_lambda.x[]];
	lambdapos.y = (*sbm[k]).y[(int)index_lambda.x[]];
	lambdapos.z = 0.;
	
#if dimension == 3
	lambdapos.z = (*sbm[k]).z[(int)index_lambda.x[]];
#endif
#if TRANSLATION
	(*qU[k]).x += -DLM_w.x[];
	(*qU[k]).y += -DLM_w.y[];
	(*qU[k]).z += -DLM_w.z[];
#endif

#if ROTATION
	/* Modify temporerly the particle center position for periodic boundary condition */
	foreach_dimension()
	  (*gci[k]).center.x += DLM_periodic_shift.x[];
		
	(*qw[k]).x += -(DLM_w.z[]*(lambdapos.y - (*gci[k]).center.y) - DLM_w.y[]*(lambdapos.z - (*gci[k]).center.z)); // -(w_z*r_y - w_y*r_z)
	(*qw[k]).y += -(DLM_w.x[]*(lambdapos.z - (*gci[k]).center.z) - DLM_w.z[]*(lambdapos.x - (*gci[k]).center.x)); // -(w_x*r_z - w_z*r_x)
	(*qw[k]).z += -(DLM_w.y[]*(lambdapos.x - (*gci[k]).center.x) - DLM_w.x[]*(lambdapos.y - (*gci[k]).center.y)); // -(w_y*r_x - w_x*r_y)

	foreach_dimension()
	  (*gci[k]).center.x -= DLM_periodic_shift.x[];
	
#endif
      }
#endif
    }
  }
#endif
    
#if DLM_Moving_particle
#if _MPI
  for (int k = 0; k < NPARTICLES; k++) {
#if TRANSLATION
    mpi_all_reduce ((*qU[k]).x, MPI_DOUBLE, MPI_SUM);
    mpi_all_reduce ((*qU[k]).y, MPI_DOUBLE, MPI_SUM);
    mpi_all_reduce ((*qU[k]).z, MPI_DOUBLE, MPI_SUM);
#endif
#if ROTATION
    mpi_all_reduce ((*qw[k]).x, MPI_DOUBLE, MPI_SUM);
    mpi_all_reduce ((*qw[k]).y, MPI_DOUBLE, MPI_SUM);
    mpi_all_reduce ((*qw[k]).z, MPI_DOUBLE, MPI_SUM);
#endif
  }
#endif
#endif

    /* (2) Invert L*t = qu with L=rho*dV/dt*I */
    foreach() {
      tu.x[] = qu.x[]*dt/(rho_f*dlmfd_dv());
      tu.y[] = qu.y[]*dt/(rho_f*dlmfd_dv());
      tu.z[] = qu.z[]*dt/(rho_f*dlmfd_dv());
    }

#if DLM_Moving_particle
    for (int k = 0; k < NPARTICLES; k++) {
#if TRANSLATION
      /* (*tU[k]).x = ((*qU[k]).x*dt)/((1. - (rho_f/p[k].rho_s))*p[k].M); */
      /* (*tU[k]).y = ((*qU[k]).y*dt)/((1. - (rho_f/p[k].rho_s))*p[k].M); */
      /* (*tU[k]).z = ((*qU[k]).z*dt)/((1. - (rho_f/p[k].rho_s))*p[k].M); */

      /* Modification for the explicit added term  */
      (*tU[k]).x = ((*qU[k]).x*dt)/( p[k].DLMFD_couplingfactor * p[k].M );
      (*tU[k]).y = ((*qU[k]).y*dt)/( p[k].DLMFD_couplingfactor * p[k].M );
      (*tU[k]).z = ((*qU[k]).z*dt)/( p[k].DLMFD_couplingfactor * p[k].M );
#endif
#if ROTATION
      /* The inertia tensor is */
      /*  Ixx -Ixy -Ixz */
      /* -Iyx  Iyy -Iyz */
      /* -Izx -Izy  Izz */ 
      /* with */
      /* Ip[0] = Ixx */
      /* Ip[1] = Iyy */
      /* Ip[2] = Izz */
      /* Ip[3] = Ixy */
      /* Ip[4] = Ixz */
      /* Ip[5] = Iyz */

      /* Allocate temporly a matrix for Ip */
      double ** Imat = malloc(3 * sizeof(double*));
      for (int ii = 0; ii < 3; ii++) {
	Imat[ii] = malloc(3 * sizeof(double));
      }

      Imat[0][0] = p[k].Ip[0];
      Imat[1][1] = p[k].Ip[1];
      Imat[2][2] = p[k].Ip[2];

      Imat[0][1] = -p[k].Ip[3];
      Imat[0][2] = -p[k].Ip[4];
      Imat[1][2] = -p[k].Ip[5];
    
      Imat[1][0] = Imat[0][1]; 
      Imat[2][0] = Imat[0][2];
      Imat[2][1] = Imat[1][2];

      double Aw[3][3];
      for (int ll = 0; ll < 3; ll++) {
	for (int nn = 0; nn < 3; nn++) {
	  /* Aw[ll][nn] = (1. - rho_f/p[k].rho_s)*Imat[ll][nn]/dt; */
	  /* Modification for the explicit added term  */
	  Aw[ll][nn] = p[k].DLMFD_couplingfactor * Imat[ll][nn] / dt;
	}
      }

      /* Takes in Aw and return Aw^-1 as Imat */
      inverse3by3matrix(Aw, Imat);
        
      (*tw[k]).x = Imat[0][0]*(*qw[k]).x + Imat[0][1]*(*qw[k]).y + Imat[0][2]*(*qw[k]).z;
      (*tw[k]).y = Imat[1][0]*(*qw[k]).x + Imat[1][1]*(*qw[k]).y + Imat[1][2]*(*qw[k]).z;
      (*tw[k]).z = Imat[2][0]*(*qw[k]).x + Imat[2][1]*(*qw[k]).y + Imat[2][2]*(*qw[k]).z;

      for (int ii = 0; ii < 3; ii++) 
	free(Imat[ii]);
      free(Imat);
      
      /* (*tw[k]).x = ((*qw[k]).x*dt)/((1. - (rho_f/p[k].rho_s))*p[k].Ip[0]); */
      /* (*tw[k]).y = ((*qw[k]).y*dt)/((1. - (rho_f/p[k].rho_s))*p[k].Ip[0]); */
      /* (*tw[k]).z = ((*qw[k]).z*dt)/((1. - (rho_f/p[k].rho_s))*p[k].Ip[0]); */
#endif
    }
#endif
         
    /* (3) Compute residual y = M*t, the residual vector y = <alpha, tu-(tU+tw^GM)>_P(t) */
    /* Sign error in eq (A.9) in J.Eng. Math, 2011 */
    /* Written -M*t, should be +M*t */
    
#if debugInterior == 0
    /* Interior points: y = M*t  */
    /* Interior points: y = M_u*tu + M_U*tU + M_w*tw */
    /* So y = <alpha, tu>_P(t) - <alpha, tU>_P(t) - <alpha, tw^r_GM>_P(t) */
    for (int k = 0; k < NPARTICLES; k++) {
      foreach_cache((*Interior[k])) {
	if ((flagfield[]  < 1) && ((int)index_lambda.y[] == k)) {
	  DLM_v.x[] = tu.x[];
	  DLM_v.y[] = tu.y[];
	  DLM_v.z[] = tu.z[]; 
#if DLM_Moving_particle
#if TRANSLATION
	  DLM_v.x[] += -(*tU[k]).x;
	  DLM_v.y[] += -(*tU[k]).y;
	  DLM_v.z[] += -(*tU[k]).z;
#endif

#if ROTATION

	  /* Modify temporerly the particle center position for periodic boundary condition */
	  foreach_dimension()
	    (*gci[k]).center.x += DLM_periodic_shift.x[];
	
	  DLM_v.x[] += -((*tw[k]).y*(z - (*gci[k]).center.z) - (*tw[k]).z*(y - (*gci[k]).center.y)); // -(t_y*r_z - t_z*r_y)
	  DLM_v.y[] += -((*tw[k]).z*(x - (*gci[k]).center.x) - (*tw[k]).x*(z - (*gci[k]).center.z)); // -(t_z*r_x - t_x*r_z)
	  DLM_v.z[] += -((*tw[k]).x*(y - (*gci[k]).center.y) - (*tw[k]).y*(x - (*gci[k]).center.x)); // -(t_x*r_y - t_y*r_x)

	  foreach_dimension()
	    (*gci[k]).center.x -= DLM_periodic_shift.x[];
	
#endif
#endif
	}
      }
    }
#endif
    
#if debugBD == 0
    /* Boundary points: y = M*t */
    /* Boundary points: y = M_u*tu + M_U*tU + M_w*tw */
    /* So y = <alpha, tu>_P(t) - <alpha, tU>_P(t) - <alpha, tw^r_GM>_P(t) */
    boundary ((scalar*){tu});

    for (int k = 0; k < NPARTICLES; k++) {
      particle * pp = &p[k];
      
      foreach_cache((*Boundary[k])) {
	if (index_lambda.x[] > -1 && ((int)index_lambda.y[] == pp->pnum)) {
	  lambdacellpos.x = x; lambdacellpos.y = y;lambdacellpos.z = 0.;
	  lambdapos.x = (*sbm[k]).x[(int)index_lambda.x[]];
	  lambdapos.y = (*sbm[k]).y[(int)index_lambda.x[]];
	  lambdapos.z = 0.;
#if dimension == 3
	  lambdacellpos.z = z;
	  lambdapos.z = (*sbm[k]).z[(int)index_lambda.x[]];
#endif
	  sum.x = 0; sum.y = 0; sum.z = 0;
	  double testweight = 0.;

	  foreach_dimension()
	    ppshift.x = DLM_periodic_shift.x[];

	  foreach_neighbor() {
	    if (level == depth()) {
	      weightcellpos.x = x; weightcellpos.y = y; weightcellpos.z = 0.;
#if dimension == 3
	      weightcellpos.z = z;
#endif
	     
	      weight = reversed_weight (pp, weightcellpos, lambdacellpos, lambdapos, Delta, ppshift, 1., 0.);
	      testweight += weight;
	      	      	      
	      sum.x += weight*tu.x[];
	      sum.y += weight*tu.y[];
	      sum.z += weight*tu.z[];
	    }
	  }


	  if ((testweight < 1. - 0.000001) || (testweight > 1. + 0.000001))
		printf("testweight = %f\n",testweight);
	  
  
	  DLM_v.x[] = sum.x;
	  DLM_v.y[] = sum.y;
	  DLM_v.z[] = sum.z; 

#if DLM_Moving_particle
#if TRANSLATION

	  DLM_v.x[] +=  -(*tU[k]).x;
	  DLM_v.y[] +=  -(*tU[k]).y;
	  DLM_v.z[] +=  -(*tU[k]).z;
#endif
#if ROTATION

	  /* Modify temporerly the particle center position for periodic boundary condition */
	  foreach_dimension()
	    (*gci[k]).center.x += DLM_periodic_shift.x[];
	
	  DLM_v.x[] += -((*tw[k]).y*(lambdapos.z - (*gci[k]).center.z) - (*tw[k]).z*(lambdapos.y - (*gci[k]).center.y)); // -(t_y*r_z - t_z*r_y)
	  DLM_v.y[] += -((*tw[k]).z*(lambdapos.x - (*gci[k]).center.x) - (*tw[k]).x*(lambdapos.z - (*gci[k]).center.z)); // -(t_z*r_x - t_x*r_z)
	  DLM_v.z[] += -((*tw[k]).x*(lambdapos.y - (*gci[k]).center.y) - (*tw[k]).y*(lambdapos.x - (*gci[k]).center.x)); // -(t_x*r_y - t_y*r_x)

	  foreach_dimension()
	    (*gci[k]).center.x -= DLM_periodic_shift.x[];
	
#endif
#endif
	}
      }
    }
#endif
     
    /* (4) Compute alpha = r^(k-1)*r^(k-1) / DLM_w.y */
    DLM_nr2_km1 = DLM_nr2;
    DLM_wv = 0.;

    foreach() {
      DLM_wv += DLM_w.x[]*DLM_v.x[] + DLM_w.y[]*DLM_v.y[] + DLM_w.z[]*DLM_v.z[];
    }
    
#if _MPI
    mpi_all_reduce (DLM_wv, MPI_DOUBLE, MPI_SUM);
#endif

    DLM_alpha = DLM_nr2_km1/DLM_wv;
        
    /* (5) Update lambda = lambda - alpha*DLM_w and r = r - alpha*y */
    foreach() {
      DLM_lambda.x[] -= DLM_alpha*DLM_w.x[];
      DLM_lambda.y[] -= DLM_alpha*DLM_w.y[];
      DLM_lambda.z[] -= DLM_alpha*DLM_w.z[];

      DLM_r.x[] -= DLM_alpha*DLM_v.x[];
      DLM_r.y[] -= DLM_alpha*DLM_v.y[];
      DLM_r.z[] -= DLM_alpha*DLM_v.z[];

      /* (7) Update u = u + alpha*t */
      
      u.x[] += DLM_alpha*tu.x[];
      u.y[] += DLM_alpha*tu.y[];
      u.z[] += DLM_alpha*tu.z[];
    }

#if DLM_Moving_particle
    for (int k = 0; k < NPARTICLES; k++) {
#if TRANSLATION
      (*U[k]).x += DLM_alpha*(*tU[k]).x;
      (*U[k]).y += DLM_alpha*(*tU[k]).y;
      (*U[k]).z += DLM_alpha*(*tU[k]).z;
#endif
#if ROTATION
      (*w[k]).x += DLM_alpha*(*tw[k]).x;
      (*w[k]).y += DLM_alpha*(*tw[k]).y;
      (*w[k]).z += DLM_alpha*(*tw[k]).z;
#endif
    }
#endif
   
    /* (8) Compute beta = nr2^k / nr2^(k-1) */
    DLM_nr2 = 0.;
    
    /* (6) compute norme ||r^k||^2 =  DLM_nr2 */
    foreach() {
      DLM_nr2 += sq(DLM_r.x[]) + sq(DLM_r.y[]) + sq(DLM_r.z[]);
    }
    
#if _MPI
    mpi_all_reduce (DLM_nr2, MPI_DOUBLE, MPI_SUM);
#endif
    
    DLM_beta = DLM_nr2/DLM_nr2_km1;

    /* (9) Update DLM_w = r + beta*DLM_w */
    /* Sign error in eq (A.15) of J.Eng. Math, 2011 */
    /* Written -beta*DLM_w, should be +beta*DLM_w */
    foreach() {
      DLM_w.x[] = DLM_r.x[] + DLM_beta*DLM_w.x[];
      DLM_w.y[] = DLM_r.y[] + DLM_beta*DLM_w.y[];
      DLM_w.z[] = DLM_r.z[] + DLM_beta*DLM_w.z[];
    }
    /* // Output nr */
    /* fprintf (stderr,"Iteration %d   ||u-u_imposed|| = %12.10e\n",ki,sqrt(DLM_nr2)); */
  }

  if (i == 0)
    fprintf (converge, "# Time Iteration \t Uzawa Iterations \t ||u-u_imposed||\n");
  fprintf (converge,"%d \t %d \t %12.10e\n", i, ki, sqrt(DLM_nr2));
  fflush(converge);

#if DLM_Moving_particle

#if TRANSLATION
  fprintf (stderr,"DLMFD: particle-0's velocity on thread %d is (%20.18f, %20.18f, %20.18f)\n",pid(), (*U[0]).x, (*U[0]).y, (*U[0]).z);
  coord cc = (*gci[0]).center;
  fprintf (stderr,"DLMFD: particle-0's position on thread %d is (%20.18f, %20.18f, %20.18f)\n",pid(), cc.x, cc.y, cc.z);
#endif
  
#if ROTATION
  fprintf (stderr,"DLMFD: particle-0's rotation on thread %d is (%20.18f, %20.18f, %20.18f)\n",pid(), (*w[0]).x, (*w[0]).y, (*w[0]).z);
#endif
  
#endif
  
  /* Compute the explicit term here */
#if DLM_alpha_coupling 

  boundary ((scalar*) {DLM_lambda});
  coord nn = {0., 0., 0.};

  foreach_dimension() {
    nn.x = 0.;
    norm n = normf (DLM_lambda.x);
    nn.x = n.max;
  }
    fprintf (stderr, "# Lambda explicit term: max of each component (%g, %g, %g)\n",  nn.x, nn.y, nn.z); // log

    for (int k = 0; k < NPARTICLES; k++) {
# if debugBD == 0
    particle * pp = &p[k];
    
    foreach_cache ((*Boundary[k])) {
      sum.x = 0.; sum.y = 0.; sum.z = 0.;
      
      weightcellpos.x = x; weightcellpos.y = y;
#if dimension == 3
      weightcellpos.z = z;
#endif
      foreach_neighbor() {
	if (((int)index_lambda.x[] > -1) && (level == depth()) && is_leaf(cell) && ((int)index_lambda.y[]) == k) {
	  lambdacellpos.x = x;
	  lambdacellpos.y = y;
	  lambdapos.x = (*sbm[k]).x[(int)index_lambda.x[]];
	  lambdapos.y = (*sbm[k]).y[(int)index_lambda.x[]];
	
#if dimension == 3
	  lambdacellpos.z = z;
	  lambdapos.z = (*sbm[k]).z[(int)index_lambda.x[]];
#endif

	  foreach_dimension()
	    ppshift.x = DLM_periodic_shift.x[];

	  weight = reversed_weight (pp, weightcellpos, lambdacellpos, lambdapos, Delta, ppshift,1.,1.);
	  
	  foreach_dimension()
	    sum.x += weight*DLM_lambda.x[];
	}
      }

      foreach_dimension() { 
	DLM_explicit.x[] += sum.x;
      }
    }
#endif
    
#if debugInterior == 0
#if POISEUILLE
  /* For the poiseuille flow, the explicit forcing term is too sharp and needs to be "smoothed" */
#define smooth 0.5
#else
#define smooth 1.
#endif
    foreach_cache ((*Interior[k])) {
      if ((flagfield[]  < 1) && ((int)index_lambda.y[] == k))
	foreach_dimension() {
	  DLM_explicit.x[] = smooth*DLM_lambda.x[];
	}
    }
#endif
  }

  boundary((scalar *) {DLM_explicit});
#endif

  /* Added by Arthur Ghigo, October 10, 2019 */
  boundary ((scalar *) {u});
  
  /* Timers and Timings */
  /* give 1 as number of cell so that timer_timing does not compute the total number of cells*/
  /* the total number of cell is used to compute the speed of the solver which does not make sens for a single step here, we will compute it globally at the end of the run*/
  timing dlmfd_timing = timer_timing (dlmfd_timer, 1, 1, mpitimings);

  dlmfd_globaltiming.cpu += dlmfd_timing.cpu;
  dlmfd_globaltiming.real += dlmfd_timing.real;
  dlmfd_globaltiming.min += dlmfd_timing.min;
  dlmfd_globaltiming.avg += dlmfd_timing.avg;
  dlmfd_globaltiming.max += dlmfd_timing.max;
}
