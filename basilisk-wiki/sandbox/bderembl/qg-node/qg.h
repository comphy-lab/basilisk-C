/**
Quasi-Geostrophic code for Basilisk

This version is for a single layer but will eventually become a multi layer model

Solves


$$
\partial_t q + \nabla (u q ) + \beta v =  - \frac{h_{EK}*f0}{h_n} q  + nu \nabla^2 q + F
$$
with $q$ the vorticity
$$
q = \nabla^2 \psi
$$
$\psi$ the stream function
$$
u = -\frac{\partial \psi}{\partial y}
$$
and 
$$
v = \frac{\partial \psi}{\partial x}
$$

$h_{EK}$ is the thickness of the bottom Ekman layer, $H$ the thickness of the
lowermost layer and $f0$ is the Coriolis parameter. $\nu$ is the harmonic
viscosity.

The forcing $F$ is for now hard coded to
$$
F = \frac{3 \tau_0 \pi}{2 h_1 L} \sin(2*\pi*y/L0)*sin(\pi*y/L0);
$$
which correspond to a wind stress curl that will force a double gyre
configuration (cyclonic at the north, anticyclonic at the south).

The advective term $\nabla (uq)$ is coded with the Arakawa Jacobian which
conserves energy and enstrophy (see Cushman Roisin and Beckers 2011)


All parameters can be changed via a "params.in" file which should be in the
location as the executable. exemple of params.in:

N  = 64
L0 = 100

tau0 = 1e-3
nu = 5
beta = 0.5
dh   = [1.0]
sbc = 0

# timestepping
DT    = 5.e-2
tend  = 200.
dtout = 10
CFL   = 0.2
TOLERANCE = 1e-5


sbc is a parametre that controls the type of boundary conditions:
sbc = -1  -> periodic BC (experimental)
sbc = 0   -> free slip BC
sbc = big number (>10) -> no slip BC

dtout is the frequency for writing outputs
tend the total time of the simulation

For this single layer version: $q$ and $zeta$ are the same field (not true for
multiple layers)


This version uses a prototype poisson solver for nodal fields written by Antoon
Van Hoft.  Hence all fields are defined at cell vertices. (still experimental)

*/

#include "predictor-corrector.h"
#include "nodal-poisson.h"

/**
   User defined constants
 */

double f0 = 1.0;
double beta = 0.;
double hEkb = 0.;
double tau0 = 0.;
double nu = 0.;
double sbc = 0.;
double tend = 100.;
double dtout = 1.; 
double dh[1] = {1.};

char dpath[80]; // name of output dir
/** 
    Fields
*/

vertex scalar psi[]; 
vertex scalar q[]; 
vertex scalar zeta[]; 
vertex scalar q_forcing[]; 


/**
   Internal variables
 */
double bc_fac = 0;
scalar * evolving = NULL;
mgstats mgpsi;


/**
   Energy diagnostics
*/

double dtdiag = -1; // non zero

/**
   Define useful operators:

   laplacian(psi) = \nabla^2 psi

   jacobian(psi, q) = dpsi/dx.dq/dy - dpsi/dy.dq/dx

   beta_effect(psi) = beta dpsi/dx = beta v
*/

#define laplacian(p) (p[1] + p[-1] + p[0,1] + p[0,-1] - 4*p[])/(sq(Delta))

#define jacobian(p,q) ((( p[1, 0 ]-p[-1, 0])*(q[0, 1]-q[ 0 ,-1]) \
                        +(p[0 ,-1]-p[ 0 ,1])*(q[1, 0]-q[-1, 0 ]) \
                        + p[1, 0 ]*( q[1,1 ] - q[1,-1 ])         \
                        - p[-1, 0]*( q[-1,1] - q[-1,-1])         \
                        - p[ 0 ,1]*( q[1,1 ] - q[-1,1 ])         \
                        + p[0 ,-1]*( q[1,-1] - q[-1,-1])         \
                        + q[ 0 ,1]*( p[1,1 ] - p[-1,1 ])         \
                        - q[0 ,-1]*( p[1,-1] - p[-1,-1])         \
                        - q[1, 0 ]*( p[1,1 ] - p[1,-1 ])         \
                        + q[-1, 0]*( p[-1,1] - p[-1,-1]))        \
                       /(12.*Delta*Delta))

#define beta_effect(p) (beta*(p[1] - p[-1])/(2*Delta))

/** 
    Dynamical core
*/


void set_bc()
{
  // derived const
  bc_fac = sbc/((0.5*sbc + 1)); 

  psi[left]   = 0;
  psi[right]  = 0;
  psi[bottom] = 0;
  psi[top]    = 0;

  zeta[left]   = bc_fac/sq(Delta)*(psi[1] - psi[]);  
  zeta[right]  = bc_fac/sq(Delta)*(psi[] - psi[1]);  
  zeta[bottom] = bc_fac/sq(Delta)*(psi[0,1] - psi[]);
  zeta[top]    = bc_fac/sq(Delta)*(psi[] - psi[0,1]);

  q[left]   = bc_fac/sq(Delta)*(psi[1] - psi[]);   
  q[right]  = bc_fac/sq(Delta)*(psi[] - psi[1]);  
  q[bottom] = bc_fac/sq(Delta)*(psi[0,1] - psi[]);
  q[top]    = bc_fac/sq(Delta)*(psi[] - psi[0,1]);

}
/**
   Invert the poisson equation $\nabla^2 \psi = q$
*/

trace
void invertq(vertex scalar psi, vertex scalar q)
{

  mgpsi = vpoisson(psi, q);
  boundary({psi});
  // need to reset the values of the BC (has to do with vertices?)
  set_bc();
  boundary({q});
}

trace
void comp_del2(scalar psi, scalar zeta, double add, double fac)
{
  foreach_inner_vertex()
      zeta[] = add*zeta[] + fac*laplacian(psi);
  
  boundary({zeta});
}

trace
void comp_q(scalar psi, scalar q)
{
  comp_del2  (psi, q, 0., 1.);
  /* comp_stretch(psi, q, 1., 1.); */
  boundary({q});
}
/**
   Advection
*/

trace
double advection_pv(scalar zeta, scalar q, scalar psi, scalar dqdt, double dtmax)
{
  foreach_inner_vertex()
    dqdt[] += -jacobian(psi, zeta) - beta_effect(psi);

  // compute dtmax (ajusted from timestep.h)
  static double previous = 0.;
  dtmax /= CFL;
  foreach_face(reduction(min:dtmax)){
      double u = (psi[0,1] - psi[])/Delta;
      if (u != 0.) {
        double dt = Delta/fabs(u);
        if (dt < dtmax) dtmax = dt;
      }
  }

  dtmax *= CFL;
  if (dtmax > previous)
    dtmax = (previous + 0.1*dtmax)/1.1;
  previous = dtmax;
  return dtmax;
}

trace
void dissip  (scalar zeta, scalar dqdt)
{
  comp_del2(zeta, dqdt, 1., nu);
}

/**
   Bottom  Ekman Friction
*/

trace
void ekman_friction  (scalar zeta, scalar dqdt)
{
  foreach_inner_vertex()
    dqdt[] -= hEkb*f0/(2*dh[0])*zeta[];
}

/**
   surface forcing
*/
trace
void surface_forcing  (scalar dqdt)
{
  foreach_inner_vertex()
    dqdt[] -= q_forcing[];
}


static void advance_qg (scalar * output, scalar * input,
                        scalar * updates, double dt)
{

  vertex scalar qi = input[0];
  vertex scalar qo = output[0];
  vertex scalar dq = updates[0];

  foreach_inner_vertex()
    qo[] = qi[] + dq[]*dt;
//  boundary(output);
}


double update_qg (scalar * evolving, scalar * updates, double dtmax)
{

  // TODO check order loop
  foreach_vertex()
    for (scalar s in updates)
        s[] = 0.;

  vertex scalar q = evolving[0];
  vertex scalar dqdt = updates[0];

  invertq(psi, q);
  comp_del2(psi, zeta, 0., 1.0);
  dtmax = advection_pv(zeta, q, psi, dqdt, dtmax);

  dissip(zeta, dqdt);
  ekman_friction(zeta, dqdt);
  surface_forcing(dqdt);
  
  return dtmax;
}


/**
   Energy diagnostics
*/

event write_1d_diag (t=0; t <= tend+1e-10; t += dtdiag){
  if (dtdiag > 0){
    FILE * fp;
    char name[80];
    sprintf (name,"%sdiag_1d.dat", dpath);

    if (i == 0){
      // create file
      if (pid() == 0) {
        if ((fp = fopen(name, "a"))) {
          fprintf(fp, "# time, ke, dissipation, forcing\n");
          fclose(fp);
        }
      }
    }
    
/**
   Using a foreach loop because foreach_vertex will count subdomain boundary
   points twice in MPI. Make sure the field is zero at the boundary
 */
    double ke = 0;
    double d_ke = 0;
    double f_ke = 0;

    foreach(reduction(+:ke) reduction(+:d_ke) reduction(+:f_ke)){
      ke -= 0.5*psi[]*laplacian(psi)*sq(Delta);
      d_ke -= nu*psi[]*laplacian(q)*sq(Delta);
      f_ke -= psi[]*q_forcing[]*sq(Delta);
    }

    if (pid() == 0) {
      if ((fp = fopen(name, "a"))) {
        fprintf(fp, "%e, %e, %e, %e\n", t, ke, d_ke, f_ke);
        fclose(fp);
      }
    }
  } // dtdiag > 0
}



void set_vars()
{

  psi.restriction = restriction_vert;
  psi.prolongation = refine_vert;

  zeta.restriction = restriction_vert;
  zeta.prolongation = refine_vert;

  q.restriction = restriction_vert;
  q.prolongation = refine_vert;


  reset ({q, psi, zeta}, 0.);

//  if (nl > 1)
//    dhc = malloc ((nl-1)*sizeof(double));


  /**
     We overload the default 'advance' and 'update' functions of the
     predictor-corrector scheme and (TODO) setup the prolongation and
     restriction methods on trees. */

  evolving = list_copy({q});
  advance = advance_qg;
  update = update_qg;
  fprintf(stdout,"ok\n");
}


event defaults (i = 0){
  set_bc();
  set_vars();
}

/**
The event below will happen after all the other initial events to take
into account user-defined field initialisations. */

void set_const() {
  comp_q(psi,q); // last part of init: invert PV


  // set surface forcing
  reset ({q_forcing}, 0.);

  foreach_vertex()
    q_forcing[] = tau0/dh[0]*3/2*pi/L0*sin(2*pi*y/L0)*sin(pi*y/L0);

  boundary (all);
}

event init (i = 0) {  
  set_const(); 
}

/**
## Cleanup */

void trash_vars(){
//  delete ({q, psi, zeta});
  free(evolving);
//  free(dhc);

}

event cleanup (t = end, last)
{
  trash_vars();
}
