/**This file aims to introduce functions for surface in the multilayer solver. 
It is composed of four parts. The first for basic characterisation functions (height, area, tangents), the second for integral surface tension implementation, the third for surfactant adsorption on the interface and surface advection-diffusion. And the last one, for solvent evaporation. 

In progress : next step is to add the possibility of a second dimension and axisymetry. Then it will remains, among other, to add projection functions and to improve area interpolation with circles.
*/

/**
# Basic functions for surface characterisation
A surface excess concentration field is created. To ensure mass conservation and independance toward the area $\mathcal{A}$, the surface field is:
$$
surface = M = \Gamma * \mathcal{A}
$$
*/
scalar surface[]; 

/**
## Total height function
This function calculate the total height function.
$$
H(x) = \sum_{layers} h(x)
$$
*/
scalar H[];
void heights (scalar* hl, scalar H)
{
  foreach(){
    H[] = zb[];
    for (scalar h in hl)
      H[] += h[];
  }
  boundary({H});
}

#define contact_angle(theta, Delta, layers) \
  val(_s) + Delta/tan(theta)/layers
/**
Nota Bene : the division by "layers" could be removed. It is added in case the boundary over the h of hl is needed.*/

/**
## Surface area
This function calculate the surface area on the faces and in a cell compare to the reference lenghth $\Delta$. This function will be later improved with a better approximation using curvature and circles.

$$
\mathcal{A} = \frac{1}{\Delta} * \sqrt{\Delta^2 + \delta h^2} = \sqrt{1 + \frac{\delta h^2}{\Delta^2}}
$$
*/

double area_f (Point point, scalar H)
{
  double hx = (H[] - H[-1])/Delta;
  return sqrt(1. + sq(hx));
}

double area_c (Point point, scalar H)
{
  double area_f0 = area_f(point, H);
  double area_f1 = area_f(neighborp(1), H);
  return (area_f0 + area_f1)*0.5;
}


/**
## Tangent
These functions compute the tangents of the interface on the faces and in a cell.
A specific structure is created to add the vertical component (which is out the grid dimension because of the layers).
*/

typedef struct Vect Vect;
struct Vect {
  coord xy;
  double v;
};

struct Vect normalize_Q (Vect T) {
  double nT = sq(T.v);
  foreach_dimension()
    nT += sq(T.xy.x);
  nT = sqrt(nT);
  foreach_dimension()
    T.xy.x = (nT != 0. ? T.xy.x/nT : 0.);
  T.v = (nT != 0. ? T.v/nT : 0.);
  return T;
}

struct Vect tangent_left (Point point, scalar H)
{
  Vect T;
  T.xy.x = 1.;
  T.v = (H[] - H[-1])/Delta;
  T = normalize_Q(T);
  return T;
}

struct Vect tangent_center (Point point, scalar H)
{
  Vect T0 = tangent_left (point, H);
  Vect T1 = tangent_left (neighborp(1), H);
  Vect T;
  foreach_dimension()
    T.xy.x = (T0.xy.x + T1.xy.x);
  T.v = (T0.v + T1.v);
  T = normalize_Q(T);
  return T;
}

/**
## Curvature
To compute the curvature, we estimate the derivatives of the height
functions. We then compute the curvature as
$$
\kappa = \frac{h_{xx}}{(1 + h_x^2)^{3/2}}
$$
in two dimensions. */

static double kappa_comp (Point point, scalar H)
{
  double hx = (H[1] - H[-1])/(2.*Delta);
  double hxx = (H[1] + H[-1] - 2.*H[])/sq(Delta);
  return hxx/pow(1. + sq(hx), 3/2.);
}

static void curvature (scalar H, scalar kappa)
{
  foreach()
    kappa[] = kappa_comp(point, H);
}


/**
#Integral formulation for surface tension
The surface tension field gam is introduced. It must be initialized to add surface tension. It can be a function, especially of surfactant surface concentration.
*/
scalar gam[];
   
/**
## Stability condition

The surface tension scheme is time-explicit so the maximum timestep is
the oscillation period of the smallest capillary wave.
$$
T = \sqrt{\frac{\rho_{m}\Delta_{min}^3}{\pi\sigma}}
$$
with $\rho_m=(\rho_1+\rho_2)/2.$ and $\rho_1$, $\rho_2$ the densities
on either side of the interface. 
*/

event stability (i++) {

  double rhom = 1.;

  /**
  We then consider each VOF interface with an associated value of
  $\sigma$ different from zero and set the maximum timestep. 
  */

  foreach(){
    if (gam[]>0) {  
      double dt = sqrt (rhom*cube(Delta)/(pi*gam[]));
      if (dt < dtmax)
	      dtmax = dt;
    }         
  }
}

/**
## Surface tension acceleration

The calculation of the acceleration is done by this event, overloaded
from [its definition](layered/hydro.h) in the multilayer Navier--Stokes solver. 
*/

event pressure (i++)
{   
  /**
  The total height H is actualised and the potential S is computed from the tangents :
  $$
  S = \frac{\gamma T}{\Delta}
  $$
  */
  heights(hl, H);
  scalar Sx[], Sz[];

  foreach(){
    Vect T = tangent_center(point, H);
    Sx[] = gam[]/Delta*T.xy.x;
  }
  boundary ({Sx});
      
  foreach_face(){
    Vect T = tangent_left(point, H);
    double gam_m = gam[];
    Sz[] = gam_m/Delta*T.v;
  }
  
  /**
  The horyzontal gradient of the tensor S is computed and added respectively to the  acceleration and the face velocities of the upper layer. 
  $$
  a^{**}_{f,l} = a^{\star}_{f,l} - \nabla S_x =  a^{\star}_l - \frac{S_{x, x - 1/2} - S_{x, x + 1/2}} {\Delta*\mathcal{A}}
  $$
  $$
  u^{**}_{f,l} = u^{\star}_{f,l} - \Delta t \frac{S_{x, x - 1/2} - S_{x, x + 1/2}} {\Delta*\mathcal{A}}
  $$
  $$
  w^{**}_l = w^{\star}_l - \frac{\Delta t}{h^{n}_l} \frac{S_{y, x - 1/2} - S_{y, x + 1/2}}{\Delta * \mathcal{A}}
  $$

  */  
  face vector uf = ufl[nl-1];
  face vector a = al[nl-1];
  foreach_face() {
    double ax = (Sx[] - Sx[-1])/area_f(point, H);
    uf.x[] += dt * ax;
	  a.x[] += ax;
  }
  boundary ((scalar *)ufl);
  boundary ((scalar *)al);

  scalar w = wl[nl-1];
  scalar h = hl[nl-1];
  foreach() {
    double az = (Sz[1] - Sz[])/area_c(point, H);
    w[] += dt*az/h[];
  }
  boundary (wl);
}



/**
# Surfactants in Multilayer

## Solute adsorption
Surfactant adsorption follows an adsorption/desorption kinetic, which can be modelled using Langmuir isotherm :
$$
\frac{d \Gamma}{dt} = \omega\, (\sigma\, c\, (1-\Gamma) - \Gamma)
\quad \textrm{then} \quad
\frac{dM}{dt} = \omega\, (\sigma\, c\, (\mathcal{A}-M) - M)
$$

The parameters $\omega$ and $\sigma$ are introduced. The adsorption rate $\omega$ must be initialized with a non-zero value to add adsorption. The adsorption term field is also declared. 
*/

double omega = 0.;
double sigma = 1.;
scalar d_ads_t[];

/**
This function calculate the adsorption term.
It will be later improved by adding other adsorption formulation, especially Frumkin adsorption. + implicit !
*/
scalar adsorption_term (scalar H, scalar c, scalar surface)
{
  scalar ads[];
  foreach(){
    ads[] = omega * (sigma * c[] * (area_c(point, H) - surface[]) - surface[]);
  }
  return ads;
}

/**
To ensure stability, the timestep must be maximized. The event set_dtmax is called to ensure that:
$$
dt<\frac{1}{5 * \omega}
\quad\textrm{and}\quad
dt<\frac{c * h}{10 \, \partial_t \Gamma}
$$
*/
event stability (i++, last) 
{
  if (omega > 0) {
    double dt_max = 1/(5.*omega);
    
    heights(hl, H);
    scalar h = hl[nl-1];
    scalar c = cl[nl-1];
    d_ads_t = adsorption_term(H, c, surface);

    foreach()
      if (d_ads_t[]>0.){
        dt_max = min(dt_max, c[]*h[]/d_ads_t[]/10.);
      }
    dtmax = min(dtmax, dt_max);
  }
}

/**
The event is called during the event viscous_term. Here, the adsorption term d_ads_t is computed a second time... This should be modified once the order decided.
*/
event viscous_term (i++, last)
{
  heights(hl, H);
  scalar h = hl[nl-1];
  scalar c = cl[nl-1];
  d_ads_t = adsorption_term(H, c, surface);
  foreach(){
    surface[] += d_ads_t[]*dt;
    c[] -= d_ads_t[]*dt/h[];
  }
  boundary({surface, c})
}

/**
## Advection of surfactant over the surface
*/
event advection_term (i++,last)
{
  vector uf = ufl[nl-1];
  scalar flux[];
  foreach_face()
    flux[] = uf.x[] * (uf.x[]>0 ? surface[] : surface[-1]);
  boundary({flux});
  foreach()
    surface[] += dt*(flux[1] - flux[]);
  boundary({surface});
}
/**
## Diffusion of surfactant over the surface
The function is called during the event diffusive_term, just after the viscous_term event.
*/

event diffusive_term (i++,last)
{  
  if (D_hor > 0.){
    scalar h_one[], area[];
    foreach() {
      h_one[] = 1.;
      area[] = area_c (point, H);
      surface[] = surface[]/area[];
    }
    boundary ({surface, h_one});
    
    horizontal_diffusion (surface, dt, h_one, area, true);
    
    foreach()
      surface[] = surface[]*area[];  
    boundary ({surface});
  }
}

/**
# Evaporation
Evaporation is controlled by the velocity $v_\mathcal{E}$. It must be initialized to add evaporation and can be a function. 
*/
scalar v_e[];

/**
This function calculate the evaporation term.
$$
\frac{dh}{dt} = v_\mathcal{E}\,\mathcal{A}
$$
As the solute does not evaporate, the concentration field must be corrected such as: 
$$
c_{old} * h_{old} = c * h
$$
*/

void evaporation (scalar h, scalar H, scalar c, double dt)
{
  foreach() {
    double h_old = h[];
      h[] -= v_e[] * area_c(point, H) * dt;
      c[] = c[]*h_old/h[];
  }
}

/**
Evaporation function is called during an event evaporation_term.
*/
  
event evaporation_term (i++,last)
{
  heights(hl, H);
  evaporation(hl[nl-1], H, cl[nl-1], dt);
}
