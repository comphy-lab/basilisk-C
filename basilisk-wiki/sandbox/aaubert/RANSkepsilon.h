/**
   Implementation of the standard k-epsilon model */

/**
The expression for the turbulent viscosity is 
$$
\mu_t=\rho C_\mu \frac{k^2}{\epsilon}
$$
with $k$ and $\epsilon$ following the equation (using Einstein notation)
$$
\frac{\partial k}{\partial t}+u_j\frac{\partial k}{\partial x_j}=\frac{\partial}{\partial x_j}\left(\left(\mu+\frac{\mu}{\sigma_k}\right)\frac{\partial k}{\partial x_j}\right)+P_k-\rho \epsilon
$$
$$
\frac{\partial \epsilon}{\partial t}+u_j\frac{\partial \epsilon}{\partial x_j}=\frac{\partial}{\partial x_j}\left(\left(\mu+\frac{\mu}{\sigma_\epsilon}\right)\frac{\partial \epsilon}{\partial x_j}\right)+C_{1\epsilon}\frac{\epsilon}{k}P_k-C_{2\epsilon}\rho\frac{\epsilon^2}{k}
$$

where $P_k=\mu_tS^2,~S=\sqrt{2S_{ij}S_{ij}}$  with $\bold{S}=\frac{1}{2}\left(\nabla \bold{u}+\nabla\bold{u}^t\right)$ the deformation tensor

At a wall, we have $k=0$ and $\epsilon=\infty$ so we use a wall model to obtain the value for $k$ and $\epsilon$. Given the isotropy of the mesh, we will never have a resolved mesh close to the wall (meaning the minimum $y^+$ at the wall will be higher than $30$) so we use an high Reynolds model. We choose the one from [Ruffin and al](https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=20582b7bb5665aed87a3180b19af8ec3a2410bed). It gives the value at the first point above the wall
$$
k=\frac{u_\tau^2}{\sqrt{c_\mu}}
$$
$$\mu_t=\kappa e^{-\kappa\bar{B}}\left(e^{\kappa u^+}-1-\kappa u^+-\frac{\left(\kappa u^+\right)^2}{2}\right)
$$

With this we can obtain the value of $\epsilon$ close to a wall
$$
\epsilon=\rho C_\mu \frac{k^2}{\mu_t}
$$

Close to the solid, we have
$$
\left(\mu+\mu_t\right)\frac{\partial u_t}{\partial y}=cst=\rho u_\tau^2
$$
To obtain the correct flux of velocity, we linearise the velocity profile and so we have a slip velocity in the solid that reads
$$
u_t=u_{t,IP}-d_{IP}\left.\frac{\partial u_t}{\partial y}\right|_{IP}
$$
with $\left.\frac{\partial u_t}{\partial y}\right|_{IP}=\rho\frac{u_\tau ^2}{\mu+\mu_t}$ and the image point (IP) located at a distance $d_{IP}$ from the wall.

The value for the constants is
$$
c_{\mu}=0.09, \sigma_k=1, \sigma_\epsilon=1.3, C_{1\epsilon}=1.44, C_{2\epsilon}=1.92, \kappa=0.41, \bar{B}=5.033
$$


*/
#include "bcg.h"   //needed for the advection part
#include "poisson.h" //needed for the diffusion part


#if !defined(WALL_FUNCTION)      //use non linearize version for the wall function if not state otherwise
#define WALL_FUNCTION 0
#endif



scalar k[];   
scalar epsilon[];

scalar mutke[];  //turbulent viscosity
face vector mut[];    //total viscosity

face vector accel[];

scalar omega[];

scalar sourcek[];   //source term
scalar sourceepsilon[];
scalar production[];

scalar molvis[];  //molecular viscosity

double Deltamin, d_IP;   //distance for the image point

double leading; //position of the leading edge
double trailing;  //position of the trailing edge

double C_1, C_2, C_mu, sigma_k, sigma_epsilon;  //constant for the k-epsilon model


event defaults(i=0) {
  mu=mut;
  a=accel;    //fixme : should not be necessary : create a bug accessing value of a constant scalar if removed
  leading=X0;
  trailing=X0+L0;

  C_1=1.44;
  C_2=1.92;
  C_mu=0.09;
  sigma_k=1.0;
  sigma_epsilon=1.3;
  
#if TREE
#if EMBED
  foreach_dimension();
  for (scalar s in {k,epsilon,sourcek,sourceepsilon}) {
    s.restriction = restriction_embed_linear;
    s.refine = s.prolongation = refine_embed_linear;
    s.depends = list_add (s.depends, cs);
  }
#endif // EMBED
#endif // TREE
  
  foreach() {
    dimensional (k[] == Delta*Delta/dt/dt);
    dimensional (epsilon[]= Delta*Delta/dt/dt/dt);
    sourcek[]=0.;
    sourceepsilon[]=0.;
    production[]=0.;
  }
}

#if EMBED

#undef center_gradient_x
#define center_gradient_x(a) (fs.x[] && fs.x[1] ? (a[1] - a[-1])/(2.*Delta) : \
			     fs.x[1] ? (a[1] - a[])/Delta :		    \
			     fs.x[]  ? (a[] - a[-1])/Delta : 0.)

#undef center_gradient_y
#define center_gradient_y(a) (fs.y[] && fs.y[0,1] ? (a[0,1] - a[0,-1])/(2.*Delta) : \
			      fs.y[0,1] ? (a[0,1] - a[])/Delta :        \
			      fs.y[]  ? (a[] - a[0,-1])/Delta : 0.)

#undef center_gradient_z
#define center_gradient_z(a) (fs.z[] && fs.z[0,0,1] ? (a[0,0,1] - a[0,0,-1])/(2.*Delta) : \
			      fs.z[0,0,1] ? (a[0,0,1] - a[])/Delta :	\
			      fs.z[]  ? (a[] - a[0,0,-1])/Delta : 0.)

#else

#undef center_gradient_x
#define center_gradient_x(a) ((a[1] - a[-1])/(2.*Delta))
			     

#undef center_gradient_y
#define center_gradient_y(a) ((a[0,1] - a[0,-1])/(2.*Delta))

#undef center_gradient_z
#define center_gradient_z(a) ((a[0,0,1] - a[0,0,-1])/(2.*Delta))

#endif


/**
   We first advect the viscosity */

void advection_source() {
  foreach_face() {
    uf.x[]=fm.x[]*face_value(u.x,0);
  }

  scalar sourcekadv[];
  scalar sourceepsilonadv[];
  foreach() {
    sourcekadv[]=sourcek[];
    sourceepsilonadv[]=sourceepsilon[];
  }

  advection((scalar *) {k,epsilon},uf,dt,(scalar *) {sourcekadv,sourceepsilonadv});

}


/**
   We need to add a correction to the viscosity computed to account for the term in the right hand side */

void correction_k_epsilon(double dt) {
  foreach() {
    k[]+=dt*sourcek[];
    epsilon[]+=dt*sourceepsilon[];
  }
}


event reaction_diffusion(i++) {


  advection_source();

  scalar mutke2[];
  foreach() {
    mutke2[]=rho[]*C_mu*sq(k[])/epsilon[];
  }
    
  face vector dissk[];
  face vector dissepsilon[];
  foreach_face() {
    dissk.x[]=fm.x[]*face_value(rho,0)*face_value(molvis,0)+fm.x[]*face_value(mutke,0)/sigma_k;    dissepsilon.x[]=fm.x[]*face_value(rho,0)*face_value(molvis,0)+fm.x[]*face_value(mutke,0)/sigma_epsilon;
  }

  //Implicit viscosity
  scalar lambda[];
  scalar bk[];
  scalar bepsilon[];
  foreach() {
    lambda[]=-rho[]/dt;
    bk[]=-k[]*rho[]/dt;
    bepsilon[]=-epsilon[]*rho[]/dt;
  }
  poisson (k,bk,dissk,lambda);
  poisson (epsilon,bepsilon,dissepsilon,lambda);
    

  trash({sourcek,sourceepsilon,production});
    
  double chi;

  double sum;
  double Pk;
  double d;              //distance to the nearest wall
  
  foreach() {
      sum=0.;
      foreach_dimension() {
        sum+=center_gradient_y(u.x)*(center_gradient_y(u.x)+center_gradient_x(u.y))+2.*sq(center_gradient_x(u.x));
#if dimension==3
	      sum+=(u.x[0,0,1]-u.x[0,0,-1])/(2.*Delta)*((u.x[0,0,1]-u.x[0,0,-1])/(2.*Delta)-(u.z[1]-u.z[-1])/(2.*Delta));
#endif      
      }

    Pk=mutke[]*sum/rho[];

      //limiter
      Pk=min(Pk,10.*epsilon[]);
      
      sourcek[]=Pk-epsilon[];
      sourceepsilon[]=C_1*epsilon[]/k[]*Pk-C_2*sq(epsilon[])/k[];
      production[]=Pk;

      double d;
#if WALL
      d=distance_to_wall(x,y);
#else
      d=HUGE;
#endif     

      if ((d<d_IP)&&(x>leading)&&(x<trailing)) {
	      sourcek[]=0.;
        sourceepsilon[]=0.;
      }
  }

   
  foreach() {
#if WALL
      d=distance_to_wall(x,y);
#else
      d=HUGE;
#endif     
      if ((d>d_IP)||(x<leading)||(x>trailing)) {
        epsilon[]=(epsilon[]+dt*C_1*epsilon[]/k[]*production[]+dt*C_2*sq(epsilon[])/sq(k[])*(dt*production[]+k[]))/(1.+2.*dt*C_2*epsilon[]/k[]+sq(dt)*C_2*sq(epsilon[])/sq(k[]));
    k[]=k[]+dt*production[]-dt*epsilon[];
      }
  }



  foreach() {
      mutke[]=rho[]*C_mu*sq(k[])/epsilon[];
  }

  foreach_face() { 
    mut.x[]=fm.x[]*face_value(rho,0)*face_value(molvis,0)+fm.x[]*face_value(mutke,0);
  }
}

#if WALL_FUNCTION
#include "wallmodelkepsilon.h"
#endif

#if 1
event acceleration (i++) {
    face vector av = a;
    scalar rhok[];
    rhok.gradient=k.gradient;
    foreach() {
      rhok[]=rho[]*k[];
    }
    foreach_face() {
      av.x[]+=-2./3.*face_gradient_x(rhok,0)/rho[];  //add term link to kinetic energy as an acceleration
  }
}
#endif