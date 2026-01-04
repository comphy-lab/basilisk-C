/**

---
title: Van der Waals diffuse interface model with Mac Cormack scheme
...
# Equations and code

This solver implements a scheme for a diffuse interface model. The intent is to use 
it for the meso-scale modelling of various physical problems involving interfaces,
such as contact lines, thin sheet break up etc.

## Basics
The Van der Waals phase field equation can be written as in
[Nadiga, B. T., \& Zaleski, S. (1995). Investigations of a two-phase fluid model. arXiv preprint comp-gas/9511003.](https://arxiv.org/abs/comp-gas/9511003) hereafter NZ95. 
The model reads
$$
\partial_t \rho + \partial_i(\rho u_i)  =  0
$$
and
$$
  \partial_t (\rho u_i ) + \partial_j ( \rho u_i u_j )  =  \partial_j \sigma_{ij} 
$$
with $\rho$ the density, $\mathbf{u} = (u_i)_i$ the
velocity. We will note $\mathbf{q} = \rho \mathbf{u}$ the momentum. 
The solver will start with the basic basilisk framework
*/
#define BGHOSTS 2
#include "run.h"
#include "timestep-vdw.h"
#include "viscosity.h"
#include "weno2.h"

/**
For implicit integration of the viscous stress
*/
mgstats mgu;

scalar normGradRho[];


/**
The variables are declared as basilisk *fields*
*/

scalar rho[];
vector q[],u[];

/**

We declare a function that will give the pressure, the user should
define this function in his code

*/

double P0( double x);

/** 
$\sigma_{ij}$ the stress tensor given by 
$$
\sigma_{ij} = - p_0 \delta_{ij} + \sigma_{ij}^{(v)}  + \sigma_{ij}^{(K)} 
$$
where $\sigma_{ij}^{(v)}$ has the Newtonian form with vanishing compression viscosity
$$
\sigma_{ij}^{(v)} = \mu ( \partial_i u_j + \partial_j u_i - (2/3) (\nabla \cdot \mathbf{u})\delta_{ij} )
$$
where $\mu$ is the viscosity. We shall thus need 
*/
face vector mu[];
/**
and
$$
\sigma_{ij}^{(K)} = \lambda \left[ \left( \frac 12 \vert \nabla \rho \vert^2 + \rho \nabla^2 \rho \right) \delta_{ij} - \partial_i \rho \partial_j \rho \right]
$$
is the Korteweg stress tensor, which can be rewritten
$$
\sigma_{ij}^{(K)} = \lambda (\rho \nabla^2 \rho \, \delta_{ij}  - T_{ij} )
$$
(notice the sign difference with NZ95)
where the tensor $T$ is expressed in 3D as
$$
T = \left( \begin{array}{ccc }
  (\partial_x \rho)^2/2 - (\partial_y \rho)^2/2 - (\partial_z \rho)^2/2 & \partial_x \rho \, \partial_y \rho                                    & \partial_x \rho \, \partial_z \rho \\
  \partial_x \rho \, \partial_y \rho                                    & (\partial_y \rho)^2/2 - (\partial_x \rho)^2/2 - (\partial_z \rho)^2/2 & \partial_y \rho \, \partial_z \rho \\
  \partial_x \rho \, \partial_z \rho                                    & \partial_y \rho \, \partial_z \rho                                    & (\partial_z \rho)^2/2 - (\partial_x \rho)^2/2 - (\partial_y \rho)^2/2
\end{array} \right)
$$
thus we define */

(const) double lambda;

/**
## Mac Cormack scheme
We shall use the predictor corrector scheme in NZ95 written as
$$
\partial_t \mathbf{f} + \partial_x\mathbf {F}_x[\mathbf{f}] + \partial_y\mathbf {F}_y[\mathbf {f}] + + \partial_z\mathbf {F}_z[\mathbf{f}] = 0,
$$
where
$$
\mathbf{f} = \left(\begin{array}{c}
\rho \\
\rho u_x  \\
\rho u_y  \\
\rho u_z  
\end{array}\right) 
$$
and $\mathbf{F}_i$ has *functional dependence* on $\mathbf{f}$ through
$$
\mathbf{F}_x[\mathbf{f}] = \mathbf{F}_x(\mathbf{f}, \partial_x \mathbf{f}, \partial_y \mathbf{f}, \partial_z \mathbf{f}, \nabla^2 \rho)
$$
(minor typo in NZ95: $\nabla \rightarrow \nabla^2$) and
$$
\mathbf{F}_y[\mathbf{f}] = \mathbf{F}_y(\mathbf{f}, \partial_y \mathbf{f}, \partial_x \mathbf{f}, \partial_z \mathbf{f}, \nabla^2 \rho) \\
\mathbf{F}_z[\mathbf{f}] = \mathbf{F}_z(\mathbf{f}, \partial_z \mathbf{f}, \partial_x \mathbf{f}, \partial_y \mathbf{f}, \nabla^2 \rho).
$$
One has
$$
\mathbf{F}_x = \left(\begin{array}{c}
\rho u_x\\
\rho u_x^2 - \sigma_{xx}  \\
\rho u_x u_y - \sigma_{xy}  \\
\rho u_x u_z - \sigma_{xz}  
\end{array}\right) 
$$
and $\mathbf{F}_y$, $\mathbf{F}_{z}$ are obtained analogously by doing the proper permutation.
(We expect basilisk to perform this exchange automatically).\\

Explicitly,
$$
\mathbf{F}_x =
\left(\begin{array}{c}
\rho u_x\\
\rho u_x^2 + p_0(\rho) \\
\rho u_x u_y \\
\rho u_x u_z 
\end{array} \right)
-
\mu \left(\begin{array}{c}
0\\
2 \partial_x u_x - (2/3)(\partial_x u_x + \partial_y u_y + \partial_z u_z)  \\
 \partial_x u_y + \partial_y u_x \\
 \partial_x u_z + \partial_z u_x 
\end{array} \right)
+ \lambda 
\left(\begin{array}{c}
  0\\
  T_{xx} -  \rho \nabla^2 \rho \\
  T_{xy} \\
  T_{xz} 
\end{array} \right).
$$
To discretize the equation we define a time step $\tau$, a spatial step $h$ and forward, backwards and centered difference operators 
$$
\partial^+_x v = \frac{v(x+h) - v(x)}h
$$
$$
\partial^-_x v = \frac{v(x) - v(x-h)}h
$$
$$
\partial^c_x v = \frac{v(x+h) - v(x-h)}{2h}
$$
The basilisk expressions are [obvious](/Basilisk C#stencils).
In addition we will need the centered Laplacian $\nabla^{2,c} = \partial_x^+\partial_x^-  + \partial_y^+\partial_y^- + \partial_z^+\partial_z^-$ that is defined and written as in [here](/Basilisk C#stencils).
The predictor step is (see NZ95)

$$
\mathbf{f}^* = \mathbf{f}^n - \tau \partial_x^-  \mathbf{F}_x(\mathbf{f}, \partial_x^+ \mathbf{f}, \partial_y^c \mathbf{f}, \partial_z^c \mathbf{f}, \nabla^{2,c} \rho)
- \tau \partial_y^-  \mathbf{F}_y(\mathbf{f}, \partial_y^+ \mathbf{f}, \partial_x^c \mathbf{f}, \partial_z^c \mathbf{f}, \nabla^{2,c} \rho)
- \tau \partial_z^-  \mathbf{F}_z(\mathbf{f}, \partial_z^+ \mathbf{f}, \partial_x^c \mathbf{f}, \partial_y^c \mathbf{f}, \nabla^{2,c} \rho),
$$

where we have written $\mathbf{f}$ for $\mathbf{f}^n$ for simplicity and the corrector step is

$$
\mathbf{f}^{n+1} = \frac12 (\mathbf{f}^n + \mathbf{f}^*)
-  \frac12\tau \partial_x^+  \mathbf{F}_x(\mathbf{f^*}, \partial_x^- \mathbf{f^*}, \partial_y^c \mathbf{f^*}, \partial_z^c \mathbf{f^*}, \nabla^{2,c} \rho^*)
-  \frac12\tau \partial_y^+  \mathbf{F}_y(\mathbf{f^*}, \partial_y^- \mathbf{f^*}, \partial_x^c \mathbf{f^*}, \partial_z^c \mathbf{f^*}, \nabla^{2,c} \rho^*)
-  \frac12\tau \partial_z^+  \mathbf{F}_z(\mathbf{f^*}, \partial_y^- \mathbf{f^*}, \partial_x^c \mathbf{f^*}, \partial_y^c \mathbf{f^*}, \nabla^{2,c} \rho^*).
$$

## Initial conditions 

What follows is shamelessly pumped from [the centered Navier-Stokes solver](/src/navier-stokes/centered.h)
*/

event defaults (i = 0)
{
  CFL = 0.1;
}


/**
After user initialisation, we initialise the fluid properties. 
*/

double dtmax;


event init (i = 0)
{
  boundary (all); 
  /**
  We update fluid properties. */

  event ("properties");

  /**
  We set the initial timestep (this is useful only when restoring from
  a previous run). */

  dtmax = DT;
  event ("stability");
}


/**
## Time integration

The timestep for this iteration is controlled by the CFL condition,
applied to the centered velocity field $\mathbf{u}$; and the
timing of upcoming events. */

event set_dtmax (i++,last) dtmax = DT;

event stability (i++,last) {
  double mydtmax = HUGE;
  foreach(reduction(min:mydtmax)){
    if ( lambda != 0.) {
      double dt = 0.6*sq(Delta)/lambda;
      if (dt < mydtmax) mydtmax = dt;
    }
  }
  dtmax = min(mydtmax,dtmax);
  foreach()
    foreach_dimension()
      u.x[] = q.x[]/rho[];
  dt = dtnext (timestep (u, q, dtmax));
}

/**
## Mac Cormack predictor_corrector event
*/
event predictor_corrector(i++,last)
{
/** ### Laplacian 
 We take the Laplacian of $\rho$ without the $1/\Delta^2$ factor, the division by the
 cell size $\Delta$ will be added below.
 We declare it as an automatic variable:
 */
  scalar laprho[],rhop[];
  vector qp[];
  trash({rhop});
/** ### Predictor step \& Laplacian 

 We take the Laplacian of $\rho$ without the $1/\Delta^2$ factor, the division by the
 cell size $\Delta$ will be added below.
 We declare it as an automatic variable.

 We also compute the predictor for the density:
$$
\rho^* = \rho - \tau \partial_i^- (\rho u_i)
$$
*/ 
  boundary({rho});
  restriction({rho});
  foreach(){
    laprho[] = rho[1] + rho[-1] + rho[0,1] + rho[0,-1] + rho[0,0,1] + rho[0,0,-1]- 6*rho[]; 
    rhop[] = rho[];
    foreach_dimension(){
      rhop[]   -= (dt/Delta)*(q.x[]  - q.x[-1,0]);
      u.x[]     = q.x[]/rho[];
    }
  }
  boundary({laprho,rhop});
  foreach(){
    normGradRho[] = 1.e-15;
    foreach_dimension(){
      normGradRho[] += u.x[]> 0. ? 
         sq(WENOdiff_x(point,rhop,-1)) : 
         sq(WENOdiff_x(point,rhop,1));
    }
    normGradRho[] = sqrt(normGradRho[]);
  }
  boundary({normGradRho,rhop});

  boundary(( scalar *) {u});
  /**
 Predictor for velocity
 $$
q_x^* = q_x - \tau \{  \partial_x^- (\rho u_x^2 + p_0 ) 
-  (2/3) \partial_x^- [\mu( 2 \partial_x^+ u_x - \partial_y^c u_y - \partial_z^c u_z )]
- (\lambda/2) \partial_x^- [ 2\rho \nabla^{2,c} \rho - (\partial_x^+ \rho)^2  +  (\partial_y^c \rho)^2 +  (\partial_z^c \rho)^2 ]
$$
$$
+ \partial_y^- (\rho u_x u_y)
-  \partial_y^- [\mu( \partial_y^+ u_x - \partial_x^c u_y )]
+  \lambda  \partial_y^-   (\partial_y^+ \rho \, \partial_x^c \rho)
$$
$$
+  \partial_z^- (\rho u_x u_z)
-  \partial_z^- [\mu( \partial_z^+ u_x - \partial_x^c u_z )]
+  \lambda  \partial_z^-   (\partial_z^+ \rho \, \partial_x^c \rho) \}
$$
        
which   after minor manipulation is

$$
u_x^* = u_x - \tau \{  \partial_x^- (\rho u_x^2 + p_0 ) + \partial_y^- (\rho u_x u_y) + \partial_z^- (\rho u_x u_z)
$$
$$
-  (2/3) \partial_x^- [\mu( 2 \partial_x^+ u_x - \partial_y^c u_y - \partial_z^c u_z )] -  \partial_y^- [\mu( \partial_y^+ u_x + \partial_x^c u_y )] -  \partial_z^- [\mu( \partial_z^+ u_x + \partial_x^c u_z )]
$$ 
This is a code (not compiled) for the last line
*/
#if naught
 - (2/3)*(  mu[ 0,0] * (2*(u.x[1,0]-u.x[ 0,0])/Delta - 0.5*(u.y[ 0,1]-u.y[ 0,-1]) - 0.5*(u.z[ 0,0,1]-u.z[ 0,0,-1]))/Delta 
          - mu[-1,0] * (2*(u.x[0,0]-u.x[-1,0])/Delta - 0.5*(u.y[-1,1]-u.y[-1,-1]) - 0.5*(u.z[-1,0,1]-u.z[-1,0,-1]))/Delta 
         )/Delta  
 - (   mu[ 0,0] * ( (u.x[0,1] - u.x[0, 0])/Delta + 0.5*(u.y[1, 0] - u.y[-1, 0])/Delta) 
     - mu[0,-1] * ( (u.x[0,0] - u.x[0,-1])/Delta + 0.5*(u.y[1,-1] - u.y[-1,-1])/Delta)
- (   mu[ 0,0,0]  * ( (u.x[0,0,1] - u.x[0,0, 0])/Delta + 0.5*(u.z[1,0, 0] - u.z[-1,0, 0])/Delta)
     - mu[0,0,-1] * ( (u.x[0,0,0] - u.x[0,0,-1])/Delta + 0.5*(u.z[1,0,-1] - u.z[-1,0,-1])/Delta)
   )/Delta   
#endif
/**
$$
+ (\lambda/2) [ \partial_x^- ( (\partial_x^+ \rho)^2  -  (\partial_y^c \rho)^2 -  (\partial_z^c \rho)^2 - 2\rho \nabla^{2,c} \rho ) + 2 \partial_y^- (\partial_y^+ \rho \partial_x^c \rho) + 2 \partial_z^- (\partial_z^+ \rho \partial_x^c \rho)] \}
$$
Here is again a code (not compiled) for the above line
*/
#if naught
 + 0.5*lambda*( sq((rho[1,0]-rho[ 0,0])/Delta) - 0.25*sq((rho[ 0,1]-rho[0 ,-1])/Delta) - 0.25*sq((rho[ 0,0,1]-rho[0 ,0,-1])/Delta) - 2*rho[ 0,0]*laprho[ 0,0] 
               -sq((rho[0,0]-rho[-1,0])/Delta) + 0.25*sq((rho[-1,1]-rho[-1,-1])/Delta) + 0.25*sq((rho[-1,0,1]-rho[-1,0,-1])/Delta) + 2*rho[-1,0]*laprho[-1,0] 
               + 2*( ((rho[0,1]-rho[ 0,0])/Delta)*(rho[1, 0]-rho[-1, 0])/(2.*Delta)
                    -((rho[0,0]-rho[0,-1])/Delta)*(rho[1,-1]-rho[-1,-1])/(2.*Delta) )/Delta
               + 2*( ((rho[0,0,1]-rho[0,0, 0])/Delta)*(rho[1,0, 0]-rho[-1,0, 0])/(2.*Delta)
                    -((rho[0,0,0]-rho[0,0,-1])/Delta)*(rho[1,0,-1]-rho[-1,0,-1])/(2.*Delta) )/Delta
              )
#endif
/**
   After some simplifications here is the actual code
   */
  trash({qp});
  foreach()
    foreach_dimension()
      qp.x[] = q.x[] - (dt/Delta)*
  (                          rho[ 0,0]*sq(u.x[ 0,0]) + P0(rho[ 0,0]) 
                           - rho[-1,0]*sq(u.x[-1,0]) - P0(rho[-1,0]) 
                           + rho[]*u.x[]*u.y[] - rho[0,-1]*u.x[0,-1]*u.y[0,-1] 
                           + rho[]*u.x[]*u.z[] - rho[0,0,-1]*u.x[0,0,-1]*u.z[0,0,-1]
                           
     
    + 0.5*lambda*( sq(rho[1,0]-rho[]    ) - 0.25*sq(rho[0,1] -rho[0,-1] ) - 0.25*sq(rho[0,0,1] -rho[0,0,-1] )   - 2*rho[ 0,0]*laprho[] 
                  -sq(rho[0,0]-rho[-1,0]) + 0.25*sq(rho[-1,1]-rho[-1,-1]) + 0.25*sq(rho[-1,0,1] -rho[-1,0,-1] ) + 2*rho[-1,0]*laprho[-1,0] 
             + (rho[0,1]-rho[0, 0]) * (rho[1, 0]-rho[-1, 0]) - (rho[]-rho[0,-1])   * (rho[1,-1]  -rho[-1,-1]) 
             + (rho[0,0,1]-rho[])   * (rho[1, 0]-rho[-1, 0]) - (rho[]-rho[0,0,-1]) * (rho[1,0,-1]-rho[-1,0,-1]) 
      )/Delta
  );
  boundary(( scalar *) {qp});
  restriction((scalar *){qp});
/** ### Corrector step 
Predicted Laplacian
*/
 foreach(){
   laprho[] = rhop[1] + rhop[-1] + rhop[0,1] + rhop[0,-1] + rhop[0,0,1] + rhop[0,0,-1]- 6*rhop[]; 
/**
Corrector for density 
$$
\rho^{n+1} = \frac 12 [ \rho^n + \rho^* ]  - \frac \tau 2  \partial_i^+ (\rho^* u^*_i)
$$
*/
   rho[] = 0.5*(rhop[]+rho[]);
   foreach_dimension(){
     rho[] -= 0.5*(dt/Delta)*(qp.x[1,0] - qp.x[]);
     u.x[] = qp.x[]/rhop[];  
   }
 }
 boundary({laprho,rho});
   
   foreach(){
    normGradRho[] = 1.e-15;
    foreach_dimension(){
      normGradRho[] += u.x[]> 0. ? 
         sq(WENOdiff_x(point,rho,-1)) : 
         sq(WENOdiff_x(point,rho,1));
    }
    normGradRho[] = sqrt(normGradRho[]);
  }

  boundary({normGradRho,rho});

  boundary(( scalar *) {u});

/**
 Corrector for velocity
 $$
q_x^{n+1} = \frac 12 (q_x^* + q_x^n )  - \frac \tau 2  \{  \partial_x^+ (\rho u_x^{* 2} + p_0^* ) + \partial_y^+ (\rho u_x^* u_y^*) + \partial_z^+ (\rho u_x^* u_z^*)
$$
$$
-  (2/3) \partial_x^+ [\mu( 2 \partial_x^- u_x^* - \partial_y^c u_y^*  - \partial_z^c u_z^* )] -  \partial_y^+ [\mu( \partial_y^- u_x^* + \partial_x^c u_y^* )] -  \partial_z^+ [\mu( \partial_z^- u_x^* + \partial_x^c u_z^* )]
$$
$$
+ (\lambda/2) [ \partial_x^+ ( (\partial_x^- \rho^*)^2  -  (\partial_y^c \rho^*)^2 -  (\partial_z^c \rho^*)^2 - 2\rho^* \nabla^{2,c} \rho^* ) 
+ 2 \partial_y^+ (\partial_y^- \rho \partial_x^c \rho^*) + 2 \partial_z^+ (\partial_z^- \rho \partial_x^c \rho^*)] \}
$$
We transcribe this, noticing that for any field f
$$ \partial_x^+ \partial_x^- f = \partial_x^- \partial_x^+ f  $$
and
$$
\partial_x^+ (\partial_x^- f)^2 = \partial_x^- (\partial_x^+ f)^2 .
$$
*/
  foreach()
    foreach_dimension()
      q.x[] = 0.5*(qp.x[] + q.x[]) - 0.5*(dt/Delta)*
  (                          rhop[1,0]*sq(u.x[1,0]) + P0(rhop[1,0]) 
                           - rhop[0,0]*sq(u.x[0,0]) - P0(rhop[0,0]) 
         + rhop[0,1]*u.x[0,1]*u.y[0,1] - rhop[0,0]*u.x[0,0]*u.y[0,0] 
         + rhop[0,0,1]*u.x[0,0,1]*u.z[0,0,1] - rhop[0,0]*u.x[0,0]*u.z[0,0] 
         
   + 0.5*lambda*( sq(rhop[1,0]-rhop[]    ) - 0.25*sq(rhop[0,1] -rhop[0,-1]  ) - 0.25*sq(rhop[0,0,1] -rhop[0,0,-1]  ) - 2*rhop[1,0]*laprho[1,0] 
                   -sq(rhop[0,0]-rhop[-1,0]) + 0.25*sq(rhop[-1,1]-rhop[-1,-1]) + 0.25*sq(rhop[-1,0,1]-rhop[-1,0,-1]) + 2*rhop[0,0]*laprho[0,0] 
             +  (rhop[0,1]-rhop[]  ) * (rhop[1,0]-rhop[-1,0]) - (rhop[]-rhop[0,-1])  *(rhop[1,-1]  -rhop[-1,-1]) 
             +  (rhop[0,0,1]-rhop[]) * (rhop[1,0]-rhop[-1,0]) - (rhop[]-rhop[0,0,-1])*(rhop[1,0,-1]-rhop[-1,0,-1]) 
             /*maybe one more line here*/
                )/Delta
  );
  boundary(( scalar *) {q});
  restriction((scalar *){q});
}

event viscous_term(i++,last){
  foreach_face(){
    mu.x[] = MU;
    if(rho[] == 0.)fprintf(stderr, "%g %g %g\n", x,y,rho[]);
  }

  boundary((scalar *){mu});
  mgu = viscosity (u, mu, rho, dt, mgu.nrelax);
}
/**
# To do
 - Work out the surface tension as a function of the equation of state. 
 - Create examples and tests as in NZ95.

# Usage
  - [bump2D-vdw-1D.c]().
  - [phasesep-1D.c]().
*/

