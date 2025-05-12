/**
# A Mac Cormack-like solver for the Phase-field Van der Waals equation

We start with basic basilisk framework
*/


#include "run.h"
#include "timestep.h"

/* 
The Van der Waals phase field equation can be written as in
[Nadiga, B. T., \& Zaleski, S. (1995). Investigations of a two-phase fluid model. arXiv preprint comp-gas/9511003.] hereafter NZ95. 
The model reads
$$
\partial_t \rho + \partial_i(\rho u_i)  =  0
$$
and
$$
  \partial_t (\rho u_i ) + \partial_j ( \rho u_i u_j )  =  \partial_j \sigma_{ij} 
$$
with $\rho$ the density, $\mathbf{u} = (u_i)_i$ the
velocity. These variables are declared as basilisk *fields*
*/

scalar rho[],rhop[],mu[];
vector u[],up[];

/* $\sigma_{ij}$ the stress tensor given by 
$$
\sigma_{ij} = - p_0 \delta_{ij} + \sigma_{ij}^{(v)}  + \sigma_{ij}^{(K)} 
$$
where $ \sigma_{ij}^{(v)} $ has the classical form with vanishing compression viscosity (see NZ95). 
We shall thus need */

scalar mu[];
(const) scalar cs2;
#define P0(x) ( cs2 * (x) )

/* 
and
$$
\sigma_{ij}^{(K)} = \lambda \left[ \left( \frac 12 \vert \nabla \rho \vert^2 + \rho \nabla^2 \rho \right) \delta_{ij} - \partial_i \rho \partial_j \rho \right]
$$
is the Korteweg stress tensor, which can be rewritten
$$
\sigma_{ij}^{(K)} = \lambda (\rho \nabla^2 \rho \, \delta_{ij}  - T_{ij} )
$$
(notice the sign difference with NZ95)
where the tensor $T$ is expressed in 2D as
$$
T = \left( \begin{array}{cc}
  (\partial_x \rho)^2/2 - (\partial_y \rho)^2/2 & \partial_x \rho \, \partial_y \rho \\
  \partial_x \rho \, \partial_y \rho              & - (\partial_x \rho)^2/2 + (\partial_y \rho)^2/2
\end{array} \right)
$$
thus we define */

(const) scalar lambda;

/*
We shall use the predictor corrector scheme in NZ95 written as
$$
\partial_t \mathbf{f} + \partial_x\mathbf {F}_x[\mathbf{f}] + \partial_y\mathbf {F}_y[\mathbf {f}] = 0,
$$
where
$$
\mathbf{f} = \left(\begin{array}{c}
\rho \\
u_x  \\
u_y  
\end{array}\right) 
$$
and $\mathbf{F}_i$ has {\em functional dependence} on $\mathbf{f}$ through
$$
\mathbf{F}_x[\mathbf{f}] = \mathbf{F}_x(\mathbf{f}, \partial_x \mathbf{f}, \partial_y \mathbf{f}, \nabla^2 \rho)
$$
(minor typo in NZ95: $\nabla \rightarrow \nabla^2$) and
$$
\mathbf{F}_y[\mathbf{f}] = \mathbf{F}_y(\mathbf{f}, \partial_y \mathbf{f}, \partial_x \mathbf{f}, \nabla^2 \rho).
$$
One has
$$
\mathbf{F}_x = \left(\begin{array}{c}
\rho u_x\\
\rho u_x^2 - \sigma_{xx}  \\
\rho u_x u_y - \sigma_{xy}  
\end{array}\right) 
$$
and $\mathbf{F}_y$ is obtained analogously by exchanging the 2nd and 3rd components and exchanging $x$ and $y$.
(We expect basilisk to perform this exchange automatically).
Explicitly,
$$
\mathbf{F}_x =
\left(\begin{array}{c}
\rho u_x\\
\rho u_x^2 + p_0(\rho) \\
\rho u_x u_y 
\end{array} \right)
-
\mu \left(\begin{array}{c}
0\\
2 \partial_x u_x - (2/3)(\partial_x u_x + \partial_y u_y)  \\
 \partial_x u_y + \partial_y u_x 
\end{array} \right)
+ \lambda 
\left(\begin{array}{c}
  0\\
  T_{xx} -  \rho \nabla^2 \rho \\
  T_{xy}
\end{array} \right).
$$
To discretize the equation we define a time step $\tau$, a spatial step $h$ and forward, backwards and explicit difference operators 
$$
\partial^+_x v = \frac{v(x+h) - v(x)}h
$$
$$
\partial^-_x v = \frac{v(x) - v(x-h)}h
$$
$$
\partial^c_x v = \frac{v(x+h) - v(x-h)}{2h}
$$
The basilisk expression shall be detailed below.
The predictor step is (see NZ95)
$$
\mathbf{f}^* = \mathbf{f}^n - \tau \partial_x^-  \mathbf{F}_x(\mathbf{f}, \partial_x^+ \mathbf{f}, \partial_y^c \mathbf{f}, \nabla^{2,c} \rho)
- \tau \partial_y^-  \mathbf{F}_y(\mathbf{f}, \partial_y^+ \mathbf{f}, \partial_x^c \mathbf{f}, \nabla^{2,c} \rho),
$$
where we have written $\mathbf{f}$ for $\mathbf{f}^n$ for simplicity and the corrector step is
$$
\mathbf{f}^{n+1} = \frac12 (\mathbf{f}^n + \mathbf{f}^*)
-  \frac12\tau \partial_x^+  \mathbf{F}_x(\mathbf{f}, \partial_x^- \mathbf{f}, \partial_y^c \mathbf{f}, \nabla^{2,c} \rho)
-  \frac12\tau \partial_y^+  \mathbf{F}_y(\mathbf{f}, \partial_y^- \mathbf{f}, \partial_x^c \mathbf{f}, \nabla^{2,c} \rho).
$$


/**
## Initial conditions */

event defaults (i = 0)
{
  CFL = 0.8;
}


/**
After user initialisation, we initialise the fluid properties. */

double dtmax;

event init (i = 0)
{
  boundary ((scalar *){u});

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
  dt = dtnext (timestep (u, dtmax));
}

/*
The detailed expression is as follows. For the predictor step
$$
u_x^* = u_x - \tau \{  \partial_x^- (\rho u_x^2 + p_0 ) 
-  (2/3) \partial_x^- [\mu( 2 \partial_x^+ u_x - \partial_y^c u_y )]
- (\lambda/2) \partial_x^- [ 2\rho \nabla^{2,c} \rho - (\partial_x^+ \rho)^2  +  (\partial_y^c \rho)^2 ]
$$
$$
+ \partial_y^- (\rho u_x u_y)
-  \partial_y^- [\mu( \partial_y^+ u_x - \partial_x^c u_y )]
+  \lambda  \partial_y^-   (\partial_y^+ \rho \, \partial_x^c \rho) \}
$$
The corrector step is
$$
\rho^{n+1} = \frac12 ( \rho^* + \rho^n) -\frac\tau 2 [ \partial_x^+ (\rho u_x) + \partial_y^+ ( \rho u_y )]
$$
and
$$
u_x^{n+1} = \frac12 (u_x^* + u_x^n) - \frac \tau 2 \{  \partial_x^+ (\rho u_x^2 + p_0 ) + \partial_y^+ (\rho u_x u_y)
$$
$$
-  (2/3) \partial_x^+ [\mu( 2 \partial_x^- u_x - \partial_y^c u_y )] -  \partial_y^+ [\mu( \partial_y^- u_x - \partial_x^c u_y )]
$$
$$
+ (\lambda/2) [ \partial_x^+ ( (\partial_x^- \rho)^2  -  (\partial_y^c \rho)^2  - 2\rho \nabla^{2,c} \rho ) + 2 \partial_y^+ (\partial_y^- \rho \partial_x^c \rho)] \}.
$$
*/
event advection_term (i++,last)
{
  foreach()
    foreach_dimension()
      {
	/*	
$$
\rho^* = \rho - \tau [ \partial_x^- (\rho u_x) + \partial_y^- ( \rho u_y )]
$$
	*/
	rhop[] = rho[] - (dt/Delta)*(2*rho[]*u.x[] - rho[1,0]*u.x[1,0] - rho[0,-1]*u.x[0,-1]);
	/* after minor manipulation
$$
u_x^* = u_x - \tau \{  \partial_x^- (\rho u_x^2 + p_0 ) + \partial_y^- (\rho u_x u_y)
$$
$$
-  (2/3) \partial_x^- [\mu( 2 \partial_x^+ u_x - \partial_y^c u_y )] -  \partial_y^- [\mu( \partial_y^+ u_x - \partial_x^c u_y )]
$$
$$
+ (\lambda/2) [ \partial_x^- ( (\partial_x^+ \rho)^2  -  (\partial_y^c \rho)^2 - 2\rho \nabla^{2,c} \rho ) + 2 \partial_y^- (\partial_y^+ \rho \partial_x^c \rho)] \}
$$	
	*/
	up.x[] = u[x] - (dt/Delta)*(rho[]*sq(u.x[]) + P0(rho[]) - rho[-1,0]*sq(u.x[-1,0]) + P0(rho[-1,0]) 
				    - (2./3.)*( mu[] * ( 2*(u.x[1,0] - u.x[])/Delta - (u.y[0,1] - u.y[0,-1])/(2*Delta))
						       
