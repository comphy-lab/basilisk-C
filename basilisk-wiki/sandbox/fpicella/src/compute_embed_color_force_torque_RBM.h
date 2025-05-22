
/**
# Compute torques and forces on _color_ cs embed + Rigid Body Motion
Using the same approach as Ghigo's [myembed.h](/sandbox/ghigo/src/myembed.h)
but accounting for non-homogeneous surface velocity.
NECESSARY for ROTATING BODIES

### What is different wrt classic [embed.h](/src/embed.h)?
Here homogeneous $\mathbf{u}$ is considered on the embedded boundary. This is Ok if you consider a rigid body in pure translation.
Under this framework, the term $\nabla \mathbf{u} \cdot \mathbf{t} = \mathbf{0}$, and some semplifications can be done in the computation of the stress tensor $\mathbf{D}$.

The hypothesis fails if the boundary $\partial\Gamma$ is in solid rotation.
*/
/**
To compute the viscous force, we need to take into account the
(Dirichlet) boundary conditions for the velocity on the
surface. We only know how to do this when computing the normal
gradient $\mathbf{\nabla}\mathbf{u}\cdot\mathbf{n}$ using the
[embed_gradient()](#embed_gradient) function. We thus
need to re-express the viscous force using only normal
derivatives of the velocity field.

If we assume that $\mathbf{u}$ is NOT constant on the boundary (i.e. RBM), then
$$
\mathbf{{\nabla}} \mathbf{u} \cdot \mathbf{t} \neq \mathbf{0}
$$
with $\mathbf{t}$ the unit tangent vector to the boundary. We
thus have the relations
$$
\mathbf{{\nabla}} \mathbf{u} = \left( \mathbf{{\nabla}} \mathbf{u}
\cdot \mathbf{n} \right) \mathbf{n} + \left( \mathbf{{\nabla}}
\mathbf{u} \cdot \mathbf{t} \right) \mathbf{t}
$$
where $\mathbf{t} = \left[ \begin{array}{c} t_x \\ t_y \end{array} \right] = \left[ \begin{array}{c} -n_y \\ n_x \end{array} \right]$, in 2D.

After some (tedious) algebra, provides the shear stress on a roto-translating boundary:
$$
\mathbf{F}_{\mu} = - \int_{\Gamma} \left(\begin{array}{c}
\mu \left[ \left( \mathbf{{\nabla}} u \cdot \mathbf{n} \right) 
(n^2_x + 1) + \left( \mathbf{{\nabla}} v \cdot \mathbf{n} \right) n_x
n_y \color{red}{+(\mathbf{\nabla} u \cdot \mathbf{t})(-n_x n_y) +(\mathbf{\nabla} v \cdot \mathbf{t})(-n_y^2)} \right]\\
\mu \left[ \left( \mathbf{{\nabla}} v \cdot \mathbf{n} \right) 
(n^2_y + 1) + \left( \mathbf{{\nabla}} u \cdot \mathbf{n} \right) n_x
n_y \color{red}{+(\mathbf{\nabla} v \cdot \mathbf{t})(+n_x n_y) +(\mathbf{\nabla} v \cdot \mathbf{t})(+n_x^2)} \right]
\end{array}\right)
$$
Terms in red here arise when accounting for a boundary whose velocity is not homogeneous, like in the case of a rigidly rotating particle.

To compute $\mathbf{\nabla} \mathbf{u} \cdot \mathbf{t}$, it is quite simple for a rigid rotating circular particle of radius $r_0$, centered in $x_0,y_0$ and having angular velocity $\Omega_z$.
$$
\mathbf{\nabla} \mathbf{u} \cdot \mathbf{t} = \left[ \begin{array}{c} \mathbf{\nabla} u \cdot \mathbf{t} \\ \mathbf{\nabla} v \cdot \mathbf{t} \end{array} \right] = \left[\begin{array}{c} -\Omega_z(x-x_0)/r_0 \\ -\Omega_z(y-y_0)/r_0\end{array} \right]
$$
*/



void embed_color_torque_RBM (scalar p, vector u, face vector mu, scalar color, coord c, coord pOmega, double pRADIUS, coord * Tp, coord * Tmu)
// Additional term here is OMEGA, the ANGULAR VELOCITY and the RADIUS of the particle...
{
	#if dimension == 3
	if(t==0)
		fprintf(stderr,"WARNING, RBM to be checked if rotation is on!!!\n");
	#endif

  coord Tps = {0}, Tmus = {0};
  foreach (reduction(+:Tps) reduction(+:Tmus))
    if (cs[] > 0. && cs[] < 1. && color[] > 0. && color[] < 1.) {

      coord n, b;
      double area = embed_geometry (point, &b, &n);
      area *= pow (Delta, dimension - 1);

      /**
      In addition to the quantities computed in the function
      *embed_force()*, we also compute the relative coordinates
      $\mathbf{x} - \mathbf{x}_{\Gamma}$. */

      // The coordinate x,y,z are not permuted with foreach_dimension()
      coord r = {x,y,z};
      // In case of a periodic domain, we shift the position of the center
      foreach_dimension() {
  r.x += b.x*Delta - c.x;
  if (Period.x) {
    if (fabs (r.x) > fabs (r.x + (L0)))
      r.x += (L0);
    if (fabs (r.x) > fabs (r.x - (L0)))
      r.x -= (L0);
  }
      }

      double Fn = area*embed_interpolate (point, p, b);
#if dimension == 2
      Tps.x += Fn*(r.x*n.y - r.y*n.x);
      Tps.y = Tps.x;
#else // dimension == 3      
      foreach_dimension()
  Tps.x += Fn*(r.y*n.z - r.z*n.y);
#endif // dimension

      if (constant(mu.x) != 0.) {
  double mua = 0., fa = 0.;
  foreach_dimension() {
    mua += mu.x[] + mu.x[1];
    fa  += fs.x[] + fs.x[1];
  }
  mua /= (fa + SEPS);

  coord dudn = embed_gradient (point, u, b, n);
  coord Fmus = {0};
/*
### Rigid Body Motion, RBM
If I'm accounting for pure translation, then Stephane's (and Arthur's) implementation is ok.
But what if the particle rotates at angular velocity pOmega.
After some (tedious but ok) derivation I obtain:
BLABLABLA, carnet jaune.
So, I need some additional quantities here to play with.
*/
// RBM
	//assert (dimension == 2);// for the moment, it works only in 2D...
  coord dudt = {0};// RBM, tangential component of the velocity gradient. it is non zero ONLY if the particle rotates.
	foreach_dimension()
		dudt.x = -pOmega.x*(r.x)/pRADIUS;
//	fprintf(stderr,"r.x r.y DUDT.x DUDT.y RADIUS %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e \n",r.x,r.y,dudt.x,dudt.y,pRADIUS);
// RBM

#if dimension == 2
  foreach_dimension()
    Fmus.x = -area*mua*(dudn.x*(sq (n.x) + 1.) +
            dudn.y*n.x*n.y);
#else // dimension == 3
  foreach_dimension()
    Fmus.x = -area*mua*(dudn.x*(sq (n.x) + 1.) +
            dudn.y*n.x*n.y +
            dudn.z*n.x*n.z);
#endif // dimension
// RBM
	coord SIGN = {-1.,1.,0.};
	foreach_dimension()
		Fmus.x += +area*mua*(
              dudt.x*SIGN.x*n.x*n.y + dudt.y*(SIGN.x*sq(n.y))
                        );

// RBM

#if dimension == 2
  Tmus.x += r.x*Fmus.y - r.y*Fmus.x;
  Tmus.y = Tmus.x;
#else // dimension == 3
  foreach_dimension()
    Tmus.x += r.y*Fmus.z - r.z*Fmus.y;
#endif // dimension
      }
    }

  *Tp = Tps; *Tmu = Tmus;
}

/**
The function *embed_color_force_RBM* computes the force on the colored
embedded boundaries, using the user defined color scalar *color*
as well as Rigid Body Motion of a solid particle.
. */
void embed_color_force_RBM (scalar p, vector u, face vector mu, scalar color, coord c, coord pOmega, double pRADIUS, coord * Fp, coord * Fmu)
// Additional term here is OMEGA, the ANGULAR VELOCITY and the RADIUS of the particle...
{
	#if dimension == 3
	if(t==0)
		fprintf(stderr,"WARNING, RBM to be checked if rotation is on!!!\n");
	#endif
  coord Fps = {0}, Fmus = {0};
  foreach (reduction(+:Fps) reduction(+:Fmus))
    if (cs[] > 0. && cs[] < 1. && color[] > 0. && color[] < 1.) {

      coord n, b;
      double area = embed_geometry (point, &b, &n);
      area *= pow (Delta, dimension - 1);

      /**
      In addition to the quantities computed in the function
      *embed_force()*, we also compute the relative coordinates
      $\mathbf{x} - \mathbf{x}_{\Gamma}$. */

      // The coordinate x,y,z are not permuted with foreach_dimension()
      coord r = {x,y,z};
      // In case of a periodic domain, we shift the position of the center
      foreach_dimension() {
  r.x += b.x*Delta - c.x;
  if (Period.x) {
    if (fabs (r.x) > fabs (r.x + (L0)))
      r.x += (L0);
    if (fabs (r.x) > fabs (r.x - (L0)))
      r.x -= (L0);
  }
      }

      double Fn = area*embed_interpolate (point, p, b);
      foreach_dimension()
  Fps.x += Fn*n.x;

//#if dimension == 2
//      Tps.x += Fn*(r.x*n.y - r.y*n.x);
//      Tps.y = Tps.x;
//#else // dimension == 3      
//      foreach_dimension()
//  Tps.x += Fn*(r.y*n.z - r.z*n.y);
//#endif // dimension

      if (constant(mu.x) != 0.) {
  double mua = 0., fa = 0.;
  foreach_dimension() {
    mua += mu.x[] + mu.x[1];
    fa  += fs.x[] + fs.x[1];
  }
  mua /= (fa + SEPS);

  coord dudn = embed_gradient (point, u, b, n);
  //coord Fmus = {0};
/*
### Rigid Body Motion, RBM
If I'm accounting for pure translation, then Stephane's (and Arthur's) implementation is ok.
But what if the particle rotates at angular velocity pOmega.
After some (tedious but ok) derivation I obtain:
BLABLABLA, carnet jaune.
So, I need some additional quantities here to play with.
*/
// RBM
	//assert (dimension == 2);// for the moment, it works only in 2D...
  coord dudt = {0};// RBM, tangential component of the velocity gradient. it is non zero ONLY if the particle rotates.
	foreach_dimension()
		dudt.x = -pOmega.x*(r.x)/pRADIUS;
//	fprintf(stderr,"r.x r.y DUDT.x DUDT.y RADIUS %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e \n",r.x,r.y,dudt.x,dudt.y,pRADIUS);
// RBM

#if dimension == 2
  foreach_dimension()
    Fmus.x += -area*mua*(dudn.x*(sq (n.x) + 1.) +
            dudn.y*n.x*n.y);
#else // dimension == 3
  foreach_dimension()
    Fmus.x += -area*mua*(dudn.x*(sq (n.x) + 1.) +
            dudn.y*n.x*n.y +
            dudn.z*n.x*n.z);
#endif // dimension
// RBM
	coord SIGN = {-1.,1.,0.};
	foreach_dimension()
		Fmus.x += +area*mua*(
              dudt.x*SIGN.x*n.x*n.y + dudt.y*(SIGN.x*sq(n.y))
                        );

// RBM
      }
    }

  *Fp = Fps; *Fmu = Fmus;
}
