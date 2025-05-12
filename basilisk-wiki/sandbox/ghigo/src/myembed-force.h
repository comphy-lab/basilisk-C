/**
# Force and torque compution with embedded boundaries */

trace
double embed_interpolate_3D (Point point, scalar s, coord p)
{
  int i = sign(p.x), j = sign(p.y);
#if dimension == 2
  if (cs[i] && cs[0,j] && cs[i,j] )
    // bilinear interpolation when all neighbors are defined
    return ((s[]*(1. - fabs(p.x)) + s[i]*fabs(p.x))*(1. - fabs(p.y)) +
        (s[0,j]*(1. - fabs(p.x)) + s[i,j]*fabs(p.x))*fabs(p.y));
#else // dimension == 3
  int k = sign(p.z);
  if (cs[i,0,0] && cs[0,j,0] && cs[i,j,0] &&
      cs[0,0,k] && cs[i,0,k] && cs[0,j,k] && cs[i,j,k] ) {
    double val_0, val_k;
    // bilinear interpolation in x-y-planes when all neighbors are defined
    val_0 = (s[0,0,0]*(1. - fabs(p.x)) + s[i,0,0]*fabs(p.x))*(1. - fabs(p.y)) +
      (s[0,j,0]*(1. - fabs(p.x)) + s[i,j,0]*fabs(p.x))*fabs(p.y);
    val_k = (s[0,0,k]*(1. - fabs(p.x)) + s[i,0,k]*fabs(p.x))*(1. - fabs(p.y)) +
      (s[0,j,k]*(1. - fabs(p.x)) + s[i,j,k]*fabs(p.x))*fabs(p.y);
    // trilinear interpolation when all neighbors are defined
    return (val_0*(1. - fabs(p.z)) + val_k*fabs(p.z));
  }
#endif
  else {
    // linear interpolation with gradients biased toward the
    // cells which are defined
    double val = s[];
    foreach_dimension() {
      int i = sign(p.x);
      if (cs[i] )
    val += fabs(p.x)*(s[i] - s[]);
      else if (cs[-i] )
    val += fabs(p.x)*(s[] - s[-i]);
    }
    return val;
  }
}

trace
void embed_force_3D (scalar p, vector u, face vector mu, coord * Fp, coord * Fmu)
{
  // fixme: this could be simplified considerably if reduction worked on vectors
  double Fp_x = 0., Fp_y = 0., Fp_z = 0., Fmu_x = 0., Fmu_y = 0., Fmu_z = 0.;
  foreach (reduction(+:Fp_x) reduction(+:Fp_y) reduction(+:Fp_z)
	   reduction(+:Fmu_x) reduction(+:Fmu_y) reduction(+:Fmu_z)) {
    if (cs[] > 0. && cs[] < 1.) {

      /**
      To compute the pressure force, we first get the coordinates of
      the barycentre of the embedded fragment, its area and normal,
      and then interpolate the pressure field on the surface. */
      
      coord n, b;
      double area = embed_geometry (point, &b, &n);
      area *= pow (Delta, dimension - 1);

      double Fn = area*embed_interpolate_3D (point, p, b);
      foreach_dimension()
	Fp_x += Fn*n.x;
          
      /**
      To compute the viscous force, we first need to retrieve the
      local value of the viscosity (ideally at the barycentre of the
      embedded fragment). This is not completely trivial since it is
      defined on the faces of the cell. We use a
      surface-fraction-weighted average value. */

      if (constant(mu.x) != 0.) {
	double mua = 0., fa = 0.;
	foreach_dimension() {
	  mua += mu.x[] + mu.x[1];
	  fa  += fs.x[] + fs.x[1];
	}
	mua /= fa;

	/**
	To compute the viscous force, we need to take into account the
	(Dirichlet) boundary conditions for the velocity on the
	surface. We only know how to do this when computing the normal
	gradient $\mathbf{\nabla}\mathbf{u}\cdot\mathbf{n}$ using the
	[embed_gradient()](#embed_gradient) function. We thus
	need to re-express the viscous force using only normal
	derivatives of the velocity field. */

	/**
	If we assume that $\mathbf{u}$ is constant on the boundary, then:
	$$
	\mathbf{{\nabla}} \mathbf{u} \cdot \mathbf{t}= \mathbf{0},
	$$
	with $\mathbf{t}$ the unit tangent vector to the boundary. We
	thus have the relations
	$$
	\mathbf{{\nabla}} \mathbf{u} = \left( \mathbf{{\nabla}} \mathbf{u}
	\cdot \mathbf{n} \right) \mathbf{n} + \left( \mathbf{{\nabla}}
	\mathbf{u} \cdot \mathbf{t} \right) \mathbf{t} = \left(
	\mathbf{{\nabla}} \mathbf{u} \cdot \mathbf{n} \right) \mathbf{n}
	$$
	$$
	\mathbf{D}= \frac{1}{2}  \left( \mathbf{{\nabla}} \mathbf{u} +
	\mathbf{{\nabla}}^T \mathbf{u} \right) = \frac{1}{2} 
	\left(\begin{array}{cc}
	2 \left( \mathbf{{\nabla}} u \cdot \mathbf{n} \right) n_x & \left(
	\mathbf{{\nabla}} u \cdot \mathbf{n} \right) n_y + \left(
	\mathbf{{\nabla}} v \cdot \mathbf{n} \right) n_x & \left(
	\mathbf{{\nabla}} u \cdot \mathbf{n} \right) n_z + \left(
	\mathbf{{\nabla}} w \cdot \mathbf{n} \right) n_x\\
	\left( \mathbf{{\nabla}} u \cdot \mathbf{n} \right) n_y + \left(
	\mathbf{{\nabla}} v \cdot \mathbf{n} \right) n_x & 2 \left(
	\mathbf{{\nabla}} v \cdot \mathbf{n} \right) n_y & \left(
	\mathbf{{\nabla}} v \cdot \mathbf{n} \right) n_z + \left(
	\mathbf{{\nabla}} w \cdot \mathbf{n} \right) n_y\\	
	\left( \mathbf{{\nabla}} u \cdot \mathbf{n} \right) n_z + \left(
	\mathbf{{\nabla}} w \cdot \mathbf{n} \right) n_x & \left(
	\mathbf{{\nabla}} v \cdot \mathbf{n} \right) n_z + \left(
	\mathbf{{\nabla}} w \cdot \mathbf{n} \right) n_y & 2 \left(
	\mathbf{{\nabla}} w \cdot \mathbf{n} \right) n_z\\
	\end{array}\right)
	$$
	$$
	\mathbf{F}_{\mu} = - \int_{\Gamma} \left(\begin{array}{c}
	\left[2 \mu \left( \mathbf{{\nabla}} u \cdot \mathbf{n} \right)
	n_x \right] n_x + \mu \left[ \left( \mathbf{{\nabla}} u \cdot \mathbf{n}
	\right) n_y + \left( \mathbf{{\nabla}} v \cdot \mathbf{n} \right) n_x
	\right] n_y + \mu \left[ \left( \mathbf{{\nabla}} u \cdot \mathbf{n}
	\right) n_z + \left( \mathbf{{\nabla}} w \cdot \mathbf{n} \right) n_x
	\right] n_z\\
	\left[2 \mu \left( \mathbf{{\nabla}} v \cdot \mathbf{n} \right)
	n_y \right] n_y + \mu \left[ \left( \mathbf{{\nabla}} u \cdot \mathbf{n}
	\right) n_y + \left( \mathbf{{\nabla}} v \cdot \mathbf{n} \right) n_x
	\right] n_x + \mu \left[ \left( \mathbf{{\nabla}} v \cdot \mathbf{n}
	\right) n_z + \left( \mathbf{{\nabla}} w \cdot \mathbf{n} \right) n_y
	\right] n_z\\
	\left[2 \mu \left( \mathbf{{\nabla}} w \cdot \mathbf{n} \right)
	n_z \right] n_z + \mu \left[ \left( \mathbf{{\nabla}} u \cdot \mathbf{n}
	\right) n_z + \left( \mathbf{{\nabla}} w \cdot \mathbf{n} \right) n_x
	\right] n_x + \mu \left[ \left( \mathbf{{\nabla}} v \cdot \mathbf{n}
	\right) n_z + \left( \mathbf{{\nabla}} w \cdot \mathbf{n} \right) n_y
	\right] n_y\\
	\end{array}\right)
	$$
	$$
	\mathbf{F}_{\mu} = - \int_{\Gamma} \left(\begin{array}{c}
	\mu \left[ \left( \mathbf{{\nabla}} u \cdot \mathbf{n} \right) 
	(n^2_x + 1) + \left( \mathbf{{\nabla}} v \cdot \mathbf{n} \right) n_x
	n_y + \left( \mathbf{{\nabla}} w \cdot \mathbf{n} \right) n_x
	n_z \right]\\
	\mu \left[ \left( \mathbf{{\nabla}} v \cdot \mathbf{n} \right) 
	(n^2_y + 1) + \left( \mathbf{{\nabla}} u \cdot \mathbf{n} \right) n_x
	n_y + \left( \mathbf{{\nabla}} w \cdot \mathbf{n} \right) n_y
	n_z \right]\\
	\mu \left[ \left( \mathbf{{\nabla}} w \cdot \mathbf{n} \right) 
	(n^2_z + 1) + \left( \mathbf{{\nabla}} u \cdot \mathbf{n} \right) n_x
	n_z + \left( \mathbf{{\nabla}} v \cdot \mathbf{n} \right) n_y
	n_z \right]
	\end{array}\right)
	$$
	*/

	coord dudn = embed_gradient (point, u, b, n);
#if dimension == 2
	foreach_dimension()
	  Fmu_x -= area*mua*(dudn.x*(sq (n.x) + 1.) +
			     dudn.y*n.x*n.y);
#else // dimension == 3
	foreach_dimension()
	  Fmu_x -= area*mua*(dudn.x*(sq (n.x) + 1.) +
			     dudn.y*n.x*n.y +
			     dudn.z*n.x*n.z);
#endif // dimension
      }
    }
  }
  foreach_dimension() {
    Fp->x = Fp_x;
    Fmu->x = Fmu_x;
  }
}

void embed_color_force_3D (scalar p, vector u, face vector mu, scalar color, coord * Fp, coord * Fmu)
{
  // fixme: this could be simplified considerably if reduction worked on vectors
  double Fp_x = 0., Fp_y = 0., Fp_z = 0., Fmu_x = 0., Fmu_y = 0., Fmu_z = 0.;
  foreach (reduction(+:Fp_x) reduction(+:Fp_y) reduction(+:Fp_z)
	   reduction(+:Fmu_x) reduction(+:Fmu_y) reduction(+:Fmu_z)) {
    if (cs[] > 0. && cs[] < 1. && color[] > 0. && color[] < 1.) {

      coord n, b;
      double area = embed_geometry (point, &b, &n);
      area *= pow (Delta, dimension - 1);

      double Fn = area*embed_interpolate_3D (point, p, b);
      foreach_dimension()
	Fp_x += Fn*n.x;
          
      if (constant(mu.x) != 0.) {
	double mua = 0., fa = 0.;
	foreach_dimension() {
	  mua += mu.x[] + mu.x[1];
	  fa  += fs.x[] + fs.x[1];
	}
	mua /= fa;

	coord dudn = embed_gradient (point, u, b, n);
#if dimension == 2
	foreach_dimension()
	  Fmu_x -= area*mua*(dudn.x*(sq (n.x) + 1.) +
			     dudn.y*n.x*n.y);
#else // dimension == 3
	foreach_dimension()
	  Fmu_x -= area*mua*(dudn.x*(sq (n.x) + 1.) +
			     dudn.y*n.x*n.y +
			     dudn.z*n.x*n.z);
#endif // dimension
      }
    }
  }
  foreach_dimension() {
    Fp->x  = Fp_x;
    Fmu->x = Fmu_x;
  }
}

/**
The torque exerted by the fluid on the solid can be written
$$
\mathbf{T}_{\Gamma} = - \int_{\partial \Gamma}
(\mathbf{x} - \mathbf{x_p})\times( - p\mathbf{I} +
2 \mu \mathbf{D}) \cdot \mathbf{n}d \partial \Gamma
$$
with $\Gamma$ the solid boundary and $\mathbf{x_p}$ its barycenter.
It can be further decomposed into a pressure (i.e. "form") torque
$$
\mathbf{T}_p = \int_{\partial \Gamma} (\mathbf{x} - \mathbf{x_p})\times
(p \mathbf{n}d) \partial \Gamma
$$
and a viscous torque
$$
\mathbf{T}_{\mu} = - \int_{\partial \Gamma} 
(\mathbf{x} - \mathbf{x_p})\times
(2 \mu \mathbf{D} \cdot \mathbf{n}d) \partial \Gamma.
$$
These two vectors are computed by the *embed_torque()* function.
*/

trace
void embed_torque_3D (scalar p, vector u, face vector mu, coord c, coord * Tp, coord * Tmu)
{
  double Tp_x = 0., Tp_y = 0., Tp_z = 0., Tmu_x = 0., Tmu_y = 0., Tmu_z = 0.;
  foreach (reduction(+:Tp_x) reduction(+:Tp_y) reduction(+:Tp_z)
	   reduction(+:Tmu_x) reduction(+:Tmu_y) reduction(+:Tmu_z))
    if (cs[] > 0. && cs[] < 1.) {

      /**
      To compute the pressure torque, we first get the coordinates of
      the barycentre of the embedded fragment, its area and normal,
      and then interpolate the pressure field on the surface. Finally,
      we compute its relative coordinates $\mathbf{x} -
      \mathbf{xp}$. */
      
      coord n, b;
      double area = embed_geometry (point, &b, &n);
      area *= pow (Delta, dimension - 1);

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
      
      double Fn = area*embed_interpolate_3D (point, p, b);
#if dimension == 2
      Tp_z += Fn*(r.x*n.y - r.y*n.x);
#else // dimension == 3      
      foreach_dimension()
	Tp_x += Fn*(r.y*n.z - r.z*n.y);
#endif // dimension

      /**
      To compute the viscous torque, we first need to retrieve the
      local value of the viscosity (ideally at the barycentre of the
      embedded fragment). This is not completely trivial since it is
      defined on the faces of the cell. We use a
      surface-fraction-weighted average value. */
      
      if (constant(mu.x) != 0.) {
	double mua = 0., fa = 0.;
	foreach_dimension() {
	  mua += mu.x[] + mu.x[1];
	  fa  += fs.x[] + fs.x[1];
	}
	mua /= fa;

	/**
	To compute the viscous torque, we need to take into account the
	(Dirichlet) boundary conditions for the velocity on the
	surface. We only know how to do this when computing the normal
	gradient $\mathbf{\nabla}\mathbf{u}\cdot\mathbf{n}$ using the
	[embed_gradient()](#embed_gradient) function. We thus
	need to re-express the viscous force using only normal
	derivatives of the velocity field. */

	coord dudn = embed_gradient (point, u, b, n);
#if dimension == 2
	double Fmu_x = 0., Fmu_y = 0.;	
	foreach_dimension()
	  Fmu_x = -area*mua*(dudn.x*(sq (n.x) + 1.) +
			     dudn.y*n.x*n.y);
#else // dimension == 3
	double Fmu_x = 0., Fmu_y = 0., Fmu_z = 0.;	
	foreach_dimension()
	  Fmu_x = -area*mua*(dudn.x*(sq (n.x) + 1.) +
			     dudn.y*n.x*n.y +
			     dudn.z*n.x*n.z);
#endif // dimension

#if dimension == 2
	Tmu_z += r.x*Fmu_y - r.y*Fmu_x;
#else // dimension == 3
	foreach_dimension()
	  Tmu_x += r.y*Fmu_z - r.z*Fmu_y;
#endif // dimension
      }
    }
#if dimension == 2
  double T_p = Tp_z, T_mu = Tmu_z;
  foreach_dimension() {
    Tp->x = T_p;
    Tmu->x = T_mu;
  }
#else // dimension == 3
  foreach_dimension() {
    Tp->x = Tp_x;
    Tmu->x = Tmu_x;
  }
#endif
}

void embed_color_torque_3D (scalar p, vector u, face vector mu, scalar color, coord c, coord * Tp, coord * Tmu)
{
  double Tp_x = 0., Tp_y = 0., Tp_z = 0., Tmu_x = 0., Tmu_y = 0., Tmu_z = 0.;
  foreach (reduction(+:Tp_x) reduction(+:Tp_y) reduction(+:Tp_z)
	   reduction(+:Tmu_x) reduction(+:Tmu_y) reduction(+:Tmu_z))
    if (cs[] > 0. && cs[] < 1. && color[] > 0. && color[] < 1.) {
      
      coord n, b;
      double area = embed_geometry (point, &b, &n);
      area *= pow (Delta, dimension - 1);

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
      
      double Fn = area*embed_interpolate_3D (point, p, b);
#if dimension == 2
      Tp_z += Fn*(r.x*n.y - r.y*n.x);
#else // dimension == 3      
      foreach_dimension()
	Tp_x += Fn*(r.y*n.z - r.z*n.y);
#endif // dimension
      
      if (constant(mu.x) != 0.) {
	double mua = 0., fa = 0.;
	foreach_dimension() {
	  mua += mu.x[] + mu.x[1];
	  fa  += fs.x[] + fs.x[1];
	}
	mua /= fa;

	coord dudn = embed_gradient (point, u, b, n);
#if dimension == 2
	double Fmu_x = 0., Fmu_y = 0.;	
	foreach_dimension()
	  Fmu_x = -area*mua*(dudn.x*(sq (n.x) + 1.) +
			     dudn.y*n.x*n.y);
#else // dimension == 3
	double Fmu_x = 0., Fmu_y = 0., Fmu_z = 0.;	
	foreach_dimension()
	  Fmu_x = -area*mua*(dudn.x*(sq (n.x) + 1.) +
			     dudn.y*n.x*n.y +
			     dudn.z*n.x*n.z);
#endif // dimension

#if dimension == 2
	Tmu_z += r.x*Fmu_y - r.y*Fmu_x;
#else // dimension == 3
	foreach_dimension()
	  Tmu_x += r.y*Fmu_z - r.z*Fmu_y;
#endif // dimension
      }
    }
#if dimension == 2
  double T_p = Tp_z, T_mu = Tmu_z;
  foreach_dimension() {
    Tp->x = T_p;
    Tmu->x = T_mu;
  }
#else // dimension == 3
  foreach_dimension() {
    Tp->x = Tp_x;
    Tmu->x = Tmu_x;
  }
#endif
}
