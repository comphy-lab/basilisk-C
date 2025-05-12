/**
# Electrohydrodynamic stresses

The EHD force density, $\mathbf{f}_e$, can be computed as the
divergence of the Maxwell stress tensor $\mathbf{M}$,
$$
M_{ij} = \varepsilon (E_i E_j - \frac{E^2}{2}\delta_{ij})
$$
where $E_i$ is the $i$-component of the electric field,
$\mathbf{E}=-\nabla \phi$ and $\delta_{ij}$ is the Kronecker delta.

We need to add the corresponding acceleration to the 
[Navier--Stokes solver](/src/navier-stokes/centered.h).

If the acceleration vector *a* (defined by the Navier--Stokes solver)
is constant, we make it variable. */

event defaults (i = 0) {
  if (is_constant (a.x))
    a = new face vector;
}

/**
The logical place of these macros is in *common.h* as others. However,
meanwhile... */

#ifndef EMBED
#define face_avg_gradient_t1_x(a,i) \
  ((a[1,i-1] + a[1,i] - a[-1,i-1] - a[-1,i])/(4.*Delta)) 
#define face_avg_gradient_t2_x(a,i) \
  ((a[1,0,i-1] + a[1,0,i] - a[-1,0,i-1] - a[-1,0,i])/(4.*Delta)) 

#define face_avg_gradient_t1_y(a,i) \
  ((a[i-1,1] + a[i,1] - a[i-1,-1] - a[i,-1])/(4.*Delta)) 
#define face_avg_gradient_t2_y(a,i) \
  ((a[0,1,i-1] + a[0,1,i] - a[0,-1,i-1] - a[0,-1,i])/(4.*Delta))

#define face_avg_gradient_t1_z(a,i) \
  ((a[i-1,0,1] + a[i,0,1] - a[i-1,0,-1] - a[i,0,-1])/(4.*Delta)) 
#define face_avg_gradient_t2_z(a,i) \
  ((a[0,i-1,1] + a[0,i,1] - a[i-1,0,-1] - a[0,i,-1])/(4.*Delta)) 
#endif

/**
If *COEF* is true the bilinear interpolation is used not only on the
embed segment but for the faces of the mixed cells. */

#if EMBED

#define COEF 1

#if COEF

static void electric_force_mixed (Point point, scalar cs, scalar phi,
				  face vector perm, coord * force)
{
  assert(cs[] > 0. && cs[] < 1.); //Only mixed cells

  double area = 0.;
  coord o, n, of, E, np;
  pseudo_t A1;
  pseudo_v A2;
  

  area = coef_bilinear_corrected_embed_gradient (point, phi, &A1, &A2, &o, &n);
  foreach_dimension()
    np.x = sign(-n.x);

  /**
  We compute first the momentum flux through the embed segment.*/
  
  double En = 0., E2 = 0., fa = 0., perma = 0.;
  foreach_dimension() {
    E.x = A1.x.x/Delta;
    En += E.x*n.x;
    E2 +=  sq(E.x);
    perma += perm.x[] + perm.x[1];
    fa += fm.x[] + fm.x[1];
  }
  perma /= fa;
  
  foreach_dimension()
    force->x = perma*(En*E.x-E2/2.*n.x)*area;

  for(int ii = 0; ii < 2; ii++)
    foreach_dimension () {
    
      of.x = -o.x + (2*ii-1)/2.;
      of.y = -o.y + np.y*(1.-fs.x[ii])/2.;
    
      foreach_dimension()
	E.x = (A1.x.x + A1.x.y*of.y
#if dimension > 2
	       + A1.x.z*of.z
#endif
	       + 2.*A2.x*of.x)/Delta;

      force->x += (2.*ii-1)*perm.x[ii]/2.*(sq(E.x)
#if dimension > 2
					   -sq(E.z)
#endif
					   -sq(E.y));
      force->y += (2.*ii-1)*perm.x[ii]*E.x*E.y;
#if dimension > 2
      force->z += (2.*ii-1)*perm.x[ii]*E.x*E.z;
#endif
    }
}
#endif
#endif

/**
We overload the 
[acceleration event](/src/navier-stokes/centered.h#acceleration-term) 
of the Navier--Stokes solver to add the electrohydrodynamics acceleration
term. */

event acceleration (i++) {
  vector f[];
  foreach_dimension() {
    face vector Mx[];
    foreach_face(x)
      Mx.x[] = epsilon.x[]/2.*(sq(face_gradient_x (phi,0))  
                               - sq(face_avg_gradient_t1_y(phi, 0))
#if dimension > 2			       
                               - sq(face_avg_gradient_t1_z(phi, 0))
#endif
			       );
    foreach_face(y)
      Mx.y[] = epsilon.y[]*face_gradient_y (phi, 0)
      *face_avg_gradient_t1_x (phi, 0);

#if dimension > 2
    foreach_face(z)
      Mx.z[] = epsilon.z[]*(face_gradient_z (phi, 0))*
                          *face_avg_gradient_t1_x (phi, 0);
#endif    
    
    /**
    The electric force is the divergence of the Maxwell stress tensor
    $\mathbf{M}$. */

    foreach() {
      double d = 0.;
      foreach_dimension()
	d += Mx.x[1,0] - Mx.x[];
      f.x[] = d/Delta;
    }
  }

  /**
  If case of embed solid the electrical stresses on the solid must be taken
  into account.*/
  
#if EMBED || AXI
  foreach() {
#if EMBED
    if(cs[] > 0 && cs[] < 1.) {
#if COEF      
      coord F;
      electric_force_mixed (point, cs, phi, epsilon, &F);
      foreach_dimension() 
      	f.x[] = F.x/Delta;
#else
      double area = 0.;
      coord grad, n;
      area = bilinear_corrected_embed_gradient (point, phi, &grad, &n);
      double E2 = 0., En = 0., epsilona = 0., fa = 0.;
      foreach_dimension() {
	En += grad.x*n.x;
	E2 += sq(grad.x);
	epsilona += epsilon.x[] + epsilon.x[1];
	fa  += fm.x[] + fm.x[1];
      }
      epsilona /= fa;
      foreach_dimension()
      	f.x[] += epsilona*(En*grad.x-E2/2.*n.x)*area/Delta;
#endif
    }
#endif
  
#if AXI
    double E2 = 0., epsilona = 0., fa = 0.;
    foreach_dimension() {
      E2 += sq(center_gradient(phi));
      epsilona += epsilon.x[] + epsilon.x[1];
      fa += fm.x[] + fm.x[1];
    }
    epsilona /= fa;
    f.y[] += epsilona*E2/2.;
#endif
  }
#endif
    
  /**
  To get the acceleration from the force we need to multiply by the
  specific volume $\alpha$. */

  face vector av = a;
#if EMBED  
  foreach_face()
    av.x[] += (cs[] > 0. && cs[-1] > 0.) ?
    (f.x[]/rho[] + f.x[-1]/rho[-1])/2. :
    (f.x[] + f.x[-1])/(rho[] + rho[-1] + SEPS);
#else  
  foreach_face()
    av.x[] += (f.x[]/rho[] + f.x[-1]/rho[-1])/2.;
#endif
}
