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
We overload the 
[acceleration event](/src/navier-stokes/centered.h#acceleration-term) 
of the Navier--Stokes solver to add the electrohydrodynamics acceleration
term. */

event acceleration (i++) {
  assert (dimension <= 2); // not 3D yet
  vector f[];
  foreach_dimension() {
    face vector Mx[];
    foreach_face(x)
      Mx.x[] = epsilon.x[]/2.*(sq(phi[] - phi[-1,0]) - 
                               sq(phi[0,1] - phi[0,-1] + 
                                  phi[-1,1] - phi[-1,-1])/16.)/sq(Delta);
    foreach_face(y)
      Mx.y[] = epsilon.y[]*(phi[] - phi[0,-1])*
      (phi[1,0] - phi[-1,0] + 
       phi[1,-1] - phi[-1,-1])/sq(2.*Delta);

    /**
    The electric force is the divergence of the Maxwell stress tensor
    $\mathbf{M}$. */

    foreach()
      f.x[] = (Mx.x[1,0] - Mx.x[] + Mx.y[0,1] - Mx.y[])/(Delta*cm[]);
  }

  /**
  If [axisymmetric cylindrical coordinates](/src/axi.h) are used an 
  additional term must be added. */

#if AXI
  foreach()
    f.y[] += (sq(phi[1,0] - phi[-1,0]) +
	      sq(phi[0,1] - phi[0,-1]))/(8.*cm[]*sq(Delta))
    *(epsilon.x[]/fm.x[] + epsilon.y[]/fm.y[] +
      epsilon.x[1,0]/fm.x[1,0] + epsilon.y[0,1]/fm.y[0,1])/4.;
#endif

  /**
  To get the acceleration from the force we need to multiply by the
  specific volume $\alpha$. */

  face vector av = a;
  foreach_face()
    av.x[] += alpha.x[]/fm.x[]*(f.x[] + f.x[-1])/2.;
}
