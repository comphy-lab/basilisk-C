/**
# Ohmic conduction flux of charged species

This function computes the fluxes due to ohmic conduction appearing in
the [Nernst--Planck
equation](http://en.wikipedia.org/wiki/Nernst%E2%80%93Planck_equation). The
species charge concentrations are then updated using the explicit
scheme
$$
c^{n+1}_i = c^n_i +\Delta t \, \nabla \cdot( K_i c^n_i \nabla \phi^n)
$$ 
where $c_i$ is the volume density of the $i$-specie, $K_i$ its volume
electric conductivity and $\phi$ the electric potential. */

extern scalar phi;

void ohmic_flux (scalar * c,         // A list of the species concentration...
		 int * z,            // ... and their corresponding valences
		 double dt,
		 vector * K = NULL)  // electric mobility (default the valence)
{
  /**
  If the volume conductivity is not provided it is set to the value of
  the valence. */
  
  if (!K) { // fixme: this does not work yet
    int i = 0;
    for (scalar s in c) {
      const face vector kc[] = {z[i], z[i]}; i++;
      K = vectors_append (K, kc); // fixme: K should be freed eventually
    }
  }

  scalar s;
  (const) face vector k;
  for (s, k in c, K) {

    /**
    The fluxes of each specie through each face due to ohmic transport
    are */

    face vector f[];
    foreach_face()
      f.x[] = k.x[]*(s[] + s[-1])*(phi[] - phi[-1])/(2.*Delta);

    /**
    The specie concentration is updated using the net amount of that
    specie leaving/entering each cell through the face in the interval
    $dt$ */

    foreach()
      foreach_dimension()
        s[] += dt*(f.x[1] - f.x[])/Delta;
  }
}
