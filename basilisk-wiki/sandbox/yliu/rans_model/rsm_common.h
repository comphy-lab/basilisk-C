#ifndef RSM_COMMON_H
#define RSM_COMMON_H

face vector mu_t[]; // turbulent viscosity
scalar nu_t[]; // Cell Centered diffusivity

// Abstract interface to implement
void correct_nut();
void set_rans_init();
void rsm_model_output(scalar**);

// Utility functions
// note: u is weighted by fm
void t_rsm_centered_gradient (scalar p, vector g)
{

  /**
  We first compute a face field $\mathbf{g}_f$. */

  face vector gf[];
  foreach_face()
    gf.x[] = fm.x[]*(p[] - p[-1])/Delta;

  /**
  We average these face values to obtain the centered gradient field. */

  trash ({g});
  foreach()
    foreach_dimension()
      g.x[] = (gf.x[] + gf.x[1])/(fm.x[] + fm.x[1] + SEPS);
}

void t_rsm_centered_gradientv (vector p, tensor g){
	foreach_dimension()
		t_rsm_centered_gradient (p.x, g.x);
}

// overload viscosity event
event viscous_term(i++) {
	// add nu_t to viscosity term
	foreach_face() {
		double mu_f = MU;
		double mut_f = 0.5 * RHO * (nu_t[]+nu_t[-1]);

		mu_t.x[] = mu_f + mut_f;
	}
	boundary((scalar *){mu_t});

	// add to dynamic viscosity
	face vector muv = mu;
	foreach_face()
	  muv.x[] = fm.x[] * mu_t.x[];
}

#endif