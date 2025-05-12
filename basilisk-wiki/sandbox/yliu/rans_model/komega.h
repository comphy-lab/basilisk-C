#ifndef KOMEGA_H
#define KOMEGA_H
#include "rsm_common.h"
#include "diffusion.h"

scalar rhok[], rhoe[];

struct KOmegaConstants {
	double sigma_k;
	double sigma_w;
	double beta_star;
	double gamma;
	double beta;
	double Clim;
	double Clim_star;
	double alpha_b_star;
	double kMin_;
	double omegaMin_;
	double nutMin_;
	double nutMax_;

	double k_0;
	double omega_0;
	double pe;
};

struct KOmegaConstants kTurbConstants = {
.sigma_k = 0.6,
.sigma_w = 0.5,
.beta_star = 0.09,
.gamma = 13.0/25.0,
.beta = 0.0708,
.Clim = 7.0/8.0,
.Clim_star = 7.0/8.0/0.3,
.alpha_b_star = 1.36,
.kMin_ = 1e-15,
.omegaMin_ = 1e-15,
.nutMin_ = 1e-15,
.nutMax_ = 1e3,

.k_0 = 1e-4,
.omega_0 = 420.0
};

event defaults(i=0){
#if TREE
	rhok.refine  = refine_linear;
	rhoe.refine  = refine_linear;

	rhok.restriction = restriction_volume_average;
	rhoe.restriction = restriction_volume_average;

	rhok.gradient = minmod2;
	rhoe.gradient = minmod2;

	nu_t.refine  = refine_linear;
	nu_t.coarsen = nu_t.restriction = restriction_volume_average;

	foreach_dimension(){
		mu_t.x.coarsen=no_restriction;
	}
#endif
}

void set_rans_init()
{
	foreach() {
		rhok[] = kTurbConstants.k_0;
		rhoe[] = kTurbConstants.omega_0;
		nu_t[] = kTurbConstants.k_0/kTurbConstants.omega_0;
	}
	boundary((scalar*){rhok, rhoe, nu_t});
}

event init (i=0,last) {
if (!restore ("restart"))
	set_rans_init();
}

void bound_komega()
{
	foreach()
	{
	  rhok[] = max(rhok[], kTurbConstants.kMin_);
	  rhoe[] = max(rhoe[], kTurbConstants.omegaMin_);
	}
	boundary({rhok, rhoe});
}

void correct_nut_bounded()
{
	tensor gradU[];
	tensor S[],tG[];
	foreach()
	{
		gradU.x.x[] = (u.x[1, 0, 0]-u.x[-1, 0, 0])/2.0/Delta;
	    gradU.x.y[] = (u.x[0, 1, 0]-u.x[0, -1, 0])/2.0/Delta;
	    gradU.y.x[] = (u.y[1, 0, 0]-u.y[-1, 0, 0])/2.0/Delta;
	    gradU.y.y[] = (u.y[0, 1, 0]-u.y[0, -1, 0])/2.0/Delta;

		// Compute stress tensor S_{ij}
		double S_tr = (gradU.x.x[] + gradU.y.y[])/3.; // remove trace; divergence free
		S.x.x[] = gradU.x.x[];
		S.y.y[] = gradU.y.y[];
		S.x.y[] = 0.5*(gradU.x.y[] + gradU.y.x[]);
		tG.x.x[] = S.x.x[] - S_tr;
		tG.y.y[] = S.y.y[] - S_tr;
		tG.x.y[] = S.x.y[];

		// update |\hat{S}| (to bound omega)
		double Gnorm = sqrt(2*(tG.x.x[]*tG.x.x[] + tG.y.y[]*tG.y.y[] + 2*(tG.x.y[]*tG.x.y[])));
		
		// Wilcox 2006
		double omega_star = kTurbConstants.Clim_star * Gnorm;
		double rhoeBounded = max(rhoe[], omega_star);

		nu_t[] = rhok[]/(rhoeBounded+1e-15);

	}
	boundary({nu_t});
}


void compute_komega_srcs(scalar Sk, scalar Somega, scalar Spk, scalar Spomega)
{
	tensor gradU[]; // grad(u)
	tensor S[],W[],tG[];
	vector dk[], de[];

	// Compute gradU and S_{ij}
	foreach(){
	    gradU.x.x[] = (u.x[1, 0, 0]-u.x[-1, 0, 0])/2.0/Delta;
	    gradU.x.y[] = (u.x[0, 1, 0]-u.x[0, -1, 0])/2.0/Delta;
	    gradU.y.x[] = (u.y[1, 0, 0]-u.y[-1, 0, 0])/2.0/Delta;
	    gradU.y.y[] = (u.y[0, 1, 0]-u.y[0, -1, 0])/2.0/Delta;

		// Compute stress tensor S_{ij}, W_{ij}, G_{ij}
		double S_tr = (gradU.x.x[] + gradU.y.y[])/3.; // remove trace; divergence free
		S.x.x[] = gradU.x.x[];
		S.y.y[] = gradU.y.y[];
		S.x.y[] = 0.5*(gradU.x.y[] + gradU.y.x[]);
		W.x.y[] = 0.5*(gradU.x.y[] - gradU.y.x[]);
		tG.x.x[] = S.x.x[] - S_tr;
		tG.y.y[] = S.y.y[] - S_tr;
		tG.x.y[] = S.x.y[];

		// update |S| (to bound omega)
		double Gnorm = sqrt(2*(tG.x.x[]*tG.x.x[] + tG.y.y[]*tG.y.y[] + 2*(tG.x.y[]*tG.x.y[])));
		
		double alphaRho = RHO;

		// compute G_{ij}duidxj
		double GbyGradU = 2.0 * ( \
				tG.x.x[] * gradU.x.x[] + tG.x.y[] * gradU.x.y[] + \
				tG.x.y[] * gradU.y.x[] + tG.y.y[] * gradU.y.y[] );
			
		// restrict unbounded growth of tke in near potential flow.(Larsen 2018)
		// Bound oemga by stability criterion
		double omega_star = 1e-15 + kTurbConstants.Clim_star * Gnorm; // Wilcox 2006
		double rhoeBounded = max(rhoe[], omega_star);

		// Production of k, P = 2.0 \nu_t S_{ij} \frac{\partial_i}{u_j}
		double P = alphaRho * GbyGradU * nu_t[];

		// Dissipation of k, E = \beta_\star rho k \omega
		double E  = alphaRho * kTurbConstants.beta_star * rhoe[];

		// Production of omega, P_\omega = \gamma \omega / \k  P
		double PO = alphaRho * kTurbConstants.gamma * (rhoe[] / rhoeBounded) *  GbyGradU;

		// Dissipation of omega, E_\omega = \beta rho \omega^2
		double EO = alphaRho * kTurbConstants.beta * rhoe[];

		// Secondary generation of omega, PP_\omega = \sigma_d / \omega \frac{\partial_i}{k} \frac{\partial_i}{\omega}
		// compute dki_dxi domegai_dxi
		// notice the treatment to \nabla\omega/\omega (assume \omega is bounded to a positive value properly) to avoid a explicit division to \omega which could lead to numerical blow-up.
	    	dk.x[] = (rhok[1, 0, 0]-rhok[-1, 0, 0])/2.0/Delta;
	    	dk.y[] = (rhok[0, 1, 0]-rhok[0, -1, 0])/2.0/Delta;
	    	de.x[] = (rhoe[1, 0, 0]-rhoe[-1, 0, 0])/2.0/Delta;
	    	de.y[] = (rhoe[0, 1, 0]-rhoe[0, -1, 0])/2.0/Delta;
		double dkdw = dk.x[]*de.x[]+dk.y[]*de.y[];
		double sigma_d = dkdw > 0.0? 0.125 : 0.0;
		double PPO = rhoe[] > kTurbConstants.omegaMin_ ? alphaRho * sigma_d * dkdw / rhoe[] : 0.0;

		// Update source terms
		Sk[] = P;
		Spk[]=-E;
		Somega[]  = PO + PPO;
		Spomega[] =-EO;

	}
	boundary({Sk, Somega, Spk, Spomega});
}

void compute_komega_mueff(face vector Dk, face vector Domega)
{
	foreach_face() {
		double mu_f = MU;
		double mut_f = 0.5 * RHO * (nu_t[]+nu_t[-1]);
		Dk.x[] = fm.x[] * (mu_f + kTurbConstants.sigma_k*mut_f);
		Domega.x[] = fm.x[] * (mu_f + kTurbConstants.sigma_w*mut_f);
	}
	boundary({Dk,Domega});
}

event solve_rans_equations (i++) {
	// Solve k and omega equations
	face vector D_k[], D_omega[];
	compute_komega_mueff(D_k, D_omega);
	scalar src_k[], src_omega[], sp_k[], sp_omega[];
	compute_komega_srcs(src_k, src_omega, sp_k, sp_omega);

	// Advection
	advection ((scalar *){rhok, rhoe}, uf, dt);
	bound_komega();

	// Diffusion
	diffusion (rhoe, dt, D_omega, beta=sp_omega, r=src_omega);
	diffusion (rhok, dt, D_k    , beta=sp_k    , r=src_k);

	// Compute nu_t
	bound_komega();
	correct_nut_bounded();
}

#if TREE
// mesh adaptation happens before solving viscous term
event adapt (i++) {
	correct_nut_bounded();
}
#endif

#endif