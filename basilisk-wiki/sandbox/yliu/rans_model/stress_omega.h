#ifndef STRESS_OMEGA_H
#define STRESS_OMEGA_H
#include "rsm_common.h"
#include "diffusion.h"

/* Implements the Full Reynolds stress model, stress-omega (Wilcox 2006) */
/* Author Y. Liu */
scalar rhok[], rhoe[];
scalar Rxx[], Rxy[], Rxz[], Ryy[], Ryz[], Rzz[];

struct StressOmegaConstants {
    double C1;
    double alpha_v;
    double beta_v;
    double gamma_v;

    double alpha;
    double beta;
    double beta0;
	double alpha_b_star;

	double sigma_k;
	double sigma_w;
	double beta_star;
	double gamma;


	double kMin_;
	double omegaMin_;
	double nutMin_;
	double nutMax_;

	double k_0;
	double omega_0;
	double pe;
};

#define C2 (10.0/19.0)
struct StressOmegaConstants kTurbConstants = {
.C1 = 9.0/5.0,
.alpha_v = (8.0 + C2)/11.0,
.beta_v = (8.0*C2 - 2.0)/11.0,
.gamma_v = (60.0*C2 - 4.0)/55.0,

.alpha = 13.0/25.0,
.beta = 0.0708,
.beta0 = 0.0708,
.alpha_b_star = 1.36,

.sigma_k = 0.6,
.sigma_w = 0.5,
.beta_star = 0.09,

.kMin_ = 1e-15,
.omegaMin_ = 1e-15,
.nutMin_ = 1e-15,
.nutMax_ = 1e3,

.k_0 = 1e-4,
.omega_0 = 420.0
};

event defaults(i=0){
#if TREE
	for (scalar s in {rhok, rhoe, Rxx, Rxy, Rxz, Ryy, Ryz, Rzz}) {
		s.refine  = refine_linear;
		s.restriction = restriction_volume_average;
		s.gradient = minmod2;
	}
	nu_t.refine  = refine_linear;

	nu_t.coarsen=no_restriction;
	foreach_dimension(){
		mu_t.x.coarsen=no_restriction;
	}
#endif
}

void set_rans_init()
{
	foreach() {
        Rxx[] = -2.0/3.0*kTurbConstants.k_0;
        Ryy[] = -2.0/3.0*kTurbConstants.k_0;
        Rxy[] = 0.0;

		Rxz[] = 0.0;
		Ryz[] = 0.0;
		Rzz[] = -2.0/3.0*kTurbConstants.k_0;

		rhoe[] = kTurbConstants.omega_0;
	}
	boundary((scalar*){Rxx, Rxy, Rxz, Ryy, Ryz, Rzz, rhoe});
}

event init (i=0,last) {
if (!restore ("restart")) {
	set_rans_init();
	correct_nut();
	}
}

void correct_nut()
{
	
	foreach()
	{
		// Realizability condition
		Rxx[] = min(Rxx[], 0);
		Ryy[] = min(Ryy[], 0);
		Rzz[] = min(Rzz[], 0);
		Rxy[] = min(max(Rxy[], -sqrt(Rxx[]*Ryy[])), sqrt(Rxx[]*Ryy[]));
		#if dimension > 2
		Ryz[] = min(max(Ryz[], -sqrt(Rzz[]*Ryy[])), sqrt(Rzz[]*Ryy[]));
		Rxz[] = min(max(Rxz[], -sqrt(Rxx[]*Rzz[])), sqrt(Rxx[]*Rzz[]));
		#endif
		rhoe[] = max(rhoe[], kTurbConstants.omegaMin_);

		rhok[] = -0.5*(Rxx[] + Ryy[] + Rzz[]);
		nu_t[] = rhok[]/(rhoe[] + 1e-15);
    }
	boundary({rhok, rhoe, nu_t});
}

#if dimension == 2
void compute_stress_omega_srcs(
	scalar Sk01, scalar Sk02, scalar Sk03,
	scalar Sk04,
	scalar Spk01, scalar Spk02, scalar Spk03,
	scalar Spk04,
scalar Somega, scalar Spomega
)
{
	tensor gradU[]; // grad(u)
	tensor S[], W[], tG[];
	vector dk[], de[];

	t_rsm_centered_gradientv(u, gradU);
	t_rsm_centered_gradient(rhok, dk);
	t_rsm_centered_gradient(rhoe, de);

	foreach(){
		// compute dui_dxi
		double divU = gradU.x.x[] + gradU.y.y[] + 0.0;

		// Compute stress tensor S_{ij}
		double S_tr = (gradU.x.x[] + gradU.y.y[] + 0.0)/3.; // remove trace; divergence free
		S.x.x[] = gradU.x.x[];
		S.y.y[] = gradU.y.y[];
		S.x.y[] = 0.5*(gradU.x.y[] + gradU.y.x[]);

		W.x.y[] = 0.5*(gradU.x.y[] - gradU.y.x[]);

		tG.x.x[] = S.x.x[] - S_tr;
		tG.y.y[] = S.y.y[] - S_tr;
		tG.x.y[] = S.x.y[];
		
		double alphaRho = RHO;

        // Compute R production and transportation
		double Pxx, Pxy, Pyy, Pzz;
		double Yxx, Yxy, Yyy, Yzz;
		double Dxx, Dxy, Dyy, Dzz;
        Pxx = (Rxx[]*gradU.x.x[] + Rxy[]*gradU.x.y[])*2;
        Pxy =  Rxx[]*gradU.y.x[] + Rxy[]*gradU.y.y[]
			 + Rxy[]*gradU.x.x[] + Ryy[]*gradU.x.y[];

        Pyy = (Rxy[]*gradU.y.x[] + Ryy[]*gradU.y.y[])*2;

        Pzz = 0.0;

        Dxx = (Rxx[]*gradU.x.x[] + Rxy[]*gradU.y.x[])*2;
        Dxy =  Rxx[]*gradU.x.y[] + Rxy[]*gradU.y.y[]
			 + Rxy[]*gradU.x.x[] + Ryy[]*gradU.y.x[];

        Dyy = (Rxy[]*gradU.x.y[] + Ryy[]*gradU.y.y[])*2;

        Dzz = 0.0;

        double P = 0.5*(Pxx + Pyy + Pzz);

        Yxx = kTurbConstants.beta_star*kTurbConstants.C1*rhoe[]*(2.0/3.0*rhok[]) - kTurbConstants.alpha_v*(Pxx - 2.0/3.0*P)
                    - kTurbConstants.beta_v*(Dxx - 2.0/3.0*P) - kTurbConstants.gamma_v*rhok[]*(tG.x.x[]);
        Yxy =-kTurbConstants.alpha_v*(Pxy) - kTurbConstants.beta_v*(Dxy) - kTurbConstants.gamma_v*rhok[]*(S.x.y[]);

        Yyy = kTurbConstants.beta_star*kTurbConstants.C1*rhoe[]*(2.0/3.0*rhok[]) - kTurbConstants.alpha_v*(Pyy - 2.0/3.0*P)
                    - kTurbConstants.beta_v*(Dyy - 2.0/3.0*P) - kTurbConstants.gamma_v*rhok[]*(tG.y.y[]);

        Yzz = kTurbConstants.beta_star*kTurbConstants.C1*rhoe[]*(0 + 2.0/3.0*rhok[]) - kTurbConstants.alpha_v*(Pzz - 2.0/3.0*P)
                    - kTurbConstants.beta_v*(Dzz - 2.0/3.0*P);

		double E = 2.0/3.0*kTurbConstants.beta_star*rhoe[]*rhok[];
		double Ecorr = kTurbConstants.beta_star*kTurbConstants.C1*rhoe[];

        // Compute the dissipation term
        // compute R_{ij}duidxj
		double RbyGradU = Rxx[] * gradU.x.x[] + Rxy[] * gradU.x.y[] + \
				          Rxy[] * gradU.y.x[] + Ryy[] * gradU.y.y[];

        double PO = kTurbConstants.alpha*RbyGradU / (rhok[] + 1e-15) * rhoe[];
        double EO = kTurbConstants.beta*rhoe[];

		// Secondary generation of omega, PP_\omega = \sigma_d / \omega \frac{\partial_i}{k} \frac{\partial_i}{\omega}
		// compute dki_dxi domegai_dxi
		// notice the treatment to \nabla\omega/\omega (assume \omega is bounded to a positive value properly) to avoid a explicit division to \omega which could lead to numerical blow-up.
        double dkdw = dk.x[]*de.x[]+dk.y[]*de.y[];
		double sigma_d = dkdw > 0.0? 0.125 : 0.0;
		double PPO = rhoe[] > kTurbConstants.omegaMin_ ? alphaRho * sigma_d * dkdw / rhoe[] : 0.0;

		// Update source terms
		Sk01[] = -Pxx-Yxx + E;
        Sk02[] = -Pxy-Yxy;
        Sk03[] = -Pyy-Yyy + E;
        Sk04[] = -Pzz-Yzz + E;

		Spk01[] = -Ecorr;
        Spk02[] = -Ecorr;
        Spk03[] = -Ecorr;
        Spk04[] = -Ecorr;

		Somega[]  = PO + PPO;
		Spomega[] = -EO;

	}
	boundary({Sk01, Sk02, Sk03, Sk04,
	Spk01, Spk02, Spk03, Spk04, Somega, Spomega});
}
#else
void compute_stress_omega_srcs(
	scalar Sk01, scalar Sk02, scalar Sk03,
	scalar Sk04, scalar Sk05, scalar Sk06,
	scalar Spk01, scalar Spk02, scalar Spk03,
	scalar Spk04, scalar Spk05, scalar Spk06,
scalar Somega, scalar Spomega
)
{
	tensor gradU[]; // grad(u)
	tensor S[], W[], tG[];
	vector dk[], de[];

	t_rsm_centered_gradientv(u, gradU);
	t_rsm_centered_gradient(rhok, dk);
	t_rsm_centered_gradient(rhoe, de);

	foreach(){
		// compute dui_dxi
		double divU = gradU.x.x[] + gradU.y.y[] + gradU.z.z[];

		// Compute stress tensor S_{ij}
		double S_tr = (gradU.x.x[] + gradU.y.y[] + gradU.z.z[])/3.; // remove trace; divergence free
		S.x.x[] = gradU.x.x[];
		S.y.y[] = gradU.y.y[];
		S.z.z[] = gradU.z.z[];
		S.x.y[] = 0.5*(gradU.x.y[] + gradU.y.x[]);
		S.x.z[] = 0.5*(gradU.x.z[] + gradU.z.x[]);
		S.y.z[] = 0.5*(gradU.y.z[] + gradU.z.y[]);

		W.x.y[] = 0.5*(gradU.x.y[] - gradU.y.x[]);
		W.x.z[] = 0.5*(gradU.x.z[] - gradU.z.x[]);
		W.y.z[] = 0.5*(gradU.y.z[] - gradU.z.y[]);

		tG.x.x[] = S.x.x[] - S_tr;
		tG.y.y[] = S.y.y[] - S_tr;
		tG.z.z[] = S.z.z[] - S_tr;
		tG.x.y[] = S.x.y[];
		tG.x.z[] = S.x.z[];
		tG.y.z[] = S.y.z[];
		
		double alphaRho = RHO;

        // Compute R production and transportation
		double Pxx, Pxy, Pxz, Pyy, Pyz, Pzz;
		double Yxx, Yxy, Yxz, Yyy, Yyz, Yzz;
		double Dxx, Dxy, Dxz, Dyy, Dyz, Dzz;
        Pxx = (Rxx[]*gradU.x.x[] + Rxy[]*gradU.x.y[] + Rxz[]*gradU.x.z[])*2;
        Pxy =  Rxx[]*gradU.y.x[] + Rxy[]*gradU.y.y[] + Rxz[]*gradU.y.z[]
			 + Rxy[]*gradU.x.x[] + Ryy[]*gradU.x.y[] + Ryz[]*gradU.x.z[];
        Pxz =  Rxx[]*gradU.z.x[] + Rxy[]*gradU.z.y[] + Rxz[]*gradU.z.z[]
			 + Rxz[]*gradU.x.x[] + Ryz[]*gradU.x.y[] + Rzz[]*gradU.x.z[];
        Pyy = (Rxy[]*gradU.y.x[] + Ryy[]*gradU.y.y[] + Ryz[]*gradU.y.z[])*2;
        Pyz =  Rxy[]*gradU.z.x[] + Ryy[]*gradU.z.y[] + Ryz[]*gradU.z.z[]
			 + Rxz[]*gradU.y.x[] + Ryz[]*gradU.y.y[] + Rzz[]*gradU.y.z[];
        Pzz = (Rxz[]*gradU.z.x[] + Ryz[]*gradU.z.y[] + Rzz[]*gradU.z.z[])*2;

        Dxx = (Rxx[]*gradU.x.x[] + Rxy[]*gradU.y.x[] + Rxz[]*gradU.z.x[])*2;
        Dxy =  Rxx[]*gradU.x.y[] + Rxy[]*gradU.y.y[] + Rxz[]*gradU.z.y[]
			 + Rxy[]*gradU.x.x[] + Ryy[]*gradU.y.x[] + Ryz[]*gradU.z.x[];
        Dxz =  Rxx[]*gradU.x.z[] + Rxy[]*gradU.y.z[] + Rxz[]*gradU.z.z[]
			 + Rxz[]*gradU.x.x[] + Ryz[]*gradU.y.x[] + Rzz[]*gradU.z.x[];
        Dyy = (Rxy[]*gradU.x.y[] + Ryy[]*gradU.y.y[] + Ryz[]*gradU.z.y[])*2;
        Dyz =  Rxy[]*gradU.x.z[] + Ryy[]*gradU.y.z[] + Ryz[]*gradU.z.z[]
			 + Rxz[]*gradU.x.y[] + Ryz[]*gradU.y.y[] + Rzz[]*gradU.z.y[];
        Dzz = (Rxz[]*gradU.x.z[] + Ryz[]*gradU.y.z[] + Rzz[]*gradU.z.z[])*2;

        double P = 0.5*(Pxx + Pyy + Pzz);

        Yxx = kTurbConstants.beta_star*kTurbConstants.C1*rhoe[]*(2.0/3.0*rhok[]) - kTurbConstants.alpha_v*(Pxx - 2.0/3.0*P)
                    - kTurbConstants.beta_v*(Dxx - 2.0/3.0*P) - kTurbConstants.gamma_v*rhok[]*(tG.x.x[]);
        Yxy =-kTurbConstants.alpha_v*(Pxy) - kTurbConstants.beta_v*(Dxy) - kTurbConstants.gamma_v*rhok[]*(S.x.y[]);
        Yxz =-kTurbConstants.alpha_v*(Pxz) - kTurbConstants.beta_v*(Dxz) - kTurbConstants.gamma_v*rhok[]*(S.x.z[]);
        Yyy = kTurbConstants.beta_star*kTurbConstants.C1*rhoe[]*(2.0/3.0*rhok[]) - kTurbConstants.alpha_v*(Pyy - 2.0/3.0*P)
                    - kTurbConstants.beta_v*(Dyy - 2.0/3.0*P) - kTurbConstants.gamma_v*rhok[]*(tG.y.y[]);
        Yyz =-kTurbConstants.alpha_v*(Pyz) - kTurbConstants.beta_v*(Dyz) - kTurbConstants.gamma_v*rhok[]*(S.y.z[]);
        Yzz = kTurbConstants.beta_star*kTurbConstants.C1*rhoe[]*(2.0/3.0*rhok[]) - kTurbConstants.alpha_v*(Pzz - 2.0/3.0*P)
                    - kTurbConstants.beta_v*(Dzz - 2.0/3.0*P) - kTurbConstants.gamma_v*rhok[]*(tG.z.z[]);

		double E = 2.0/3.0*kTurbConstants.beta_star*rhoe[]*rhok[];
		double Ecorr = kTurbConstants.beta_star*kTurbConstants.C1*rhoe[];

        // Compute the dissipation term
        // compute R_{ij}duidxj
		double RbyGradU = Rxx[] * gradU.x.x[] + Rxy[] * gradU.x.y[] + Rxz[] * gradU.x.z[] + \
				          Rxy[] * gradU.y.x[] + Ryy[] * gradU.y.y[] + Ryz[] * gradU.y.z[] + \
                          Rxz[] * gradU.z.x[] + Ryz[] * gradU.z.y[] + Rzz[] * gradU.z.z[] ;

        double PO = kTurbConstants.alpha*RbyGradU / (rhok[] + 1e-15) * rhoe[];
        double EO = kTurbConstants.beta*rhoe[];

		// Secondary generation of omega, PP_\omega = \sigma_d / \omega \frac{\partial_i}{k} \frac{\partial_i}{\omega}
		// compute dki_dxi domegai_dxi
		// notice the treatment to \nabla\omega/\omega (assume \omega is bounded to a positive value properly) to avoid a explicit division to \omega which could lead to numerical blow-up.
        double dkdw = dk.x[]*de.x[]+dk.y[]*de.y[]+dk.z[]*de.z[];
		double sigma_d = dkdw > 0.0? 0.125 : 0.0;
		double PPO = rhoe[] > kTurbConstants.omegaMin_ ? alphaRho * sigma_d * dkdw / rhoe[] : 0.0;

		// Update source terms
		Sk01[] = -Pxx-Yxx + E;
        Sk02[] = -Pxy-Yxy;
        Sk03[] = -Pxz-Yxz;
        Sk04[] = -Pyy-Yyy + E;
        Sk05[] = -Pyz-Yyz;
        Sk06[] = -Pzz-Yzz + E;

		Spk01[] = -Ecorr;
        Spk02[] = -Ecorr;
        Spk03[] = -Ecorr;
        Spk04[] = -Ecorr;
        Spk05[] = -Ecorr;
        Spk06[] = -Ecorr;

		Somega[]  = PO + PPO;
		Spomega[] = -EO;

	}
	boundary({Sk01, Sk02, Sk03, Sk04, Sk05, Sk06,
	Spk01, Spk02, Spk03, Spk04, Spk05, Spk06, Somega, Spomega});
}
#endif


void compute_komega_mueff(face vector Dk, face vector Domega)
{
	foreach_face() {
		double mu_f = MU;
		double mut_f = 0.5 * (nu_t[]+nu_t[-1]);
		Dk.x[] = fm.x[] * (mu_f + kTurbConstants.sigma_k*mut_f);
		Domega.x[] = fm.x[] * (mu_f + kTurbConstants.sigma_w*mut_f);
	}
	boundary({Dk,Domega});
}

#if dimension == 2
event turbulence_correction (i++) {
	// incorporate turbulence dissipation in R and omega equations
	face vector D_k[], D_omega[];
	compute_komega_mueff(D_k, D_omega);
	scalar Sk01[], Sk02[], Sk03[], Sk04[];
	scalar Spk01[], Spk02[], Spk03[], Spk04[];
    scalar src_omega[], sp_omega[];
	compute_stress_omega_srcs(Sk01, Sk02, Sk03, Sk04,
	Spk01, Spk02, Spk03, Spk04, src_omega, sp_omega);

	// Advection
	advection ((scalar *){Rxx, Rxy, Ryy, Rzz, rhoe}, uf, dt);
		
	// Diffusion
	diffusion (rhoe, dt, D_omega, beta=sp_omega, r=src_omega);

    scalar tmpR, tmpSU, tmpSP;
	scalar *Rlist = (scalar *){Rxx, Rxy, Ryy, Rzz};
	scalar *SUlist = (scalar *){Sk01, Sk02, Sk03, Sk04};
	scalar *SPlist = (scalar *){Spk01, Spk02, Spk03, Spk04};

    for (tmpR, tmpSU, tmpSP in Rlist, SUlist, SPlist)
        diffusion (tmpR, dt, D_k, beta=tmpSP , r=tmpSU);

	// Compute nu_t
	correct_nut();
}
#else
event turbulence_correction (i++) {
	// incorporate turbulence dissipation in R and omega equations
	face vector D_k[], D_omega[];
	compute_komega_mueff(D_k, D_omega);
	scalar Sk01[], Sk02[], Sk03[], Sk04[], Sk05[], Sk06[];
	scalar Spk01[], Spk02[], Spk03[], Spk04[], Spk05[], Spk06[];
    scalar src_omega[], sp_omega[];
	compute_stress_omega_srcs(Sk01, Sk02, Sk03, Sk04, Sk05, Sk06,
	Spk01, Spk02, Spk03, Spk04, Spk05, Spk06, src_omega, sp_omega);

	// Advection
	advection ((scalar *){Rxx, Rxy, Rxz, Ryy, Ryz, Rzz, rhoe}, uf, dt);
		
	// Diffusion
	mgOmega = diffusion (rhoe, dt, D_omega, beta=sp_omega, r=src_omega);

    scalar tmpR, tmpSU, tmpSP;
	scalar *Rlist = (scalar *){Rxx, Rxy, Rxz, Ryy, Ryz, Rzz};
	scalar *SUlist = (scalar *){Sk01, Sk02, Sk03, Sk04, Sk05, Sk06};
	scalar *SPlist = (scalar *){Spk01, Spk02, Spk03, Spk04, Spk05, Spk06};

    for (tmpR, tmpSU, tmpSP in Rlist, SUlist, SPlist)
        diffusion (tmpR, dt, D_k, beta=tmpSP , r=tmpSU);

	// Compute nu_t
	correct_nut();
}
#endif

#if TREE
// mesh adaptation happens before solving viscous term
event adapt (i++) {
	correct_nut();
}
#endif

#endif

