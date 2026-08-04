#include "poisson.h"
#include <complex.h>

double complex om;//Complex growth rate (Guess)

double rho1, rho2;
double mu1, mu2;
double g, h;

# define rho(f) (clamp(f,0.,1.)*(rho1 - rho2) + rho2)
# define mu(f)  (clamp(f,0.,1.)*(mu1 - mu2) + mu2)


struct PoissonSystem {
  scalar REPsi, RErhsPsi;
  scalar IMPsi, IMrhsPsi;
  scalar rho;
  (const) face vector alpha1;
  (const) face vector alpha2;
  (const) scalar lambda1;
  (const) scalar lambda2;
  (const) scalar rho0;//+
  double tolerance;
  int nrelax, minlevel;
  scalar * res;
};


/*  Computation of the growth rate using the Kinematik condition  */
double complex get_omega (double RePsi0, double IMPsi0, double k) {

	double complex omega_sq=( RePsi0 + I*IMPsi0)*sq(k);
	double omega_R = sqrt(  (  sqrt(  sq(creal(omega_sq)) + sq(cimag(omega_sq))  ) + creal(omega_sq)   )/2  );
	double omega_I = sqrt(  (  sqrt(  sq(creal(omega_sq)) + sq(cimag(omega_sq))  ) - creal(omega_sq)   )/2  );

	return omega_R + I*omega_I;
}

static void relax_system (scalar * al, scalar * bl, int l, void * data)
{
  scalar REPsi   = al[0], RErhsPsi   = bl[0];
  scalar IMPsi   = al[1], IMrhsPsi   = bl[1];
  scalar REOmega = al[2], RErhsOmega = bl[2];
  scalar IMOmega = al[3], IMrhsOmega = bl[3];

  struct PoissonSystem * p = (struct PoissonSystem *) data;
  scalar rho = p->rho;
  (const) face vector alpha1 = p->alpha1;
  (const) face vector alpha2 = p->alpha2;
  (const) scalar lambda1 = p->lambda1;
  (const) scalar lambda2 = p->lambda2;
  (const) scalar rho0 = p->rho0;//+

  double complex numPsi, PsiR, PsiL, Psi, Psirhs;
  double complex numOmega, OmegaR, OmegaL, Omega, Omegarhs, Omegaspe;

#if JACOBI
  scalar cREPsi[],cIMPsi[], cREOmega[],cIMOmega[];
#else
  scalar cREPsi = REPsi, cIMPsi=IMPsi, cREOmega = REOmega, cIMOmega = IMOmega;
#endif

  


  foreach_level_or_leaf (l) {
    Psi    = REPsi[] + I*IMPsi[];
    Psirhs = RErhsPsi[] + I*IMrhsPsi[];

    Omega    = REOmega[] + I*IMOmega[];
    Omegarhs = RErhsOmega[] + I*IMrhsOmega[];

    numPsi = lambda2[]*Omega*sq(Delta) -sq(Delta)*Psirhs;// + lambda1[]*sq(Delta)*Omega;//Dernier terme de couplage attention

    double denPsi = -lambda1[]*sq(Delta);

    numOmega = -sq(Delta)*Omegarhs;// + Omegarhs*Psi*Delta;//A verifier si le signe est niquel pour le terme gradient(2nd terme)
    double complex denOmega = -lambda2[]*sq(Delta) - I*om*rho0[]*sq(Delta);//Premiere Erreur /sq(k) (Je vais faire une petite modif sachant que lambda1[] = - rho0[]/sq(k))
    

    foreach_dimension() {
	    
      PsiR = REPsi[1] + I*IMPsi[1] ;
      PsiL = REPsi[-1] + I*IMPsi[-1] ; 

      OmegaR = REOmega[1] + I*IMOmega[1] ;
      OmegaL = REOmega[-1] + I*IMOmega[-1] ; 
        


      Omegaspe = (rho[1] - rho[-1])/Delta/2;
      //Omegaspe = (alpha1.x[1] - alpha1.x[])/Delta;

      numPsi += alpha1.x[1]*PsiR + alpha1.x[]*PsiL + alpha2.x[1]*(OmegaR - Omega) - alpha2.x[]*(Omega - OmegaL);
      denPsi += alpha1.x[1] + alpha1.x[];

      numOmega += alpha2.x[1]*OmegaR + alpha2.x[]*OmegaL+ Omegaspe*(PsiR -PsiL) * Delta/2; 
      denOmega += alpha2.x[1] + alpha2.x[] ;

    }


   
      Psi = numPsi/denPsi;
      Omega = numOmega/denOmega;
 

      cREPsi[] = creal(Psi);
      cIMPsi[] = cimag(Psi);

      cREOmega[] = creal(Omega);
      cIMOmega[] = cimag(Omega);

  }

#if JACOBI
  foreach_level_or_leaf (l) {
    REPsi[] = (REPsi[] + 2.*cREPsi[])/3.;
    IMPsi[] = (IMPsi[] + 2.*cIMPsi[])/3.;
    REomega[] = (REOmega[] + 2.*cREOmega[])/3.;
    IMomega[] = (IMOmega[] + 2.*cIMOmega[])/3.;
  }
#endif

}

/**
The equivalent residual function is obtained in a similar way in the
case of a Cartesian grid, however the case of the tree mesh
requires more careful consideration... */

static double residual_system (scalar * al, scalar * bl, scalar * resl, void * data)
{
  scalar REPsi = al[0],   RErhsPsi = bl[0],   REres = resl[0];
  scalar IMPsi = al[1],   IMrhsPsi = bl[1],   IMres = resl[1];
  scalar REOmega = al[2], RErhsOmega = bl[2], REresOmega = resl[2];
  scalar IMOmega = al[3], IMrhsOmega = bl[3], IMresOmega = resl[3];


  struct PoissonSystem * p = (struct PoissonSystem *) data;
  scalar rho = p->rho;
  (const) face vector alpha1 = p->alpha1;
  (const) face vector alpha2 = p->alpha2;
  (const) scalar lambda1 = p->lambda1;
  (const) scalar lambda2 = p->lambda2;
  (const) scalar rho0 = p->rho0;//+

  double maxres = 0.;
  iter=0;
  double complex resPsi, Psi , PsiR ,PsiL , Psirhs;
  double complex resOmega, Omega, OmegaR , OmegaL, Omegarhs,Omegaspe;
  
#if TREE
  /* conservative coarse/fine discretisation (2nd order) */
  face vector REgradPsi[] ,IMgradPsi[];
  face vector REgradOmega[] ,IMgradOmega[];

  foreach_face() {
    REgradPsi.x[] = alpha1.x[]*face_gradient_x (REPsi, 0);
    IMgradPsi.x[] = alpha1.x[]*face_gradient_x (IMPsi, 0);

    REgradOmega.x[] = alpha2.x[]*face_gradient_x (REOmega, 0);
    IMgradOmega.x[] = alpha2.x[]*face_gradient_x (IMOmega, 0);
  }
  foreach (reduction(max:maxres), nowarning) {
    Psi    = REPsi[] + I*IMPsi[];
    Psirhs = RErhsPsi[] + I*IMrhsPsi[];

    Omega    = REOmega[] + I*IMOmega[];
    Omegarhs = RErhsOmega[] + I*IMrhsOmega[];

    

    resPsi = Psirhs - lambda1[]*Psi - lambda2[] * Omega;

    resOmega = Omegarhs - lambda2[]*Omega - I*om*rho0[]*Omega; // Deuxieme Erreur il manque toujours la /sq(k) (Meme modif que pour relax system)


    foreach_dimension() {

      PsiR = REPsi[1] + I*IMPsi[1] ;
      PsiL = REPsi[-1] + I*IMPsi[-1] ; 

      OmegaR = REOmega[1] + I*IMOmega[1] ;
      OmegaL = REOmega[-1] + I*IMOmega[-1] ; 
      

      Omegaspe = (rho[1] - rho[-1])/Delta/2;
      //Omegaspe = (alpha1.x[1] - alpha1.x[])/Delta;


      resPsi -= (REgradPsi.x[1] - REgradPsi.x[] + I*(IMgradPsi.x[1] - IMgradPsi.x[]))/Delta + (REgradOmega.x[1] - REgradOmega.x[] + I*(IMgradOmega.x[1] - IMgradOmega.x[]))/Delta;

      resOmega -= (REgradOmega.x[1] - REgradOmega.x[] + I*(IMgradOmega.x[1] - IMgradOmega.x[]))/Delta + Omegaspe*(PsiR - PsiL)/2/Delta;
    }

    REres[] = creal(resPsi);
    IMres[] = cimag(resPsi);

    REresOmega[] = creal(resOmega);
    IMresOmega[] = cimag(resOmega);



  }
    
    if (cabs (resPsi) > maxres)
      maxres = cabs (resPsi);
    
    if (cabs (resOmega) > maxres)
      maxres = cabs (resOmega);
#else // !TREE
  /* "naive" discretisation (only 1st order on trees) */
  foreach (reduction(max:maxres), nowarning) {
    Psi    = REPsi[] + I*IMPsi[];
    Psirhs = RErhsPsi[] + I*IMrhsPsi[];

    Omega    = REOmega[] + I*IMOmega[];
    Omegarhs = RErhsOmega[] + I*IMrhsOmega[];

    resPsi = Psirhs - lambda1[]*Psi - lambda2[]*Omega;
    resOmega = Omegarhs - lambda2[]*Omega + I*om*lambda1[]*Omega;// - Omegarhs * Psi /Delta;

    foreach_dimension() {

      
      PsiR = REPsi[1] + I*IMPsi[1] ;
      PsiL = REPsi[-1] + I*IMPsi[-1] ; 

      OmegaR = REOmega[1] + I*IMOmega[1] ;
      OmegaL = REOmega[-1] + I*IMOmega[-1] ; 

      Omegaspe = (rho[1] - rho[-1])/Delta/2;

      resPsi += (alpha1.x[0]*face_gradient_x (REPsi, 0) -
		alpha1.x[1]*face_gradient_x (REPsi, 1))/Delta
		+  I*(alpha1.x[0]*face_gradient_x (IMPsi, 0) -
		alpha1.x[1]*face_gradient_x (IMPsi, 1))/Delta + (alpha2.x[0]*face_gradient_x (REOmega, 0) -
		alpha2.x[1]*face_gradient_x (REOmega, 1))/Delta
		+  I*(alpha2.x[0]*face_gradient_x (IMOmega, 0) -
		alpha2.x[1]*face_gradient_x (IMOmega, 1))/Delta;  

      resOmega += (alpha2.x[0]*face_gradient_x (REOmega, 0) -
		alpha2.x[1]*face_gradient_x (REOmega, 1))/Delta
		+  I*(alpha2.x[0]*face_gradient_x (IMOmega, 0) -
		alpha2.x[1]*face_gradient_x (IMOmega, 1))/Delta -  Omegaspe*(PsiR - PsiL)/2/Delta;//Omegarhs * //+ Omegarhs * PsiL /Delta; 

    }

    REres[] = creal(resPsi);
    IMres[] = cimag(resPsi);

    REresOmega[] = creal(resOmega);
    IMresOmega[] = cimag(resOmega);



    if (cabs(resPsi) > maxres)
      maxres = cabs (resPsi);

    if (cabs (resOmega) > maxres)
      maxres = cabs (resOmega);

  }
#endif // !TREE
  return maxres;
}

double complex solve_stability_RT (
		scalar REPsi, scalar IMPsi, 
		scalar REOmega, scalar IMOmega, 
		scalar f,
		double rho1, double rho2, double sigma, double g, double k,
		double tolerance = 0.,
		int nrelax = 4,
		int minlevel = 0,
		scalar * res = NULL,
		double (* flux) (Point, scalar, vector, double *) = NULL)
{


	scalar rho0[];

	face vector rhof[];
	face vector muf[];

	scalar lambda1[];
	scalar lambda2[];

	scalar RErhsOmega[] ;
	scalar IMrhsOmega[] ;
	scalar RErhsPsi[];
  	scalar IMrhsPsi[];

	om = om/sq(k);

	/* Definition of the parametres of the PDE */
	foreach(){
	        rho0[] = rho(f[]);
		lambda1[] = -rho(f[])*sq(k);
		lambda2[] = -mu(f[])*sq(k);
		double ftop = f[0,1];
		double fbot = f[0,-1];
		RErhsOmega[] = (rho(ftop)-rho(fbot))/Delta/2;
		IMrhsOmega[] =0;
    		RErhsPsi[] =((rho(ftop)-rho(fbot))/Delta/2) * (1. - sigma*sq(k)/(rho1-rho2))*g*h;//Nouveau dimentionnement *g*h ;
    		IMrhsPsi[] = 0;
	}


	foreach_face(){
		double ff = (f[-1] + f[])/2.;
		rhof.x[] = rho(ff);
		muf.x[]  = mu(ff);
	}

	restriction ({rho0, rhof, muf, lambda1, lambda2});

	double defaultol = TOLERANCE;
	if (tolerance)
		TOLERANCE = tolerance;


	struct PoissonSystem p = {REPsi, RErhsPsi, IMPsi, IMrhsPsi, rho0, rhof, muf, lambda1, lambda2,rho0, tolerance, nrelax, minlevel, res };//Warning utilisation de rho0 2 fois, trouver comment regler ça

#if EMBED
	if (!flux && a.boundary[embed] != symmetry)
		p.embed_flux = embed_flux;
	else
		p.embed_flux = flux;
#endif // EMBED
	mgstats s = mg_solve ({REPsi,IMPsi,REOmega, IMOmega}, {RErhsPsi,IMrhsPsi,RErhsOmega, IMrhsOmega}, residual_system, relax_system, &p,
			nrelax, res, max(1, minlevel));
	if (tolerance)
		TOLERANCE = defaultol;

	om= get_omega (interpolate (REPsi, 0, 0), interpolate (IMPsi, 0, 0), k);

	return om;

}