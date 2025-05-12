#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"


/** Physicals parameters */
double Ranf;
#if nanofluid
double Prnf;
double Phi;
double dp;
double Lenf, St, tau;
#define H Phi
#else
double Pr = Nuf/Kappaf;
//double Pr = 0.71;
#endif

int MAXLEVEL;

/** init temperature field(tracer.h + diffusion.h) */
#if nanofluid
scalar Np[], T[], CLT[], rt[], rbrow[], rthermo[], rc[];
scalar * tracers = {T,Np};
mgstats mgT;
mgstats mgC;
#else
scalar T[], CLT[];
scalar * tracers = {T};
mgstats mgT;
#endif

/** We need a new field to define the acceleration, viscosity, temperature diffusion, phi diffusion. */
face vector av[], muc[], D[];
#if nanofluid
face vector Dc[];
#endif

#if nanofluid
#define mu(T) (((2.414e-5*pow(10,247.8/((T)-140.))) + (Rhop*VBB(T)*sq(dp))/(72.*CC(T)*deltaa))/MUnf)
#define VBB(T) (sqrt((18*kb*(T))/(M_PI*Rhop*dp))*1./dp)
#define deltaa (pow(M_PI/(6*Phi),1./3)*dp)
#define CC(T) (((-0.000001133*dp*1e9-0.000002771)*Phi+(0.00000009*dp*1e9-0.000000393))/(2.414e-5*pow(10,247.8/((T)-140.))))
#endif

/** ## Set viscosity  */
event properties (i++) {
  foreach_face(){
#if nanofluid
    double ff = face_value(T,0);
    muc.x[] = fm.x[]*Prnf/sqrt(Ranf)*(mu(ff+T0));
#else
    muc.x[] = fm.x[]*Pr/sqrt(Ranf);
#endif
  }
  boundary ((scalar*){muc});
}

/** 
## Set acceleration 
*/
event acceleration (i++) {
  coord grav = {0, 1};
  foreach_face(){
#if nanofluid
    av.x[] += grav.x*fm.x[]*Prnf*(T[] + T[0,-1])/2.;
#else
    av.x[] += grav.x*fm.x[]*Pr*(T[] + T[0,-1])/2.;
#endif
  }
  boundary ((scalar*){av});
}

/** 
## Set temperature diffusion + nanoparticle diffusion (Np)
*/
#if nanofluid
#define D4ORDNPx ((Np[-1,0] - Np[1,0])/(2*Delta)) //difference finis ordre 2 centre derive premiere
#define D4ORDNPy ((Np[0,-1] - Np[0,1])/(2*Delta))
#define D4ORDTx ((T[-1,0]-T[1,0])/(2*Delta))
#define D4ORDTy ((T[0,-1]-T[0,1])/(2*Delta))

#define D4ORD2Tx ((T[1,0] - 2*T[] + T[-1,0])/(sq(Delta)))
#define D4ORD2Ty ((T[0,1] - 2*T[] + T[0,-1])/(sq(Delta)))


//#define Prf (Nuf/Kappaf)
#define RhofT(T) (2446. -20.674*(T) + 0.11576*sq(T) -3.12895e-4*cube(T) + 4.0505e-7*pow(T,4) -2.0546e-10*pow(T,5))
#define kfT(T) (-0.76761 + 7.535211e-3*(T) -0.98249e-5*sq(T))
#define CpfT(T) (exp((8.29041-0.012557*(T))/(1.-1.52373e-3*(T))))
#define PrfT(T) (((2.414e-5*pow(10,247.8/((T)-140.)))/RhofT(T))/(kfT(T)/(RhofT(T)*CpfT(T))))

#define Re(T) ((2*RhofT(T)*kb*(T))/(M_PI*sq((2.414e-5*pow(10,247.8/((T)-140.))))*dp))
#define k(T) (((1+4.4*pow(Re(T),0.4)*pow(PrfT(T),0.66)*pow((T)/Tfr,10)*pow(kp/kfT(T),0.03)*pow(Phi,0.66))*kfT(T))/knf)
#endif

event tracer_diffusion (i++) {
  foreach_face(){
#if nanofluid
    double ff = face_value(T,0);
    D.x[] = fm.x[]*1./sqrt(Ranf)*k(ff+T0);
    if(y!=Y0 && y!=fabs(Y0)) //conservation des nanoparticules
      Dc.x[] = fm.x[]*1./(sqrt(Ranf)*Lenf);
    else
      Dc.x[] = 0.;
#else
    D.x[] = fm.x[]*1./sqrt(Ranf);
#endif
  }
  boundary ((scalar*){D});
#if nanofluid
  boundary ((scalar*){Dc});
/**
  The energy equation for the nanofluid :
$$
        \partial_{t} \widetilde{\theta} + \widetilde{u}.\widetilde{\nabla} \widetilde{\theta} = \underbrace{\frac{\widetilde{\nabla}^{2}\widetilde{\theta}}{\sqrt{Ra_{nf}}}}_{\text{Thermal diffusion D[]}} + \underbrace{\frac{\tau \widetilde{\nabla} \Phi . \widetilde{\nabla} \widetilde{\theta}}{Le_{nf}\sqrt{Ra_{nf}}}}_{\text{Brownian motion rbrow[]}} + \underbrace{\frac{\tau \widetilde{\nabla} \widetilde{\theta} . \widetilde{ \nabla} \widetilde{\theta}}{S_{t}\sqrt{Ra_{nf}}}}_{\text{thermophoresis rthermo[]}}
$$
*/
  foreach(){
    //if(x!=X0 && x!=fabs(X0)){
      rt [] = tau/(Lenf*sqrt(Ranf))*(D4ORDNPx * D4ORDTx + D4ORDNPy * D4ORDTy)
	+ (tau*St)/(sqrt(Ranf))*(D4ORDTx * D4ORDTx + D4ORDTy* D4ORDTy);
      rbrow[] = tau/(Lenf*sqrt(Ranf))*(D4ORDNPx * D4ORDTx + D4ORDNPy * D4ORDTy);
      rthermo[] = (tau*St)/(sqrt(Ranf))*(D4ORDTx * D4ORDTx + D4ORDTy* D4ORDTy);
      //}
    //else
      //rt[] =0.;
  }
  /**
  Continuity equation for nanoparticles :
$$
        \partial_{t}\Phi + \widetilde{u}. \widetilde{\nabla} \Phi = \frac{1}{\sqrt{Ra_{nf}}}[\underbrace{\frac{\widetilde{\nabla}^{2}\Phi}{Le_{nf}}}_\text{nanoparticle diffusion Dc[]} + \underbrace{S_{t} \widetilde{\nabla}^{2} \widetilde{\theta}}_\text{rb[]}]
$$
*/
  foreach(){
    if(y!=Y0 && y!=fabs(Y0)) //conservation des nanoparticules
      rc[] = (St/(sqrt(Ranf)))*(D4ORD2Tx + D4ORD2Ty);
    else
      rc[] = 0.;
  }
  //double dtmaxT = sq(L0/pow(2,MAXLEVEL))/min(tau/(sqrt(Ranf)*Lenf),(tau*St/(sqrt(Ranf))));
  //dt=dtnext(dtmaxT);
  mgT = diffusion (T, dt, D, r = rt);
  boundary({T});
  //double dtmaxC = sq(L0/pow(2,MAXLEVEL))/min(1./(sqrt(Ranf)*Lenf),(St/(sqrt(Ranf))));
  //dt=dtnext(dtmaxC);
  mgC = diffusion(Np, dt, Dc, r = rc);
  //mgC = diffusion(Np, dt, Dc);
  boundary({Np});
  char name[256];
  sprintf(name, "stats_%2.2e.dat", Ranf);
  static FILE * fc = fopen (name, "w");
  if (i == 0)
    fprintf(fc, "#t dt mgT.i mgT.nrelax mgC.i mgC.nrelax\n");
  fprintf (fc, "%f %g %d %d %d %d\n", t, dt, mgT.i, mgT.nrelax, mgC.i, mgC.nrelax);
  fflush(fc);
#else
  mgT = diffusion (T, dt, D);
  boundary({T});
  char name[256];
  sprintf(name, "stats_%2.2e.dat", Ranf);
  static FILE * fc = fopen (name, "w");
  if (i == 0)
    fprintf(fc, "#t dt mgT.i mgT.nrelax\n");
  fprintf (fc, "%f %g %d %d\n", t, dt, mgT.i, mgT.nrelax);
  fflush(fc);
#endif
}