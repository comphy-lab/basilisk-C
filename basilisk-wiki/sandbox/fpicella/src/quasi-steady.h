#define QUASISTEADY 1
/**
# Unsteady to steady Navier-Stokes solver.
Overload end_timestep event so to iterate up to convergence to a quasi-steady state.
It is designed to work on top of navier-stokes/centered.h
 */
event end_timestep(i++,last)
{
/** 
Extra parameters so to force the system to converge at EACH iteration towards a STEADY solution.
I just re-call some steps from the navier-stokes/centered.h solver, a given number of times (limit iSTOKES_max)
so to get up a given convergence ratio (iSTOKES_tol) */
//  int iSTOKES_max = 100;
//  double iSTOKES_tol = 1e-5;
  scalar un[]; // field where I store the previous solution...
  int iSTOKES, iSTOKES_STOP=iSTOKES_max;
  double du_STOP;
  for (iSTOKES = 0; iSTOKES<iSTOKES_max; iSTOKES++){
    event("advection_term");
    event("viscous_term");
    event("acceleration");
    event("projection");
    event("HOOK");// an additional hook, so that I can play around with. FP, 20250316 11h11, with LoloGege
    double avg = normf(u.x).avg, du = change (u.x, un)/(avg + SEPS);
//    fprintf(stderr,"i= %d, t = %f, iSTOKES = %d, du = %f \n",i,t,iSTOKES,du);
    if (du<iSTOKES_tol){
      iSTOKES_STOP = iSTOKES;
      iSTOKES=iSTOKES_max;
    }
    du_STOP = du;
  }
  if(iSTOKES_STOP < iSTOKES_max)
  fprintf(stderr,"-> -> -> STOKES CONVERGED, i= %06d, t = %6.5e, iSTOKES_STOP = %04d, du_STOP = %6.5e \n",i,t,iSTOKES_STOP,du_STOP);
  if(iSTOKES_STOP ==iSTOKES_max)
  fprintf(stderr,"STOKES UNCONVERGED <-<-<-, i= %06d, t = %6.5e, iSTOKES_STOP = %04d, du_STOP = %6.5e \n",i,t,iSTOKES_STOP,du_STOP);
	event ("when_steady"); // hook to play with...
}
event when_steady(i++){
	tnext += DT_when_steady - dt;
}
