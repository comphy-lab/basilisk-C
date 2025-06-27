#include "simple_discretization.h"

#include "LS_recons.h"
#include "phase_change_velocity.h"
#if CURVE_LS
 #include "LS_curvature.h"
#else
#include "curvature.h"
#endif

struct LSspeed{
  scalar dist;
  double L_H;
  scalar cs;
  face vector fs;
  scalar TS;
  scalar TL;
  double T_eq;
  vector vpc;
  vector vpcf;
  double lambda1;
  double lambda2;
  double epsK;
  double epsV;
  double eps4;
  double deltat;
  int itrecons;
  double tolrecons;
  double NB_width;
};

void LS_speed(struct LSspeed p){

  scalar dist      = p.dist;
  double L_H       = p.L_H;
  scalar cs        = p.cs;
  face vector fs   = p.fs;
  scalar TS        = p.TS;
  scalar TL        = p.TL;
  double T_eq      = p.T_eq;
  vector vpc       = p.vpc;
  vector vpcf      = p.vpcf;
  double lambda[2] = {p.lambda1, p.lambda2};
  double epsK      = p.epsK;
  double epsV      = p.epsV;
  double eps4      = p.eps4;
  double deltat    = p.deltat;
  int    itrecons  = p.itrecons;
  double tolrecons = p.tolrecons;
  double NB_width = p.NB_width;

  scalar curve[];
#if LS_perf
double start; 
double end; 
#endif
  if(NB_width==0){fprintf(stderr, "Erreur NB_width LS_speed\n" );exit(1);}

/** Curvature must be updated*/
#if CURVE_LS
  curvature_LS(dist, curve);
#else
  curvature (cs,curve);
  boundary({curve});
#endif

#if LS_perf
start = omp_get_wtime(); 
#endif 
  phase_change_velocity_LS_embed (cs, fs ,TL, TS, T_eq, vpc, L_H, 
    lambda,epsK, epsV, eps4, curve);
#if LS_perf
end = omp_get_wtime(); 
fprintf(stderr,"phase change velo %f seconds\n", end - start);
#endif 
/**
We copy this value in vpcr.
*/
  vector vpcr[];
  foreach(){
    foreach_dimension(){
      if(interfacial(point,cs))vpcr.x[] = vpc.x[];
      else vpcr.x[] = 0.;
    }
  }
  boundary((scalar * ){vpcr});
  restriction((scalar * ){vpcr});
#if dimension == 1
  scalar * speed_recons  = {vpcr.x};
#elif dimension == 2
  scalar * speed_recons  = {vpcr.x,vpcr.y};
#else
  scalar * speed_recons  = {vpcr.x,vpcr.y,vpcr.z};
#endif

/**
We reconstruct a cell-centered vpc field in the vicinity of the interface where
vpcr is the solution to the bilinear interpolation of vpc on face centroids.
*/
  double err = 0;



#if LS_perf
start = omp_get_wtime(); 
#endif

  recons_speed(dist, deltat, speed_recons,
   tolrecons, &err, 
   itrecons, 
   cs, fs,
   NB_width);
#if LS_perf
end = omp_get_wtime(); 
fprintf(stderr,"recons_speed %f seconds\n", end - start);
#endif 
 
  foreach()
    foreach_dimension(){
      if(fabs(dist[])< 0.99*NB_width)
        vpcf.x[] = vpcr.x[];
      else
        vpcf.x[] = 0.;
    }

  boundary((scalar *){vpcf});
  restriction((scalar *){vpcf});
}

event stability(i++){
  double lambda1 = lambda[0], lambda2 = lambda[1], dtmax; 
  LS_speed(
  dist,latent_heat,cs,fs,TS,TL,T_eq,
  vpc,vpcf,lambda1,lambda2,
  epsK,epsV,eps4,deltat=DT_LS,
  itrecons = itrecons,tolrecons = tolrecons,
  NB_width
  );

#if dtLS // only diffusion
  double DT3 = timestep_LS(vpcf,DT,dist,NB_width);
  tnext = t+DT3;
  dt = DT3;
#else // navier stokes
  dtmax = timestep_LS(vpcf,DT,dist,NB_width);
#endif
}