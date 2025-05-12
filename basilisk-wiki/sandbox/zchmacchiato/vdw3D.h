#define BGHOSTS 2
#define growthTimestep 0.1
#include "run.h"
#include "timestep-vdw.h"
#include "difference.h"

scalar normGradRho[];
scalar rho[],rhop[];
vector q[],qp[],u[];
scalar mu[];
scalar e0[];
scalar laprho[];
double P0( double x);

(const) double lambda;
(const) double miu;

event defaults (i = 0)
{
  CFL = 0.1;
}

double dtmax;

event init (i = 0)
{
  boundary (all);
  event ("properties");
  dtmax = DT;
  event ("stability");
}

event set_dtmax (i++,last) dtmax = DT;

event stability (i++,last) {
  double mydtmax = HUGE;
  foreach(reduction(min:mydtmax)){
    if ( miu != 0.) {

      double dt = 0.1*sq(Delta)/miu;

      if (dt < mydtmax) mydtmax = dt;
    }
  }
  dtmax = min(mydtmax,dtmax);
  foreach()
    foreach_dimension()
      u.x[] = q.x[]/rho[];
  dt = dtnext (timestep (u, q, dtmax));
}

event predictor_corrector(i++,last)
{



  trash({rhop});
    foreach()
     {
      rhop[] = rho[]; 
      foreach_dimension(){
        rhop[] -=dt/Delta*(q.x[]-q.x[-1]);
	      u.x[] = q.x[]/rho[];}
     }

      boundary({rhop});
      foreach(){
      //normGradRho[] = abs(rhop[1,0]-rhop[-1,0])/2./Delta;
      e0[]=fabs(-8./3.*rhop[]*0.95*log(1./rhop[]-1./3.)-3.*sq(rhop[])+1.1759*(rhop[]-1.46)+2.5287);
      laprho[]=(rho[1,0,0]+rho[-1,0,0]+rho[0,1,0]+rho[0,-1,0]+rho[0,0,1]+rho[0,0,-1]-6.*rho[])/sq(Delta);
      }
      boundary({normGradRho});
      boundary(( scalar *) {u});

    boundary(all);
    trash({qp});
    foreach(){
      foreach_dimension()
        qp.x[] = q.x[] - dt*         ((rho[]*sq(u.x[])-rho[-1,0]*sq(u.x[-1,0]) + 
                                       rho[]*u.x[]*u.y[]-rho[0,-1]*u.x[0,-1]*u.y[0,-1]+
                                       rho[]*u.x[]*u.z[]-rho[0,0,-1]*u.x[0,0,-1]*u.z[0,0,-1]+  
                                       P0(rho[])-P0(rho[-1,0]))/Delta  
                       -             ((0.5*lambda*(sq((rho[0,1]-rho[0,-1])/2./Delta))-
                                       0.5*lambda*(sq((rho[-1,1]-rho[-1,-1])/2./Delta)))/Delta+
                                      (0.5*lambda*(sq((rho[0,0,1]-rho[0,0,-1])/2./Delta))-
                                       0.5*lambda*(sq((rho[-1,0,1]-rho[-1,0,-1])/2./Delta)))/Delta-
                                      (0.5*lambda*(sq((rho[1,0]-rho[])/Delta))-
                                       0.5*lambda*(sq((rho[0,0]-rho[-1,0])/Delta)))/Delta+
                                      (lambda*rho[]*laprho[]-lambda*rho[-1,0]*laprho[-1,0])/Delta-
                                      (lambda*(rho[0,1]-rho[])/Delta*(rho[1,0]-rho[-1,0])/2./Delta-
                                       lambda*(rho[0,0]-rho[0,-1])/Delta*(rho[1,-1]-rho[-1,-1])/2./Delta)/Delta-
                                      (lambda*(rho[0,0,1]-rho[])/Delta*(rho[1,0,0]-rho[-1,0,0])/2./Delta-
                                       lambda*(rho[]-rho[0,0,-1])/Delta*(rho[1,0,-1]-rho[-1,0,-1])/2./Delta)/Delta) 
                       -             ((mu[]*(4.0/3.0*(u.x[1,0]-u.x[])/Delta-2.0/3.0*(u.y[0,1]-u.y[0,-1]+u.z[0,0,1]-u.z[0,0,-1])/2./Delta)-
                                       mu[]*(4.0/3.0*(u.x[0,0]-u.x[-1,0])/Delta-2.0/3.0*(u.y[-1,1]-u.y[-1,-1]+u.z[-1,0,1]-u.z[-1,0,-1])/2./Delta))/Delta+
                                      (mu[]*((u.x[0,1]-u.x[])/Delta+(u.y[1,0]-u.y[-1,0])/2./Delta)-
                                       mu[]*((u.x[0,0]-u.x[0,-1])/Delta+(u.y[1,-1]-u.y[-1,-1])/2./Delta))/Delta+
                                      (mu[]*((u.x[0,0,1]-u.x[])/Delta+(u.z[1,0,0]-u.z[-1,0,0])/2./Delta)-
                                       mu[]*((u.x[]-u.x[0,0,-1])/Delta+(u.z[1,0,-1]-u.z[-1,0,-1])/2./Delta))/Delta
                                        ));
    }                                    

  boundary(( scalar *) {qp});
  foreach()
  {
    rho[] = 0.5*(rhop[]+rho[]);
    foreach_dimension(){
	    rho[] -= 0.5*dt/Delta*(qp.x[1]-qp.x[]);
      u.x[] = qp.x[]/rhop[];
  }}
  boundary({rho});

      foreach(){
      //normGradRho[] = abs(rho[1,0]-rho[-1,0])/2./Delta;
      e0[]=fabs(-8./3.*rho[]*0.95*log(1./rho[]-1./3.)-3.*sq(rho[])+1.1759*(rho[]-1.46)+2.5287);
      laprho[]=(rhop[1,0,0]+rhop[-1,0,0]+rhop[0,1,0]+rhop[0,-1,0]+rhop[0,0,1]+rhop[0,0,-1]-6.*rhop[])/sq(Delta);
      }
  boundary({normGradRho});
  boundary(( scalar *) {u});
  boundary(all);
  foreach(){
    foreach_dimension()
      q.x[] = 0.5*(qp.x[] + q.x[]) - 0.5*dt*    ((rhop[1,0]*sq(u.x[1,0])-rhop[]*sq(u.x[])+
                                                  rhop[0,1]*u.x[0,1]*u.y[0,1]-rhop[]*u.x[]*u.y[]+
                                                  rhop[0,0,1]*u.x[0,0,1]*u.z[0,0,1]-rhop[]*u.x[]*u.z[]+
                                                  P0(rhop[1,0])-P0(rhop[]))/Delta
                                    -            ((0.5*lambda*(sq((rhop[1,1]-rhop[1,-1])/2./Delta))-
                                                    0.5*lambda*(sq((rhop[0,1]-rhop[0,-1])/2./Delta)))/Delta+
                                                  (0.5*lambda*(sq((rhop[1,0,1]-rhop[1,0,-1])/2./Delta))-
                                                    0.5*lambda*(sq((rhop[0,0,1]-rhop[0,0,-1])/2./Delta)))/Delta-
                                                  (0.5*lambda*(sq((rhop[1,0]-rhop[0,0])/Delta))-
                                                    0.5*lambda*(sq((rhop[]-rhop[-1,0])/Delta)))/Delta +
                                                  (lambda*rhop[1,0]*laprho[1,0]-lambda*rhop[]*laprho[])/Delta-
                                                  (lambda*(rhop[0,1]-rhop[0,0])/Delta*(rhop[1,1]-rhop[-1,1])/2./Delta-
                                                   lambda*(rhop[]-rhop[0,-1])/Delta*(rhop[1,0]-rhop[-1,0])/2./Delta)/Delta-
                                                  (lambda*(rhop[0,0,1]-rhop[])/Delta*(rhop[1,0,1]-rhop[-1,0,1])/2./Delta-
                                                   lambda*(rhop[]-rhop[0,0,-1])/Delta*(rhop[1,0,0]-rhop[-1,0,0])/2./Delta)/Delta)
                                    -            ((mu[]*(4.0/3.0*(u.x[1,0]-u.x[0,0])/Delta-2.0/3.0*(u.y[1,1]-u.y[1,-1]+u.z[1,0,1]-u.z[1,0,-1])/2./Delta)-
                                                    mu[]*(4.0/3.0*(u.x[]-u.x[-1,0])/Delta-2.0/3.0*(u.y[0,1]-u.y[0,-1]+u.z[0,0,1]-u.z[0,0,-1])/2./Delta))/Delta+    
                                                  (mu[]*((u.x[0,1]-u.x[0,0])/Delta+(u.y[1,1]-u.y[-1,1])/2./Delta)-
                                                   mu[]*((u.x[]-u.x[0,-1])/Delta+(u.y[1,0]-u.y[-1,0])/2./Delta))/Delta+
                                                  (mu[]*((u.x[0,0,1]-u.x[])/Delta+(u.z[1,0,1]-u.z[-1,0,1])/2./Delta)-
                                                   mu[]*((u.x[]-u.x[0,0,-1])/Delta+(u.z[1,0,0]-u.z[-1,0,0])/2./Delta))/Delta));

}
  boundary(( scalar *) {q});
  restriction((scalar *){q});

}
	