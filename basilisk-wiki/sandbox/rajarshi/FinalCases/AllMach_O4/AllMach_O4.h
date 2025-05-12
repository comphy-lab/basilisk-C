/**
# An all mach flow solver
*/

#include "run.h"
#include "timestep.h"
#include "viscosity_O4.h"
#include "O5_Flux.h"

vector q[];
scalar p[];
face vector uf[];
(const) face vector alpha = unityf, mu = zerof;
(const) scalar rho = unity;

scalar ps[];
(const) scalar rhoc2 = zeroc;
(const) face vector a = zerof;

vector g[];



event defaults (i=0) {

  if(alpha.x.i = unityf.x.i)
      alpha = fm;

}

event init (i=0) {
  
  boundary({q,p,rho});
  event("properties");

  tensor grad_q[];
/*
  for(scalar s in {grad_q})
     s.gradient = zero;
*/
  gradients({q},{grad_q});
  foreach_face(){
      if(q.x[] >= 0)
          uf.x[] = alpha.x[]*WENO5_Reconstruction_Left_x  (point,q.x,grad_q.x.x[-2,0]);
      else
          uf.x[] = alpha.x[]*WENO5_Reconstruction_Right_x (point,q.x,grad_q.x.x[ 1,0]);
  }
  boundary((scalar *){uf});

}

double dtmax;

event set_dtmax(i++,last) dtmax = DT;

event stability (i++,last) {
  dt = dtnext (timestep (uf,dtmax));
}

event vof (i++,last);

event tracer_advection(i++,last){
  advection ((scalar *){q},uf,dt,(scalar *){g});  
}

event properties (i++,last){

  boundary ({rho,rhoc2,ps,alpha,mu});
  if (!is_constant(a.x)){
     face vector af = a;
     foreach_face()
        af.x[] = 0;
  }

}

event acceleration (i++,last){
  
  boundary ((scalar *){a});

}


mgstats mgp,mgu;

event pressure (i++,last){

  if (constant(mu.x) != 0.) {
     foreach()
        foreach_dimension()
           q.x[] = (q.x[] + dt*g.x[])/rho[];
     boundary((scalar *){q});
     mgu = viscosity (q,mu,rho,dt,mgu.nrelax);
     foreach()
        foreach_dimension()
           q.x[] = q.x[]*rho[] - dt*g.x[];
     boundary((scalar *){q});
     }

  tensor grad_q[];
/*
  for(scalar s in {grad_q})
     s.gradient = zero;
*/
  gradients({q},{grad_q});

  foreach_face(){
      if(q.x[] >= 0)
          uf.x[] = alpha.x[]*WENO5_Reconstruction_Left_x  (point,q.x,grad_q.x.x[-2,0]) + dt*fm.x[]*a.x[];
      else
          uf.x[] = alpha.x[]*WENO5_Reconstruction_Right_x (point,q.x,grad_q.x.x[ 1,0]) + dt*fm.x[]*a.x[];
  }
  boundary((scalar *){uf});

  scalar lambda = rhoc2, rhs = ps;
  foreach() {
     if(constant(lambda) == 0.)
        rhs[] = 0.;
     else{
        lambda[] = -cm[]/(sq(dt)*rhoc2[]);
        rhs[] = lambda[]*ps[];
        }
     double div = 0.;
     foreach_dimension()
       div += uf.x[1]-uf.x[];
     rhs[] += div/(dt*Delta);
   }
   boundary ({lambda,rhs});

   mgp = poisson (p,rhs,alpha,lambda,tolerance = TOLERANCE/sq(dt));

   face vector gf[];
   foreach_face() {
      double dp = alpha.x[]*(p[-2] - 15.*p[-1] + 15.*p[] - p[1])/(12.*Delta);
      uf.x[] -= dt*dp;
      gf.x[] = a.x[] - dp/fm.x[];
   }
   boundary_flux({gf});

   foreach()
     foreach_dimension() {
        g.x[] = rho[]*(-gf.x[-1] + 13.*gf.x[] + 13.*gf.x[1] - gf.x[2])/24.;
        q.x[] += dt*g.x[];
     }
     boundary ((scalar *){q,g,uf});
}

#if TREE
event adapt (i++,last){
  event ("properties");
}
#endif