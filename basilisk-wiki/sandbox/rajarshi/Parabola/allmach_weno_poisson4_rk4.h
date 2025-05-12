/**
# Navier Stokes All Mach Formulation - RK4 + WENO-5 + O4-(Viscosity & Projection) 
*/

#include "run.h"
#include "timestep.h"
#define BGHOSTS 2
#include "O5NS_Flux.h"
#include "viscosity_O4.h"

vector q[];
scalar p[];
face vector uf[];
face vector alpha[],mu[];
scalar rho[];

scalar ps[];
(const) scalar rhoc2 = zeroc;
(const) face vector a = zerof;

vector g[];

event defaults (i = 0) {

  foreach_face(){
    alpha.x[] = 1.;
    mu.x[] = 0.;
  }
  boundary((scalar *){alpha,mu});

  foreach()
    rho[] = 1.;
  boundary({rho});
  
  Prolongation_Weight_Initialization();

}

double dtmax;

extern void properties (scalar hp, scalar rhoc2p, scalar psp, vector alphap, vector mup);
extern void acceleration (face vector alphap, face vector ap);

event init (i = 0) {

  boundary ({q,p,rho});
  properties(rho,rhoc2,ps,alpha,mu);

  tensor gq[];
  foreach()
     foreach_dimension()
        gq.x.x[] = (q.x[1]-q.x[-1])/(2.*Delta);
  boundary((scalar *){gq});
 
  trash({uf});
  foreach_face(){
      if(q.x[] >= 0)
          uf.x[] = fm.x[]*alpha.x[]*weno5_left_x  (point,q.x,gq.x.x[-2]);
      else
          uf.x[] = fm.x[]*alpha.x[]*weno5_right_x (point,q.x,gq.x.x[ 1]);
  }
  boundary ((scalar *){uf});

}

event set_dtmax (i++,last) dtmax = DT;

event stability (i++,last) {
  dt = dtnext (timestep (uf, dtmax));
}

mgstats mgp,mgu;

static void LOperator (scalar * fields, double t, scalar * slopefields){

  scalar * temp_fields = list_clone(fields);

  scalar i,j;
  for (i,j in fields,temp_fields){
       if(!is_constant(i) && i.face)
           foreach_face()
              j[] = i[];
       else
           foreach()
              j[] = i[];
  }
  boundary(temp_fields);    

  vector q_temp,g_temp;
  scalar p_temp,ps_temp;
  face vector uf_temp;
  scalar rho_temp, rhoc2_temp;
  face vector a_temp,alpha_temp,mu_temp;

  q_temp     = vector(temp_fields[0]);
  g_temp     = vector(temp_fields[1]);
  p_temp     = temp_fields[2];
  ps_temp    = temp_fields[3];
  uf_temp    = vector(temp_fields[4]);
  rho_temp   = temp_fields[5];
  rhoc2_temp = temp_fields[6];
  alpha_temp = vector(temp_fields[7]);
  mu_temp    = vector(temp_fields[8]);
  a_temp     = vector(temp_fields[9]);

  advection ( (scalar *){q_temp}, uf_temp, dt);
  advection ( (scalar *){rho_temp}, uf_temp, dt);

  properties(rho_temp,rhoc2_temp,ps_temp,alpha_temp,mu_temp);
  acceleration(alpha_temp,a_temp);

  foreach()
    printf("\n %g %g",x,rho_temp[]);
 
/* 
  if (constant(mu_temp.x)!=0.){
       foreach()
          foreach_dimension()
             q_temp.x[] = (q_temp.x[] + dt*g_temp.x[])/rho_temp[];
       boundary ((scalar *){q_temp});
       mgu = viscosity (q_temp,mu_temp,rho_temp,dt,mgu.nrelax);
       foreach()
          foreach_dimension()
             q_temp.x[] = q_temp.x[]*rho_temp[] - dt*g_temp.x[];
       boundary((scalar *){q_temp});
     }
*/

  tensor gq_temp[];
  foreach()
     foreach_dimension()
        gq_temp.x.x[] = (q_temp.x[1]-q_temp.x[-1])/(2.*Delta);
  boundary((scalar *){gq_temp});

 // foreach()
   // printf("\n %g %g",x,q_temp.x[]);

  foreach_face(){
      if(q_temp.x[] >= 0)
          uf_temp.x[] = fm.x[]*alpha_temp.x[]*weno5_left_x  (point,q_temp.x,gq_temp.x.x[-2]);
      else
          uf_temp.x[] = fm.x[]*alpha_temp.x[]*weno5_right_x (point,q_temp.x,gq_temp.x.x[ 1]);
  }
  boundary((scalar *){uf_temp});

  scalar lambda = rhoc2_temp, rhs = ps_temp;
  foreach() {
    if (constant(lambda) == 0.)
      rhs[] = 0.;
    else {
      lambda[] = - cm[]/(sq(dt)*rhoc2_temp[]);
      rhs[] = lambda[]*ps_temp[];
    }
 
    double div = 0.;
    foreach_dimension()
      div += uf_temp.x[1] - uf_temp.x[];
    rhs[] += div/(dt*Delta);
  }
  boundary ({lambda, rhs});
   
  mgp = poisson (p_temp, rhs, alpha_temp, lambda);

  face vector gf[];
  foreach_face() {
    double dp = alpha_temp.x[]*(p_temp[-2] - 15.*p_temp[-1] + 15.*p_temp[] - p_temp[1])/(12.*Delta);
    uf_temp.x[] -= dt*dp;
    gf.x[] = a_temp.x[] - dp/fm.x[];
  }
  boundary((scalar *){gf});
 
  foreach()
    foreach_dimension() {
      g_temp.x[]  = rho_temp[]*( -1.*gf.x[-1] + 13.*gf.x[0] + 13.*gf.x[1] - gf.x[2] )/24.; 
      q_temp.x[] += dt*g_temp.x[];
    }
  boundary ((scalar *){q_temp,g_temp,uf_temp});
  
  foreach(){
     scalar i,j,k;
     for (i,j,k in fields,temp_fields, slopefields)
          k[] = (j[]-i[])/dt;
  }    

  boundary(slopefields);
  delete(temp_fields);
  free(temp_fields);

}

#include "runge-kutta.h"
event TimeMarching (i++,last)
{
  runge_kutta({q,g,p,ps,uf,rho,rhoc2,alpha,mu,a},t,dt,LOperator,4);
}

/**
After mesh adaptation fluid properties need to be updated. */

#if TREE
event adapt (i++,last) {
  event ("properties");
}
#endif
