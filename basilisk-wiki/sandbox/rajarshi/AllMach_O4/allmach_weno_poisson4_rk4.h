/**
# Navier Stokes All Mach Formulation - RK4 + WENO-5 + O4-(Viscosity & Projection) 
*/

#include "run.h"
#include "timestep.h"
#define BGHOSTS 2
#include "O5_Flux.h"
#include "viscosity_O4.h"

vector q[];
scalar p[];
scalar pf[];
face vector uf[];
(const) face vector alpha = unityf, mu = zerof;
(const) scalar rho = unity;

scalar ps[];
(const) scalar rhoc2 = zeroc;
(const) face vector a = zerof;

vector g[];

event defaults (i = 0) {

  CFL=0.95;
  p.nodump = true;
  if(alpha.x.i == unityf.x.i){
     alpha = fm;
     rho = cm;
    }
  else if(!(is_constant(alpha.x))){
     face vector alphav = alpha;
     foreach_face()
        alphav.x[] = fm.x[];
     boundary((scalar *){alpha});
    }
  Prolongation_Weight_Initialization();

}

double dtmax;

event init (i = 0) {

  boundary ({q,p,rho});

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
  
  event ("properties");
  dtmax = DT;
  event ("stability");

}

event set_dtmax (i++,last) dtmax = DT;

event stability (i++,last) {
  dt = dtnext (timestep (uf, dtmax));
}

event properties (i++,last)
{
  boundary ({rho, rhoc2, ps, alpha, mu}); 
}

mgstats mgp,mgu;

static void LOperator (scalar * fields, double t, scalar * slopefields){

  scalar * temp_fields = list_clone(fields);  
  foreach(){
    scalar i,j;
    for(i,j in fields,temp_fields)
       j[] = i[];
  }
  boundary (temp_fields);

  vector q_temp,g_temp;
  scalar p_temp;
  face vector uf_temp;

  q_temp.x = temp_fields[0];
  q_temp.y = temp_fields[1];
  g_temp.x = temp_fields[2];
  g_temp.y = temp_fields[3];
  p_temp = temp_fields[4];
  uf_temp.x = temp_fields[5];
  uf_temp.y = temp_fields[6];

  advection ( (scalar *){q_temp}, uf_temp, dt);
  
  if (constant(mu.x)!=0.){
       foreach()
          foreach_dimension()
             q_temp.x[] = (q_temp.x[] + dt*g_temp.x[])/rho[];
       boundary ((scalar *){q_temp});
       mgu = viscosity (q_temp,mu,rho,dt,mgu.nrelax);
       foreach()
          foreach_dimension()
             q_temp.x[] = q_temp.x[]*rho[] - dt*g_temp.x[];
       boundary((scalar *){q_temp});
     }

  tensor gq_temp[];
  foreach()
     foreach_dimension()
        gq_temp.x.x[] = (q_temp.x[1]-q_temp.x[-1])/(2.*Delta);
  boundary((scalar *){gq_temp});
  foreach_face(){
      if(q_temp.x[] >= 0)
          uf_temp.x[] = fm.x[]*alpha.x[]*weno5_left_x  (point,q_temp.x,gq_temp.x.x[-2]);
      else
          uf_temp.x[] = fm.x[]*alpha.x[]*weno5_right_x (point,q_temp.x,gq_temp.x.x[ 1]);
  }
  boundary((scalar *){uf_temp});

  scalar lambda = rhoc2, rhs = ps;
  foreach() {
    if (constant(lambda) == 0.)
      rhs[] = 0.;
    else {
      lambda[] = - cm[]/(sq(dt)*rhoc2[]);
      rhs[] = lambda[]*ps[];
    }
 
    double div = 0.;
    foreach_dimension()
      div += uf_temp.x[1] - uf_temp.x[];
    rhs[] += div/(dt*Delta);
  }
  boundary ({lambda, rhs});
   
  mgp = poisson (p_temp, rhs, alpha, lambda);

  face vector gf[];
  foreach_face() {
    double dp = alpha.x[]*(p_temp[-2] - 15.*p_temp[-1] + 15.*p_temp[] - p_temp[1])/(12.*Delta);
    uf_temp.x[] -= dt*dp;
    gf.x[] = -dp/fm.x[];
  }
  boundary((scalar *){gf});
 
  foreach()
    foreach_dimension() {
      g_temp.x[]  = rho[]*( -1.*gf.x[-1] + 13.*gf.x[0] + 13.*gf.x[1] - gf.x[2] )/24.; 
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
event TimeMarching (i++,last) {
  //runge_kutta ({q,g,p,uf},t,dt,LOperator,4);
  runge_kutta ({q.x,q.y,g.x,g.y,p,uf.x,uf.y},t,dt,LOperator,4);
}

/**
After mesh adaptation fluid properties need to be updated. */

#if TREE
event adapt (i++,last) {
  event ("properties");
}
#endif