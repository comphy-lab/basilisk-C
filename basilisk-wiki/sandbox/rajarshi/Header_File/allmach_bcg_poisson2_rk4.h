/**
## All mach header file - RK1
*/

#include "run.h"
#include "timestep.h"
#include "bcg.h"
#include "viscosity.h"

vector q[];
scalar p[];
face vector uf[];
face vector alpha, mu;
scalar rho;

scalar ps[];
scalar rhoc2;
face vector a;

vector g[];

event defaults (i=0);

double dtmax;
extern void properties (scalar hp, scalar rhoc2p, scalar psp, vector alphap, vector mup);
extern void acceleration (face vector alphap, face vector ap);

event init (i = 0) {
  boundary ({q,p,rho});
  properties(rho,rhoc2,ps,alpha,mu);
  foreach_face()
    uf.x[] = alpha.x[]*(q.x[] + q.x[-1])/2.;
  boundary ((scalar *){uf});
}

event set_dtmax (i++,last) dtmax = DT;

event stability (i++,last) {
  dt = dtnext (timestep (uf, dtmax));
}

mgstats mgp, mgu;

static void LOperator (scalar * fields, double t, scalar * slopefields){

  scalar * temp_fields = list_clone(fields);

  scalar * fields_c = NULL;
  vector * fields_f = NULL;
  scalar * temp_fields_c = NULL;
  vector * temp_fields_f = NULL;
  scalar * slopefields_c = NULL;
  vector * slopefields_f = NULL;

  scalar i,j,k;
  for (i,j,k in fields,temp_fields,slopefields){
      if(i.face){
         fields_f = vectors_add(fields_f, i.v);
         temp_fields_f = vectors_add(temp_fields_f, j.v);
         slopefields_f = vectors_add(slopefields_f, k.v);
       } 
      else{
         fields_c = list_add(fields_c, i);
         temp_fields_c = list_add(temp_fields_c, j);
         slopefields_c = list_add(slopefields_c, k);
       }
    }


  vector l,m,n;

  for (l,m in fields_f,temp_fields_f)
       foreach_face()
          m.x[] = l.x[];

  for (i,j in fields_c,temp_fields_c){
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

  advection ( (scalar *){q_temp}, uf_temp, dt, (scalar *){g_temp});
  advection ( (scalar *){rho_temp}, uf_temp, dt);

  properties(rho_temp,rhoc2_temp,ps_temp,alpha_temp,mu_temp);
  acceleration(alpha_temp,a_temp);

/*  
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
*/

  trash({uf_temp});
  foreach_face()
    uf_temp.x[] = alpha_temp.x[]*(q_temp.x[] + q_temp.x[-1])/2. + dt*fm.x[]*a_temp.x[];
  boundary ((scalar *){uf_temp});

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
   
  mgp = poisson (p_temp, rhs, alpha_temp, lambda, tolerance = TOLERANCE/sq(dt));

  face vector gf[];
  foreach_face() {
    double dp = alpha_temp.x[]*(p_temp[] - p_temp[-1])/Delta;
    uf_temp.x[] -= dt*dp;
    gf.x[] = a_temp.x[] - dp/fm.x[];
  }
  boundary_flux ({gf});
 
  foreach()
    foreach_dimension() {
      g_temp.x[] = rho_temp[]*(gf.x[] + gf.x[1])/2.;
      q_temp.x[] += dt*g_temp.x[];
    }
  boundary ((scalar *){q_temp,g_temp,uf_temp});

  for (l,m,n in fields_f,temp_fields_f,slopefields_f)
       foreach_face()
          n.x[] = (m.x[]-l.x[])/dt;

  for (i,j,k in fields_c,temp_fields_c,slopefields_c)
       foreach()
          k[] = (j[]-i[])/dt;
v v v v v v v
  }    

  boundary(slopefields);
=============

  boundary(slopefields); 

  free(fields_c);
  free(fields_f);
  free(temp_fields_c);
  free(temp_fields_f);
  free(slopefields_c);
  free(slopefields_f);
*************

  boundary(slopefields); 

  free(fields_c);
  free(fields_f);
  free(temp_fields_c);
  free(temp_fields_f);
  free(slopefields_c);
  free(slopefields_f);
^ ^ ^ ^ ^ ^ ^
  delete(temp_fields);
  free(temp_fields);

}

#include "runge-kutta-new.h"
event TimeMarching (i++,last)
{
  runge_kutta({q,g,p,ps,uf,rho,rhoc2,alpha,mu,a},t,dt,LOperator,1);
}