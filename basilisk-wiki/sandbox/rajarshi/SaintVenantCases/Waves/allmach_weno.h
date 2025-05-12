/**
# Navier Stokes All Mach Formulation - RK4 + WENO5-Advection + O4-(Viscosity & Projection) 
*/

#include "run.h"
#include "timestep.h"
#include "Weno_flux.h"
#include "poisson_O4.h"

vector q[];
scalar p[];
face vector uf[];
face vector alpha,mu;
scalar rho;

scalar ps[];
scalar rhoc2;
face vector a;

vector g[];

event defaults (i = 0);

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
          uf.x[] = alpha.x[]*weno5_left_x  (point,q.x,gq.x.x[-2]);
      else
          uf.x[] = alpha.x[]*weno5_right_x (point,q.x,gq.x.x[ 1]);
  }
  boundary ((scalar *){uf});

}

event set_dtmax (i++,last) dtmax = 0.1;

event stability (i++,last) {
  dt = dtnext (timestep (uf, dtmax));
}

mgstats mgp,mgu;

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

  scalar p_temp, ps_temp, rho_temp, rhoc2_temp ;
  vector q_temp,g_temp;
  face vector uf_temp, a_temp, alpha_temp, mu_temp;

  p_temp     = temp_fields[0];
  ps_temp    = temp_fields[1];
  rho_temp   = temp_fields[2];
  rhoc2_temp = temp_fields[3];
  q_temp     = vector(temp_fields[4]);
  g_temp     = vector(temp_fields[4 + 1*dimension]);
  uf_temp    = vector(temp_fields[4 + 2*dimension]);
  alpha_temp = vector(temp_fields[4 + 3*dimension]);
  mu_temp    = vector(temp_fields[4 + 4*dimension]);
  a_temp     = vector(temp_fields[4 + 5*dimension]);

  advection ( (scalar *){q_temp, rho_temp}, uf_temp, dt);
  properties(rho_temp,rhoc2_temp,ps_temp,alpha_temp,mu_temp);
  acceleration(alpha_temp,a_temp);

  tensor gq_temp[];
  foreach()
     foreach_dimension()
        gq_temp.x.x[] = (q_temp.x[1]-q_temp.x[-1])/(2.*Delta);
  boundary((scalar *){gq_temp});

  foreach_face(){
      if(q_temp.x[] >= 0)
          uf_temp.x[] = alpha_temp.x[]*weno5_left_x  (point,q_temp.x,gq_temp.x.x[-2]) + dt*fm.x[]*a_temp.x[];
      else
          uf_temp.x[] = alpha_temp.x[]*weno5_right_x (point,q_temp.x,gq_temp.x.x[ 1]) + dt*fm.x[]*a_temp.x[];
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
 
  vector gf_manipulator[];
  foreach()
    foreach_dimension()
        gf_manipulator.x[] = gf.x[] + gf.x[1];
  boundary ( ( scalar * ) {gf_manipulator} ); 

  foreach()
    foreach_dimension() {
      g_temp.x[]  = rho_temp[]*( -1.*(gf_manipulator.x[-1] - gf.x[]) + 13.*gf.x[] + 13.*gf.x[1] -1.*(gf_manipulator.x[1] - gf.x[1]) )/24.; 
      q_temp.x[] += dt*g_temp.x[];
    }
  boundary ((scalar *){q_temp,g_temp,uf_temp});

  
  for (l,m,n in fields_f,temp_fields_f,slopefields_f)
       foreach_face()
          n.x[] = (m.x[]-l.x[])/dt;

  for (i,j,k in fields_c,temp_fields_c,slopefields_c)
       foreach()
          k[] = (j[]-i[])/dt;

  boundary(slopefields); 

  free(fields_c);
  free(fields_f);
  free(temp_fields_c);
  free(temp_fields_f);
  free(slopefields_c);
  free(slopefields_f);
  delete(temp_fields);
  free(temp_fields);

}

#include "runge-kutta.h"
event TimeMarching (i++,last)
{
  runge_kutta({p,ps,rho,rhoc2,q,g,uf,alpha,mu,a},t,dt,LOperator,1);
}

/**
After mesh adaptation fluid properties need to be updated. */

#if TREE
event adapt (i++,last) {
  event ("properties");
}
#endif
