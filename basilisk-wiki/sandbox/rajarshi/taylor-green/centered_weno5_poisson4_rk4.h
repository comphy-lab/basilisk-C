/**
# Navier Stokes Centered Formulation - RK4 + WENO-5 + O4-(Viscosity & Projection) 
*/

#include "run.h"
#include "timestep.h"
#include "O5NS_Flux.h"
#include "viscosity_O4.h"
#define BGHOSTS 2

scalar p[];
vector u[],g[];
face vector uf[];

(const) face vector mu = zerof, a = zerof, alpha = unityf;
(const) scalar rho = unity;
mgstats mgp, mgu;
bool stokes = false;

event defaults (i=0){

  CFL = 0.9;
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
  #if TREE
     uf.x.refine = refine_face_solenoidal;              // need to make order 4 ?
  #endif

}

double dtmax;

event init (i=0){
 
  boundary((scalar *){u});

  tensor gu[];
  foreach()
     foreach_dimension()
        gu.x.x[] = (u.x[1]-u.x[-1])/(2.*Delta);
  boundary((scalar *){gu});
 
  trash({uf});
  foreach_face(){
      if(u.x[] >= 0)
          uf.x[] = fm.x[]*weno5_left_x  (point,u.x,gu.x.x[-2]);
      else
          uf.x[] = fm.x[]*weno5_right_x (point,u.x,gu.x.x[ 1]);
  }
  boundary_flux({uf});

  event ("properties");
  dtmax = DT;
  event ("stability");

}

event set_dtmax (i++,last) dtmax = DT;

event stability (i++,last){
  dt = dtnext(timestep(uf,dtmax));
}

event properties (i++,last){
  boundary({alpha,mu,rho});
}

static void correction (vector u_correction, vector g_correction, double dt_correction){
  foreach()
    foreach_dimension()
      u_correction.x[] += dt_correction*g_correction.x[];
  boundary((scalar *){u_correction});
}

static void LOperator (scalar * fields, double t, scalar * slopefields){

   scalar * temp_fields = list_clone(fields);  
   foreach(){
      scalar i,j;
      for (i,j in fields,temp_fields)
         j[] = i[];
   }
   boundary (temp_fields);

   vector u_temp,g_temp;
   scalar p_temp;
   face vector uf_temp;

   u_temp.x = temp_fields[0];
   u_temp.y = temp_fields[1];
   g_temp.x = temp_fields[2];
   g_temp.y = temp_fields[3];
   p_temp = temp_fields[4];
   uf_temp.x = temp_fields[5];
   uf_temp.y = temp_fields[6]; 

   if(!(stokes))
       advection ((scalar *){u_temp}, uf_temp, dt);
  
   if(constant(mu.x)!=0){
      correction(u_temp,g_temp,dt);
      mgu = viscosity(u_temp,mu,rho,dt,mgu.nrelax);
      correction(u_temp,g_temp,-dt);
     }
 
   trash({uf_temp});

   tensor gu_temp[];
   foreach()
      foreach_dimension()
         gu_temp.x.x[] = (u_temp.x[1]-u_temp.x[-1])/(2.*Delta);
   boundary((scalar *){gu_temp});
   foreach_face(){
      if(u_temp.x[] >= 0)
          uf_temp.x[] = fm.x[]*weno5_left_x  (point,u_temp.x,gu_temp.x.x[-2]);
      else
          uf_temp.x[] = fm.x[]*weno5_right_x (point,u_temp.x,gu_temp.x.x[ 1]);
   }
   boundary_flux({uf_temp});

   mgp = project (uf_temp,p_temp,alpha,dt,mgp.nrelax);
   face vector gf[];
   foreach_face()
      gf.x[] = -1.*(alpha.x[]/fm.x[])*(p_temp[-2] - 15.*p_temp[-1] + 15.*p_temp[] - p_temp[1])/(12.*Delta);
   boundary_flux ({gf});
   boundary((scalar *){gf});
   trash({g_temp});
   foreach()
     foreach_dimension()
       g_temp.x[] = ( -1.*gf.x[-1] + 13.*gf.x[0] + 13.*gf.x[1] - gf.x[2] )/24.;    //(gf.x[] + gf.x[1])/2.;
   boundary((scalar *){g_temp});
   correction (u_temp,g_temp,dt);    

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
  runge_kutta ({u.x,u.y,g.x,g.y,p,uf.x,uf.y},t,dt,LOperator,4);
}

#if TREE
event adapt (i++,last) {
  event ("properties");
}
#endif