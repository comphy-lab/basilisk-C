/**
# Navier Stokes Centered Formulation - EulerTimeMarching + WENO-5 + O4-(Viscosity & Projection) 
*/

#include "run.h"
#include "timestep.h"
#include "O5_Flux.h"
#include "viscosity_O4.h"
#define BGHOStS 2

scalar p[];
vector u[],g[];
face vector uf[];

(const) face vector mu = zerof, a = zerof, alpha = unityf;
(const) scalar rho = unity;
mgstats mgp, mgu;
bool stokes = false;

event defaults (i=0){

  CFL = 0.2;
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

  #if TREE
     uf.x.refine = refine_face_solenoidal;              // need to make order 4 ?
  #endif

}

double dtmax;

event init (i=0){
 
  boundary((scalar *){u});
  trash({uf});
  foreach_face()
    uf.x[] = fm.x[]*(-1.*u.x[-2] + 7.*u.x[-1] + 7.*u.x[] - u.x[1])/12.;
  boundary((scalar *){uf});
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

static void Update (vector u_temp, scalar p_temp, face vector uf_temp, vector g_temp,  double dt_temp, vector slope){

   vector u_in[];
   foreach()
     foreach_dimension()
        u_in.x[] = u_temp.x[];
   boundary((scalar *){u_in});

   if(!(stokes))
       advection ((scalar *){u_temp}, uf_temp, dt_temp);
  
   if(constant(mu.x)!=0){
      correction(u_temp,g_temp,dt_temp);
      mgu = viscosity(u_temp,mu,rho,dt_temp,mgu.nrelax);
      correction(u_temp,g_temp,-dt_temp);
     }
 
   trash({uf_temp});
   foreach_face()
      uf_temp.x[] = fm.x[]*(-1.*u_temp.x[-2] + 7.*u_temp.x[-1] + 7.*u_temp.x[] - u_temp.x[1])/12.;
   boundary_flux({uf_temp});

   mgp = project (uf_temp,p_temp,alpha,dt_temp,mgp.nrelax);
   face vector gf[];
   foreach_face()
      gf.x[] = -1.*(alpha.x[]/fm.x[])*(p_temp[-2] - 15.*p_temp[-1] + 15.*p_temp[] - p_temp[1])/(12.*Delta);
   boundary_flux ({gf});
   trash({g_temp});
   foreach()
     foreach_dimension()
       g_temp.x[] = ( -1.*gf.x[-1] + 13.*gf.x[0] + 13.*gf.x[1] - gf.x[2] )/24.;
   boundary((scalar *){g_temp});
   correction (u_temp,g_temp,dt_temp);    

   foreach()
      foreach_dimension()
         slope.x[] = (u_temp.x[]-u_in.x[])/dt_temp;
   boundary((scalar *){slope});

}

event EulerFT (i++,last){
  vector k[];
  Update(u,p,uf,g,dt,k);
}

#if TREE
event adapt (i++,last) {
  event ("properties");
}
#endif