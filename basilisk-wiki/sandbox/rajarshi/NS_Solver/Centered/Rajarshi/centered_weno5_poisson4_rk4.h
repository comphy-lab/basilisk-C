/**
# Navier Stokes Centered Formulation - RK4 + WENO-5 + O4-(Viscosity & Projection) 
*/

#include "run.h"
#include "timestep.h"
#include "O5_Flux.h"
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

static void Update (vector u_temp, scalar p_temp, face vector uf_temp, vector g_temp,  double dt_temp, vector su, vector sg, face vector suf, scalar sp){

   vector u_in[],g_in[];
   scalar p_in[];
   face vector uf_in[];

   foreach(){
     p_in[] = p_temp[];
     foreach_dimension(){
        u_in.x[] = u_temp.x[];
        g_in.x[] = g_temp.x[];
     }
   }
   foreach_face()
     uf_in.x[] = uf_temp.x[];

   boundary({p_in});
   boundary((scalar *){u_in,g_in,uf_in});

   if(!(stokes))
       advection ((scalar *){u_temp}, uf_temp, dt_temp);
  
   if(constant(mu.x)!=0){
      correction(u_temp,g_temp,dt_temp);
      mgu = viscosity(u_temp,mu,rho,dt_temp,mgu.nrelax);
      correction(u_temp,g_temp,-dt_temp);
     }
 
   trash({uf_temp});
   tensor gu_temp[];
   for (scalar s in {gu_temp})
      s.gradient = zero;
   gradients({u_temp},{gu_temp});

   foreach_face(){
      if(u_temp.x[] >= 0)
          uf_temp.x[] = fm.x[]*WENO5_Reconstruction_Left_x (point,u_temp.x,gu_temp.x.x[-2]);
      else
          uf_temp.x[] = fm.x[]*WENO5_Reconstruction_Right_x (point,u_temp.x,gu_temp.x.x[1]);
   }
   boundary_flux({uf_temp});

   mgp = project (uf_temp,p_temp,alpha,dt_temp,mgp.nrelax);
   face vector gf[];
   foreach_face()
      gf.x[] = -1.*(alpha.x[]/fm.x[])*(p_temp[-2] - 15.*p_temp[-1] + 15.*p_temp[] - p_temp[1])/(12.*Delta);
   boundary_flux ({gf});
   boundary((scalar *){gf});
   trash({g_temp});
   foreach()
     foreach_dimension()
       g_temp.x[] = ( -1.*gf.x[-1] + 13.*gf.x[0] + 13.*gf.x[1] - gf.x[2] )/24.;
   boundary((scalar *){g_temp});
   correction (u_temp,g_temp,dt_temp);    

   foreach(){
     sp[] = (p_temp[] - p_in[])/dt_temp;  
     foreach_dimension(){
        su.x[] = (u_temp.x[] - u_in.x[])/dt_temp;
        sg.x[] = (g_temp.x[] - g_in.x[])/dt_temp;
     }
   }
   foreach_face()
      suf.x[] = (uf_temp.x[] - uf_in.x[])/dt_temp;

   boundary({sp});
   boundary((scalar *){su,sg,suf});
}

#include "runge-kutta.h"
event TimeMarching (i++,last) {
  RungeKutta4();
}

#if TREE
event adapt (i++,last) {
  event ("properties");
}
#endif