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


vector u1[],u2[],u3[],u4[],g1[],g2[],g3[],g4[];
scalar p1[],p2[],p3[],p4[];
face vector uf1[],uf2[],uf3[],uf4[];
vector k1[],k2[],k3[],k4[];

event RungeKutta4 (i++,last) {

   // Compute slope k1 - Time Integration [ t , t + Delta_t/2 ] 
   
   foreach()
      foreach_dimension(){
         u1.x[] = u.x[];
         g1.x[] = g.x[];
      }

   foreach_face() 
      uf1.x[] = uf.x[]; 

   boundary((scalar *){u1,g1,uf1});
   Update (u1,p1,uf1,g1,dt/2.,k1);   
   
   
   // Compute slope k2 - Time Integration [ t + Delta_t/2, t + Delta_t ]
   
   foreach()
      foreach_dimension(){
         u2.x[] = u1.x[];
         g2.x[] = g1.x[];
      }

   foreach_face() 
      uf2.x[] = uf1.x[]; 

   boundary((scalar *){u2,g2,uf2}); 
   Update (u2,p2,uf2,g2,dt/2.,k2);


   // Compute slope k3 - Time Integration [ t + Delta_t/2, t + Delta_t ]
   
   foreach()
      foreach_dimension()
         u3.x[] = u.x[] + (dt/2.)*k2.x[];
   boundary((scalar *){u3});   

   trash({uf3});
   foreach_face()
      uf3.x[] = fm.x[]*(-1.*u3.x[-2] + 7.*u3.x[-1] + 7.*u3.x[] - u3.x[1])/12.;
   boundary_flux({uf3});

   mgp = project (uf3,p3,alpha,dt/2.,mgp.nrelax);
   face vector gf[];
   foreach_face()
      gf.x[] = -1.*(alpha.x[]/fm.x[])*(p3[-2] - 15.*p3[-1] + 15.*p3[] - p3[1])/(12.*Delta);
   boundary_flux ({gf});
   trash({g3});
   foreach()
     foreach_dimension()
       g3.x[] = ( -1.*gf.x[-1] + 13.*gf.x[0] + 13.*gf.x[1] - gf.x[2] )/24.;
   boundary((scalar *){g3});
   correction (u3,g3,dt/2.);

   Update (u3,p3,uf3,g3,dt/2.,k3);

   // Compute slope k4 - Time Integration [ t + Delta_t, t + 2.*Delta_t ]

   foreach()
      foreach_dimension()
         u4.x[] = u.x[] + dt*k3.x[];
   boundary((scalar *){u4});   

   trash({uf4});
   foreach_face()
      uf4.x[] = fm.x[]*(-1.*u4.x[-2] + 7.*u4.x[-1] + 7.*u4.x[] - u4.x[1])/12.;
   boundary_flux({uf4});

   mgp = project (uf4,p4,alpha,dt,mgp.nrelax);
   trash({gf});
   foreach_face()
      gf.x[] = -1.*(alpha.x[]/fm.x[])*(p4[-2] - 15.*p4[-1] + 15.*p4[] - p4[1])/(12.*Delta);
   boundary_flux ({gf});
   trash({g4});
   foreach()
     foreach_dimension()
       g4.x[] = ( -1.*gf.x[-1] + 13.*gf.x[0] + 13.*gf.x[1] - gf.x[2] )/24.;
   boundary((scalar *){g4});
   correction (u4,g4,dt);

   Update (u4,p4,uf4,g4,dt,k4);


   // Time Marching -->   Y_(n+1) = Y_(n) + (Delta_t / 6) * (k1 + 2*k2 + 2*k3 + k4)

   foreach()
      foreach_dimension()
         u.x[] += (dt/6.)*( k1.x[] + 2.*k2.x[] + 2.*k3.x[] + k4.x[]);
   boundary((scalar *){u});

   trash({uf});
   foreach_face()
      uf.x[] = fm.x[]*(-1.*u.x[-2] + 7.*u.x[-1] + 7.*u.x[] - u.x[1])/12.;
   boundary_flux({uf});

   mgp = project (uf,p,alpha,dt,mgp.nrelax);
   trash({gf});
   foreach_face()
      gf.x[] = -1.*(alpha.x[]/fm.x[])*(p[-2] - 15.*p[-1] + 15.*p[] - p[1])/(12.*Delta);
   boundary_flux ({gf});
   trash({g});
   foreach()
     foreach_dimension()
       g.x[] = ( -1.*gf.x[-1] + 13.*gf.x[0] + 13.*gf.x[1] - gf.x[2] )/24.;
   boundary((scalar *){g});
   correction (u,g,dt);      
    
}

#if TREE
event adapt (i++,last) {
  event ("properties");
}
#endif