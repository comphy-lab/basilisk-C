/**
# Incompressible Navier--Stokes solver (centered formulation)
*/

#include "run.h"
#include "timestep.h"
#include "bcg.h"
#include "viscosity.h"

scalar p[];
vector u[],g[];
face vector uf[];
scalar pf[];

(const) face vector mu = zerof, a = zerof, alpha = unityf;
(const) scalar rho = unity;
mgstats mgp, mgu, mgpf;
bool stokes = false;

event defaults (i=0){

  CFL = 0.95;
  p.nodump = pf.nodump = true;
  if (alpha.x.i == unityf.x.i) {
    alpha = fm;
    rho = cm;
  }
  else if (!is_constant(alpha.x)) {
    face vector alphav = alpha;
    foreach_face()
      alphav.x[] = fm.x[];
    boundary ((scalar *){alpha});
  }

#if TREE
  uf.x.refine = refine_face_solenoidal;
#endif
}

double dtmax;

event init (i = 0)
{
  boundary ((scalar *){u});
  trash ({uf});
  foreach_face()
    uf.x[] = fm.x[]*(u.x[] + u.x[-1])/2.;
  boundary ((scalar *){uf});

  event ("properties");
  dtmax = DT;
  event ("stability");
}

event set_dtmax (i++,last) dtmax = DT;

event stability (i++,last) {
  dt = dtnext (timestep (uf, dtmax));
}

event properties (i++,last) {
  boundary ({alpha, mu, rho});
}

void prediction(vector u_temporary, face vector uf_temporary, double dt_temporary)
{
  vector du;
  foreach_dimension() {
    scalar s = new scalar;
    du.x = s;
  }

  if (u.x.gradient)
    foreach()
      foreach_dimension()
        du.x[] = u.x.gradient (u_temporary.x[-1], u_temporary.x[], u_temporary.x[1])/Delta;
  else
    foreach()
      foreach_dimension()
        du.x[] = (u_temporary.x[1] - u_temporary.x[-1])/(2.*Delta);
  boundary ((scalar *){du});

  trash ({uf});
  foreach_face() {
    double un = dt_temporary*(u_temporary.x[] + u_temporary.x[-1])/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    uf_temporary.x[] = u_temporary.x[i] + (g.x[] + g.x[-1])*dt_temporary/4. + s*(1. - s*un)*du.x[i]*Delta/2.;
    #if dimension > 1
      double fyy = u_temporary.y[i] < 0. ? u_temporary.x[i,1] - u_temporary.x[i] : u_temporary.x[i] - u_temporary.x[i,-1];
      uf_temporary.x[] -= dt_temporary*u_temporary.y[i]*fyy/(2.*Delta);
    #endif
    #if dimension > 2
      double fzz = u_temporary.z[i] < 0. ? u_temporary.x[i,0,1] - u_temporary.x[i] : u_temporary.x[i] - u_temporary.x[i,0,-1];
      uf_temporary.x[] -= dt_temporary*u_temporary.z[i]*fzz/(2.*Delta);
    #endif
    uf_temporary.x[] *= fm.x[];
  }
  boundary ((scalar *){uf_temporary});

  delete ((scalar *){du});
}

static void correction (vector u_c, vector g_c, double dt_c)
{
  foreach()
    foreach_dimension()
      u_c.x[] += dt_c*g_c.x[];
  boundary ((scalar *){u_c});  
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

   if(!stokes) {
      prediction(u_temp,uf_temp,dt_temp);
      mgpf = project (uf_temp, pf, alpha, dt_temp/2., mgpf.nrelax);
      advection ((scalar *){u_temp}, uf_temp, dt_temp, (scalar *){g_temp});
     }
  
   if(constant(mu.x)!=0){
      correction(u_temp,g_temp,dt_temp);
      mgu = viscosity(u_temp,mu,rho,dt_temp,mgu.nrelax);
      correction(u_temp,g_temp,-dt_temp);
     }
   trash({uf_temp});
   foreach_face()
      uf_temp.x[] = fm.x[]*(u_temp.x[] + u_temp.x[-1])/2.;
   boundary_flux({uf_temp});

   mgp = project (uf_temp,p_temp,alpha,dt_temp,mgp.nrelax);
   face vector gf[];
   foreach_face()
      gf.x[] = -1.*(alpha.x[]/fm.x[])*(p_temp[] - p_temp[-1])/(Delta);
   boundary_flux ({gf});
   trash({g_temp});
   foreach()
     foreach_dimension()
       g_temp.x[] = (gf.x[] + gf.x[1])/2.;
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