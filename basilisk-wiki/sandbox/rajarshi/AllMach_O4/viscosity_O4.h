/**
# Header file required for fourth order implicit viscosity solver
*/

#include "poisson_O4.h"

struct Viscosity {
  
   vector u;
   face vector mu;
   scalar rho;
   double dt;
   int nrelax;
   scalar * res;
};

#if dimension == 1
#define lambda ((coord){1.})
#elif dimension == 2
#define lambda ((coord){1.,1.})
#elif dimension == 3
#define lambda ((coord){1.,1.,1.})
#endif

static void relax_viscosity (scalar *a, scalar *b, int l, void * data){

   struct Viscosity * p = (struct Viscosity *) data;
   (const) face vector mu = p->mu;
   (const) scalar rho = p->rho;
   double dt = p->dt;
   vector u = vector(a[0]), r = vector(b[0]);

   #if Jacobi
      vector w[];
   #else
      vector w = u;
   #endif

   foreach_level_or_leaf (l) {
      foreach_dimension()
         w.x[] =   (dt/rho[]* (   2.*mu.x[1]*(u.x[-1] + 15.*u.x[1] - u.x[2])/12.
                                 - 2.*mu.x[0]*(u.x[-2] - 15.*u.x[-1] - u.x[1])/12.
         #if dimension > 1
                   + mu.y[0,1]*(    (u.x[0,-1]+15.*u.x[0,1]-u.x[0,2])/(12.)
 
                                      + (   1.*(-u.y[-2,2] + 7.*u.y[-2,1] + 7.*u.y[-2,0] - u.y[-2,-1]) 
                                          - 8.*(-u.y[-1,2] + 7.*u.y[-1,1] + 7.*u.y[-1,0] - u.y[-1,-1])  
                                          + 8.*(-u.y[1,2]  + 7.*u.y[1,1]  + 7.*u.y[1,0]  - u.y[1,-1])
                                          - 1.*(-u.y[2,2]  + 7.*u.y[2,1]  + 7.*u.y[2,0]  - u.y[2,-1])
                                        )/(12.*12.)
                               )

                   - mu.y[0,0]*(    (u.x[0,-2]-15.*u.x[0,-1]-u.x[0,1])/(12.)
 
                                      + (   1.*(-u.y[-2,1] + 7.*u.y[-2,0] + 7.*u.y[-2,-1] - u.y[-2,-2]) 
                                          - 8.*(-u.y[-1,1] + 7.*u.y[-1,0] + 7.*u.y[-1,-1] - u.y[-1,-2])  
                                          + 8.*(-u.y[1,1]  + 7.*u.y[1,0]  + 7.*u.y[1,-1]  - u.y[1,-2])
                                          - 1.*(-u.y[2,1]  + 7.*u.y[2,0]  + 7.*u.y[2,-1]  - u.y[2,-2])
                                        )/(12.*12.)
                               )
         #endif
                               )  + r.x[]*sq(Delta) ) / (lambda.x*sq(Delta)

                + (dt/rho[])*( (5./4.)*(2.*mu.x[1] + 2.*mu.x[0]
         #if dimension > 1
                + mu.y[0,1] + mu.y[0,0]
         #endif
                 )));
   }

   #if Jacobi
      foreach_level_or_leaf (l)
         foreach_dimension()
             u.x[] = (u.x[] + 2.*w.x[])/3.;
   #endif
 
   #if TRASH
      vector u1[];
      foreach_level_or_leaf (l)
         foreach_dimension()
            u1.x[] = u.x[];
      trash ({u});
      foreach_level_or_leaf (l)
         foreach_dimension()
            u.x[] = u1.x[];
   #endif

}

static double residual_viscosity (scalar * a, scalar * b, scalar * resl, void * data){

   struct Viscosity * p = (struct Viscosity *) data;
   (const) face vector mu = p->mu;
   (const) scalar rho = p->rho;
   double dt = p->dt;
   vector u = vector(a[0]), r = vector(b[0]), res = vector(resl[0]);
   double maxres = 0.;

   #if TREE
      foreach_dimension() {
          face vector taux[];
          foreach_face (x)
              taux.x[] = 2.*mu.x[]*(u.x[-2]-15.*u.x[-1]+15.*u.x[]-u.x[1])/(12.*Delta);
          #if dimension > 1
          foreach_face (y)
              taux.y[] =    mu.y[]*( (u.x[0,-2]-15.*u.x[0,-1]+15.*u.x[0,0]-u.x[0,1])/(12.*Delta) 
                       + (   1.*(-u.y[-2,1] + 7.*u.y[-2,0] + 7.*u.y[-2,-1] - u.y[-2,-2]) 
                           - 8.*(-u.y[-1,1] + 7.*u.y[-1,0] + 7.*u.y[-1,-1] - u.y[-1,-2])  
                           + 8.*(-u.y[1,1]  + 7.*u.y[1,0]  + 7.*u.y[1,-1]  - u.y[1,-2])
                           - 1.*(-u.y[2,1]  + 7.*u.y[2,0]  + 7.*u.y[2,-1]  - u.y[2,-2])
                         )/(12.*12.*Delta) );
          #endif            
          boundary_flux({taux});
          foreach (reduction(max:maxres)) {
             double d = 0.;
             foreach_dimension()
                d += taux.x[1] - taux.x[];
             res.x[] = r.x[] - lambda.x*u.x[] + (dt/rho[])*(d/Delta);
             if (fabs(res.x[]) > max)
                maxres = fabs(res.x[]);
          }
       }   

   #else
      foreach(reduction(max:maxres))
         foreach_dimension(){
             res.x[] = r.x[] - lambda.x*u.x[] + (dt/rho[])*(
                                                           2.*mu.x[1,0]*(u.x[-1] - 15.*u.x[0]  + 15.*u.x[1] -u.x[2])/12. 
                                                         - 2.*mu.x[0,0]*(u.x[-2] - 15.*u.x[-1] + 15.*u.x[]  -u.x[1])/12.
             #if dimension > 1
                       + mu.y[0,1]*(    (u.x[0,-1]-15.*u.x[0,0]+15.*u.x[0,1]-u.x[0,2])/(12.)
 
                                      + (   1.*(-u.y[-2,2] + 7.*u.y[-2,1] + 7.*u.y[-2,0] - u.y[-2,-1]) 
                                          - 8.*(-u.y[-1,2] + 7.*u.y[-1,1] + 7.*u.y[-1,0] - u.y[-1,-1])  
                                          + 8.*(-u.y[1,2]  + 7.*u.y[1,1]  + 7.*u.y[1,0]  - u.y[1,-1])
                                          - 1.*(-u.y[2,2]  + 7.*u.y[2,1]  + 7.*u.y[2,0]  - u.y[2,-1])
                                        )/(12.*12.)
                                   )

                       - mu.y[0,0]*(    (u.x[0,-2]-15.*u.x[0,-1]+15.*u.x[0,0]-u.x[0,1])/(12.)
 
                                      + (   1.*(-u.y[-2,1] + 7.*u.y[-2,0] + 7.*u.y[-2,-1] - u.y[-2,-2]) 
                                          - 8.*(-u.y[-1,1] + 7.*u.y[-1,0] + 7.*u.y[-1,-1] - u.y[-1,-2])  
                                          + 8.*(-u.y[1,1]  + 7.*u.y[1,0]  + 7.*u.y[1,-1]  - u.y[1,-2])
                                          - 1.*(-u.y[2,1]  + 7.*u.y[2,0]  + 7.*u.y[2,-1]  - u.y[2,-2])
                                        )/(12.*12.)
                                   )
             #endif
                       )/sq(Delta);
             if(fabs(res.x[]) > maxres)
                maxres = fabs (res.x[]);
         }
   #endif

   boundary (resl);
   return (maxres);

}   

#undef lambda

mgstats viscosity (struct Viscosity p) {

  vector u = p.u, r[];
  foreach()
     foreach_dimension()
        r.x[] = u.x[];
  face vector mu = p.mu;
  scalar rho = p.rho;
  restriction ({mu,rho});
  return mg_solve ( (scalar *){u}, (scalar *){r}, residual_viscosity, relax_viscosity, &p, p.nrelax, p.res);

}

mgstats viscosity_explicit (struct Viscosity p){

  vector u = p.u, r[];
  mgstats mg = {0};
  mg.resb = residual_viscosity ((scalar *){u}, (scalar *){u}, (scalar *){r}, &p);
  foreach()
     foreach_dimension()
        u.x[] += r.x[];
  boundary ((scalar *){u});
  return mg;

} 