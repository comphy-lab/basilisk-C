#include "poisson.h"
#include "penalization.h"

#if AXI
#define lambda ((coord){1., 1. + dt/rho[]*(mu.x[] + mu.x[1] + \
                                mu.y[] + mu.y[0,1])/2./sq(y)})
#else // not AXI
    #if dimension == 1
        #define lambda ((coord){1.})
    #elif dimension == 2
        #define lambda ((coord){1.,1.})
    #elif dimension == 3
        #define lambda ((coord){1.,1.,1.})
    #endif
#endif

#if dimension == 1
    #define scalar_a_by_b(a, b) (a.x[]*b.x[])
#elif dimension == 2
    #define scalar_a_by_b(a, b) (a.x[]*b.x[] + a.y[]*b.y[])
#else // dimension == 3
    #define scalar_a_by_b(a, b) (a.x[]*b.x[] + a.y[]*b.y[] + a.z[]*b.z[])
#endif

#if dimension == 1
    #define m_scalar_a_by_b(a, b) (a.x[]*(b.x[1] - b.x[-1]))
#elif dimension == 2
    #define m_scalar_a_by_b(a, b) (a.x[]*(b.x[1] - b.x[-1]) + a.y[]*(b.y[1] - b.y[-1]))
#else // dimension == 3
    #define m_scalar_a_by_b(a, b) (a.x[]*(b.x[1] - b.x[-1]) + a.y[]*(b.y[1] - b.y[-1]) + a.z[]*(b.z[1] - b.z[-1]))
#endif

struct Viscosity {
    vector u;
    face vector mu;
    scalar rho;
    double dt;
    int nrelax;
    scalar * res;
};

static void relax_viscosity (scalar*a, scalar*b, int l, void*data)
{
    struct Viscosity*p = (struct Viscosity *) data;
    (const) face vector mu = p->mu;
    (const) scalar rho = p->rho;
    double dt = p->dt;
    vector u = vector(a[0]), r = vector(b[0]);
    #if JACOBI
        vector w[];
    #else
        vector w = u;
    #endif
    foreach_level_or_leaf (l) {
        foreach_dimension(){
            w.x[] = (frhs*(dt/rho[])*(2.*mu.x[1]*u.x[1] + 2.*mu.x[]*u.x[-1]
                #if dimension > 1
                    + mu.y[0,1]*(u.x[0,1] +
                        (u.y[1,0] + u.y[1,1])/4. -
                        (u.y[-1,0] + u.y[-1,1])/4.)
                    - mu.y[]*(- u.x[0,-1] +
                         (u.y[1,-1] + u.y[1,0])/4. -
                         (u.y[-1,-1] + u.y[-1,0])/4.)
                #endif
                #if dimension > 2
                    + mu.z[0,0,1]*(u.x[0,0,1] +
                          (u.z[1,0,0] + u.z[1,0,1])/4. -
                          (u.z[-1,0,0] + u.z[-1,0,1])/4.)
                    - mu.z[]*(- u.x[0,0,-1] +
                         (u.z[1,0,-1] + u.z[1,0,0])/4. -
                         (u.z[-1,0,-1] + u.z[-1,0,0])/4.)
                #endif
                    )
                    PLUS_NUMERATOR_BRINKMAN
                    + r.x[]*sq(Delta))/(sq(Delta)*lambda.x
                    PLUS_DENOMINATOR_BRINKMAN
                    + frhs*dt/rho[]*(2.*mu.x[1] + 2.*mu.x[]
                #if dimension > 1
                    + mu.y[0,1] + mu.y[]
                #endif
                #if dimension > 2
                    + mu.z[0,0,1] + mu.z[]
                #endif
                    ) );
        }
    }

    #if JACOBI
        foreach_level_or_leaf (l)
            foreach_dimension()
                u.x[] = (u.x[] + 2.*w.x[])/3.;
    #endif

    #if TRASH
        vector u1[];
        foreach_level_or_leaf (l)  foreach_dimension() u1.x[] = u.x[];
        trash ({u});
        foreach_level_or_leaf (l) foreach_dimension()  u.x[] = u1.x[];
    #endif
}

static double residual_viscosity (scalar * a, scalar * b, scalar * resl, void * data) {
    struct Viscosity * p = (struct Viscosity *) data;
    (const) face vector mu = p->mu;
    (const) scalar rho = p->rho;
    double dt = p->dt;
    vector u = vector(a[0]), r = vector(b[0]), res = vector(resl[0]);
    double maxres = 0, d = 0;
    #ifdef CALC_GRAD_U_TAU
        calc_target_U(u, target_U, n_sol);
    #endif
    foreach_dimension() {
        face vector taux[];
        foreach_face(x)
        taux.x[] = 2.*mu.x[]*(u.x[] - u.x[-1])/Delta;
        #if dimension > 1
        foreach_face(y)
              taux.y[] = mu.y[]*(u.x[] - u.x[0,-1] +
                     (u.y[1,-1] + u.y[1,0])/4. -
                     (u.y[-1,-1] + u.y[-1,0])/4.)/Delta;
        #endif
        #if dimension > 2
        foreach_face(z)
              taux.z[] = mu.z[]*(u.x[] - u.x[0,0,-1] +
                     (u.z[1,0,-1] + u.z[1,0,0])/4. -
                     (u.z[-1,0,-1] + u.z[-1,0,0])/4.)/Delta;
        #endif
        boundary_flux ({taux});
        foreach (reduction(max:maxres)) {
            d = 0;
            foreach_dimension() d += taux.x[1] - taux.x[];
            res.x[] = r.x[]
                      PLUS_BRINKMAN_RHS
                      - lambda.x*u.x[] + frhs*dt*d/(Delta*rho[]);
            if (fabs (res.x[]) > maxres) maxres = fabs (res.x[]);
        }
    }
    boundary (resl);
    #ifdef DEBUG_BRINKMAN_PENALIZATION
        foreach (){
            foreach_dimension() {
                dbp.x[] = PLUS_BRINKMAN_RHS/dt;
                total_rhs.x[] = (res.x[] - r.x[] + lambda.x*u.x[])/dt;
            }
        }
    #endif
//    fprintf(ferr, "maxres=%g\n", maxres);
    return maxres;
}

#undef lambda

trace
mgstats viscosity (struct Viscosity p){
    #if AXI
        fprintf(ferr, "Not correct solution for AXI in viscosity-weugene.h");
        exit(9);
    #endif
    #if BRINKMAN_PENALIZATION_ERROR_BC
        fprintf(ferr, "Check BRINKMAN_PENALIZATION_NO_SLIP and BRINKMAN_PENALIZATION_SLIP in viscosity-weugene.h");
        exit(1);
    #endif
    vector u = p.u, r[];
    foreach() foreach_dimension(){ r.x[] = u.x[];}
    face vector mu = p.mu;
    scalar rho = p.rho;
    restriction ({mu,rho});
    return mg_solve ((scalar *){u}, (scalar *){r}, residual_viscosity, relax_viscosity, &p, p.nrelax, p.res);
}
