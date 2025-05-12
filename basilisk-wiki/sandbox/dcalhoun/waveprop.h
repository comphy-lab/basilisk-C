/**

<script type="text/x-mathjax-config">
MathJax.Hub.Config({
                tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']],
                    processEscapes: true},
            menuSettings: {
                    zoom: "Click",
               mpContext: true,
         mpMouse: true},
         errorSettings: {
             message: ["[Math Error]"] }
         });
</script>

<script type="text/javascript"
    src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML">
</script>


# Wave propagation algorithm

This implements the wave propagation algorithm (WPA) for solving systems of hyperbolic problems in either conservation form

$$
\mathbf q_t + \mathbf f_1(\mathbf q)_x + \mathbf f_2(\mathbf q)_y = 0
$$

or in non-conservative form

$$
\mathbf q_t + 
A_1(\mathbf q) \mathbf q_x + 
A_2(\mathbf q) \mathbf q_y = 0
$$
where $\mathbf q \in \mathbf R^m$ is a vector of state variables, $\mathbf f_1(\mathbf q), f_2(\mathbf q)\in \mathcal R^m$ are flux functions, and $A_1, A_2 \in \mathcal R^{m \times m}$ are  matrices.  For hyperbolic problems, we must have that the any linear combination of $A_1(\mathbf q)$ and $A_2(\mathbf q)$ is diagonizable with real eigenvalues. 

A complete description of the multi-dimensional wave propagation algorithm can be found in [LeVeque, 1997](le:1997) and [LeVeque, 2002](le:2002).

*/

#include "utils.h"   

/**

## Global arrays and parameters

### State variables

State variables, defined by the user in example files, are a list of Basilisk scalars. Components of velocity vectors are stored as individual scalars. 

*/

scalar *statevars = NULL;

/**

### Auxiliary variables

Many hyperbolic problems require auxiliary variables that store metric terms,  bathymetry, material properties, and so on.  These variables are stored in the ${\tt auxvars}$ variable and is typically set up by individual solvers. 

*/

scalar *auxvars = NULL;

/**
### Numerical parameters

The WPA can be run in either conservation or non-conservative form.  The solver will set the variable below, depending on the form of the problem and corresponding Riemann solver. 

Additionally, the F-wave approach (see [Bale et al., 2003](a-le-mi-ro:2003))  may be used and can offer higher order accuracy in some cases.  

The numerical order of the scheme can be specified as first order upwind ${\tt order}=1$ or  second order ${\tt order}=2$.
*/

bool conservation_law = true;   
bool use_fwaves = false;
int order = 2;               

/**
The WPA algorithm requires an estimate of the initial time step.
*/

bool dt_fixed = false;     
double dt_initial;
double dt;


/**
The WPA expects two layers of ghost cells for each stencil. 
*/

#define BGHOSTS 2

/** 
To create finer cells from coarse cells, we refine using linear interpolation.  To create a coarse grid cell from four finer grid cells, we use volume averaging.
*/

event defaults (i = 0) {
#if TREE
    for (scalar q in statevars)  {
        q.refine = q.prolongation = refine_linear;
        q.restriction = restriction_volume_average;
    }
#endif
}


/** 
It is sometimes useful to keep track of which plot frame is being created. This allows us to create and name an output file for later plotting.
*/

int Frame = 0; 


/**
These routine will set any specific user parameters.   We put this here so that the user can override any values set in the defaults event defined by the 
solvers.
*/

event wpa_parameters(i=0) {
    CFL = 0.9;

    extern bool limit_each_wave;    
    extern int* mthlim;
    extern int mwaves;
    extern int wave_limiter;

    if (!limit_each_wave) {
        if (mthlim == NULL) {
            fprintf(stderr,"waveprop.h : mthlim is a NULL pointer;  " \
                    "Set to valid memory location.\n");
            exit(0);
        }
        for(int mw=0; mw < mwaves; mw++)
            mthlim[mw] = wave_limiter;
    }
}

/**
Initialize global arrays at the beginning of a run.
*/
event init (i = 0) {
    boundary (statevars);
    if (auxvars != NULL){
        boundary (auxvars);        
    }
}

/**
Free up global arrays at the end of a run. 
*/

event cleanup(i=end, last) {
    free(statevars);
    if (auxvars != NULL) {
        free(auxvars);        
    }
}

/** 
At the beginning of a run, local arrays for accumulating fluxes and fluctuations are allocated.  
*/

trace
void wpa_initialize(vector **wpa_fm, vector **wpa_fp, vector **wpa_flux)  {
    *wpa_fm = NULL;
    *wpa_fp = NULL;
    *wpa_flux = NULL;

    int meqn = list_len(statevars);

    if (conservation_law)
        for(int m = 0; m < meqn; m++) {
            vector f = new face vector;
            *wpa_flux = vectors_append(*wpa_flux,f);            
        }
    else
        for(int m = 0; m < meqn; m++) {
            vector fmv = new face vector;
            *wpa_fm = vectors_append(*wpa_fm,fmv);

            vector fpv = new face vector;
            *wpa_fp = vectors_append(*wpa_fp,fpv);
        }
}

/** 
At the end of a run, this routine is called to clean up local variables used to store fluxes and fluctuations.
*/

void wpa_cleanup(vector** wpa_fm, vector** wpa_fp, vector** wpa_flux) {
    if (conservation_law) {
        delete((scalar*) *wpa_flux);
        free(*wpa_flux);   
    }
    else {
        delete((scalar*) *wpa_fp);
        free(*wpa_fp);

        delete((scalar*) *wpa_fm);
        free(*wpa_fm);        
    }                
}

attribute {
    int limiter;
}


/**
## Riemann solver

Riemann solvers are defined in solver files.  The Riemann solver for the WPA takes as input arguments the state and auxiliary variables from two horizontally or vertically adjacent cells and returns waves, speeds and fluctuations at the common cell interface. 
*/

void wpa_rpn2(int dir, int meqn, int mwaves, 
              double *ql, double *qr, 
              double *auxl, double *auxr, 
              double *waves, double *speeds, 
              double *amdq, double *apdq, double *flux);




/**
## Wave limiters

The simplest way to set limiters is to specify a single limiter for all waves.  To do this, the 
boolean variable ${\tt limit\_each\_wave}$ must be set to false, and a single numeric value (see below) will be used for each wave.   This is the default option for the wave propagation algorithm.  The default limiter is the generalized minmod limiter.  
*/

bool limit_each_wave = false;
int wave_limiter = 6;    /* Generalized minmod limiter */

/**
Individual waves can be limited separately by setting entries in the array ${\tt mthlim}$.  Again, this is usually managed by the solver file, and so the user does not explicitly set $\tt mthlim$, but will set a proxy variable defined by the solver. 

The variable ${\tt mwaves}$ is set by the solver file.   Often, we will have $\tt mwaves$ equal to the number of state variables $\tt meqn$.  But there are common cases where the number of waves is less  than the number of state variables.  For example, a tracer transport problem may involve several tracers all moving with a single wave speed.  

For this reason, there is no correct default value for $\tt mwaves$.   Similarly, the number of components in the array $\tt mthlim$ cannot also be known here, but is specified by individual solvers.
*/

int mwaves = 0;    /* Set by the solver */
int *mthlim = NULL;

/** 
The numeric values below can be used to specify wave limiters.
*/

double wpa_limiter(double a, double b,int mlim) {
    double wlimiter, r;
    r = b/a;

    switch (mlim) {
        double c,th;
        case 0:
            wlimiter = 1;
            break;

        /* Minmod */
        case 1:
            wlimiter = max(0., min(1., r));
            break;

        /* Superbee */
        case 2:
            wlimiter = max(0., max(min(1., 2.*r), min(2., r)));
            break;

        /* Van Leer */
        case 3:

            wlimiter = (r + fabs(r)) / (1. + fabs(r));
            break;

        /* Monotonized - Centered */
        case 4:
            c = (1. + r)/2.;
            wlimiter = max(0., min(c, min(2., 2.*r)));
            break;

        /* Sweby */
        case 5:
            wlimiter = max(0., max(min(r,1.5),min(1.5*r,1.)));
            break;

        /* Generalized minmod (default) */
        case 6:
            th = 1.3;
            wlimiter = max(0., min(th*r,min((1+r)/2.,th)));
            break;

        default:
            printf("Invalid limiter method.\n");
            exit(0);

    }
    return wlimiter; 
}

/**
## Advancing the solution

This is the main routine that advances the solution from one time step to the next

*/


double wpa_advance(double dt, double *cflmax) {
    double dtnew = DT;
    double cflmax_local = 0;

    boundary(statevars);
    if (auxvars != NULL){
        boundary (auxvars);        
    }

    vector *wpa_fp, *wpa_fm, *wpa_flux;
    wpa_initialize(&wpa_fm,&wpa_fp, &wpa_flux);

    int meqn = list_len(statevars);
    int maux = list_len(auxvars); 

/**

### Multidimensional sweeping

The basic WPA algorithm implemented here is a dimensionally split algorithm in which the solution is first updated using Riemann solvers in the x direction, and then updated a second time  using Riemann data from sweeps in the y-direction.

*/

    int dir = 0;
    foreach_dimension() {
        foreach_face(x, reduction(min:dtnew) reduction(max:cflmax_local)) {
            double amdq[meqn];
            double apdq[meqn];

            double speeds[mwaves];

            double flux[meqn];
            double waves[meqn*mwaves];    /* waves at interface I */

            double ql[meqn];
            double qr[meqn];

            double auxl[maux];
            double auxr[maux];

            int m = 0;        
            for (scalar q in statevars) {
                qr[m] = q[0];
                ql[m] = q[-1];
                m++;
            }

            m = 0;
            for (scalar a in auxvars) {
                auxr[m] = a[0];
                auxl[m] = a[-1];                    
                m++;
            }

            /* Check s.v.x.i == s.i   --> is x component */

/**

For conservation laws, we modify the original algorithms as implemented in Clawpack. This  requirement ensures conservation, but requires that the user supply the value of the flux at state variables.

*/
            for(int m = 0; m < meqn; m++)
                flux[m] = 0;

/**

### First order terms

Call the Riemann solver to get waves, speeds and fluctuations at the face between two adjacent cell values.

*/

            wpa_rpn2(dir,meqn, mwaves, ql, qr, auxl, auxr, waves, speeds, amdq, apdq, flux); 

/**

We use the results of the Riemann solve to compute a time step to be used in the next update (not this one). 

*/

            for(int mw = 0; mw < mwaves; mw++) {
                /* Compute next step */
                double s = fabs(speeds[mw]);

                if (s == 0)
                    continue;  /* speed of 0 does not constrain the time step */

                double dt_local = CFL*Delta/s;
                if (dt_local < dtnew)
                    dtnew = dt_local;  
                                  
                /* Test current CFL */
                double cfl = dt*s/Delta;
                if (cfl > cflmax_local)
                    cflmax_local = cfl;
            } 

/**

Fluctuations are accumulated in arrays ${\mathcal F}$ (for a conservation law) or  ${\mathcal F^-}$ and ${\tt \mathcal F^+}$ (for problems in non-conservative form).  These arrays are initialized with first order fluctuations, compute in the Riemann sovler. 
*/

            if (conservation_law) {
                m = 0;
                for (vector f in wpa_flux) {
                    f.x[] = flux[m] - apdq[m];   
                    m++;
                }
            }
            else {
                vector fm;
                vector fp;
                m = 0;
                for (fp, fm in wpa_fp, wpa_fm) {
                    fm.x[] = amdq[m];                
                    fp.x[] = -apdq[m]; 
                    m++;
                }                
            }

/**
### Second order corrections

The second order terms are computed using first order waves and speeds, and added to accumulation terms ${\tt \mathcal F^-}$ and ${\tt \mathcal F^+}$.

*/

            if (order == 2) {
                double wavesl[meqn*mwaves];   
                double wavesr[meqn*mwaves];               
                double s[mwaves]; 
                double ap[meqn];
                double am[meqn];
                double flx[meqn];

                double qm2[meqn];
                double qm1[meqn];
                double q0[meqn];
                double qp1[meqn];

                double auxm2[maux];
                double auxm1[maux];
                double aux0[maux];
                double auxp1[maux];

                int m = 0;
                for (scalar q in statevars) {
                    qm2[m] = q[-2];
                    qm1[m] = q[-1];
                    q0[m]  = q[0];
                    qp1[m] = q[1];

                    m++;
                }

                m = 0;
                for (scalar a in auxvars) {
                    auxm2[m] = a[-2];
                    auxm1[m] = a[-1];
                    aux0[m]  = a[0];
                    auxp1[m] = a[1];
                    m++;                    
                }

/**

We test to see if we are applying wave limiters.  If none of the waves requires any limiting, we can skip two additional Riemann solves required to 
limit the waves in place.  
*/

                bool use_limiting = false;
                for(int mw = 0; mw < mwaves; mw++)
                    use_limiting = use_limiting || mthlim[mw] > 0;

                double wlimiter[mwaves];

/**
Entries in the array ${\tt wlimiter}$ are scalar values in $[0,1]$ which will be used to multiply components of each wave.  If we are not using limiting, this value is set to 1.  
*/

                if (!use_limiting) {
                    for(int mw = 0; mw < mwaves; mw++)
                        wlimiter[mw] = 1;                                    
                }
                else {

/**
At least one wave requires limiting, and so we must make two additional calls to Riemann solvers so we can compare waves from neighboring cell faces.  
*/

                    wpa_rpn2(dir,meqn, mwaves, qm2, qm1, auxm2, auxm1, wavesl, 
                             s, am, ap,flx);                           
                    wpa_rpn2(dir,meqn, mwaves, q0,  qp1, aux0, auxp1, wavesr, 
                             s, am, ap,flx);

/**
Limiters are applied to each wave independently.  
*/

                    for(int mw = 0; mw < mwaves; mw++) {
                        double w2 = 0;
                        double dr = 0;  
                        double dl = 0;
                        for (int m = 0; m < meqn; m++) {
                            double w = waves[m + mw*meqn];
                            double wr = wavesr[m + mw*meqn];
                            double wl = wavesl[m + mw*meqn];
                            w2 += w*w;
                            dr += w*wr;
                            dl += w*wl;
                        }

                        if (w2 == 0 || mthlim[mw] == 0)
                            wlimiter[mw] = 1;
                        else {                
                            if (speeds[mw] >= 0)
                                wlimiter[mw] = wpa_limiter(w2,dl,mthlim[mw]);
                            else if (speeds[mw] < 0)
                                wlimiter[mw] = wpa_limiter(w2,dr,mthlim[mw]);
                        }
                    }
                }

/**

We use the limiting terms to compute second order corrections.  These corrections will be accumulated in the flux or fluctuations arrays $\mathcal F$, $\mathcal F^-$ or $\mathcal F^+$.  

*/

                double dtdx = dt/Delta;
                double cqxx[meqn];
                for(int m = 0; m < meqn; m++)
                    cqxx[m] = 0;

                for(int mw = 0; mw < mwaves; mw++) {
                    double abs_sign;
                    if (use_fwaves)
                        abs_sign = sign(speeds[mw]);
                    else
                        abs_sign = fabs(speeds[mw]);

                    double cq = abs_sign*(1 - fabs(speeds[mw])*dtdx)*
                                wlimiter[mw];
                    for(int m = 0; m < meqn; m++)
                        cqxx[m] += cq*waves[m+mw*meqn];
                }

                if (conservation_law) {
                    m = 0;
                    for (vector f in wpa_flux) {
                        f.x[]  += 0.5*cqxx[m];  
                        m++;
                    }  
                }
                else {
                    vector fm;
                    vector fp;
                    m = 0;
                    for (fp,fm in wpa_fp, wpa_fm) {
                        fm.x[] += 0.5*cqxx[m];
                        fp.x[] += 0.5*cqxx[m];   
                        m++;
                    }  
                }
            }   
        }

        if (conservation_law)
            boundary_flux(wpa_flux);
        else {
            boundary_flux(wpa_fp);
            boundary_flux(wpa_fm);            
        }

/**

### Solution update

We update the solution after each dimensional sweep.  Because the sweep in the y direction will use data from the sweep in the x direction, this scheme is stable up to CFL = 1. 

*/


        foreach() {
            double dtdx = dt/Delta;

            if (conservation_law) {            
                vector f;
                scalar q;
                for (f,q in wpa_flux, statevars) 
                    q[] += -dtdx*(f.x[1] - f.x[]);
            }
            else {
                vector fp;
                vector fm;
                scalar q;
                for (fp, fm, q in wpa_fp, wpa_fm, statevars) 
                    q[] += -dtdx*(fm.x[1] - fp.x[]);  
            }
        }
        dir++;
    } 

    wpa_cleanup(&wpa_fm, &wpa_fp, &wpa_flux);

    *cflmax = cflmax_local;
    return dtnew;
}

/**

## Main time stepping loop

*/


void run() {
    t = 0.;
    iter = 0;
    init_grid (N);  

    double cflmax = 0;
    double dtnew;
    //double dt_curr;

    perf.nc = perf.tnc = 0;
    perf.gt = timer_start();
    while (events (true)) {
        dtnew = wpa_advance(dt, &cflmax);
        if (cflmax > 1) {
            /* CFL using time step dt */
            printf("CFL exceeds 1; cflmax = %g\n",cflmax);
            //exit(0);
        }    
#if 0        
        dt_curr = dt;
        t += dt_curr; /* Use step just taken */        

        if (cflmax > 1) 
        {
            /* CFL using time step dt */
            printf("CFL exceeds 1; cflmax = %g\n",cflmax);
            exit(0);
        }    
        if (dt_fixed)
        {
            dt = dt_initial;  
        }    
        else
        {
            dt = dtnext(dtnew);            
            //t = tnext;  This is t += dtnew;  we want t += dt (above). 
        }
#endif
        dt = dtnext(dtnew);
        t = tnext;
        iter = inext;
        fprintf(stdout,"Step %6d : dt = %8.4e; cflmax = %12.6f; Time = %12.4f\n",
                iter, dt, cflmax,t);
        fflush(stdout);
    }
    timer_print (perf.gt, iter, perf.tnc);

    free_grid();
}


/**

## References

~~~bib
@article{le:1997,
    Author = {LeVeque, Randall J.},
    Date-Added = {2013-12-28 22:42:46 +0000},
    Date-Modified = {2016-04-28 20:51:11 +0000},
    Doi = {10.1006/jcph.1996.5603},
    Journal = {J. Comput. Phys.},
    Month = {Mar},
    Number = {2},
    Pages = {327-353},
    Publisher = {Elsevier -- Academic Press},
    Title = {Wave Propagation Algorithms for Multidimensional Hyperbolic Systems},
    Volume = {131},
    Year = {1997},
    DOI = {http://dx.doi.org/10.1006/jcph.1996.5603}}

@book{le:2002,
    Author = {Randall J. LeVeque},
    Date-Added = {2013-12-26 20:56:55 +0000},
    Date-Modified = {2013-12-27 02:43:06 +0000},
    Publisher = {Cambridge University Press},
    DOI = {10.1017/CBO9780511791253.004},
    Title = {Finite volume methods for hyperbolic problems},
    Year = 2002}

@article{ba-le-mi-ro:2003,
    Author = {Bale, Derek and LeVeque, Randall J. and Mitran, Sorin and Rossmanith, James},
    Date-Added = {2014-07-20 19:52:47 +0000},
    Date-Modified = {2014-07-20 19:53:57 +0000},
    Doi = {10.1137/S106482750139738X},
    Eprint = {http://dx.doi.org/10.1137/S106482750139738X},
    Journal = {SIAM J. Sci. Comput.},
    Number = {3},
    Pages = {955-978},
    Title = {{A Wave Propagation Method for Conservation Laws and Balance Laws with Spatially Varying Flux Functions}},
    Url = {http://dx.doi.org/10.1137/S106482750139738X},
    Volume = {24},
    Year = {2003}}

~~~
*/



