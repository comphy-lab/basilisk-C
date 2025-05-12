/**
A multilayer compatible convection diffusion solver. 
We start by defining the diffusion scheme and a small tolerance function (To avoid dividing by zero). 

**Schemes**

 - scheme 0 is the first order upwind scheme. Fast, but first order and best suited for high Peclet numbers.
 
 - scheme 1 is the van Albada TVD scheme (from *A comparative study of computational methods in cosmic gas dynamics* by van Albada et al. 1982. This scheme is second order in space, and stable at all Peclet numbers. 
*/

#define scheme 1
#define EPS_ML_TRACER 1e-8

/**
Here we define a lot of package specific variables.
Of these, only the ml_tracer struct is intended to be used by the user. 
*/
typedef int Ml_Tracer;

typedef struct {
  double kappa_c;
  double tol;
  double rho_c;
  scalar c;
  double (*bf)(Point point, double t, Ml_Tracer tracer);
  double (*Sf)(Point point, double t, Ml_Tracer tracer);
} ml_tracer;

typedef struct {
  double aP;
  double aW;
  double aE;
  double aT;
  double aB;
  #if dimension > 1
    double aN;
    double aS;
  #endif
  double Sp;
  double Su;
  double div;
} TDMA_fluxes;

/**
Here we perform some book-keeping where we initialize a list of all tracers the user initializes. Each tracer gets a unique ID corresponding to its index in the tracer list array ```tl```. 
*/
typedef ml_tracer * ml_tracers;
ml_tracers tl = NULL;
long unsigned int num_tracers = 0;

/**
 The vertical velocity w is not defined in the hydrostatic solver. 
 Hence, we define the vertical velocity and then set it to zero.
*/
#if NH
#else
  scalar w;
#endif

event make_w(i=0){
  #if NH
  #else
    w = new scalar[nl];
    foreach(){
      foreach_layer(){
        w[] = 0.0;
      }
    }
  #endif
}

event clean_w(t=end){
  #if NH
  #else
    delete((scalar *) {w});
  #endif
}

/** 
 Here starts the first important function of the tracer library, namely the init_tracer function. This function takes four arguments, namely

- ```scalar c``` , a user defined scalar field which acts as the initial condition of the tracer. 
- ```double kappa_c``` is the diffusion coefficient of the tracer. 
- ```double rho_c``` is the density of the tracer. Currently, the relative density between the tracer and the fluid in the multilayer solver is not considered. 
- ```double tol``` is the numerical tolerance. The algorithm is an iterative solver, and this sets the cutoff point for the solver, as the iterations will terminate when the residual in terms of the $L^2$-norm is less than this tolerance. 

The functions bf0 and Sf0 are the default boundary condition and source functions, and can be overwritten by a user using the set_tracer_BC and set_source_function functions. The boundary and source functions need to take three arguments, namely the spatial coordinate ```Point point```, the time ```double t``` and the tracer index ```Ml_Tracer tracer```. 
*/
double bf0(Point point,  double t, Ml_Tracer tracer){
  return 0.0;
}
double Sf0(Point point,  double t, Ml_Tracer tracer){
  return 0.0;
}
Ml_Tracer init_tracer(scalar c, double kappa_c, double rho_c, double tol){
  ml_tracer tracer;
  tracer.kappa_c = kappa_c; 
  tracer.tol = tol;
  tracer.rho_c = rho_c;
  tracer.c = c;
  //boundary and source functions
  tracer.bf = &bf0;
  tracer.Sf = &Sf0;
  num_tracers++;
  tl = realloc (tl, num_tracers*sizeof(ml_tracer));
  tl[num_tracers-1] = tracer;

  return num_tracers-1;
}

Ml_Tracer set_tracer_BC(Ml_Tracer _i_tracer,double (*bf)(Point point, double t, Ml_Tracer tracer)){
  tl[_i_tracer].bf = bf;
  return _i_tracer;
}

Ml_Tracer set_source_function(Ml_Tracer _i_tracer, double (*Sf)(Point point, double t, Ml_Tracer tracer)){
  tl[_i_tracer].Sf = Sf;
  return _i_tracer;
}

double _tracer_boundary_function(Ml_Tracer _i_tracer, Point point,  double t){
  ml_tracer tracer = tl[_i_tracer];
  return tracer.bf(point, t, _i_tracer);
}

double _tracer_source_function(Ml_Tracer _i_tracer, Point point, double t){
  ml_tracer tracer = tl[_i_tracer];
  return tracer.Sf(point, t, _i_tracer);
}

/**
Here we define some useful (and some potentially useful) shortcuts, and the foreach_ml_tracer() iterator, which allow a user to iterate over all initialized tracers.
*/

#define tml() tl[_i_tracer]

@def TRACER_VARIABLES
  double kappa_c = tml().kappa_c; NOT_UNUSED(kappa_c);
  double tol     = tml().tol; NOT_UNUSED(tol);
  scalar c       = tml().c; NOT_UNUSED(c);
  double rho_c   = tml().rho_c; NOT_UNUSED(rho_c);
@

@def foreach_ml_tracer()
OMP_PARALLEL() {
  int _i_tracer;
  for (_i_tracer = 0; _i_tracer < num_tracers; _i_tracer++){
    TRACER_VARIABLES

@
@def end_foreach_ml_tracer()
  }
}
@

/**
This package solves the convection diffusion equations using the [Thomas algorithm](https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm), and here we define the row reduction function. 
*/

void TDMA(double * a_w, double * a_c, double * a_e, double * b, double * c, int n){
  double  a_d[n], rhs[n];
  int _i;
  double d;
  for (_i = 0; _i < n; _i++){
    a_d[_i] = a_c[_i];
    rhs[_i] = b[_i];
  }
  
  for (_i = 1; _i < n; _i++){
    d  = a_w[_i]/ a_d[_i-1];
    a_d[_i] -= d*a_e[_i];
    rhs[_i] -= d*rhs[_i-1];
  }

  for (_i = n - 1; _i > 0; _i--) {
    d  = a_e[_i-1]/a_d[_i];
    rhs[_i-1] -= d*rhs[_i];
  }

  for (_i = 0; _i < n; _i++) {
    c[_i] = rhs[_i]/a_d[_i];
  }
}

/**
TVD schemes use a TVD-limiter function. More functionality can easily be added here to support more schemes
*/

#if scheme == 0
  // Upwind
  double TVD_limiter(double r){
    return 0.0;
  }
#elif scheme == 1
  // Van Albada
  double TVD_limiter(double r){
    return (r + sq(r))/(1 + sq(r));
  }
#endif

/**
This (way too big) function computes the tracer fluxes entering each cell.
Currently, this is done manually for each face, and can in all likelihood be simplified dramatically. 
*/
TDMA_fluxes compute_fluxes(Point point, int l, int _i_tracer, scalar c_0){
  double rho_c = tml().rho_c;
  double kappa_c = tml().kappa_c;
  double aW = 0, aE = 0, aT = 0, aB = 0;
  double Su = 0.0;
  double Sp = 0.0;
  double aP = 1;
  double div = 0;
  #if dimension > 1
    double aN = 0, aS = 0;
  #endif
  if (h[] > dry) {
    // Implement with matrix entries from Versteeg
    double Dh = kappa_c/(Delta);
    double Dv = kappa_c/h[0,0,l];
    double Fw = 0, Fe = 0, Fb = 0, Ft = 0;
    Su = _tracer_source_function(_i_tracer, point, t);
    #if dimension > 1
      double Fn = 0, Fs = 0;
    #endif
    // West side
    if (point.i > GHOSTS) {
      Fw = rho_c*u.x[-1,0,l];
      aW = (Dh +  max(Fw,0))/Delta;
    } else {
      Fw = rho_c*u.x[0,0,l];
      aW = 0.0;
      Sp -= (2*Dh + max(Fw,0))/Delta;
      //Add BC effects 
      Su += (2*Dh + max(Fw,0))*_tracer_boundary_function(_i_tracer, point ,t)/Delta;
    }

    // East side
    if (point.i < point.n + GHOSTS -1){
      Fe = rho_c*u.x[1,0,l];
      aE = (Dh + max(-Fe,0))/Delta;
    } else {
      Fe = rho_c*u.x[0,0,l];
      aE = 0.0;
      Sp -= (2*Dh + max(-Fe,0))/Delta;
      Su += (2*Dh + max(-Fe,0))*_tracer_boundary_function(_i_tracer, point ,t)/Delta;
    }
    #if dimension > 1
      // North side
      if (point.j < point.n + GHOSTS- 1){
        Fn = rho_c*u.y[0,1,l];
        aN = (Dh + max(-Fn,0))/Delta;
      } else {
        Fn = rho_c*u.y[0,0,l];
        aN = 0.0;
        Sp -= (2*Dh + max(-Fn,0))/Delta;
        Su += (2*Dh + max(-Fn,0))*_tracer_boundary_function(_i_tracer, point ,t)/Delta;
      }

      // South side
      if (point.j > GHOSTS){
        Fs = rho_c*u.y[0,-1,l];
        aS = (Dh + max(Fs,0))/Delta;
      } else {
        Fs = rho_c*u.y[0,0,l];
        aS = 0.0;
        Sp -= (2*Dh + max(Fs,0))/Delta;
        Su += (2*Dh + max(Fs,0))*_tracer_boundary_function(_i_tracer, point ,t)/Delta;
      }
    #endif
    // Top side
    if (l < nl - 1 && nl > 1){
      Ft = rho_c*w[0,0,l+1];
      aT = (Dv + max(Ft,0))/h[0,0,l];
    } else {
      Ft = 0.0;
      aT = 0.0;
      Sp -= (2*Dv + Ft)/h[0,0,l];
      Su += (2*Dv + Ft)*c_0[0,0,l]/(h[0,0,l]);
    }

    // Bottom side
    if (l > 0 && nl > 1) {
      Fb = rho_c*w[0,0,l-1];
      aB = (Dv + Fb)/h[0,0,l];
    } else {
      aB = 0.0;
      Fb = 0.0;
      Sp -= (2*Dv + Fb)/h[0,0,l];
      Su += (2*Dv + Fb)*c_0[0,0,l]/h[0,0,l];
    }
    div = Fe + Ft - Fw - Fb;
    #if dimension > 1
      div += Fn - Fs;
    #endif
  }
  TDMA_fluxes cf;
  cf.aP = aP; 
  cf.aW = aW;
  cf.aE = aE;
  cf.aT = aT;
  cf.aB = aB;
  cf.div = div;
  #if dimension > 1
    cf.aN = aN;
    cf.aS = aS;    
  #endif
  cf.Sp = Sp;
  cf.Su = Su;
  return cf;
}

/**
Finally, this is where the deferred correction term in the TVD algorithm is computed. In the upwind case, this function rturns 0. Similarily to the function above, this code can probably be simplified drstically using e.g. foreach_dimension()
*/
double TVD_deferred_correction(Point point, scalar c, int l, int _i_tracer){
  double Sdc = 0;
  
  #if scheme == 0
  #else
    double rho_c = tml().rho_c;
    double Fe = point.i < point.n + GHOSTS ? rho_c*u.x[ 1,0,l] : rho_c*u.x[0,0,l];
    double Fw = point.i > GHOSTS ? rho_c*u.x[ -1,0,l] : rho_c*u.x[0,0,l];
    double Fb = l > 0 ? rho_c*w[0,0,l-1] : 0;
    double Ft = l < nl - 1 ? rho_c*w[0,0,l+1] : 0;
    double re, rw, rt, rb;
    #if dimension > 1
      double rn, rs;
      double Fn = rho_c*u.y[0,1,l], Fs = rho_c*u.y[0,-1,l];
      int alpha_n = Fn > 0 ? 1 : 0;
      int alpha_s = Fs > 0 ? 1 : 0;
    #endif
    int alpha_w = Fw > 0 ? 1 : 0;
    int alpha_e = Fe > 0 ? 1 : 0; 
    int alpha_t = Ft > 0 ? 1 : 0;
    int alpha_b = Fb > 0 ? 1 : 0;
    double c0 = 2*_tracer_boundary_function(_i_tracer, point ,t) - c[0,0,l];

    if (alpha_e > 0) {
      if (point.i > GHOSTS && point.i < point.n + GHOSTS -1) {
        re = (c[0,0,l] - c[-1, 0,l])/(c[1,0,l] - c[0,0,l] + EPS_ML_TRACER);
      } else if (point.i == GHOSTS){
        re = (c[0,0,l] - c0)/(c[1,0,l] - c[0,0,l] + EPS_ML_TRACER);
      } else {
        re =  0.0;
      }
      // Do not need to check bounds, as TVD_limiter is zero at point.i == point.n + GHOSTS
    } else {
      if (point.i >= GHOSTS && point.i < point.n + GHOSTS - 2) {
        re = (c[1,0,l] - c[2,0,l])/(c[0,0,l] - c[1,0,l] + EPS_ML_TRACER);
      } else if (point.i == point.n + GHOSTS -2){
        re = (c[1,0,l] - c0)/(c[0,0,l] - c[1,0,l] + EPS_ML_TRACER);
      } else {
        re =  0.0;
      }
      // Do not need to check bounds, as TVD_limiter is zero at point.i == point.n + GHOSTS
    }
    Sdc -= 0.5*fabs(Fe)*TVD_limiter(re)*(c[1,0,l] - c[0,0,l]);
    if (alpha_w > 0) {
      if (point.i > GHOSTS +  1 && point.i < point.n + GHOSTS) { 
        rw = (c[-1,0,l] - c[-2, 0,l])/(c[0,0,l] - c[-1,0,l] + EPS_ML_TRACER);
      } else if (point.i == GHOSTS + 1){
        rw = (c[-1,0,l] - c0)/(c[0,0,l] - c[-1,0,l] + EPS_ML_TRACER);
      } else {
        rw =  0.0;
      }
    } else {
      if (point.i > GHOSTS && point.i < point.n + GHOSTS - 1) { 
        rw = (c[0,0,l] - c[1,0,l])/(c[-1,0,l] - c[0,0,l] + EPS_ML_TRACER);
      } else if (point.i == point.n + GHOSTS -1){
        rw = (c[0,0,l] - c0)/(c[-1,0,l] - c[0,0,l] + EPS_ML_TRACER);
      } else {
        rw =  0.0;
      }
    }
    Sdc += 0.5*fabs(Fw)*TVD_limiter(rw)*(c[0,0,l] - c[-1,0,l]);
    #if dimension > 1
      if (alpha_n > 0) {
        if (point.j > GHOSTS && point.j < point.n + GHOSTS -1) { 
          rn = (c[0,0,l] - c[0, -1,l])/(c[0,1,l] - c[0,0,l] + EPS_ML_TRACER);
        } else if (point.j == GHOSTS){
          rn = (c[0,0,l] - c0)/(c[0,-1,l] - c[0,0,l] + EPS_ML_TRACER);
        } else {
          rn =  0.0;
        }
      } else {
        if (point.j >= GHOSTS && point.j < point.n + GHOSTS - 2) { 
          rn = (c[0,1,l] - c[0,2,l])/(c[0,0,l] - c[0,1,l] + EPS_ML_TRACER);
        } else if (point.j == point.n + GHOSTS - 2){
          rn = (c[0,1,l] - c0)/(c[0,0,l] - c[0,1,l] + EPS_ML_TRACER);
        } else {
          rn =  0.0;
        }
      }
      Sdc -= 0.5*fabs(Fn)*TVD_limiter(rn)*(c[0,1,l] - c[0,0,l]);
      if (alpha_s > 0) {
        if (point.j > GHOSTS +  1 && point.j < point.n + GHOSTS) { 
          rs = (c[0,-1,l] - c[0, -2,l])/(c[0,0,l] - c[0,-1,l]+ EPS_ML_TRACER);
        } else if (point.j == GHOSTS + 1){
          rs = (c[0,-1,l] - c0)/(c[0,0,l] - c[0,-1,l] + EPS_ML_TRACER);
        } else {
          rs =  0.0;
        }
      } else {
        if (point.j > GHOSTS && point.j < point.n + GHOSTS - 1) { 
          rs = (c[0,0,l] - c[0,1,l])/(c[0,-1,l] - c[0,0,l] + EPS_ML_TRACER);
        } else if (point.j == point.n + GHOSTS -1){
          rs = (c[0,0,l] - c0)/(c[0,-1,l] - c[0,0,l] + EPS_ML_TRACER);
        } else {
          rs =  0.0;
        }
      }
      Sdc += 0.5*fabs(Fs)*TVD_limiter(rs)*(c[0,0,l] - c[0,-1,l]);
    #endif
    if (alpha_t > 0) {
      if (l > 0 && l < nl - 1) { 
        rt = (c[0,0,l] - c[0,0,l-1])/(c[0,0,l-1] - c[0,0,l] + EPS_ML_TRACER);
      } else if (l == nl){
        rt = (c[0,0,l])/(c[0,0,l-1] - c[0,0,l] + EPS_ML_TRACER);
      } else {
        rt =  0.0;
      }
    } else {
      if (l >= 0 && l < nl - 2) { 
        rt = (c[0,0,l+1] - c[0,0,l+2])/(c[0,0,l] - c[0,0,l+1] + EPS_ML_TRACER);
      } else if (l == nl-2){
        rt = (c[0,0,l+1] - c0)/(c[0,0,l] - c[0,0,l+1] + EPS_ML_TRACER);
      } else {
        rt =  0.0;
      }
    }
    if (l < nl - 1)
      Sdc -= 0.5*fabs(Ft)*TVD_limiter(rt)*(c[0,0,l+1] - c[0,0,l]);
    if (alpha_b > 0) {
      if (l > 1 && l < nl) { 
        rb = (c[0,0,l-1] - c[0, 0,l-2])/(c[0,0,l] - c[0,0,l-1] + EPS_ML_TRACER);
      } else if (l == 1){
        rb = (c[0,0,l-1] - c0)/(c[0,0,l] - c[0,0,l-1] + EPS_ML_TRACER);
      } else {
        rb =  0.0;
      }

    } else {
      if (l > 0 && l < nl - 2) { 
        rb = (c[0,0,l] - c[0,0,l+1])/(c[0,0,l-1] - c[0,0,l] + EPS_ML_TRACER);
      } else if (l == nl - 2){
        rb = (c[0,0,l] - c0)/(c[0,0,l-1] - c[0,0,l] + EPS_ML_TRACER);
      } else {
        rb =  0.0;
      }
    }
    if (l > 1)
      Sdc += 0.5*fabs(Fb)*TVD_limiter(rb)*(c[0,0,l] - c[0,0,l-1]);
    #endif
  return Sdc;
}

/**
 Cleanup and importation of the transport-function
*/
#if dimension == 1
  #include "mltracer2D.h"
#else
  #include "mltracer3D.h"
#endif


event free_tracers (t = end){
  free (tl);
  tl = NULL;
}