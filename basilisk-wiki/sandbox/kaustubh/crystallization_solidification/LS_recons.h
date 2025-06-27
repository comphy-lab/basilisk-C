/**
# Reconstruction of a velocity off an interface
*/  
#include "alex_functions.h"
#include "weno5.h"
/**
This include is to use my_minmod
*/

void buildCache(scalar dist, double NB_width, Cache * c) {
  foreach (serial, noauto)
    if (fabs(dist[]) < NB_width && !interfacial (point, cs))
      cache_append (c, point, 0);
}

// void myFE_WENO5(scalar f, scalar fi, vector u, double dt, scalar dist, double
//   maxd, Cache H, Cache H2){

//   vector upfluxp[], upfluxm[];
//   foreach(){ // WENO5 requires 2 more cells
//     foreach_dimension(){
//       upfluxp.x[] = (fi[1] - fi[])/Delta;
//       upfluxm.x[] = (fi[] - fi[-1])/Delta;
//     }
//   }
//   boundary((scalar *){upfluxm,upfluxp});
//   restriction((scalar *){upfluxm,upfluxp});

//   foreach_cache(H){
//     foreach_dimension()
//       f[] -= dt*u.x[]*
//       (u.x[] >= 0 ? WENO5_x(point, upfluxm.x, -1) : WENO5_x(point,upfluxp.x, 1));
//   }
//   boundary({f});
//   restriction({f});
// }


void myFE(scalar f, scalar fi, vector ndist, double deltat, scalar dist, 
  double maxd, Cache H){
/**
Simple Forward Euler scheme.
*/

  foreach_cache(H){
    foreach_dimension(){
      double graplus  = WENOdiff_x(point, fi,-1);
      double graminus = WENOdiff_x(point, fi,1);
      f[] -= deltat * ( max(ndist.x[],0.)*graplus
        + min(ndist.x[],0.)*graminus);
    }
  }
  boundary({f});
  restriction({f});
}

/**
Runge Kutta 2 advection scheme
*/
void myRK2(scalar f, vector ndist, double deltat, scalar dist, double maxd,
  Cache H){
  scalar f1[],f2[],f3[];

  foreach(){
    f1[] = f[];
    f2[] = f[];
  }
  boundary({f1,f2});
  restriction({f1,f2});

  myFE(f2, f1, ndist, deltat, dist, maxd, H);
  
  foreach(){
    f3[] = f2[];
  }
  boundary({f3});
  restriction({f3});

  myFE(f3, f2, ndist, deltat, dist, maxd, H);
  foreach(){
    f[] = (f3[] + f1[])*0.5;
  }
  boundary({f});
  restriction({f});
}

/**
Runge Kutta 3 advection scheme
*/


void myRK3(scalar f, vector ndist, double deltat, scalar dist, double maxd,
  Cache H){
  // Cache H, Cache H2){
  static const double coeff[2][2] = {
    {0.75, 0.25},
    {1./3.,2./3.}
  };
  scalar f1[],f2[];

  foreach(){
    f1[] = f[];
  }
  boundary({f1});
  restriction({f1});

  myFE(f1, f, ndist, deltat, dist, maxd, H);
  // myFE_WENO5(f1, f, ndist, deltat, dist, maxd, H , H2);
  
  foreach(){
    f2[] = f1[];
  }
  boundary({f2});
  restriction({f2});

  myFE(f2, f1, ndist, deltat, dist, maxd, H);
  // myFE_WENO5(f2, f1, ndist, deltat, dist, maxd, H , H2);
  foreach(){
    f1[] = coeff[0][0]*f[]+coeff[0][1]*f2[];
    f2[] = f1[];
  }
  boundary({f1,f2});
  restriction({f1,f2});
  myFE(f2, f1, ndist, deltat, dist, maxd, H);
  // myFE_WENO5(f2, f1, ndist, deltat, dist, maxd, H , H2);
  foreach(){
    f[] = coeff[1][0]*f[] + coeff[1][1]*f2[];
  }
  boundary({f});
  restriction({f});

}



void buildNdist(scalar dist, vector ndist, Cache H){

  foreach()
    foreach_dimension()
      ndist.x[] = 0.;

  foreach_cache(H){

    double sum=1e-15;
    coord d3;
    coord d4;
    foreach_dimension(){
      d3.x = (dist[1,0,0]-dist[-1,0,0])/2.;
      d4.x = dist[-2,0,0]/12. - 2./3*dist[-1,0,0] 
          + 2./3*dist[1,0,0] - dist[2,0,0]/12.;
    }

    double d = max(sum,min(norm2(d3),norm2(d4)));
    foreach_dimension()
        ndist.x[] = d3.x;
    if(d == norm2(d4)){
      foreach_dimension()
        ndist.x[] = d4.x;
    }

    foreach_dimension(){
      ndist.x[] = sign2(dist[])*ndist.x[]/(d+1.e-15);
    }
  }

  boundary((scalar *) {ndist});
  restriction((scalar *) {ndist});
}

void corrInterf(scalar f, scalar f2, scalar fi, scalar cs, face vector fs,
  Cache cut_cells){
  scalar error[];
  foreach_cache(cut_cells){
    coord n, p;
    embed_geometry(point, &p, &n);

/**
Because our level-set function of cell-centered, one has to choose
the interpolation stencil according to the position of the face
centroid.
*/

    coord p_interp = p;

/**
We do a quadratic interpolation with a correction using the interpolation
coefficient of the central cell in the 3x3 stencil.
*/
#if dimension==2
    double f_temp = mybiquadratic(point , f2, p_interp,0);
#else // dimension ==3
    double f_temp = mytriquadratic(point , f2, p_interp);
#endif
    error[]  = fi[] - f_temp;
    foreach_dimension(){
      error[] /=(1-sq(p.x));
    }
    f[] += error[];
  }
  boundary({f});
  restriction({f});
}


/**
Calculation of the interpolation between the reconstructed continuous phase
change velocity and the velocity originally calculated on the interface
centroids.
*/
double interp_error(scalar f, scalar fi, Cache cut_cells){
  double s2 = 0;
  foreach_cache (cut_cells, reduction(max:s2)){
    coord n, p;
    embed_geometry(point, &p, &n);
    int Stencil[2];
    InterpStencil(p, Stencil);
    coord p_interp = p;
    
#if dimension  == 2
    double f_temp = fabs(fi[] - mybiquadratic(point , f, p_interp,0));
#else
    double f_temp = fabs(fi[] - mytriquadratic(point , f, p_interp));
#endif
    s2 = max(s2,f_temp);
  }
  return s2;
}

struct LS_recons {
  scalar dist;
  double deltat;
  scalar * LS_speed;
  double tolerance;
  double * err;
  int nb_iter;
  scalar cs;
  face vector fs;
  double NB_width;
};

void recons_speed(struct LS_recons p){
  scalar dist       = p.dist;
  double deltat     = p.deltat;
  scalar * LS_speed = p.LS_speed;
  double tolerance  = p.tolerance; // default tolerance is 1.e-8
  double * err      = p.err;
  int nb_iter       = p.nb_iter; // default is 35
  scalar cs         = p.cs;
  face vector fs    = p.fs;
  double NB_width   = p.NB_width;

#if LS_perf
double start; 
double end; 
#endif
/**
The phase change field $v_{pc}$ is only defined in interfacial cells, this
function modifies the velocity modifies the velocity field in the vicinity of
the interface with the following set of equations:

$$
\left\{\begin{array}{c}
\dfrac{\partial v_{pc}}{\partial t}  + S(\phi) \frac{\nabla \phi}{|\nabla\phi|}. \nabla v_{pc}= 0\\
\mathcal{L}{v_{pc}}(x_{\Gamma}) = \left. v_{pc}\right|_{\Gamma} 
\end{array}\right.
$$

For the first equation, we use the method described by [Peng et al., 1999](#peng_pde-based_1999) 

First, we calculate : 
$$ 
S(\phi) \frac{\nabla \phi}{|\nabla \phi|}
$$

using a centered approximation (second or fourth-order accurate).
*/

  *err = 0.;
  if(tolerance == 0){
    tolerance = 1.e-5;
  }
  if(NB_width == 0){
    NB_width = 4.*L0/(1 << (grid->depth));
  }

  double second_tolerance = 1.e-2;

  vector ndist[];
  scalar tag[];

  Cache H = {0};
  Cache cut_cells = {0};
  foreach(serial, noauto)
    if (cs[]*(1. - cs[]) != 0.)
      cache_append (&cut_cells, point, 0);
  cache_shrink (&cut_cells);

  buildCache(dist,NB_width,&H);
  cache_shrink(&H);

  buildNdist(dist,ndist,H);
  boundary((scalar *){ndist});
  restriction((scalar *){ndist});


/**
Then, do the advection a few times to extend the velocity from the 
surface along the normal to the surface.

*/ 

  int ii = 0;
  if(nb_iter < 5){
    nb_iter = 5;
  }

/**
We iterate on the `scalar * LS_speed  = vpc` pointer, therefore a
`foreach_dimension`
wouldn't work.
*/

/**
TODO: Double residual calculation.
Add a mask for 2-3 cells next to the interface, and do a
residual calculation, once it has converged, apply the correction step in the
interfacial cells, and calculate a residual for the interpolation error (second
residual is already calculated).
**/
  scalar myres[];
  for (scalar f in LS_speed){

    double velo_res = 0., s2 = 0.;
    scalar fi[];


    foreach_cache(cut_cells, reduction(max:velo_res)){
      velo_res  = max(velo_res,fabs(f[]));
      fi[] = f[];
    }
    // boundary({fi});

    velo_res *= second_tolerance;
    // fprintf(stderr, "##LS_recons velo_res %g\n", velo_res);
  /**
  First, we advect the velocity without modifying the velocity in interfacial
  cells.*/
    for (ii=0; ii<=nb_iter; ii++){
#if LS_perf
      start = omp_get_wtime(); 
#endif 
/**
RK3 advection scheme. Collocated variables...
*/

      foreach_cache(H)
        myres[] = f[];
      // boundary({myres});
      // restriction({myres});

      myRK3(f,ndist,deltat,dist, NB_width, H);
      // myRK3(f,ndist,deltat,dist, NB_width, H,H2);
      bool stop  = 0;
      double residual= 0.;
      foreach_cache(H, reduction(max:residual)){
        residual = max(residual,fabs(f[]-myres[]));
      }
      // fprintf(stderr, "##LS_recons residual %g\n", residual);
  /**

This second part is a bit tricky since our system has been initialized with the
velocity calculated on the face centroids $\left. v_{pc} \right|_{\Gamma}$ (see
Figure 1, in red). Thus, while we advect the velocity off the interface, we use
an iterative method to correct the cell-centered values upc (in blue) so that their
interpolation matches the correct interfacial boundary condition
$\left. v_{pc}\right|_{\Gamma}$.

![Fig. 1: Initial face values (red) and cell-centered values (blue) interfacial cells](interp_abstract.svg)
 
When we have reached convergence, we check that our interpolation re-gives the
interfacial velocity, if not, we correct the values in interfacaial cells.

   */
      if(residual<velo_res){
        stop = 1;
        scalar f2[];
        foreach(){
          f2[] = f[];
        }
        boundary({f2});
        restriction({f2});
        corrInterf(f,f2,fi,cs,fs,cut_cells);
        s2 = interp_error(f, fi, cut_cells);
        // fprintf(stderr, "### CORR %g\n", s2);
      }
      if(cut_cells.n==0){
        *err = 0;
        break;
      }
#if DEBUG > 1
      if(stop){
        fprintf(stderr, "%g %g\n", s2, tolerance);
      }
#endif
      if(stop){
        if(s2 < tolerance){
          *err = s2;
#if DEBUG
          fprintf(stderr, "### CONVERGED %g %g\n", s2, tolerance);
#endif
          break;
        }
        else{
#if DEBUG
          fprintf(stderr, "### NOT CONVERGED %g %g\n", s2, tolerance);
#endif
        }
      }
#if LS_perf
      end = omp_get_wtime(); 
      fprintf(stderr,"LS_recons %f seconds\n", end - start);
#endif 
    }
    *err = max(*err,s2);

  }
  free(H.p);
  free(cut_cells.p);
}
/**
##TODO 
Add two test cases : 
- planar interface, check that velocity is constant.
- circular interface, radial velocity of constant magnitude.

Control the order of convergence, speed of convergence.
AND THEN try to optimize, but not too early...


## References

~~~bib

@article{peng_pde-based_1999,
  title = {A {PDE}-Based Fast Local Level Set Method},
  volume = {155},
  issn = {0021-9991},
  url = {http://www.sciencedirect.com/science/article/pii/S0021999199963453},
  doi = {10.1006/jcph.1999.6345},
  pages = {410--438},
  number = {2},
  journaltitle = {Journal of Computational Physics},
  author = {Peng, Danping and Merriman, Barry and Osher, Stanley and Zhao, Hongkai and Kang, Myungjoo},
  urldate = {2019-09-09},
  date = {1999-11}
}
~~~
*/
