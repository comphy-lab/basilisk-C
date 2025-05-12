/**
Rayleigh-Taylor instability using the Semushin method for
advecting the interface. A one meter wide, 
four meters high rectangular domain is discretized using a 64 × 256 grid.
 */

#define Nl 256
#define length 4

// #define ADAPT
// #define NORM
#define SEMU2VOF

// #define OUT_X
// #define OUT_HDF5
#define OUT_GNUPLOT

#define IT (20)
#define TEST "Rayleigh"

/**
 The volume fraction is stored in scalar field `f` which is listed as
 an *interface* for the VOF solver. We do not advect any tracer with
 the default (diffusive) advection scheme of the advection solver. */

scalar f[], cf[];
scalar * interfaces = {f}, * tracers = NULL;
scalar tmp1[];

face vector alphav[], muv[], av[];
scalar rhov[];

#include "navier-stokes/centered.h"
#include "lchirco/semushin.h"
#include "tension.h"


#define rho1 1.225
#define rho2 0.1694
#define mu1  .00313
#define mu2  .00313
#define rho(f) (clamp(f,0,1)*(rho1 - rho2) + rho2)
#define mu(f) (1./(clamp(f,0,1)*(1./(mu1) - 1./(mu2)) + 1./(mu2)))
#define sigmaf 1e-9


int count_it = 0;
const char *testName = TEST;

int main() {
  TOLERANCE = 1e-4;
  size(length);
  origin (0., 0.);
  init_grid (Nl);
  a = av;
  mu = muv;
  alpha = alphav;
  rho = rhov;
  f.sigma = sigmaf;
  f.tracers = {cf};
  DT = 1.e-2;
  run();
}


event properties (i++) {
  #if TREE
  f.prolongation = refine_bilinear;
  boundary ({f});
  #endif
  
  foreach_face() {
    double ff = (f[] + f[-1])/2.;
    alphav.x[] = fm.x[]/rho(ff);
    muv.x[] = fm.x[]*mu(ff);
  }
  foreach()
    rhov[] = cm[]*rho(f[]);
  
  #if TREE
  f.prolongation = fraction_refine;
  boundary ({f});
  #endif
}

/**
 We define the levelset function $\phi$ on each vertex of the grid and
 compute the corresponding volume fraction field $f$ and face fractions
 $s\_tmp$. Then, we initialize the position of the interface points 
 (or Semushin fractions) on every edge. We assume $s = 0$ if the marker point
 coincides with the left or bottom vertices and $s = 1$ for the right
 and top ones. */


event init (i = 0) {
  
  mask (x > 1. ? right :
        x < 0. ? left : none);

  vertex scalar phi[];
  foreach_vertex(){
    phi[] = - ( 2.001 - y - 0.05*sin(pi*x)) ;
  }
  boundary ({phi});
  fractions (phi, f, s_tmp);
  boundary ((scalar *){s_tmp});
  
  /** 
   In `full` cells (volume fraction $f=1$) there is no interface and then the
   Semushin fraction is zero. To obtain the Semushin interface position we have
   to `rotate` the face fractions $s\_tmp$, since they are computed with respect
   to the normal, while in Semushin method they always refer to fixed vertexes
   (e.g. left and bottom). */
  
  foreach_face(y){
    if (exist(s_tmp.y[])) {    
      s.y[] = s_tmp.y[];
//       if (x < xcenter)  
//         s.y[] = 1. - s.y[];
    }
  }
  foreach_face(x){
    if (exist(s_tmp.x[])) { 
      s.x[] = s_tmp.x[];
//       if (y < ycenter)  
        s.x[] = 1. - s.x[];
    }
  }
  boundary ((scalar *){s});
  
/** 
   The boundary condtions for the scalar field used for the Semushin markers
   are set. They are not wokring correctly! */
  
  snew.t[left] = 0.;
  snew.t[right] = 0.;
  snew.n[left] = s.x[1,0];
  snew.n[right] = s.x[0,0];
  boundary ((scalar *){snew}); 

}

/** 
The vertical acceleration is added here. */

event acceleration (i++) {
  foreach_face(y)
    av.y[] -= 9.81;
  boundary ((scalar *){av});
}

/** The timestep `dt` and the velocity field are set. */
event stability (i++, first) {
  
  tTime += dt;
  
  if (pid()==0){
    printf("i: %d  Time: %1.4f  test: %s \n",i, tTime, testName);
  }
}


#ifdef ADAPT
event adapt (i++) {
  
  foreach()
    tmp1[] = 0.;
  boundary({tmp1});
  
  foreach(){
    if (exist(s.x[]) || exist (s.x[1,0]) || exist (s.x[-1,0]) 
      || exist (s.x[2,0]) || exist (s.x[-2,0])||  exist(s.y[]) 
      || exist(s.y[0,1]) ||  exist(s.y[0,-1]) || exist(s.y[0,2]) 
      || exist(s.y[0,-2])){
      tmp1[] = noise();
      }
  }
  boundary({tmp1}); 
  
  adapt_wavelet ({tmp1}, (double[]){0}, maxlevel = MAXLEV-1, 
                 minlevel = MAXLEV-2);
//   adapt_wavelet ({f}, (double[]){0.001}, maxlevel = 9, minlevel = 6);
}
#endif

/**  We can output the results in HDF5 format for visualization with Paraview.*/
#ifdef OUT_HDF5
#include "lchirco/output_hdf.h"
event calc_and_print (t = 0; i+=1) {
  
  scalar * list = {s.x,s.y,f,ff};
  vector * vlist = {u};
  
  char buf[100];
  char time[100];
  char itime[100];
  count_it++;
  
  snprintf(buf, sizeof(buf), "out_%05d.xmf", count_it);
  snprintf(itime, sizeof(itime), "%05d", count_it);
  snprintf(time, sizeof(time), "%06d", i);
  
  FILE * fp = fopen(buf, "w");
  output_xmf_h5_foreach(list, vlist, 64, fp, itime, time); 
  fclose(fp);
  
  char name[500];
  sprintf(name, "%s-%.5d-CC.vtk", "interface", i);
  output_paraview_IF(name, t, 0.0, f);   
}
#endif




