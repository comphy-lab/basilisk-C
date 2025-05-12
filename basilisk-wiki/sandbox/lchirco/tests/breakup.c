#include "axi.h"
#include "navier-stokes/centered.h"
#include "vof.h"
#include "tension.h"
#include "output_hdf.h"
#include "signature.h"

scalar f[], cf[], phii[], sign[], M[];
scalar * interfaces = {f};

const int max_change = 10; 
bool large = true;

#define MAXLEV 13
#define MINLEV 2
#define SIGN_LEV 11

#define We 2.5
#define Re 1090.
#define rd 1110.
#define rv 90.9

#define R 1
#define H 12*R
#define U0 0.1

#define rho2 1.
#define rho1 (rho2*rd)
#define mu2  (rho2*U0*R/Re)
#define mu1  mu2*rv
#define sigmaf (rho2*R*sq(U0)/We)

#define rho(f) (clamp(f,0,1)*(rho1 - rho2) + rho2)
#define mu(f) (1./(clamp(f,0,1)*(1./(mu1) - 1./(mu2)) + 1./(mu2)))

#define We_c (rho2*sq(U0)*R/sigmaf)
#define Re_c (rho2*R*U0/mu2)
#define Oh_c (mu1/sqrt(sigmaf*rho1*R))
#define rd_c (rho1/rho2)
#define rv_c (mu1/mu2)

face vector alphav[], muv[], av[];
scalar rhov[];

#define Nl 32

/**
  We set the grid dimension and position.
*/

u.n[left]  = dirichlet(U0);
u.t[left]  = dirichlet(0);
p[left]    = neumann(0);
u.n[right] = neumann(0);
p[right]   = dirichlet(0);

timer tt;
int count_it = 0;

int main (int argc, char * argv[]) {
  if (pid()==0)  
    printf("Computed Re %g We %g Oh %g rd %g rv %g \n", Re_c, We_c, Oh_c, rd_c, rv_c);
  size (H);
  origin (-3*R, 0.);  //2D axisymmetryc
  f.tracers = {cf};
  init_grid(Nl);
  a = av;
  mu = muv;
  alpha = alphav;
  rho = rhov;
  f.sigma = sigmaf;
  DT = 5e-3;
  tt = timer_start();
  run();
}

event calc_and_print (t = 0.0; t += 5; t <= 1800.){
  
  foreach(){
    phii[] = 2*f[] - 1;
    sign[] = 7;
  }
  
    /** 
   We choose at which level (`l_sign`) we want to compute
   the quadratic moments and the signature.
   Then, the indicator function `phii` is restricted at level `l_sign`.
   */ 
  
  int l_sign = SIGN_LEV;
  
  for (int ilev = depth() - 1; ilev >= l_sign; ilev--)  
    foreach_level(ilev){
      if(is_refined(cell))
      restriction_average(point, phii);
    }
  
  compute_signature_neigh_level (f, phii, sign, l_sign);
  
  /** 
   The signature `sign` is available only at the level `l_sign`. 
   We need to prolong it onto the finest grid. */
  if (pid()==0)  
    printf("time %g level used for moments %d and depth is %d \n", t, l_sign, depth()); 
  
  for (int ilev = l_sign; ilev < depth(); ilev++)  
    foreach_level(ilev){
      sign.prolongation = phii.prolongation = refine_injection;
      if(is_refined(cell)){
        sign.prolongation (point, sign);
        phii.prolongation (point, phii);
      }
    }

//  change_topology (f, sign, M, l_sign, max_change, large);
  
  /**
     Finally we output the results in HDF5 format.  */
  
  scalar * list = {sign, f};
  vector * vlist = {u};
  
//   output_ppm (sign, linear = true, file = "f.mp4", n = 200);
  char buf[100];
  char time[100];
  char itime[100];
  count_it++;
  
  snprintf(buf, sizeof(buf), "out_%05d.xmf", count_it);
  snprintf(itime, sizeof(itime), "%05d", count_it);
  snprintf(time, sizeof(time), "%06g", t);
  
  FILE * fp = fopen(buf, "w");
  output_xmf_h5_foreach(list, vlist, 64, fp, itime, time); 
  fclose(fp);
  
  char dumpname[80];
  sprintf (dumpname, "snapshot-%g", t);
  dump (dumpname);
  
}


event init (i = 0) {
  #if TREE
  scalar f1[];
  foreach()
    f1[] = (sq(x) + sq(y) + sq(z) - sq(R) <=0); //circle(0, 0, 0);
  boundary ({f1});
  astats s;
  do {
    s = adapt_wavelet ({f1}, (double[]){0.0}, MAXLEV, MINLEV, list = NULL);
    foreach()
      f1[] = (sq(x) + sq(y) + sq(z) - sq(R) <=0);
    boundary ({f1});
  } while (s.nf);
  foreach()
    f[] = (sq(x) + sq(y) + sq(z) - sq(R) <=0);
  boundary ({f});
  #else
  foreach()
    f[] = (sq(x) + sq(y) + sq(z) - sq(R) <=0);
  boundary ({f});
  #endif
  foreach()
    foreach_dimension()
      u.x[] = 0;
  boundary ({f,u,p});
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


#if TREE
event adapt (i++) {
  adapt_wavelet ({f, u}, (double[]){0.00005,0.00001,0.00001,0.00001},
                   MAXLEV, MINLEV);
}
#endif
