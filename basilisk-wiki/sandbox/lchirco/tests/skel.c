/**
 
# Thin sheets detection in a four vortexes field with no surface tension 
 
 */

// #include "grid/multigrid.h"
#include "advection.h"
#include "vof.h"
#include "output_hdf.h"
#include "signature.h"

#define Nl 128
#define length 2.
#define ADAPT

/**
   The volume fraction is stored in scalar field `f` which is listed as
   an *interface* for the VOF solver. We do not advect any tracer with
   the default (diffusive) advection scheme of the advection solver. */

scalar f[], cf[], phii[], s[];
scalar * interfaces = {f}, * tracers = NULL;

/**
   We set the grid dimension and position. */

double min_x = 0.;
double min_y = 0.;
double step = length/Nl;

int main()
{
  size(length);
  origin (min_x, min_y);
  DT = .1;
  
  /**
   The scalar field `cf` is a "vof concentration" associated with phase
   `f`. */
  
  f.tracers = {cf};
  init_grid (Nl);
  run();
}

/**
 The initial interface is a circle centered in
 (-1, -1). We use the levelset
 function `circle()` to define this interface. */

#define circle(x,y) (sq(0.3) - (sq(x - 1) + sq(y - 1)))

/**
  We define the levelset function $\phi$ on each vertex of the grid and
  compute the corresponding volume fraction field. */

event init (i = 0)
{
  fraction (f, circle(x,y));
  foreach()
    cf[] = f[];
  boundary ({cf});
  event("adapt");
}

  /**
    This event defines the velocity field.
    The velocity field is defined through a streamfunction $\psi$, defined
    on the vertices of the grid. This streamfunction gives four vortices*/

event velocity (i++) {
  
  vertex scalar psi[];
  double tmp = 0.;
  double angle = 0*pi/6.;
  foreach_vertex(){
    tmp = x;
    x = 1 + cos(angle)*(x - 1) - sin(angle)*(y - 1);
    y = 1 + sin(angle)*(tmp - 1) + cos(angle)*(y - 1);
    psi[] = - sin(x * pi) * sin(y * pi)/5*pi;
  }
  
  /**
     We can then differentiate the streamfunction to get the velocity
     components. This guarantees that the velocity field is exactly
     non-divergent. */
  
  trash ({u});
  struct { double x, y; } f = {-1., 1.};
  foreach_face()
    u.x[] = f.x*(psi[0, 1] - psi[])/Delta;
  
  boundary ((scalar *){u});
}

/** 
   This event computes thin sheets. 
   At the moment it works only on the coarsest grid. 
*/ 

event calc_and_print (t = 0; t += 0.06){
  
  foreach()
    phii[] = 2*f[] - 1;
  
/** 
   We determine at which level (`l_sign`) we want to compute
   the quadratic moments and the signature.
   Then, the indicator function `phii` is restricted at level `l_sign`.
*/ 
  
  double R = 0.02;  
  double target_delta = R*2./5.;
  
  int l_sign = 7 /*find_moments_level(f, length, target_delta)*/;
  
  for (int ilev = depth() - 1; ilev >= l_sign; ilev--)  
    foreach_level(ilev){
      if(is_refined(cell))
      restriction_average(point, phii);
    }
  
  compute_signature_neigh_level (f, phii, s, l_sign);
  
  /** 
   The signature `s` is available only at the level `l_sign`. 
   We need to prolong it onto the finest grid. */
  
  printf("\n level used for moments %d and depth is %d \n", l_sign, depth()); 
  
  for (int ilev = l_sign; ilev < depth(); ilev++)  
    foreach_level(ilev){
      s.prolongation = phii.prolongation = refine_injection;
      if(is_refined(cell)){
        s.prolongation (point, s);
        phii.prolongation (point, phii);
      }
    }
  
//    change_topology (f, s, Nl, min_x, min_y, step);
  
  /**
     Finally we output the results in HDF5 format.  */
  
  scalar * list = {phii, s, f};
  vector * vlist = {u};
  
  char buf[15];
  char itime[15];
  
  snprintf(buf, sizeof(buf), "out_%05d.xmf", i);

  snprintf(itime, sizeof(itime), "%05d", i);
  FILE * fp = fopen(buf, "w");
  output_xmf_h5_foreach(list, vlist, 64, fp, itime); 
  fclose(fp);
  
  char name[500];
  sprintf(name, "%s-%.5d-CC.vtk", "out", i);
  output_paraview_IF(name, t, 0.0, f);
  
//   char dumpname[80];
//   sprintf (dumpname, "snapshot-%g", t);
//   dump (dumpname);
  
}

event end(t = 0.6){
  return 0;
}

#ifdef ADAPT
event adapt (t = 0; t += 0.05) {
  int maxLev = 9, minLev = 7;
  adapt_wavelet ({f}, (double[]){1e-4}, maxlevel = maxLev,
                 minlevel = minLev);
}
#endif




