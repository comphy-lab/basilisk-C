/**
# 3D Single Vortex with the EBIT method */

#define ADAPT 1

scalar f[];
scalar * interfaces = {f}, * tracers = NULL;

#include "grid/octree.h"
#include "advection-ebit.h"
event stability (i++) {} // this ensures that time step is computed before marker advection
#include "ebit-3d.h"

const char *OUTNAME = "vortex";
int level = 6, maxlvl, minlvl;
double Cfl = 0.125 [0], R = 0.15, xcenter, ycenter, zcenter, \
  ddt = 1. [0,1], EndT = 1. [0,1];
int IT;

int main(int argc, char * argv[]) {
  if (argc > 1)
    level = atoi (argv[1]);
  maxlvl = level;
  minlvl = max(level - 4, 3);

  init_grid (1 << level);

  ddt = Cfl/N *1. [0,1];

  xcenter = 0.35;
  ycenter = 0.35;
  zcenter = 0.35;

  EndT = 3.;

  IT = (int) EndT/ddt;

  run();
}


event init (i = 0) {
  vertex scalar phi[];
  foreach_vertex() {
    phi[] = sq(R) - (sq(x - xcenter) + sq(y - ycenter) + sq(z - zcenter));
  }

  init_markers (phi);

  semu2vof();
  volume0 = volume;
}

/** The timestep `dt` and the velocity field are set. */
event stability (i++, i < IT, first) {
  double u0 = 1. [1,-1];
  dt = dtnext (ddt);

  double cost = cos(pi*tTime/EndT)*u0;
  foreach() {
    double xx = x/L0, yy = y/L0, zz = z/L0;
    u.x[] = 2.*cost*sq(sin(xx*pi))*sin(2.*yy*pi)*sin(2.*zz*pi);
    u.y[] = -cost*sin(2.*xx*pi)*sq(sin(yy*pi))*sin(2.*zz*pi);
    u.z[] = -cost*sin(2.*xx*pi)*sin(2.*yy*pi)*sq(sin(zz*pi));
  }

  boundary ((scalar *){u});
  tTime += dt;
}

#if ADAPT
event adapt (i++) {
  adapt_wavelet ({mask_intf}, (double[]) {0.02}, maxlevel = maxlvl, minlevel = minlvl);
}
#endif

// event interface_out (i++, last) {
//   if (2*(i + 1) % (int) max(IT, 1) == 0){
//     int ii = 2*(i + 1)/max(IT, 1);
//     char name[80], name_fem[80], name_con[80];

//     sprintf (name, "%s%s_%d_%d_%d.dat", OUTPUTPATH, "vortex_3d_ebit", N, (int) EndT, ii);
//     sprintf (name_fem, "%s%s_%d_%d_%d_fem.dat", OUTPUTPATH, "vortex_3d_ebit", N, (int) EndT, ii);
//     sprintf (name_con, "%s%s_%d_%d_%d_con.dat", OUTPUTPATH, "vortex_3d_ebit", N, (int) EndT, ii);

//     output_facets_semushin (name);
//     output_facets_semushin (name_fem, tri = true, file_con = name_con);
//   }
// }

/**  
 We can compute the shape error $E_{shape}=\max_{i}|$`dist`$(\boldsymbol{x}_i)|$.
 For the tests where the reference solution is a 
 circle centered in $(x_c,y_c)$ and with radius $R$, we have 
 `dist`$(\boldsymbol{x}_i)=\sqrt{(x_i-x_c)^2+(y_i-y_c)^2}-R$
 */
event calc_infty_norm (t = end) {
  // This event is correct only if there is no markers on the computational boundary
  #if ADAPT
  refine (fabs(x - xcenter) <= 0.25*L0 && fabs(y - ycenter) <= 0.25*L0
  && fabs(z - zcenter) <= 0.25*L0 && level < maxlvl);
  #endif

  double l_inf = 0.;
  foreach(serial, reduction(max:l_inf)) {
    coord xo = {x - 0.5*Delta, y - 0.5*Delta, z - 0.5*Delta}, \
      xcen = {xcenter, ycenter, zcenter};
    double dist = 0.;

    foreach_dimension() {
      if (with_marker.x[] > 0.) {
        dist = sq(s.x[]*Delta + xo.x - xcen.x) + sq(xo.y - xcen.y) + sq(xo.z - xcen.z);
        dist = fabs(sqrt(dist) - R);
        if (dist > l_inf )
          l_inf = dist;
      }
    }
  }

  if (pid() == 0) {
    printf ("l-Infty norm is %e\n", l_inf);
    printf ("volume0: %e  volume: %e; Error %e\n", volume0, volume, fabs(volume0 - volume)/volume0);
  }

  // Initial shape
  vertex scalar phi[];
  scalar f0[];
  foreach()
    f0[] = f[];

  foreach_vertex() {
    phi[] = sq(R) - (sq(x - xcenter) + sq(y - ycenter) + sq(z - zcenter));
  }

  init_markers (phi);
  semu2vof();

  foreach()
    swap(double, f[], f0[]);

  double v0 = 0., v1 = 0., delta_v = 0.;
  foreach(reduction(+:v0) reduction(+:v1) reduction(+:delta_v)) {
    v0 += f0[]*dv();
    v1 += f[]*dv();
    delta_v += fabs(f[] - f0[])*dv();
  }
  if (pid() == 0) {
    printf ("vol0: %e vol: %e\n", v0, v1);
    printf ("E_m Error: %.8e, E_g Error: %.8e\n", fabs(v0 - v1)/v0, delta_v/v0);
  }
}
