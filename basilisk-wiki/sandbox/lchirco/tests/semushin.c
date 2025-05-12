#define Nl 64
#define length 1
#define Cfl (0.125)
#define ycenter (0.75)
#define R (0.15)

// #define ADAPT
#define NORM
#define SEMU2VOF

#define TRANSLA
// #define VORTEX
// #define ZALESAK

// #define OUT_X
// #define OUT_HDF5
#define OUT_GNUPLOT

#ifdef TRANSLA
#define xcenter (0.25)
#define ddt ((1./Nl)*Cfl)
#define IT (Nl/Cfl)
#define EndT 1
#define TEST "transla"
#endif

#ifdef VORTEX
#define EndT 2
#define xcenter (0.5)
#define ddt (0.0005)
#define IT (EndT/ddt)
#define TEST "vortex"
#endif

#ifdef ZALESAK
#define EndT 1
#define xcenter (0.5)
#define ddt (0.0005)
#define IT (EndT/ddt)
#define TEST "zalesak"
#endif

/**
 The volume fraction is stored in scalar field `f` which is listed as
 an *interface* for the VOF solver. We do not advect any tracer with
 the default (diffusive) advection scheme of the advection solver. */

scalar f[];
scalar * interfaces = {f}, * tracers = NULL;
vertex scalar psi[];
scalar tmp1[];

#include "advection.h"
#include "lchirco/semushin.h"


/**
 The following tests can be run. 
 TRANSLA: translation of a circle along the diagonal.
 VORTEX: the single vortex test.
 ZALESAK: advection of Zalesak's notched disk.
 Here we initialize the interface, depending on the test case chosen.
 */

#define mycircle(x,y) (sq(R) - (sq(x - xcenter) + sq(y - ycenter)))

#ifdef ZALESAK
double rectangle(double x, double y, coord center, coord size) {
  
  double P1_Plus = x - size.x/2.  - center.x;
  double P1_Minus = x + size.x/2. - center.x;
  double P1 = max (P1_Plus, -P1_Minus);
  
  double P2_Plus = y - size.y/2.  - center.y;
  double P2_Minus = y + size.y/2. - center.y;
  double P2 = max (P2_Plus, -P2_Minus);
  
  double c = max ( P1,P2 );
  return c;
}

double circle(double x, double y,  coord center, double radius) {
  double R2  =  sq(x - center.x) + sq (y - center.y);
  return ( sqrt(R2) - radius);
}

double geometry(double x, double y) {
  
  coord center_circle, center_rectangle, size_rectangle;
  center_circle.x = center_rectangle.x = 0.5;
  center_circle.y = 0.75;
  
  center_rectangle.y = 0.725;
  
  size_rectangle.x = 0.05;
  size_rectangle.y = 0.25; 
  
  double s = -circle (x, y, center_circle, 0.15);
  double r = -rectangle (x, y, center_rectangle, size_rectangle);
  
  double zalesak = difference(s,r) ;
  
  return zalesak;
}
#endif


int count_it = 0;
const char *testName = TEST;

int main() {
  size(length);
  origin (0., 0.);
  init_grid (Nl);
  run();
}

/**
 We define the levelset function $\phi$ on each vertex of the grid and
 compute the corresponding volume fraction field $f$ and face fractions
 $s\_tmp$. Then, we initialize the position of the interface points 
 (or Semushin fractions) on every edge. We assume $s = 0$ if the marker point
 coincides with the left or bottom vertices and $s = 1$ for the right
 and top ones. */


event init (i = 0) {
  foreach_dimension(){
    s.x.restriction = myrestriction_face;
    s.x.coarsen = myrestriction_face;
    s.x.prolongation = myrefine_face_x;
    s.x.refine = myrefine_face_x;
  }
  tmp1.prolongation = refine_injection;
  tmp1.refine = refine_injection;
  vertex scalar phi[];
  foreach_vertex()
    #if defined(TRANSLA) || defined(VORTEX)
    phi[] = mycircle(x,y) ;
  #endif    
  #ifdef ZALESAK
  phi[] = geometry(x,y) ;  
  #endif    
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
      if (x < xcenter)  
        s.y[] = 1. - s.y[];
    }
  }
  foreach_face(x){
    if (exist(s_tmp.x[])) { 
      s.x[] = s_tmp.x[];
      if (y < ycenter)  
        s.x[] = 1. - s.x[];
    }
  }
  boundary ((scalar *){s});
  foreach_vertex(){
    psi[] =  sq(sin(x * pi)) * sq(sin(y * pi))/pi;
  }
}

/** The timestep `dt` and the velocity field are set. */
event stability (i++, first) {
  
  dt = dtnext (ddt);
  tTime += dt;
  
  if (pid()==0){
    printf("i: %d  Time: %1.4f  test: %s \n",i, tTime, testName);
  }
  
  struct { double x, y; } f = {1., -1.};
  
  #ifdef ZALESAK  
  foreach_face(x)
    u.x[] = 2*pi*(0.5-y);
  
  foreach_face(y)
    u.y[] = 2*pi*(x-0.5);  
  #endif
  
  foreach_dimension(){
    
    #ifdef VORTEX  
    foreach_face(x)
      u.x[] = cos(pi*tTime/EndT)*f.x*(psi[0, 1] - psi[])/Delta;
    #endif  
    
    #ifdef TRANSLA  
    double reversed = 1.;
    if (i>IT/2.)
      reversed = -1.;
    foreach_face(x){
      u.x[] =  f.x * reversed ;
    }
    #endif  
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
}
#endif

/**  We can output the results in HDF5 format for visualization with Paraview.*/
#ifdef OUT_HDF5
#include "lchirco/output_hdf.h"
event calc_and_print (t = 0; i++) {
  
  scalar * list = {s.x,s.y,f};
  vector * vlist = {u};
  
  char buf[100];
  char time[100];
  char itime[100];
  count_it++;
  
  snprintf(buf, sizeof(buf), "out_%05d.xmf", count_it);
  snprintf(itime, sizeof(itime), "%05d", count_it);
  snprintf(time, sizeof(time), "%06g", tTime);
  
  FILE * fp = fopen(buf, "w");
  output_xmf_h5_foreach(list, vlist, 64, fp, itime, time); 
  fclose(fp);
}
#endif

/**  
 We can compute the shape error $E_{shape}=\max_{i}|$`dist`$(\boldsymbol{x}_i)|$.
 For the tests where the reference solution is a 
 circle centered in $(x_c,y_c)$ and with radius $R$, we have 
 `dist`$(\boldsymbol{x}_i)=\sqrt{(x_i-x_c)^2+(y_i-y_c)^2}-R$
 */
#ifdef NORM 
event calc_infty_norm (t = end) {
  
  double l_inf = 0.;
  foreach(serial){
    double xx[2], yy[2], dist = 0.;
    xx[0] = xx[1] = yy[0] = yy[1] = 0.;
    int ii = 0;
    if (exist(s.x[0,0]) || exist(s.x[1,0]) 
      || exist(s.y[0,0]) || exist(s.y[0,1])){
      if (exist(s.x[0,0])){
        yy[ii] = y - Delta/2. + s.x[0,0]*Delta;
        xx[ii] = x - Delta/2.;
        ii++;
      }
      if (exist(s.x[1,0])){
        yy[ii] = y - Delta/2. + s.x[1,0]*Delta;
        xx[ii] = x + Delta/2.;
        ii++;
      }
      if (exist(s.y[0,0])){
        xx[ii] = x - Delta/2. + s.y[0,0]*Delta;
        yy[ii] = y - Delta/2.;
        ii++;
      }
      if (exist(s.y[0,1])){
        xx[ii] = x - Delta/2. + s.y[0,1]*Delta;
        yy[ii] = y + Delta/2.;
        ii++;
      }      
      dist = fabs(sqrt(sq(xx[0] - xcenter) + sq(yy[0] - ycenter)) - R);
      if (dist > l_inf ) l_inf = dist;
      dist = fabs(sqrt(sq(xx[1] - xcenter) + sq(yy[1] - ycenter)) - R);  
      if (dist > l_inf ) l_inf = dist;
    }
  }
  printf("l-Infty norm is %f\n", l_inf);
}
#endif

/**
~~~gnuplot Translation test 
set ytics 0.125
set xtics 0.125
set xrange [0.:1]
set yrange [0.:1]
set grid
set size square
set object 1 circle at 0.25,0.75 size first 0.15 fc rgb "navy" 
plot 'transla_64.dat' w l t "64" lc 'red', NaN lc 'navy' title "Ref."
~~~
*/