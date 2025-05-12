 // #include "grid/multigrid3D.h"
#include "embed.h"
#include "navier-stokes/centered.h"
// #include "navier-stokes/double-projection.h"
#include "mixed-embed.h"
#include "view.h"
// #include "navier-stokes/perfs.h"

double R1 = 0.25 [0], R2 = 0.5 [0];
double l_val = 1000.0;
// scalar lambda[];
scalar un[];
scalar slip_length[], plate_vel[];

#define powerlaw(r) (r*(0.5+(-l_val)) + sq(0.5)*((-l_val)-0.5)/r)*(sq(0.25)/(sq(0.25)*(0.5+(-l_val))+sq(0.5)*((-l_val)-0.5)))
// #define powerlaw(r) (r*(R2+(-l_val)) + sq(R2)*((-l_val)-R2)/r)*(sq(R1)/(sq(R1)*(R2+(-l_val))+sq(R2)*((-l_val)-R2)))

int main()
{
  /**
  Space and time are dimensionless. This is necessary to be able to
  use the 'mu = fm' trick. */

  size (1.1 [0]);
  DT = HUGE [0];
  
  // origin (-L0/2., -L0/2., -L0/2.);
  origin (-L0/2., -L0/2.);
  stokes = true;
  TOLERANCE = 1e-7;
  
  // for (N = 16; N <= 64; N *= 2)
  // slip = lambda;
  lambda = slip_length;
  U0 = plate_vel;
  // for (l_val = 0.; l_val <= 1.0; l_val += 0.2){
    // for (N = 16; N <= 128; N *= 2){
        lambda = slip_length;
        U0 = plate_vel;
        run();
    // }
  // }
}


event init (t = 0) {

  /**
  The gravity vector is aligned with the channel and viscosity is
  unity. */
  
  // a[] = {0.,0.,0.};
  a[] = {0.,0.};
  mu = fm;

  /**
  The channel geometry is defined using Constructive Solid Geometry. */  

  solid (cs, fs, difference (sq(R2) - sq(x) - sq(y), sq(R1) - sq(x) - sq(y)));
  fractions_cleanup (cs, fs);
  foreach(){
    slip_length[] = l_val;
    plate_vel[] = 0.;
  }
  // foreach(){
  //   if (sq(x) + sq(y) > 0.14){
  //     lambda[] = l_val;
  //   }
  // }
    

  /**
  The boundary condition is zero velocity on the embedded boundaries. */
// dirichlet (x*x + y*y > 0.14 ? nodata : -y);
  u.n[embed] = (sq(x) + sq(y) > 0.14 ? dirichlet(mx.x[]) : dirichlet(-y));
  u.t[embed] = (sq(x) + sq(y) > 0.14 ? dirichlet(mx.y[]) : dirichlet(x));
  // u.r[embed] = dirichlet(0);

  // u.n[embed] = dirichlet(y);
  // u.t[embed] = dirichlet(-x);
  /**
  We use "third-order" [face flux interpolation](/src/embed.h). */

  // for (scalar s in {u})
  //   s.third = true;
  
  /**
  We initialize the reference velocity. */
  
  foreach()
    un[] = u.x[];
}


event logfile (t += 0.1; i <= 1000) {
  double du = change (u.x, un);
  if (i > 0 && du < 1e-6){
    return 1; /* stop */
  }
}


event profile (t = end)
{
  // scalar utheta[], e[];
  // foreach() {
  //   double theta = atan2(y, x), r = sqrt(x*x + y*y);
  //   if (cs[] > 0.) {
  //     utheta[] = - sin(theta)*u.x[] + cos(theta)*u.y[];
  //     e[] = utheta[] - powerlaw (r);
  //   }
  //   else
  //     e[] = p[] = utheta[] = nodata;
  // }

  // norm n = normf (e);

  // char name[50];
  // sprintf (name, "sqoppfile-Navier-lambda%g.dat", l_val);
  // static FILE * fp = fopen (name,"w");

  // fprintf (fp, "%d %.3g %.3g %.3g %d %d %d %d %d\n",
	//    N, n.avg, n.rms, n.max, i, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax);
  // dump();
  
  // draw_vof ("cs", "fs", filled = -1, fc = {1,1,1});
  // squares ("utheta", spread = -1);
  // save ("utheta.png");

  // draw_vof ("cs", "fs", filled = -1, fc = {1,1,1});
  // squares ("p", spread = -1);
  // save ("p.png");

  // draw_vof ("cs", "fs", filled = -1, fc = {1,1,1});
  // squares ("e", spread = -1);
  // save ("e.png");
  char name[360];
  sprintf (name, "slip%g", l_val);
FILE * fp = fopen(name, "a");
char name3[360];
sprintf (name3, "flow_slip%g_N%d", l_val, N);
FILE * fp1 = fopen(name3, "w");

  scalar utheta[], e[];
  foreach() {
    double theta = atan2(y, x), r = sqrt(x*x + y*y);
    if (cs[] > 0.) {
      utheta[] = - sin(theta)*u.x[] + cos(theta)*u.y[];
      e[] = utheta[] - powerlaw (r);
      fprintf (fp1, "%g %g %g\n", r, utheta[], powerlaw(r));
    }
    else
      e[] = p[] = utheta[] = nodata;
  }
  fclose(fp1);
  norm n = normf (e);

// norm n = normf (e);
// scalar utan[];
// foreach(){
//   utan[] = 0.707106781*u.x[] + 0.707106781*u.y[];
// }

fprintf (fp, "%d %.3g %.3g %.3g %d %d %d %d %d\n",
   N, n.avg, n.rms, n.max, i, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax);
fclose(fp);
  char name2[360];
    sprintf (name2, "ut_slip%g.png", l_val);
  if (N == 32){
    vector csu[];
    foreach()
      foreach_dimension()
        csu.x[] = cs[]*u.x[];
    view(width = 1024, height = 1024);
    box();
    draw_vof ("cs", "fs", lw = 3);
    vectors ("csu", scale = 0.025);
    save (name2);
  }
}