
#include "grid/octree.h"
#include "nsf4t.h"
scalar * tracers = NULL;
#include "lambda2.h"
#include "view.h"
#include "filamentsG_circ.h"
#include "bwatch.h"
#include "stokes-particles-G_dsn_nolift.h"

# define R      (sqrt(sq(x) + sq(z)))
# define radi(yp)   (R <= 0.01 ? major_rt : (sqrt(sq(x - major_rt * (x / R)) + sq(z - major_rt * (z / R)) + sq(y - yp))))


double ue = 5e-3;
int maxlevel = 11;
int segs = 999;
double major_rt = 1.5;
double d1 = 82;
double muv = 1./100.;
double rhof = 0.001;
double Gamma_LO = 49.6695, a_LO = 0.3886;

double tend = 3.69;

double rp;

double ux_interp;
double uy_interp;
double uz_interp;


FILE * Vc;


int main(int argc, char ** argv) {

  foreach_dimension()
    periodic(left);
  L0 = 100;
  origin(-L0 / 2, 0, -L0 / 2);
  G.y = -981;
  const scalar muc[] = muv;
  nu = muc;
  Vc = fopen("Uvr_solverG_Repz_03_rrel_01.txt","w");
  run();
}


coord ring1 (double t) {
  coord C;
  C.x = major_rt*cos(t);
  C.y = d1 ;
  C.z = -major_rt*sin(t);
  return C;
}


event init (t = 0) {
  refine ((radi(d1) < 12) && level < maxlevel - 2);
  refine ((radi(d1) < 6)  && level < maxlevel - 1);
  refine ((radi(d1) < 3)  && level < maxlevel);
  

  new_inertial_particles(1);

  foreach_particle() {
    p().y = 68.2355;
    p().x = 0.30000000000000004;
    p().z = 0.0;
    p().u2.x = 0.0009898649;
    p().u2.y = 0.25;
  }



 vector A[], omg[], omg2[];
 TOLERANCE = 1e-5;
 get_vor_vector (omg,  ring1, 0, 2*pi, segs, Lamb_Oseen);
 foreach_dimension()
   A.x.prolongation = refine_4th;
 foreach() {
   foreach_dimension() {
     A.x[] = 0.;
   }
 }
 foreach_dimension() {
   stats o = statsf (omg.x);
   foreach() 
     omg.x[] -= o.sum/o.volume;
 }
 boundary ((scalar*){A});
 foreach_dimension()
   poisson (A.x, omg.x);
 vector uc[];
 foreach_dimension()
   uc.x.prolongation = refine_4th;
 foreach() { 
   foreach_dimension() {
     uc.x[] = (8.*(A.z[0,1] - A.z[0,-1]) + A.z[0,-2] - A.z[0,2] -
	       8.*(A.y[0,0,1] - A.y[0,0,-1]) - A.y[0,0,-2] + A.y[0,0,2])/(12*Delta);
   }
 }
 boundary ((scalar*){uc});
 vector_to_face (uc);
 project (u, p);

}


event logger (i += 5) {
  double ke2 = 0., vd2 = 0.;
  double ke4 = 0., vd4 = 0.;
  vector uv[];
  scalar dissv[], Ev[];
  face_to_vector(uv);
  foreach(reduction (+:ke2) reduction (+:vd2)) {
    foreach_dimension() {
      ke2 += dv()*sq(uv.x[]);
      vd2 += dv()*(sq(uv.x[1]     - uv.x[-1]) +
		   sq(uv.x[0,1]   - uv.x[0,-1]) +
		   sq(uv.x[0,0,1] - uv.x[0,0,-1]))/sq(2.*Delta);
    }
  }
  foreach_dimension() {
    uv.x.prolongation = refine_vert5;
    uv.x.restriction = restriction_vert;
  }
  dissv.prolongation = refine_vert5;
  dissv.restriction = restriction_vert;
  Ev.prolongation = refine_vert5;
  Ev.restriction = restriction_vert;
  foreach() {
    foreach_dimension()
      uv.x[] = FACE_TO_VERTEX_4 (u.x);
  }
  boundary ((scalar*){uv});
  foreach() {
    Ev[] = dissv[] = 0;
    foreach_dimension() {
      Ev[] += sq(u.x[]);
      dissv[] += (sq(8*(u.x[1] - u.x[-1]) + u.x[-2] - u.x[2]) +
		 sq(8*(u.x[0,1] - u.x[0,-1]) + u.x[0,-2] - u.x[0,2]) +
		 sq(8*(u.x[0,0,1] - u.x[0,0,-1]) + u.x[0,0,-2] - u.x[0,0,2]))/sq(12*Delta);
    }
  }
  boundary ({dissv, Ev});
  foreach(reduction (+:ke4) reduction (+:vd4)) {
    coord cc = {x, y , z};
    foreach_child() {
      coord ccc = {0};
      foreach_dimension()
	     ccc.x = cc.x + child.x*Delta/sqrt(3);
      ke4 += dv()*interpolate_vertex_4 (Ev,    ccc.x, ccc.y, ccc.z);
      vd4 += dv()*interpolate_vertex_4 (dissv, ccc.x, ccc.y, ccc.z);
    }
  }
  ke2 /= 2.;
  ke4 /= 2.;
  vd2 *= muv;
  vd4 *= muv;
}


event get_previous_u(i++) {
  vector uc[];
  face_to_vector (uc);
  
  foreach_particle() {
    rp = p().u2.y;

        ux_interp = interpolate(uc.x, p().x, p().y, p().z);
        uy_interp = interpolate(uc.y, p().x, p().y, p().z);
        uz_interp = interpolate(uc.z, p().x, p().y, p().z);

    }
}


event trac_p (t += 0.01) {	

  foreach_particle()
    fprintf(stderr, "%g %g %g %g %g %g %g %g\n", t, x, y, z, p().u.x, p().u.y, p().u.z, dt);
}

event trac_vr (t += 0.05) {	

  vector uc[];
  face_to_vector (uc);
  
  double deltap = L0/(1<<maxlevel);
  
  
  for (double ly=0.40*L0; ly < 0.83*L0; ly += deltap){
    fprintf(Vc, "%g %g %g %g %g\n", t, ly, interpolate(uc.x, 0, ly, 0), interpolate(uc.y, 0, ly, 0), interpolate(uc.z, 0, ly, 0));
  }
}


event adapt (i++) {
  adapt_flow (ue, 99, 1);
  double xp, yp, zp, rp;
  foreach_particle(){
    xp = x;
    yp = y;
    zp = z;
    rp = p().u2.y;
   }
  for(int kk; kk < 4; kk++){
    refine((sq(x - xp) + sq(y - yp) + sq(z - zp) < sq(2. * rp)) & & level < maxlevel);
  }
}


event stop (t = tend) {
  return 1;
}


