The interfacial gas flow through a liquid medium is solved by this code.

#include "navier-stokes/centered.h"
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "view.h"
 
#define d 1.0
#define H1 (0.2*d)
#define H2 (0.2*d)
#define L (2048.0*d/8.0)
#define epsilon 0.117

#define Re 23.46
#define Ca 0.156
#define M 492.04
#define R 863.17
#define Q (R/M)

#define Tu 100.0

double dx;
double tEnd = Tu*(d/epsilon)*d/Q;

double aNoise;
double deformationPast = 1.0;

int maxlevel = 12;
double uemax = 0.0001;

u.n[top]  = dirichlet(0);
u.t[top]  = dirichlet(0);
u.n[bottom]  = dirichlet(0);
u.t[bottom]  = dirichlet(0);

f[left] = dirichlet (y < H1 ? 1.0 : (y > d-H2 ? 1.0 : 0.0));
u.n[left]  = dirichlet(uInlet(y));
p[right]  = dirichlet(0);
pf[right] = dirichlet(0);

u.n[right] = neumann(0);
u.t[left]  = dirichlet(0);
//u.t[right] = neumann(0);
p[left]    = neumann(0);
pf[left] = neumann(0);

The function which defines the inlet velocity is written below.
  
  double uInlet (double y){

  double muL = 1.0/Re;
  double muG = 1.0/Re/M;

  double Z1 = sq(H1)/2.0/muL-sq(d-H2)/2.0/muL+sq(d)/2.0/muL-sq(H1)/2.0/muG+sq(d-H2)/2.0/muG;
  double Z2 = -H1/muL+(d-H2)/muL-d/muL+H1/muG-(d-H2)/muG;
  double dpdx = Q/(pow(d-H2,3.0)/6.0/muG+Z1*sq(d-H2)/2.0/muG/Z2+sq(H1)*(d-H2)*(1.0/muL-1.0/muG)/2.0+Z1*H1*(d-H2)*(1.0/muL-1.0/muG)/Z2-pow(H1,3.0)/6.0/muG-Z1*sq(H1)/2.0/Z2/muG+pow(H1,3.0)*(1.0/muG-1.0/muL)/2.0-Z1*sq(H1)*(1.0/muL-1.0/muG)/Z2);

  double c2 = 0.0;
  double c5 = dpdx*Z1/Z2;
  double c1 = c5;
  double c3 = c5;
  double c4 = -dpdx*sq(d)/2.0/muL-c5*d/muL;
  double c6 = dpdx*sq(H1)/2.0/muL+c5*H1/muL-dpdx*sq(H1)/2.0/muG-c5*H1/muG;

  double U_inlet = y < H1 ? dpdx*sq(y)/2.0/muL+c1*y/muL+c2+aNoise : (y > d-H2 ? dpdx*sq(y)/2.0/muL+c3*y/muL+c4+aNoise : dpdx*sq(y)/2.0/muG+c5*y/muG+c6);

  return U_inlet;

}

The function which defines the initial fluid flow is written below.
  
  double uInitial (double y){

  double muL = 1.0/Re;
  double muG = 1.0/Re/M;

  double Z1 = sq(H1)/2.0/muL-sq(d-H2)/2.0/muL+sq(d)/2.0/muL-sq(H1)/2.0/muG+sq(d-H2)/2.0/muG;
  double Z2 = -H1/muL+(d-H2)/muL-d/muL+H1/muG-(d-H2)/muG;
  double dpdx = Q/(pow(d-H2,3.0)/6.0/muG+Z1*sq(d-H2)/2.0/muG/Z2+sq(H1)*(d-H2)*(1.0/muL-1.0/muG)/2.0+Z1*H1*(d-H2)*(1.0/muL-1.0/muG)/Z2-pow(H1,3.0)/6.0/muG-Z1*sq(H1)/2.0/Z2/muG+pow(H1,3.0)*(1.0/muG-1.0/muL)/2.0-Z1*sq(H1)*(1.0/muL-1.0/muG)/Z2);

  double c2 = 0.0;
  double c5 = dpdx*Z1/Z2;
  double c1 = c5;
  double c3 = c5;
  double c4 = -dpdx*sq(d)/2.0/muL-c5*d/muL;
  double c6 = dpdx*sq(H1)/2.0/muL+c5*H1/muL-dpdx*sq(H1)/2.0/muG-c5*H1/muG;

  double U_initial = y < H1 ? dpdx*sq(y)/2.0/muL+c1*y/muL+c2 : (y > d-H2 ? dpdx*sq(y)/2.0/muL+c3*y/muL+c4 : dpdx*sq(y)/2.0/muG+c5*y/muG+c6);

  return U_initial;

}

The below function is used to define fraction at init event.

  double Geometry (double x,double y){

  double Line_1 = y-H1;
  double Line_2 = y-d+H2;
  double Line_3 = x+L/2.;
  double Line_4 = x-L/2.;

  double Geom_y = max (Line_2, -Line_1);
  double Geom_x = max (Line_4, -Line_3);
  double Geom_f = max (Geom_y, Geom_x);

  return Geom_f;
}

Three below functions are used to employ in the refine function in order to improve the mesh generation.

  double refineGeom (double x,double y){

  double Line_1 = y;
  double Line_2 = y-d;
  double Line_3 = x+L/2.;
  double Line_4 = x-L/2.;

  double Geom_y = max (Line_2, -Line_1);
  double Geom_x = max (Line_4, -Line_3);
  double Geom_f = max (Geom_y, Geom_x);

  return Geom_f;
}

  double refineGeomD (double x,double y){

  double refineParameter = 0.05;

  double Line_1 = y-H1+refineParameter*d;
  double Line_2 = y-H1-refineParameter*d;
  double Line_3 = x+L/2.;
  double Line_4 = x-L/2.;

  double Geom_y = max (Line_2, -Line_1);
  double Geom_x = max (Line_4, -Line_3);
  double Geom_f = max (Geom_y, Geom_x);

  return Geom_f;
}

  double refineGeomU (double x,double y){

  double refineParameter = 0.05;

  double Line_1 = y-d+H2+refineParameter*d;
  double Line_2 = y-d+H2-refineParameter*d;
  double Line_3 = x+L/2.;
  double Line_4 = x-L/2.;

  double Geom_y = max (Line_2, -Line_1);
  double Geom_x = max (Line_4, -Line_3);
  double Geom_f = max (Geom_y, Geom_x);

  return Geom_f;
}

int main (int argc, char * argv[])
{
  if (argc > 1)
    maxlevel = atoi (argv[1]);
  if (argc > 2)
    uemax = atof (argv[2]);

  origin (-L/2., -L/2.+d/2.0);
  size (L);

  rho1 = 1.0;
  rho2 = 1.0/R;
  mu1  = 1.0/Re;
  mu2  = 1.0/Re/M;
  f.sigma = R/Re/sq(M)/sq(epsilon)/d/Ca;
  init_grid (2048);

  DT = 1e-2;
  TOLERANCE = 1e-6;
  NITERMAX  = 1000;

  run();
}

event init (t = 0) {

  double WN = 0.8;
  mask (y > d ? top : (y < 0.0 ? bottom : none));

  refine (refineGeom(x,y) && level < maxlevel+2);
  refine (refineGeomD(x,y) && level < maxlevel+2);
  refine (refineGeomU(x,y) && level < maxlevel+2);

//  refine ((x >= -length/2.0 || x <= length/2.0) && (refineGeomD(x,y) || refineGeomU(x,y)) && level < maxlevel);

  fraction (f, Geometry(x,y));

  foreach() {
     u.x[] = uInitial(y);
  }

  boundary ((scalar *){u});
  aNoise = WN*noise();

}

The below event is used to create white-noise in order to add to the inlet velocity. In summary, V=V0+V1 where V is the inlet velocity, V0 is the velocity which is obtained by the analytical solution when the surface tension is neglected, v1 is the white-noise which is the multiplication of amplitude and a random number between -1 and 1; v1 is calculated by the below event in each iteration.

event noiseCreate (i++) {
  
  double WN = 0.8;
  aNoise = WN*noise();
//  fprintf (ferr, "%g \n", aNoise);
}

event adapt (i++) {
  adapt_wavelet ({f,u}, (double[]){0.0001,uemax,uemax}, maxlevel, maxlevel);
}

The problem is about the below event. The bug is occurred here. When the below event is added to the code, the results will be changed. If it is called at i++, it will not change the results; however if the event is added in the following format, it will affect the resaults. Moreover, if the frequency of calling is changed, for example as i += 10 or i += 100 and something like that, the results will change again. Furthermore, it should be noted that the adding of this event in the following format modify the results and the obtained results become physically reasonable which is surprising for me because as far as I know the properties event is called several times during each iteration in navier-stokes/centered.h.

event properties (i += 2; t <= tEnd) {
  boundary ({alpha, mu, rho});
}

event movie1 (i += 10; t <=tEnd)
{
  view (fov = 19.2, quat = {0.0,0.0,0.0,0.0},
	  tx = 0.0, ty = -0.3, bg = {1,1,1},
	  width = 1600, height = 1600, sy = 150.0);
  clear();
  draw_vof ("f", lc = {1,0,0}, lw = 2);
  squares ("f", linear = true, spread = 10);
  box ();
  save ("f.mp4");
}

event snapshot (t = tEnd/50.0; t += tEnd/50.0; t <= tEnd) {
  char name[80];
  sprintf (name, "snapshot-%g", t);
  dump (file = name);
}

The below event is used to calculate the deformations.

event deformation (t += 0.1; t <= tEnd) {
  double deformation = 0.0;
  double deformation1 = 0.0;
  double deformation2 = 0.0;
  double rad1;
  double rad2;
  double length1 = 0.0;
  double length2 = 0.0;
// caculation of deformation
  foreach (reduction(+:deformation1) reduction(+:deformation2) reduction(+:length1) reduction(+:length2))
    if (f[] > 0 && f[] < 1) {
      coord p;
      coord n = mycs (point, f);
      double alpha = plane_alpha (f[], n);
      plane_area_center (n, alpha, &p);
      rad1 = y < d/2.0 ? fabs(y + Delta*p.y - H1) : 0.0;
      rad2 = y > d/2.0 ? fabs(y + Delta*p.y - d + H2) : 0.0;
      length1 = y < d/2.0 ? length1+Delta : length1;
      length2 = y > d/2.0 ? length2+Delta : length2;
      deformation1 = deformation1 + rad1*Delta;
      deformation2 = deformation2 + rad2*Delta;
    }
  deformation1 = deformation1/length1;
  deformation2 = deformation2/length2;
  deformation = deformation1+deformation2;
  double dd = fabs(deformation-deformationPast)/fabs(deformationPast);
//  if (i > 0 && dd < 1e-3)
//  return 1; /* stop */
  deformationPast = deformation;
  fprintf (stderr, "%ld %g %g %g %g %g %g %g %g %g %g \n", i, t, t*Q*epsilon/sq(d), dt, deformation1, deformation2, deformation, deformation1/d, deformation2/d, deformation/d, dd);
}

The below event is used to calculate an integral which calculates the magnitude of symmetricity of deformations in the domain.

event caculationC (t += 0.1; t <= tEnd) {
  double depP21;
  double depP22;
  double dep1;
  double dep2;
  double de1 = 0.0;
  double de2 = 0.0;
  double de;
  double nu1;
  double nu2;
  double nu;
  double length1;
  double length2;
  double Tnu = 0.0;

// caculation of C
  foreach (reduction(+:de1) reduction(+:de2))
    if (f[] > 0 && f[] < 1) {
      coord p;
      coord n = mycs (point, f);
      double alpha = plane_alpha (f[], n);
      plane_area_center (n, alpha, &p);
      depP21 = y < d/2.0 ? sq(y + Delta*p.y - H1) : 0.0;
      depP22 = y > d/2.0 ? sq(d - H2 - y - Delta*p.y) : 0.0;
      de1 = de1 + depP21*Delta;
      de2 = de2 + depP22*Delta;
    }
  de = sqrt(de1*de2);

  int jjjMax = 2048;
  double rjjjMax = 2048.0;
  int jjj;
  double step_i = L/rjjjMax;
  for (jjj = 0; jjj <= jjjMax - 1; jjj++) {
    nu1 = 0.0;
    nu2 = 0.0;
    length1 = 0.0;
    length2 = 0.0;
    foreach (reduction(+:nu1) reduction(+:nu2) reduction(+:length1) reduction(+:length2))
      if (x >= jjj*step_i - L/2.0 && x <= (jjj+1.0)*step_i - L/2.0) {
        if (f[] > 0 && f[] < 1) {
          coord p;
          coord n = mycs (point, f);
          double alpha = plane_alpha (f[], n);
          plane_area_center (n, alpha, &p);
          dep1 = y < d/2.0 ? y + Delta*p.y - H1 : 0.0;
          dep2 = y > d/2.0 ? d - H2 - y - Delta*p.y : 0.0;
          length1 = y < d/2.0 ? length1+Delta : length1;
          length2 = y > d/2.0 ? length2+Delta : length2;
          nu1 = nu1 + dep1*Delta;
          nu2 = nu2 + dep2*Delta;
        }
      }
    nu1 = nu1/length1;
    nu2 = nu2/length2;
    nu = nu1*nu2*step_i;
    Tnu = Tnu + nu;
  }
  double C = Tnu/de;
  static FILE * fp = fopen ("deTnuC", "w");
  fprintf (fp, "%g %g %g %g %g\n", t*Q*epsilon/sq(d), t, de, Tnu, C);
  fflush (fp);
}

The below event gives the interface deformations at the end of simulation.

event interface (t = end) {
  output_facets (f);
}

/*
gnuplot commands:
set xlabel 'x'
set ylabel 'y'
set key bottom Right
plot  'out.ppm' u 1:2 w l lw 2 t 'interface'
*/