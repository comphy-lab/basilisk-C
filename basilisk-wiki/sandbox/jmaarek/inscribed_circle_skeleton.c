#define BGHOSTS 2
#define QUADRATIC 1

#include "grid/quadtree.h"


#include "fractions.h"
#include "utils.h"
#include "curvature.h"
#include "vof2front.h"
#include "alex_functions.h"
#include "LS_reinit.h"
#include "display.h"

static double xc = 0.0, yc = 0.0;

double interface_area2 (scalar cs, face vector fs)
{
  double area = 0.;
  foreach (reduction(+:area))
    if (interfacial(point,cs)) {
      coord n = facet_normal (point, cs, fs), p;
      double alpha = plane_alpha (cs[], n);
      area += pow(Delta, dimension - 1)*plane_area_center (n, alpha, &p);
    }
  return area;
}

scalar f[], f1[], f2[];
vector h[];
vertex scalar phi[], phi2[];
vertex scalar magnitude[];

scalar cs[], cs2[], cs3[];
face vector fs[], fs2[], fs3[];

face vector grad_phi[];
face vector grad_phi_transverse[];
face vector grad_phi_error[];


int main(){

  f.refine = f.prolongation = fraction_refine;
  f.height = h;

  f1.refine = f1.prolongation = fraction_refine;
  f2.refine = f2.prolongation = fraction_refine;

  L0 = 5.;
  X0 = Y0 = Z0 = -L0/2;
  for (int ii = 4; ii<=8; ii++){
  init_grid(1 << ii);

  fraction(f1, pow(sq(x - xc) + sq(y - yc), 0.5) - 2.0);
  fraction(f2, pow(sq(x - xc) + sq(y - yc), 0.5) - 1.0);

  foreach_vertex()
    phi2[] = pow(sq(x - xc) + sq(y - yc), 0.5) - 1.5;
  boundary({phi2});

  foreach()
    f[] = f2[] - f1[];
  boundary({f});

//seed function : second order
  vof2dist(f, phi);
  fractions(phi, cs, fs);
  double area = interface_area2(cs,fs);

//propogate seed to full domain with LS_reinit : second order
  LS_reinit(phi, dt = 0.3 * L0/(1 << grid->maxdepth), it_max = 10/3*(1 << grid->maxdepth));

  fractions(phi, cs2, fs2);
  double area2 = interface_area2(cs2,fs2);

  double error = 0.0;
  double nn = 0.0;

  foreach_vertex(){
    error += fabs(fabs(phi[] - 0.5) - fabs(phi2[]));
    nn++;}


double error_dx = 0.0;
double thickness_error = 0.0;
double nn_2 = 0.0;
double grad_error = 0.0;

fractions(phi2, cs3, fs3);


foreach_face()
  if (fabs(x) < 2.3 && fabs(y) < 2.3){ //need to do an approximation of +x,y,z boundaries, develop tagging function for thin structures
    if (((phi[] > phi[-1]) && (phi[] >= phi[1]) && (phi[-1] < phi[1])) || ((phi[1] >= phi[0]) && (phi[1] > phi[2]) && (phi[2] < phi[0]))){ //determines if the ridge crosses the face
        grad_phi.x[] = fabs(-phi[-1] + phi[] + phi[1] - phi[2])/2/Delta;}
    else{
        grad_phi.x[] = fabs(phi[1] - phi[0])/Delta;}
    //grad_phi_error.x[] = fabs(fabs(x/max(sqrt(sq(x)+sq(y)),1e-16)) - fabs(grad_phi.x[]));
  }
boundary_flux ({grad_phi});


foreach_face(){
  grad_phi_transverse.x[] = (grad_phi.y[-1,0] + grad_phi.y[-1,1] + grad_phi.y[0,0] + grad_phi.y[0,1])/4;
}
boundary_flux ({grad_phi_transverse});

foreach_face(){
double summ = max(sqrt(sq(grad_phi.x[]) + sq(grad_phi_transverse.x[])),1e-16);
grad_phi.x[]/= summ;
}
boundary_flux ({grad_phi});


foreach_face(x){
  grad_phi_error.x[] = fabs(fabs(x/max(sqrt(sq(x)+sq(y)),1e-16)) - fabs(grad_phi.x[]));
  grad_phi_error.y[] = fabs(fabs((y-Delta/2)/max(sqrt(sq(x+Delta/2)+sq(y-Delta/2)),1e-16)) - fabs(grad_phi.y[]));
}

boundary_flux ({grad_phi_error});

foreach_face()
  if (fs3.x[] > 0.0 && fs3.x[] < 1.0){
    grad_error += grad_phi_error.x[];
    double dx = ((phi[0]-phi[1])/grad_phi.x[] - 1)/-2;
    error_dx += fabs(sqrt(sq(x + (-0.5+dx)*Delta)+ sq(y)) - 1.5);
    thickness_error += fabs(0.5 - (phi[] + fabs(grad_phi.x[])*dx*Delta));
    nn_2++;
  }

fprintf(fout, "%d %12e %12e %12e %12e %12e %12e\n", ii, error/nn, fabs(area - 6*M_PI), fabs(area2 - 6*M_PI), error_dx/nn_2, thickness_error/nn_2, grad_error/nn_2);
}

  /*fields_stats();

  fputc ('\n', stderr);
  display_url (stderr);
  fputc ('\n', stderr);

  display ("box();");
  while (1) {
    if (display_poll (-1))
      display_update (INT_MAX);
  }
  display_destroy();*/
}
