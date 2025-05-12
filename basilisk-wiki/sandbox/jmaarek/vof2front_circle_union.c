#include "vof2front.h"
#include "display.h"

/*
Visual demonstration of propsosed VOF2Front method. Attempts to reconstruct an interface created from the unsion of two circles. The method is computed on a quadtree where the grid size on the interface is allowed to vary
*/


scalar f[];
vertex scalar phi_analytical[], phi_reconstruct[];
vector h[];

double init_circle(double x, double y, double a, double r){
  if((a-x)/sqrt(sq(a-x)+sq(y))>= a/r &&
     (a+x)/sqrt(sq(a+x)+sq(y))>= a/r){
    return -min(sqrt(sq(x) + sq(y + sqrt(sq(r) - sq(a))) ),
               sqrt(sq(x) + sq(y - sqrt(sq(r) - sq(a))) )) ;
  }
  else
  return min(sqrt(sq(x+a)+sq(y))-r,
   sqrt(sq(x-a)+sq(y))-r) ;
}


int main() {

  phi_analytical.restriction = restriction_vertex;
  phi_analytical.refine = phi_analytical.prolongation =  prolongation_vertex;

  phi_reconstruct.restriction = restriction_vertex;
  phi_reconstruct.refine = phi_reconstruct.prolongation =   prolongation_vertex;

  f.refine = f.prolongation = fraction_refine;
  f.height = h;


  size(3);
  origin(-L0/2,-L0/2);

  init_grid (1 << 9);

  foreach_vertex(){
    phi_analytical[] = max(init_circle(x,y,-0.4,0.5), init_circle(x,y,0.4,0.5));
  }
  boundary({phi_analytical});
  restriction({phi_analytical});

  fractions(phi_analytical, f);

  adapt_wavelet ({f}, (double[]) {1e-3}, 9);
  adapt_wavelet ({f}, (double[]) {1e-3}, 9);
  adapt_wavelet ({f}, (double[]) {1e-3}, 9);
  adapt_wavelet ({f}, (double[]) {1e-3}, 9);
  adapt_wavelet ({f}, (double[]) {1e-3}, 9);

  vof2front(f, phi_reconstruct);

  fields_stats();

  fputc ('\n', stderr);
  display_url (stderr);
  fputc ('\n', stderr);

  display ("box();");
  while (1) {
    if (display_poll (-1))
      display_update (INT_MAX);
  }
  display_destroy();

  }
