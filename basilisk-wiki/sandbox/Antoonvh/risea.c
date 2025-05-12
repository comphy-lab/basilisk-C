#include "nsf4t.h"
#include "view.h"
scalar b[], * tracers = {b};
scalar zeros[]; // Just zeros

double zeta = 0.00025;

int main(int argc, char ** argv) {
  foreach_dimension()
    periodic (left);
  int l = 9;
  if (argc >= 2) 
    l = atoi(argv[1]);
  L0 = 10;
  X0 = Y0 = Z0 = -L0/2;
  const scalar muz[] = 0.001;
  kappa = nu = muz;
  N = 1 << l;
  a.x = zeros;
  a.y = b;
  for (zeta = 0.016; zeta >= 0.00025; zeta/= 2)
    run();
}

event init (t = 0) {
  TOLERANCE = 1e-5;
  do {
    foreach_vert() 
      b[] = exp (-sq((x - exp(1)/10.)) - sq(y + 1) - sq(z - pi/10.)) ;
    boundary ({b});
  } while (adapt_list({b}, (double[]){zeta*20}, 99, 1).nc > 100);
  boundary ({b});
}

event mov (t += 0.1) {
  output_ppm (b, file = "b.mp4", n = 300);
  scalar lev[];
  foreach()
    lev[] = level;
  output_ppm (lev, file = "l.mp4", n = 300);
  
}
event adapt (i++) {
  adapt_list ({b, u}, (double[]){zeta*20., zeta, zeta},99, 1);
}

#define OFFSET(l)  (((1 << ((2*l)+2)) - 1)/3)
#define _O (-GHOSTS)
#define C_IND(i,j,l) ((i + _O) + (1 << l)*(j + _O))
#define INDEX (OFFSET(level - 1) + C_IND(point.i, point.j, level))

event dumper (t = {1,3,5}) {
  double e = 0, em = -1;
  char fname[99];
  sprintf (fname, "%d-risedata", (int)(t + .1));
  FILE * fp = fopen (fname, "rb");
  foreach(reduction(+:e) reduction(max:em)) {
    fseek (fp, INDEX*sizeof(double), SEEK_SET);
    double bi;
    fread (&bi, 1, sizeof(double), fp);
    double el = fabs(b[] - bi);
    e += dv()*el;
    if (el > em)
      em = el;
    b[] = bi;
  }
  printf ("%d %d %g %g %g %g\n", (int)(t + .1), depth(), sqrt(grid->tn), sqrt((perf.tnc)/(i )), e, em);
} 

event drawer (t = {1,3,5}) {
  view (fov = 15, ty = -1/30., samples = 4, width = 600, height = 750);
  cells();
  squares("b", map = cool_warm, min = -.8, max = .8);
  char iname[99];
  translate (z = 0.1) {
    if (fabs(t - 1) < 0.1)
      isoline ("b", .5, lw = 10, lc = {31./255., 119./255., 180./255.});
    if (fabs(t - 3) < 0.1)
      isoline ("b", .5, lw = 10, lc = {255/255., 127./255., 14/255.});
    if (fabs(t - 5) < 0.1)
      isoline ("b", .5, lw = 10, lc = {44./255., 160./255., 44/255.});
  }
  sprintf (iname, "%d-drawa.png", (int)(t + .1));
  save (iname);
} 
