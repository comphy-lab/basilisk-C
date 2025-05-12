/**
# Vortex Ejection from a mode 3 instability

According to [Kizner et
al. (2013)](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/instabilities-of-the-flow-around-a-cylinder-and-emission-of-vortex-dipoles/2327C2CFA76059D27461A2BC6A09F146),
a flow around a no-slip cylinder with radius $R$ maybe unstable and
eject three dipolar vortex pairs. We study the flow using embedded
boundaries and the Navier-Stokes solver. Furthermore, we will `view`
our results.
*/
#include "embed.h"
#include "navier-stokes/centered.h"
#include "view.h"
/**
The maximum resolution is set to $\Delta_{min}=R/100$. This allows to
run on the sandbox server.
*/
int maxlevel = 12;
double Re = 30000., ue = 0.003;
face vector muc[];

int main(){
  init_grid(64);
  L0 = 40;
  mu = muc;
  X0 = Y0 = -L0/2;
  run();
}

event properties (i++){
  foreach_face()
    muc.x[] = fm.x[]/Re;
  boundary ((scalar*){muc});
}
/**
   The cylinder is defined and the flow field is initialized c.f.  Kizner
   et al. with an $m=3$ perturbation.
*/
#define RAD (pow(sq(x) + sq(y), 0.5))
#define THETA(M) (M*asin(x/RAD))
#define RADP(P, M) ((1 + P*sin(THETA(M)))/(pow(1 + 0.5*sq(P), 0.5)))
double a1 = 1.5, b1 = 2.25;
double P = 0.005, m = 3;

event init (t = 0){
  double gamma = (sq(a1) - 1.)/(sq(b1) - sq(a1));
  refine (RAD < b1  && level < (maxlevel - 1));
  refine (RAD < 1.05 && RAD > 0.95  && level < maxlevel);
  vertex scalar phi[];
  foreach_vertex()
    phi[] = RAD - 1;
  fractions (phi, cs, fs);
  foreach(){
    double r = RAD;
    double r1 = RADP(P,m)*r;
    double vr;
    if (r1 > 0.9 && r1 < a1)
      vr = r1 - 1./r1;
    else if (r1 >= a1 && r1 <= b1)
      vr = -gamma*r1 + ((1 + gamma)*sq(a1) - 1.)/r1;
    else // (0.9 > r || r > b)
      vr = 0;
    u.x[] = cm[]*0.5*vr*r*-y/(sq(x) + sq(y));
    u.y[] = cm[]*0.5*vr*r* x/(sq(x) + sq(y));
  }
  /**
     The boundary conditions on the embedded boundary are set:
  */
  u.t[embed] = dirichlet (0.);
  u.n[embed] = dirichlet (0.);
  boundary (all);
}
/**
   The grid is adaptedively refined and coarsened to properly represent
the boundary and the evolution of the flow field. We set
$\zeta_{u,v}\approx U_{max}/150$.
*/
event adapt (i++)
  adapt_wavelet ({cs, u.x, u.y}, (double[]){0.01, ue, ue}, maxlevel);

/**
## Output

Movies are generated displaying the vorticity field and the
adaptive-grid-cell structure.

![](kizner2/kizner.mp4)

![](kizner2/kizner_cells.mp4)
*/

event movie (t += 0.5){
  scalar omega[];
  vorticity (u, omega);
  boundary ({omega});
  view (fov = 8, width = 600, height = 600);
  clear();
  draw_vof ("cs", filled = -1, fc = {1., 1., 1.});
  squares ("omega", min = -0.75, max = 0.75,
	   map = cool_warm, linear = true);
  draw_vof ("cs", "fs", lw = 2);
  save ("kizner.mp4");
  clear();
  view (fov = 5);
  cells();
  save ("kizner_cells.mp4");
  clear();
  view (fov = 2, width = 1080, height = 1080);
  draw_vof ("cs", filled = -1, fc = {0.6, 0.6, 0.6});
  squares ("omega", min = -0.75, max = 0.75,
	   map = cool_warm);
  cells();
  save ("kizner_cells_zoom.mp4");
}
/**
Also there is this movie:

![A zoom](kizner2/kizner_cells_zoom.mp4)(width="50%")

We log some solver data:
 */
event logger (t += 0.1){
  fprintf (stderr, "%d %g %d %d %d %g %g\n",
	   i, t, mgp.i, mgu.i, mgpf.i, DT, dt);
  int cn = 0;
  foreach(reduction(+:cn))
    cn++;
  if (pid() == 0){
    static FILE * fp = fopen ("data", "w");
    fprintf (fp, "%g\t%d\t%d\t%g\t%g\n", t, i, cn, dt, DT);
    fflush (fp);
  }
}
/**
Including the number of grid cells:

~~~gnuplot These numbers may be compared against the millions of cells that Kizner et al. (2013) employed.
set yr [0 : 25000]
set xlabel 'time'
set ylabel 'Cells'
set key off
plot 'data' u 1:3 w l lw 2
~~~
*/

event the_last_event (t = 100);

/**
## Flow direction partitioning

We can study the radially-averaged absolute velocity in the azimutal
and radial direction as a function of time and the distance from the
centre. Therefore, we render movies of the radial profiles using
`profile6.h`, `gnuplot` and `ffmpeg`.

![A more quantative perspective on the instability. Mind the
 logaritmic vertical axis](kizner2/mov.mp4)
*/
#include "profile6.h"
FILE * gnuplotPipe;
event init (t = 0){
  gnuplotPipe = popen ("gnuplot", "w");
  fprintf(gnuplotPipe,
	  "set term pngcairo\n"
	  "set xr [1 : 6]\n"
	  "set yr [0.001 : 1]\n"
	  "set logscale y\n"
	  "set key top right\n"
	  "set grid\n"
	  "set xlabel 'r/R'\n"
	  "set ylabel 'E'\n"
	  "set size ratio 0.5\n");
}

int frame = 0;
event profs(t += 0.5){
  scalar ur[], ua[] ,* list;
  list = {ur, ua};
  foreach(){
    double r = sqrt(sq(x) + sq(y));
    ur[] = fabs(u.x[]*x/r + u.y[]*y/r);
    ua[] = fabs((u.y[]*x/r) - (u.x[]*y/r));
  }
  profile (list, sqrt(sq(x) + sq(y)), "prof.dat"); 
  fprintf (gnuplotPipe, "set output 'plot%d.png'\n", frame);
  fprintf (gnuplotPipe, "set title 't = %d' font ',25'\n",
	   (int)(t + 0.9));
  fprintf (gnuplotPipe, "plot 'prof.dat' u 1:2 w l lw 5 lt rgb'#11BB11' t 'Radial',"
	   "'prof.dat' u 1:3 w l lw 5 lt rgb '#BB11BB' t 'Azimutal'\n");
  fflush (gnuplotPipe);
  frame++;
}

event stop (t = end){
  system("ffmpeg -r 25 -f image2 -i plot%d.png -c:v libx264 -vf format=yuv420p -y mov.mp4 && "
	 "rm -f plot*");
  return 1;
}
/**
## Reference

Kizner, Z., Makarov, V., Kamp, L., & Van Heijst,
G. (2013). *Instabilities of the flow around a cylinder and emission
of vortex dipoles*. Journal of Fluid Mechanics, 730,
419-441. doi:10.1017/jfm.2013.345
*/
 
