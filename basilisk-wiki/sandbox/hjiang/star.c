#define SURFACE 1
#define F_ERR 1.e-6
#define S_ERR 1.e-6

#include "navier-stokes/centered.h"
#include "surfactant/vof-surfactant.h"
#include "tension.h"
#include "surfactant/adapt_wavelet_leave_interface.h"
#include "view.h"

/**
### Boundary Conditions

Open (outflow) boundary conditions are imposed on every boundary, since the
star is initialized in the center of the domain. 
*/

u.n[top] = neumann(0.);
u.t[top] = neumann(0.);
p[top] = dirichlet(0.);
pf[top] = dirichlet(0.);

u.n[right] = neumann(0.);
u.t[right] = neumann(0.);
p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);

u.n[left] = neumann(0.);
u.t[left] = neumann(0.);
p[left] = dirichlet(0.);
pf[left] = dirichlet(0.);

u.n[bottom] = neumann(0.);
u.t[bottom] = neumann(0.);
p[bottom] = dirichlet(0.);
pf[bottom] = dirichlet(0.);

double density = 0.5;

scalar c[], *interfaces = {c};

double r0(double theta) { return 0.75 * (1. - 0.2 * sin(7. * theta)); }
double drdtheta(double theta) { return -0.75 * 0.2 * 7. * cos(7. * theta); }
double phi0(double x, double y)
{
  double theta = atan2(y, x);
  double r = sqrt(x * x + y * y);
  return r0(theta) - r;
}

int maxlevel = 9;
int minlevel = 4;

/**
### Model Data

Additional data useful for the simulation are the maximum level of refinement
and the initial droplet radius.
*/

int main()
{
  L0 = 4;
  origin(-L0 / 2., -L0 / 2.);
  c.sigma = 10;
  for (maxlevel = 7; maxlevel <= 9; maxlevel++) {
    init_grid(1 << maxlevel);
    run();
  }
}

#define circle(x, y, R) (sq(R) - sq(x) - sq(y))

event init(t = 0)
{
  // scalar d[];
  // foreach() d[] = phi0(x, y);
  // vertex scalar phi[];
  // foreach_vertex() phi[] = (d[] + d[-1] + d[0, -1] + d[-1, -1]) / 4.;
  // fractions(phi, c);

  fraction(c, phi0(x, y));
  foreach () {
    surfactant[] = (c[] > F_ERR && c[] < 1. - F_ERR) ? density : 0.0;
  }
  boundary({surfactant});
}

event end_timestep(i++)
{
  foreach_face()
  {
    coord o = {x, y, z};
    double r = sqrt(sq(x) + sq(y));
    uf.x[] = o.x / (r);
  }

  foreach () {
    coord o = {x, y, z};
    double r = sqrt(sq(x) + sq(y));
    foreach_dimension() u.x[] = o.x / (r);
  }
  event("stability");
}

#if TREE
event adapt(i++) { adapt_wavelet_leave_interface({surfactant}, {c}, (double[]){1.e-3}, maxlevel, minlevel = 2, 1); }
#endif

/**
## Post-Processing

The following lines of code are for post-processing purposes.
*/

double exact(double time, double theta)
{
  double r = r0(theta);
  double rp = drdtheta(theta);
  return density * sqrt(sq(r) + sq(rp)) / sqrt(sq(r + time) + sq(rp));
}

event export_error(t = {0.025, 0.125, 0.25, 0.375, 0.5})
{
  char errname[128];
  sprintf(errname, "ErrorData-%d.dat", maxlevel);

  char partname[128];
  sprintf(partname, "OutputData-%d-%.3f-pid%d.dat", maxlevel, t, pid());
  FILE *fpart = fopen(partname, "w");
  fprintf(fpart, "# theta surfactant_num surfactant_exact rel_error\n");

  double L1 = 0.;
  double Linf = 0.;
  double ltot = 0.;

  foreach (reduction(+ : L1) reduction(max : Linf) reduction(+ : ltot)) {
    if (st.Astar[] > 1.e-2 * L0 / (1 << grid->maxdepth)) {
      coord n = interface_normal(point, c), p;
      double a = plane_alpha(c[], n);
      double dl = plane_area_center_surfactant(n, a, &p) * Delta;

      double xp = x + p.x * Delta;
      double yp = y + p.y * Delta;
      double theta = atan2(yp, xp);
      if (theta < 0.) theta += 2. * M_PI;

      if (theta >= pi / 4. && theta <= 3. * pi / 4.) {
        double s_exact = exact(t, theta);
        double err = 0.;
        if (fabs(s_exact) > 1e-30) err = (surfactant[] - s_exact) / s_exact;

        fprintf(fpart, "%.12g %.12g %.12g %.12g\n", theta, surfactant[], s_exact, err);

        double aerr = fabs(err);
        L1 += aerr * dl;
        if (aerr > Linf) Linf = aerr;
        ltot += dl;
      }
    }
  }

  fclose(fpart);

  if (ltot > 0.) L1 /= ltot;

#if _MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if (pid() == 0) {
    char outname[128];
    sprintf(outname, "OutputData-%d-%.3f.dat", maxlevel, t);

    char cmd[512];
    sprintf(cmd,
      "LC_ALL=C cat OutputData-%d-%.3f-pid*.dat | "
      "grep -v '^#' | sort -g -k1,1 > %s",
      maxlevel, t, outname);
    system(cmd);

    FILE *fout = fopen(outname, "a");
    fprintf(fout, "\n# L1 %.12g\n", L1);
    fprintf(fout, "# Linf %.12g\n", Linf);
    fclose(fout);

    FILE *ferr = fopen(errname, t == 0.025 ? "w" : "a");
    if (t == 0.025) fprintf(ferr, "# t L1 Linf\n");
    fprintf(ferr, "%.12g %.12g %.12g\n", t, L1, Linf);
    fclose(ferr);
  }
}

#if 0
event snapshots(i++)
{
  //if (t > 0.01) {
    char name[80];
    // sprintf(name, "snapshots-%.3f", t);
    sprintf(name, "snapshots-%d", i);
    dump(name);
  //}
}
#endif

event movie(t += 0.01)
{
  if (maxlevel == 9) {
    clear();
    squares("(c[] > 1e-10 && c[] < 1.-1e-10) ? surfactant : nodata",
        min = 0., max = density);
    //draw_vof("c", lw = 1.5, lc = {255, 255, 255});
    box();
    save("surfactant.mp4");
  }
}

event end(t = 0.5);

/**
## Results
*/

/**
![Movie](star/surfactant.mp4)

~~~gnuplot Average vs Analytic

set datafile commentschars "#"

r0(th)  = 0.75*(1. - 0.2*sin(7.*th))
dr0(th) = -0.75*0.2*7.*cos(7.*th)
G0 = 0.5
Gs(th,t) = G0*sqrt(r0(th)**2 + dr0(th)**2)/sqrt((r0(th)+t)**2 + dr0(th)**2)

set samples 400

#do for [lev=7:9] {
do for [lev=9:9] {

    do for [tt in "0.025 0.125 0.250 0.375 0.500"] {
        outfile = sprintf("Analytic-%d-%s.dat", lev, tt)
        set table outfile
        plot [pi/4:3.*pi/4] Gs(x, real(tt))
        unset table
    }

    #set terminal qt size 800,600 enhanced font "Arial Bold,28"

    #set xlabel "{/Symbol q}" font "Arial Bold,28"
    #set ylabel "Surfactant Concentration" font "Arial Bold,28"

    #set xtics font "Arial Bold,28"
    #set ytics font "Arial Bold,28"
    #set key font "Arial Bold,20"
    unset key

    #set border lw 2
    #set grid lw 1

    set xrange [pi/4:3.*pi/4]
    set xtics ("{/Symbol p}/4" pi/4, "{/Symbol p}/2" pi/2, "3{/Symbol p}/4" 3.*pi/4)
    set yrange [0.25:0.5]

    #set size ratio 0.75
    #set tmargin 1
    #set bmargin 3
    #set lmargin 7
    #set rmargin 1

    plot \
    sprintf("OutputData-%d-0.025.dat", lev) u 1:2 w p pt 6 ps 1.5 lc rgb "black" title "Simulation, t = 0.025", \
    sprintf("Analytic-%d-0.025.dat", lev) u 1:2 w l lw 1.9 lc rgb "black" title "Analytic, t = 0.025", \
    sprintf("OutputData-%d-0.125.dat", lev) u 1:2 w p pt 6 ps 1.5 lc rgb "blue" title "Simulation, t = 0.125", \
    sprintf("Analytic-%d-0.125.dat", lev) u 1:2 w l lw 1.9 lc rgb "blue" title "Analytic, t = 0.125", \
    sprintf("OutputData-%d-0.250.dat", lev) u 1:2 w p pt 6 ps 1.5 lc rgb "red" title "Simulation, t = 0.25", \
    sprintf("Analytic-%d-0.250.dat", lev) u 1:2 w l lw 1.9 lc rgb "red" title "Analytic, t = 0.25", \
    sprintf("OutputData-%d-0.375.dat", lev) u 1:2 w p pt 6 ps 1.5 lc rgb "brown" title "Simulation, t = 0.375", \
    sprintf("Analytic-%d-0.375.dat", lev) u 1:2 w l lw 1.9 lc rgb "brown" title "Analytic, t = 0.375", \
    sprintf("OutputData-%d-0.500.dat", lev) u 1:2 w p pt 6 ps 1.5 lc rgb "purple" title "Simulation, t = 0.5", \
    sprintf("Analytic-%d-0.500.dat", lev) u 1:2 w l lw 1.9 lc rgb "purple" title "Analytic, t = 0.5"

    pause -1
}

~~~

~~~gnuplot Relative Errors

reset
set datafile commentschars "#"

set print "ErrorsAtT0.500.dat"
do for [lev=7:9] {
    stats sprintf("ErrorData-%d.dat", lev) using (abs($1 - 0.5) < 1e-12 ? $2 : 1/0) nooutput
    L1 = STATS_mean

    stats sprintf("ErrorData-%d.dat", lev) using (abs($1 - 0.5) < 1e-12 ? $3 : 1/0) nooutput
    Linf = STATS_mean

    print sprintf("%d %d %.12g %.12g", lev, 2**lev, L1, Linf)
}
unset print

f1(x) = a1*x**(-p1)
f8(x) = a8*x**(-p8)

a1 = 1.
p1 = 1.
a8 = 1.
p8 = 1.

fit f1(x) "ErrorsAtT0.500.dat" u 2:3 via a1,p1
fit f8(x) "ErrorsAtT0.500.dat" u 2:4 via a8,p8

#set terminal qt size 900,700 enhanced font "Arial Bold,28"

#set xlabel "N" font "Arial Bold,28"
#set ylabel "Error at t = 0.5" font "Arial Bold,28"

#set xtics font "Arial Bold,24"
#set ytics font "Arial Bold,24"
#set key font "Arial Bold,20"

#set border lw 2
#set grid lw 1
set logscale x 2
set logscale y
set format y "10^{%L}"

set xr [2**6:2**10]

set style line 1 lc rgb "black" lw 2 pt 6 ps 2
set style line 2 lc rgb "red"   lw 2 pt 8 ps 2
set style line 3 lc rgb "black" lw 2 dt 2
set style line 4 lc rgb "red"   lw 2 dt 2

plot \
"ErrorsAtT0.500.dat" u 2:3 w p ls 1 title "L^{1}", \
"ErrorsAtT0.500.dat" u 2:4 w p ls 2 title "L^{inf}", \
f1(x) w l ls 3 title sprintf("fit: N^{-%.2f}", p1), \
f8(x) w l ls 4 title sprintf("fit: N^{-%.2f}", p8)

~~~

*/
