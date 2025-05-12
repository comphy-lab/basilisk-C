#include "embed.h"
#include "poisson.h"
#include "view.h"

double lamb_oseen (double r, double a){
  return 1.0/(pi*sq(a)) * exp( -sq(r/a) );
}

int main(){
  L0 = 16;
  for (N = 32; N <= 2048; N *= 2) {
    origin (-L0/2., -L0/2.);
    init_grid (N);

    vertex scalar phi[];
    foreach_vertex() {
      double r = sqrt (sq(x) + sq(y));
      phi[] = L0/2. - r;
    }
    boundary ({phi});
    fractions (phi, cs, fs);
#if TREE
    cs.refine = cs.prolongation = fraction_refine;
#endif
    boundary ({cs,fs});
    restriction ({cs,fs});

    cm = cs;
    fm = fs;

    /**
    Conditions on the box boundaries are set.
    We use "third-order" [face flux interpolation](/src/embed.h).
    */

    scalar a[], b[];
    a[embed] = dirichlet_homogeneous();
    a.third = true;

    /**
    The right-hand-side is defined using the coordinates of the
    barycenter of the cut cell (xc,yc), which is calculated from the
    cell and surface fractions. */
    double bsum = 0.;
    foreach(reduction(+:bsum)) {
      a[] = cs[] > 0. ? 0. : nodata;
      double xc = x, yc = y;
      if (cs[] > 0. && cs[] < 1.) {
        coord n = facet_normal (point, cs, fs), p;
        double alpha = plane_alpha (cs[], n);
        line_center (n, alpha, cs[], &p);
        xc += p.x*Delta, yc += p.y*Delta;
      }
      double r = sqrt(sq(xc) + sq(yc));
      b[] = lamb_oseen(r, 1.0) * cs[];
      bsum += b[]*dv();
    }
    boundary ({a,b});
    fprintf (stderr, "bsum: %g\n", bsum);

    /**
    The Poisson equation is solved. */
    struct Poisson p;
    p.alpha = fs;
    p.lambda = zeroc;
    p.embed_flux = embed_flux;
    scalar res[];

    // FIXME: need to set minlevel to 4
    timer t = timer_start();
    mgstats s = poisson (a, b, alpha = fs, tolerance = 1e-10, minlevel = 4);
    double dt = timer_elapsed (t);
    printf ("%d %g %d %d\n", N, dt, s.i, s.nrelax);

    double maxp = residual ({a}, {b}, {res}, &p), maxf = 0.;
    foreach()
      if (cs[] == 1. && fabs(res[]) > maxf)
        maxf = fabs(res[]);
    fprintf (stderr, "maxres %d %.3g %.3g\n", N, maxf, maxp);

    vector u[];
    u.n[embed] = dirichlet_homogeneous();
    u.t[embed] = neumann_homogeneous();

    coord f = {-1.,1.};
    foreach()
      foreach_dimension()
        u.y[] = f.y * center_gradient(a);
    boundary ((scalar *){u});

    scalar omega[];
    omega[embed] = neumann_homogeneous();
    vorticity (u, omega);

    scalar psi[];
    psi[embed] = dirichlet_homogeneous();
    psi.third = true;

    foreach()
      psi[] = a[];
    boundary ({psi, omega});

    t = timer_start();
    s = poisson (psi, omega, alpha = fs, tolerance = 1e-10, minlevel = 4);
    dt = timer_elapsed (t);
    printf ("%d %g %d %d\n", N, dt, s.i, s.nrelax);

    maxp = residual ({psi}, {omega}, {res}, &p), maxf = 0.;
    foreach()
      if (cs[] == 1. && fabs(res[]) > maxf)
        maxf = fabs(res[]);
    fprintf (stderr, "maxres %d %.3g %.3g\n", N, maxf, maxp);

    scalar e1[], e2[];
    foreach(){
      e1[] = a[]-psi[];
      e2[] = b[]-omega[];
    }

    /**
    The solution is displayed using bview.
    */

    squares ("a", spread = -1);
    isoline ("a", n = 21);
    draw_vof ("cs", "fs");
    save ("a.png");
    clear();

    squares ("b", spread = -1);
    isoline ("b", n = 21);
    draw_vof ("cs", "fs");
    save ("b.png");
    clear();

    squares ("psi", spread = -1);
    isoline ("psi", n = 21);
    draw_vof ("cs", "fs");
    save ("psi.png");
    clear();

    squares ("omega", spread = -1);
    isoline ("omega", n = 21);
    draw_vof ("cs", "fs");
    save ("omega.png");
    clear();

    squares ("u.x", spread = -1);
    isoline ("u.x", n = 21);
    draw_vof ("cs", "fs");
    save ("u.png");
    clear();

    squares ("u.y", spread = -1);
    isoline ("u.y", n = 21);
    draw_vof ("cs", "fs");
    save ("v.png");
    clear();

    squares ("e1", spread = -1);
    isoline ("e1", n = 21);
    draw_vof ("cs", "fs");
    save ("e1.png");
    clear();

    squares ("e2", spread = -1);
    isoline ("e2", n = 21);
    draw_vof ("cs", "fs");
    save ("e2.png");
    clear();

    dump ("dump");

    stats s0 ;
    s0 = statsf (a); fprintf (stderr, "a: %g %g %g \n", s0.min, s0.max, s0.sum);
    s0 = statsf (b); fprintf (stderr, "b: %g %g %g \n", s0.min, s0.max, s0.sum);
    s0 = statsf (psi); fprintf (stderr, "psi: %g %g %g \n", s0.min, s0.max, s0.sum);
    s0 = statsf (omega); fprintf (stderr, "omega: %g %g %g \n", s0.min, s0.max, s0.sum);
    s0 = statsf (e1); fprintf (stderr, "e(psi): %g %g %g \n", s0.min, s0.max, s0.sum);
    s0 = statsf (e2); fprintf (stderr, "e(omega): %g %g %g \n", s0.min, s0.max, s0.sum);

    double ekin = 0.0, circ = 0.0, enst = 0.0, area_tot=0.0, utan=0.0;
    foreach(reduction(+:ekin), reduction(+:circ), reduction(+:enst), reduction(+:area_tot), reduction(+:utan)){
      double vort = omega[];
      double uvel = u.x[];
      double vvel = u.y[];
      double area = dv();
      if (cs[] < 1. && cs[] > 0){
        coord b, n;
        area *= embed_geometry (point, &b, &n);
        uvel = embed_interpolate (point, u.x, b);
        vvel = embed_interpolate (point, u.y, b);
        vort = 0.; // embed_vorticity (point, u, b, n); Embed vorticity works for dirichlet conditions
        utan += (-uvel * n.y + vvel * n.x) * pow(Delta, dimension - 1) * embed_geometry (point, &b, &n);
      }
      ekin += area * (sq(uvel) + sq(vvel));
      circ += area * vort;
      enst += area * sq(vort);
      area_tot += area;
    }

    s0 = statsf (omega);

    double circ2 = 0., mu_x = 0., mu_y = 0.;
    foreach(reduction(+:circ2), reduction(+:mu_x), reduction(+:mu_y)){
      coord p = {x,y,z};
      double vort = omega[];
      double area = dv();
      if (cs[] < 1. && cs[] > 0){
        coord b, n;
        area *= embed_geometry (point, &b, &n);
        vort = 0.;
        p = (coord) {p.x + b.x*Delta, p.y + b.y*Delta};
      }
      circ2 += area * vort ;
      mu_x += area * vort * p.x;
      mu_y += area * vort * p.y;
    }

    fprintf (stderr, "E: %g \t C: %g %g %g \t Z: %g \t A: %g \t max(omega): %g \t mu: %g %g \n", ekin, circ, circ2, utan, enst, area_tot, s0.max, mu_x, mu_y);
    fprintf (stderr, "\n");
  }
}
