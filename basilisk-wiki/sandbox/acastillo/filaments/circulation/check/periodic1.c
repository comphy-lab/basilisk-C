#include "poisson.h"
#include "view.h"
#include "fractions.h"

double lamb_oseen (double x, double y, double a, double b){
  return 1.0/(pi*sq(a)) * exp( -(sq(x-b)+sq(y))/sq(a) );
}


int main(){
  L0 = 16;
  for (N = 32; N <= 4096; N *= 2) {
    origin (-L0/2., -L0/2.);
    init_grid (N);

    /**
    Conditions on the box boundaries are set.
    We use "third-order" [face flux interpolation](/src/embed.h).
    */

    scalar a[], b[];
    periodic(top);
    periodic(left);

    /**
    The right-hand-side is defined using the coordinates of the
    barycenter of the cut cell (xc,yc), which is calculated from the
    cell and surface fractions. */
    double bsum = 0.;
    foreach(reduction(+:bsum)) {
      b[] = lamb_oseen(x, y, 0.1, 1.0) + lamb_oseen(x, y, 0.1, -1.0);
      bsum += b[]*dv();
    }
    fprintf (stderr, "bsum (before): %g\n", bsum);

    double omega_mean = -0.5*bsum/sq(L0);

    bsum = 0.;
    foreach(reduction(+:bsum)) {
      double r = sqrt(sq(x) + sq(y));
      a[] = -sq(r)*omega_mean/2.0;
      b[] += 2*omega_mean;
      bsum += b[]*dv();
    }
    boundary ({a,b});
    fprintf (stderr, "bsum (after): %g\n", bsum);

    /**
    The Poisson equation is solved. */

    scalar res[];
    scalar * list = {res};

    timer t = timer_start();
    mgstats s = poisson (a, b, tolerance = 1e-15, minlevel = 4, res=list);
    double dt = timer_elapsed (t);
    printf ("%d %g %d %d\n", N, dt, s.i, s.nrelax);

    vector u[];
    coord f = {-1.,1.};
    foreach()
      foreach_dimension()
        u.y[] = f.y * center_gradient(a);
    boundary ((scalar *){u});

    scalar omega[];
    vorticity (u, omega);

    scalar psi[];
    foreach()
      psi[] = a[];
    boundary ({psi, omega});

    scalar res2[];
    scalar * list2 = {res2};
    t = timer_start();
    s = poisson (psi, omega, tolerance = 1e-15, minlevel = 4, res=list2);
    dt = timer_elapsed (t);
    printf ("%d %g %d %d\n", N, dt, s.i, s.nrelax);

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
    save ("a.png");
    clear();

    squares ("b", spread = -1);
    isoline ("b", n = 21);
    save ("b.png");
    clear();

    squares ("res", spread = -1);
    save ("res.png");
    clear();

    squares ("psi", spread = -1);
    isoline ("psi", n = 21);
    save ("psi.png");
    clear();

    squares ("omega", spread = -1);
    isoline ("omega", n = 21);
    save ("omega.png");
    clear();

    squares ("res2", spread = -1);
    save ("res2.png");
    clear();

    squares ("u.x", spread = -1);
    isoline ("u.x", n = 21);
    save ("u.png");
    clear();

    squares ("u.y", spread = -1);
    isoline ("u.y", n = 21);
    save ("v.png");
    clear();

    squares ("e1", spread = -1);
    isoline ("e1", n = 21);
    save ("e1.png");
    clear();

    squares ("e2", spread = -1);
    isoline ("e2", n = 21);
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

    foreach(){
      omega[] -= 2*omega_mean;
      u.x[] += y*omega_mean;
      u.y[] -= x*omega_mean;
    }

    double ekin = 0.0, circ = 0.0, enst = 0.0, area_tot=0.0;
    foreach(reduction(+:ekin), reduction(+:circ), reduction(+:enst), reduction(+:area_tot)){
      ekin += dv() * (sq(u.x[]) + sq(u.y[]));
      circ += dv() * omega[];
      enst += dv() * sq(omega[]);
      area_tot += dv();
    }

    s0 = statsf (omega);

    double circ2[2] = {0.,0.}, mu_x[2] = {0.,0.}, mu_y[2] = {0.,0.};
    foreach(reduction(+:circ2[:2]), reduction(+:mu_x[:2]), reduction(+:mu_y[:2])){
      coord p = {x,y,z};
      int nj = x > 0 ? 0 : 1 ;
      circ2[nj] += dv() * omega[] ;
      mu_x[nj] += dv() * omega[]  * p.x;
      mu_y[nj] += dv() * omega[]  * p.y;
    }

    for (int nj = 0; nj < 2; nj++){
      mu_x[nj] /= circ2[nj];
      mu_y[nj] /= circ2[nj];
    }

    double M20[2] = {0.,0.}, M02[2] = {0.,0.}, M11[2] = {0.,0.};
    foreach(reduction(+:M20[:2]), reduction(+:M02[:2]), reduction(+:M11[:2])){
      coord p = {x,y,z};
      int nj = x > 0 ? 0 : 1 ;
      M20[nj] += dv() * omega[] * sq(p.x - mu_x[nj]);
      M02[nj] += dv() * omega[] * sq(p.y - mu_y[nj]);
      M11[nj] += dv() * omega[] * (p.x - mu_x[nj]) * (p.y - mu_y[nj]);
    }

    for (int nj = 0; nj < 2; nj++){
      M20[nj] /= circ2[nj];
      M02[nj] /= circ2[nj];
      M11[nj] /= circ2[nj];
    }

    fprintf (stderr, "E: %g \t C: %g \t Z: %g \t A: %g \t max(omega): %g \n", ekin, circ, enst, area_tot, s0.max);
    int nj = 0;
    fprintf (stdout, "Vortex 1: %d %g %g %g %g %g %g %g \n", N, circ2[nj], mu_x[nj], mu_y[nj], interpolate(omega, mu_x[nj], mu_y[nj]), sqrt(M20[nj] + M02[nj]), sqrt((M20[nj] + M02[nj] + sqrt(4*sq(M11[nj]) + sq(M20[nj] - M02[nj])))), sqrt((M20[nj] + M02[nj] - sqrt(4*sq(M11[nj]) + sq(M20[nj] - M02[nj])))));
    nj = 1;
    fprintf (stdout, "Vortex 2: %d %g %g %g %g %g %g %g \n", N, circ2[nj], mu_x[nj], mu_y[nj], interpolate(omega, mu_x[nj], mu_y[nj]), sqrt(M20[nj] + M02[nj]), sqrt((M20[nj] + M02[nj] + sqrt(4*sq(M11[nj]) + sq(M20[nj] - M02[nj])))), sqrt((M20[nj] + M02[nj] - sqrt(4*sq(M11[nj]) + sq(M20[nj] - M02[nj])))));
    fprintf (stdout, "\n");
    fprintf (stderr, "\n");

    for (int i = 0; i < 2; i++){
      for (int j = 1; j <= 4; j++){
        vertex scalar phi[];
        foreach_vertex()
          phi[] = sqrt(sq(x-mu_x[i])+sq(y)) - j;

        boundary ({phi});

        scalar c[];
        face vector s[];
        fractions (phi, c, s);

        squares ("c", spread = -1);
        save ("c.png");

        double area = 0., utan=0.;

        FILE * fp = fopen("interface.asc", "w");
        foreach (reduction(+:area), reduction(+:utan)){
          if (c[] > 0 && c[] < 1.) {
            coord n = interface_normal (point, c), p;
            double alpha = plane_alpha (c[], n);
            double dl = pow(Delta, dimension - 1)*plane_area_center (n, alpha, &p);
            normalize (&n);

            double uvel = interpolate (u.x, x+p.x*Delta, y+p.y*Delta, z + p.z*Delta );
            double vvel = interpolate (u.y, x+p.x*Delta, y+p.y*Delta, z + p.z*Delta );
            double nvel = ( uvel * n.x + vvel * n.y);
            double tvel = ( uvel * n.y - vvel * n.x);
            utan += tvel * dl;
            area += dl;

            fprintf (fp, "%g %g %g %g %g %g %g %g %g \n", x+p.x*Delta, y+p.y*Delta, n.x, n.y, uvel, vvel, nvel, tvel, dl);
          }
        }
        fclose(fp);
        fprintf (stderr, "Vortex %d Perimeter %d: %.15g Circulation: %.15g \n", i, j, area, utan);
      }
    }
  }
}
