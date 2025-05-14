#include "grid/octree.h"
#include "poisson.h"
#include "view.h"
#include "fractions.h"

double lamb_oseen (double r, double a){
  return 1.0/(pi*sq(a)) * exp( -sq(r/a) );
}

trace
void vorticity3d(vector u, vector omega){
  vector du[], dv[], dw[]; // Would be nice to use a tensor here.
  foreach(){
    du.x[] = (u.x[1] - u.x[-1])/(2.*Delta);
    dv.x[] = (u.y[1] - u.y[-1])/(2.*Delta);
    dw.x[] = (u.z[1] - u.z[-1])/(2.*Delta);

    du.y[] = (u.x[0,1] - u.x[0,-1])/(2.*Delta);
    dv.y[] = (u.y[0,1] - u.y[0,-1])/(2.*Delta);
    dw.y[] = (u.z[0,1] - u.z[0,-1])/(2.*Delta);

    du.z[] = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);
    dv.z[] = (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta);
    dw.z[] = (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta);
  }

  foreach(){
    omega.x[] = dw.y[] - dv.z[];
    omega.y[] = du.z[] - dw.x[];
    omega.z[] = dv.x[] - du.y[];
  }
  boundary ((scalar*){omega});
}


int main(){
  L0 = 16;
  for (int maxlevel = 5; maxlevel <= 9; maxlevel += 1) {
    int N = 1 << maxlevel;
    init_grid (1 << maxlevel);
    origin (-L0/2., -L0/2.,-L0/2.);
    /**
    Conditions on the box boundaries are set.
    We use "third-order" [face flux interpolation](/src/embed.h).
    */

    vector c[], b[];
    periodic(top);
    periodic(left);
    periodic(front);

    /**
    The right-hand-side is defined using the coordinates of the
    barycenter of the cut cell (xc,yc), which is calculated from the
    cell and surface fractions. */
    double bsum = 0.;
    foreach(reduction(+:bsum)) {
      double r = sqrt(sq(x) + sq(y));
      foreach_dimension(){
        c.x[] = 0;
        b.x[] = 0;
      }
      b.z[] += lamb_oseen(r, 1.0);
      bsum += b.z[]*dv();
    }
    fprintf (stderr, "bsum (before): %g\n", bsum);

    double omega_mean = -0.5*bsum/pow(L0,3);

    bsum = 0.;
    foreach(reduction(+:bsum)) {
      double r = sqrt(sq(x) + sq(y));
      c.z[] = -sq(r)*omega_mean/2.0;
      b.z[] += 2*omega_mean;
      bsum += b.z[]*dv();
    }
    boundary ((scalar *){c,b});
    fprintf (stderr, "bsum (after): %g\n", bsum);

    /**
    The Poisson equation is solved. */

    foreach_dimension(){
      timer t = timer_start();
      mgstats s = poisson (c.x, b.x, tolerance = 1e-15, minlevel = 4);
      double dt = timer_elapsed (t);
      printf ("%d %g %d %d\n", N, dt, s.i, s.nrelax);
    }

    vector u[];
    foreach(){
      u.x[] = -((c.z[0,1,0] - c.z[0,-1,0]) - (c.y[0,0,1] - c.y[0,0,-1]))/(2.*Delta);
      u.y[] = -((c.x[0,0,1] - c.x[0,0,-1]) - (c.z[1,0,0] - c.z[-1,0,0]))/(2.*Delta);
      u.z[] = -((c.y[1,0,0] - c.y[-1,0,0]) - (c.x[0,1,0] - c.x[0,-1,0]))/(2.*Delta);
    }

    vector omega[];
    vorticity3d (u, omega);

    vector psi[];
    foreach()
      foreach_dimension()
        psi.x[] = c.x[];
    boundary ((scalar *){psi, omega});

    foreach_dimension(){
      timer t = timer_start();
      mgstats s = poisson (psi.x, omega.x, tolerance = 1e-15, minlevel = 4);
      double dt = timer_elapsed (t);
      printf ("%d %g %d %d\n", N, dt, s.i, s.nrelax);
    }

    scalar e1[], e2[];
    foreach(){
      e1[] = c.z[]-psi.z[];
      e2[] = b.z[]-omega.z[];
    }

    /**
    The solution is displayed using bview.
    */

    view(camera="iso");
    squares ("c.z", spread = -1);
    save ("a.png");
    clear();

    squares ("b.z", spread = -1);
    save ("b.png");
    clear();

    squares ("psi.z", spread = -1);
    save ("psi.png");
    clear();

    squares ("omega.z", spread = -1);
    save ("omega.png");
    clear();

    squares ("u.x", spread = -1);
    save ("u.png");
    clear();

    squares ("u.y", spread = -1);
    save ("v.png");
    clear();

    squares ("u.z", spread = -1);
    save ("w.png");
    clear();

    squares ("e1", spread = -1);
    save ("e1.png");
    clear();

    squares ("e2", spread = -1);
    save ("e2.png");
    clear();

    cells();
    save ("cells.png");
    clear();

    dump ("dump");

    stats s0 ;
    s0 = statsf (c.z); fprintf (stderr, "a: %g %g %g \n", s0.min, s0.max, s0.sum);
    s0 = statsf (b.z); fprintf (stderr, "b: %g %g %g \n", s0.min, s0.max, s0.sum);
    s0 = statsf (psi.z); fprintf (stderr, "psi: %g %g %g \n", s0.min, s0.max, s0.sum);
    s0 = statsf (omega.z); fprintf (stderr, "omega: %g %g %g \n", s0.min, s0.max, s0.sum);
    s0 = statsf (e1); fprintf (stderr, "e(psi): %g %g %g \n", s0.min, s0.max, s0.sum);
    s0 = statsf (e2); fprintf (stderr, "e(omega): %g %g %g \n", s0.min, s0.max, s0.sum);

    foreach(){
      omega.z[] -= 2*omega_mean;
      u.x[] += y*omega_mean;
      u.y[] -= x*omega_mean;
    }

    double ekin = 0.0, circ = 0.0, enst = 0.0, area_tot=0.0;
    foreach(reduction(+:ekin), reduction(+:circ), reduction(+:enst), reduction(+:area_tot)){
      ekin += dv() * (sq(u.x[]) + sq(u.y[]) + sq(u.z[]));
      circ += dv() * omega.z[];
      enst += dv() * sq(omega.z[]);
      area_tot += dv();
    }

    s0 = statsf (omega.z);

    double circ2[1] = {0.}, mu_x[1] = {0.}, mu_y[1] = {0.};
    foreach(reduction(+:circ2[:1]), reduction(+:mu_x[:1]), reduction(+:mu_y[:1])){
      coord p = {x,y,z};
      int nj = 0;
      circ2[nj] += dv() * omega.z[] ;
      mu_x[nj] += dv() * omega.z[]  * p.x;
      mu_y[nj] += dv() * omega.z[]  * p.y;
    }

    for (int nj = 0; nj < 1; nj++){
      mu_x[nj] /= circ2[nj];
      mu_y[nj] /= circ2[nj];
    }

    double M20[1] = {0.}, M02[1] = {0.}, M11[1] = {0.};
    foreach(reduction(+:M20[:1]), reduction(+:M02[:1]), reduction(+:M11[:1])){
      coord p = {x,y,z};
      int nj = 0;
      M20[nj] += dv() * omega.z[] * sq(p.x - mu_x[nj]);
      M02[nj] += dv() * omega.z[] * sq(p.y - mu_y[nj]);
      M11[nj] += dv() * omega.z[] * (p.x - mu_x[nj]) * (p.y - mu_y[nj]);
    }

    for (int nj = 0; nj < 1; nj++){
      M20[nj] /= circ2[nj];
      M02[nj] /= circ2[nj];
      M11[nj] /= circ2[nj];
    }

    fprintf (stderr, "E: %g \t C: %g \t Z: %g \t A: %g \t max(omega): %g \n", ekin, circ, enst, area_tot, s0.max);
    int nj = 0;
    fprintf (stdout, "Vortex 1: %d %g %g %g %g %g %g %g \n", N, circ2[nj], mu_x[nj], mu_y[nj], interpolate(omega.z, mu_x[nj], mu_y[nj], 0), sqrt(M20[nj] + M02[nj]), sqrt((M20[nj] + M02[nj] + sqrt(4*sq(M11[nj]) + sq(M20[nj] - M02[nj])))), sqrt((M20[nj] + M02[nj] - sqrt(4*sq(M11[nj]) + sq(M20[nj] - M02[nj])))));
    fprintf (stdout, "\n");
    fprintf (stderr, "\n");

    for (int i = 0; i < 1; i++){
      for (int j = 1; j <= 6; j++){
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
            double wvel = interpolate (u.z, x+p.x*Delta, y+p.y*Delta, z + p.z*Delta );
            double nvel = ( uvel * n.x + vvel * n.y);
            double tvel = ( uvel * n.y - vvel * n.x);
            utan += tvel * dl;
            area += dl;

            fprintf (fp, "%g %g %g %g %g %g %g %g %g %g %g %g \n", x+p.x*Delta, y+p.y*Delta, z+p.z*Delta, n.x, n.y, n.z, uvel, vvel, wvel, nvel, tvel, dl);
          }
        }
        fclose(fp);
        fprintf (stderr, "Vortex %d Perimeter %d: %.15g Circulation: %.15g \n", i, j, area, utan);
      }
    }
  }
}
