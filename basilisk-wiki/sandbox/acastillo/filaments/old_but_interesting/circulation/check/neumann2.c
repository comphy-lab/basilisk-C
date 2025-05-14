#include "poisson.h"
#include "view.h"

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
    a[left]   = dirichlet(0);
    a[right]  = dirichlet(0);
    a[top]    = dirichlet(0);
    a[bottom] = dirichlet(0);

    /**
    The right-hand-side is defined using the coordinates of the
    barycenter of the cut cell (xc,yc), which is calculated from the
    cell and surface fractions. */
    double bsum = 0.;
    foreach(reduction(+:bsum)) {
      a[] = 0.;
      double xc = x, yc = y;
      b[] = lamb_oseen(xc, yc, 0.1, 1.0) - lamb_oseen(xc, yc, 0.1, -1.0);
      bsum += b[]*dv();
    }
    boundary ({a,b});
    fprintf (stderr, "bsum: %g\n", bsum);

    /**
    The Poisson equation is solved. */
    timer t = timer_start();
    mgstats s = poisson (a, b, tolerance = 1e-10, minlevel = 4);
    double dt = timer_elapsed (t);
    printf ("%d %g %d %d\n", N, dt, s.i, s.nrelax);

    vector u[];
    u.n[top]    = dirichlet_homogeneous(); u.t[top]    = neumann_homogeneous();
    u.n[bottom] = dirichlet_homogeneous(); u.t[bottom] = neumann_homogeneous();
    u.n[left]   = dirichlet_homogeneous(); u.t[left]   = neumann_homogeneous();
    u.n[right]  = dirichlet_homogeneous(); u.t[right]  = neumann_homogeneous();

    coord f = {-1.,1.};
    foreach()
      foreach_dimension()
        u.y[] = f.y * center_gradient(a);
    boundary ((scalar *){u});

    scalar omega[];
    omega[top]    = neumann_homogeneous();
    omega[bottom] = neumann_homogeneous();
    omega[left]   = neumann_homogeneous();
    omega[right]  = neumann_homogeneous();
    vorticity (u, omega);

    scalar psi[];
    psi[left]   = dirichlet(0);
    psi[right]  = dirichlet(0);
    psi[top]    = dirichlet(0);
    psi[bottom] = dirichlet(0);

    foreach()
      psi[] = a[];
    boundary ({psi, omega});

    t = timer_start();
    s = poisson (psi, omega, tolerance = 1e-10, minlevel = 4);
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

    squares ("psi", spread = -1);
    isoline ("psi", n = 21);
    save ("psi.png");
    clear();

    squares ("omega", spread = -1);
    isoline ("omega", n = 21);
    save ("omega.png");
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

    double ekin = 0.0, circ = 0.0, enst = 0.0, area_tot=0.0;
    foreach(reduction(+:ekin), reduction(+:circ), reduction(+:enst), reduction(+:area_tot)){
      double vort = omega[];
      double uvel = u.x[];
      double vvel = u.y[];
      double area = dv();
      ekin += area * (sq(uvel) + sq(vvel));
      circ += area * vort;
      enst += area * sq(vort);
      area_tot += area;
    }
    double utan1=0.0;
    foreach_boundary (left, reduction(+:utan1)){
      coord n = {-1,0};
      utan1 += (-u.x[] * n.y + u.y[] * n.x) * pow(Delta, dimension - 1);
    }
    double utan2=0.0;
    foreach_boundary (right, reduction(+:utan2)){
      coord n = {1,0};
      utan2 += (-u.x[] * n.y + u.y[] * n.x) * pow(Delta, dimension - 1);
    }
    double utan3=0.0;
    foreach_boundary (top, reduction(+:utan3)){
      coord n = {0,1};
      utan3 += (-u.x[] * n.y + u.y[] * n.x) * pow(Delta, dimension - 1);
    }
    double utan4=0.0;
    foreach_boundary (bottom, reduction(+:utan4)){
      coord n = {0,-1};
      utan4 += (-u.x[] * n.y + u.y[] * n.x) * pow(Delta, dimension - 1);
    }
    double utan= utan1 + utan2 + utan3 + utan4;

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

    fprintf (stderr, "E: %g \t C: %g %g \t Z: %g \t A: %g \t max(omega) %g \n", ekin, circ, utan, enst, area_tot, s0.max);
    int nj = 0;
    fprintf (stdout, "Vortex 1: %d %g %g %g %g %g %g %g \n", N, circ2[nj], mu_x[nj], mu_y[nj], interpolate(omega, mu_x[nj], mu_y[nj]), sqrt(M20[nj] + M02[nj]), sqrt((M20[nj] + M02[nj] + sqrt(4*sq(M11[nj]) + sq(M20[nj] - M02[nj])))), sqrt((M20[nj] + M02[nj] - sqrt(4*sq(M11[nj]) + sq(M20[nj] - M02[nj])))));
    nj = 1;
    fprintf (stdout, "Vortex 2: %d %g %g %g %g %g %g %g \n", N, circ2[nj], mu_x[nj], mu_y[nj], interpolate(omega, mu_x[nj], mu_y[nj]), sqrt(M20[nj] + M02[nj]), sqrt((M20[nj] + M02[nj] + sqrt(4*sq(M11[nj]) + sq(M20[nj] - M02[nj])))), sqrt((M20[nj] + M02[nj] - sqrt(4*sq(M11[nj]) + sq(M20[nj] - M02[nj])))));
    fprintf (stderr, "\n");
  }
}
