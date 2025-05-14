#include "poisson.h"
#include "view.h"

double lamb_oseen (double r, double a){
  return 1.0/(pi*sq(a)) * exp( -sq(r/a) );
}

double omega_mean ;
int main(){
  L0 = 16;

  for (N = 32; N <= 2048; N *= 2) {
    origin (-L0/2., -L0/2.);
    init_grid (N);

    /**
    Conditions on the box boundaries are set.
    We use "third-order" [face flux interpolation](/src/embed.h).
    */

    scalar a[], b[];

    /**
    The right-hand-side is defined using the coordinates of the
    barycenter of the cut cell (xc,yc), which is calculated from the
    cell and surface fractions. */
    double bsum = 0.;
    foreach(reduction(+:bsum)) {
      double r = sqrt(sq(x) + sq(y));
      b[] = lamb_oseen(r, 1.0);
      bsum += b[]*dv();
    }
    fprintf (stderr, "bsum (before): %g\n", bsum);

    omega_mean = -0.5*bsum/sq(L0);

    a[left]   = dirichlet(-(sq(L0/2) + sq(y))*omega_mean/2.0);
    a[right]  = dirichlet(-(sq(L0/2) + sq(y))*omega_mean/2.0);
    a[top]    = dirichlet(-(sq(x) + sq(L0/2))*omega_mean/2.0);
    a[bottom] = dirichlet(-(sq(x) + sq(L0/2))*omega_mean/2.0);

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
    u.n[top]    = dirichlet(-x*omega_mean); u.t[top]    = neumann( omega_mean);
    u.n[bottom] = dirichlet( x*omega_mean); u.t[bottom] = neumann(-omega_mean);
    u.n[left]   = dirichlet( x*omega_mean); u.t[left]   = neumann(-omega_mean);
    u.n[right]  = dirichlet(-x*omega_mean); u.t[right]  = neumann( omega_mean);

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
    psi[left]   = dirichlet(-(sq(L0/2) + sq(y))*omega_mean/2.0);
    psi[right]  = dirichlet(-(sq(L0/2) + sq(y))*omega_mean/2.0);
    psi[top]    = dirichlet(-(sq(x) + sq(L0/2))*omega_mean/2.0);
    psi[bottom] = dirichlet(-(sq(x) + sq(L0/2))*omega_mean/2.0);

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

    foreach()
      omega[] -= 2*omega_mean;

    double ekin = 0.0, circ = 0.0, enst = 0.0, area_tot=0.0;
    foreach(reduction(+:ekin), reduction(+:circ), reduction(+:enst), reduction(+:area_tot)){
      double vort = omega[];
      double uvel = u.x[] + y*omega_mean;
      double vvel = u.y[] - x*omega_mean;
      double area = dv();
      ekin += area * (sq(uvel) + sq(vvel));
      circ += area * vort;
      enst += area * sq(vort);
      area_tot += area;
    }
    double utan1=0.0;
    foreach_boundary (left, reduction(+:utan1)){
      coord n = {-1,0};
      utan1 += (-(0.5*u.x[]+0.5*u.x[ghost] + y*omega_mean) * n.y + (0.5*u.y[]+0.5*u.y[ghost] - (-L0/2)*omega_mean) * n.x) * pow(Delta, dimension - 1);
    }
    double utan2=0.0;
    foreach_boundary (right, reduction(+:utan2)){
      coord n = {1,0};
      utan2 += (-(0.5*u.x[]+0.5*u.x[ghost] + y*omega_mean) * n.y + (0.5*u.y[]+0.5*u.y[ghost] - (L0/2)*omega_mean) * n.x) * pow(Delta, dimension - 1);
    }
    double utan3=0.0;
    foreach_boundary (top, reduction(+:utan3)){
      coord n = {0,1};
      utan3 += (-(0.5*u.x[]+0.5*u.x[ghost] + (L0/2)*omega_mean) * n.y + (0.5*u.y[]+0.5*u.y[ghost] - x*omega_mean) * n.x) * pow(Delta, dimension - 1);
    }
    double utan4=0.0;
    foreach_boundary (bottom, reduction(+:utan4)){
      coord n = {0,-1};
      utan4 += (-(0.5*u.x[]+0.5*u.x[ghost] + (-L0/2)*omega_mean) * n.y + (0.5*u.y[]+0.5*u.y[ghost] - x*omega_mean) * n.x) * pow(Delta, dimension - 1);
    }
    double utan= utan1 + utan2 + utan3 + utan4;

    s0 = statsf (omega);

    double circ2 = 0., mu_x = 0., mu_y = 0.;
    foreach(reduction(+:circ2), reduction(+:mu_x), reduction(+:mu_y)){
      coord p = {x,y,z};
      circ2 += dv() * omega[] ;
      mu_x += dv() * omega[]  * p.x;
      mu_y += dv() * omega[]  * p.y;
    }

    double M20 = 0., M02 = 0., M11 = 0.;
    foreach(reduction(+:M20), reduction(+:M02), reduction(+:M11)){
      coord p = {x,y,z};
      M20 += dv() * omega[] * sq(p.x - mu_x);
      M02 += dv() * omega[] * sq(p.y - mu_y);
      M11 += dv() * omega[] * (p.x - mu_x) * (p.y - mu_y);
    }

    M20 /= circ2;
    M02 /= circ2;
    M11 /= circ2;

    fprintf (stderr, "E: %g \t C: %g %g %g \t Z: %g \t A: %g \t max(omega): %g \t mu: %g %g \t a: %g %g %g \n", ekin, circ, circ2, utan, enst, area_tot, s0.max, mu_x, mu_y, sqrt(M20 + M02), sqrt((M20 + M02 + sqrt(4*sq(M11) + sq(M20 - M02)))), sqrt((M20 + M02 - sqrt(4*sq(M11) + sq(M20 - M02)))));
    fprintf (stderr, "\n");
  }
}
