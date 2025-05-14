/**
# The dynamics of a viscous vortex dipole

The evolution of a viscous vortex dipole is simulated as in [../2d/corotating.c]()
but in 3D.

![Vorticity field](corotating/omega.png)

## Initial conditions and simulation parameters
The initial condition used in the direct numerical simulation (DNS) is chosen to
be a superposition of two Lamb–Oseen vortices with circulation $\Gamma$ located
at $(\pm b/2, 0)$ in the ($x,y$)-plane. Velocity field is initialized from the
vorticity field  which requires solving a Poisson problem for each vorticity
component.
*/

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "lambda2.h"
#include "view.h"

trace
double lamb_oseen (double x, double y, double a, double b){
  return 1.0/(pi*sq(a)) * exp( -(sq(x-b)+sq(y))/sq(a) );
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

#define sortme(x) a = x[i]; x[i] = x[j]; x[j] = a;
#define vecdist2(a,b) (sq((a).x - (b).x) + sq((a).y - (b).y) + sq((a).z - (b).z))
#define MINLEVEL 5

double omega_zero ;
int main(){
  L0 = 8;
  X0 = Y0 = Z0 = -L0/2;
  init_grid (1 << MINLEVEL);
  periodic(top);
  periodic(left);
  periodic(front);
  run();
}

event adapt (i++){
  adapt_wavelet ((scalar*){u}, (double[]){1e-4,1e-4,1e-4}, MAXLEVEL, MINLEVEL);
}

event acceleration (i++) {
  face vector av = a;
  foreach_face(x)
  	av.x[] = -2.0*omega_zero*uf.y[];

  foreach_face(y)
    av.y[] = 2.0*omega_zero*uf.x[];

  foreach_face(z)
    av.z[] = 0.0;
}

event init (t = 0){

  a = new face vector;

#ifdef _INIT
  printf( "Creating initial conditions... \n" );
  double ad = 0.1, bd=0.5, Gamma= 1.0;

  omega_zero = (-Gamma/sq(L0));

  refine  ((sqrt(sq(x) + sq(y))) < 16*bd && level < MAXLEVEL-4);
  refine  ((sqrt(sq(x) + sq(y))) < 8*bd && level < MAXLEVEL-3);
  refine  ((sqrt(sq(x) + sq(y))) < 4*bd && level < MAXLEVEL-2);
  refine  ((sqrt(sq(x) + sq(y))) < 2*bd && level < MAXLEVEL-1);
  refine  ((sqrt(sq(x-bd) + sq(y))) < 8*ad && level < MAXLEVEL);
  refine  ((sqrt(sq(x+bd) + sq(y))) < 8*ad && level < MAXLEVEL);


  vector c[], b[];
  foreach(){
    foreach_dimension(){
      b.x[] = 0.;
      c.x[] = 0.;
    }
  }

  double bsum = 0.;
  foreach(reduction(+:bsum)) {
    double r = sqrt(sq(x) + sq(y));
    b.z[] += Gamma*lamb_oseen(x,y,ad,-bd) + Gamma*lamb_oseen(x,y,ad,bd) + 2.0*omega_zero;
    c.z[] -= sq(r)*omega_zero/2.0;
    bsum += b.z[]*dv();
  }

  fprintf (stderr, "omega_zero: %.15g\n", omega_zero);
  fprintf (stderr, "bsum: %g\n", bsum);

  /** The Poisson equation is solved. */
  foreach_dimension(){
    timer t = timer_start();
    mgstats s = poisson (c.x, b.x, tolerance=1e-10, minlevel = 4);
    double dt = timer_elapsed (t);
    printf ("%d %g %d %d\n", N, dt, s.i, s.nrelax);
  }
  boundary ((scalar *){c});

  cells(); save ("cells.png"); clear();
  squares ("c.x", spread = -1); save ("a_x.png"); clear();
  squares ("c.y", spread = -1); save ("a_y.png"); clear();
  squares ("c.z", spread = -1); save ("a_z.png"); clear();

  squares ("b.x", spread = -1); save ("b_x.png"); clear();
  squares ("b.y", spread = -1); save ("b_y.png"); clear();
  squares ("b.z", spread = -1); save ("b_z.png"); clear();

  foreach(){
    u.x[] = -((c.z[0,1,0] - c.z[0,-1,0]) - (c.y[0,0,1] - c.y[0,0,-1]))/(2.*Delta);
    u.y[] = -((c.x[0,0,1] - c.x[0,0,-1]) - (c.z[1,0,0] - c.z[-1,0,0]))/(2.*Delta);
    u.z[] = -((c.y[1,0,0] - c.y[-1,0,0]) - (c.x[0,1,0] - c.x[0,-1,0]))/(2.*Delta);
  }
  boundary ((scalar *){u});

  vector omega[];
  vorticity3d (u, omega);
  boundary ((scalar *){omega});

  dump ("dump0", list = (scalar *){u});

  foreach_dimension(){
    stats s0_x = statsf (omega.x);
    fprintf (stderr, "omega: %g %g %g \n", s0_x.min, s0_x.max, s0_x.sum);
  }

  squares ("omega.x", spread = -1); cells(); save ("omega_x.png"); clear();
  squares ("omega.y", spread = -1); cells(); save ("omega_y.png"); clear();
  squares ("omega.z", spread = -1); cells(); save ("omega_z.png"); clear();

  printf( "Done! \n" );
#else
  omega_zero = -1./sq(L0);
  printf( "Loading initial conditions... \n" );
  restore("dump0");
  printf( "Done! \n" );
#endif

  FILE * fp;
  char header[180];
  sprintf(header, "[1]t\t [2]tag\t [3]Atot\t [4]Gamma\t [5]mu_x\t [6]mu_y\t [7]mu_z\t [8]M20\t [9]M02\t [10]M11\t [11]a\t [12]b\t [13]c\t [14]e\t [15]maxvor \n");
  fp = fopen("vortex_z0_l2.asc", "w"); fputs(header, fp); fclose (fp);
  fp = fopen("vortex_z0.asc", "w"); fputs(header, fp); fclose (fp);
  fp = fopen("globals.asc", "w"); fclose (fp);
}







struct VorticityCentroids {
  int nm;
  coord n;
  double _alpha;
  double omega_zero;
  double * circulation;
  double * centroids_x;
  double * centroids_y;
  double * centroids_z;
  double * Atot;
};

struct VorticityMoments {
  double * M20;
  double * M02;
  double * M11;
  double * omax;
};

trace
void scatter (double x, double y, double z){
  glPointSize(40.0f);
  glColor3f(0.9f,0.9f,0.9f);
  glBegin (GL_POINTS);
  glVertex3f(x, y, z);
  glEnd();
}

#include "tag.h"
#define dA (0.5 * sq(Delta))
trace
void vorticity_centroids_plane(scalar omega, scalar m, struct VorticityCentroids pcen){

  int nm = pcen.nm;
  double omega_zero = pcen.omega_zero, _alpha=pcen._alpha;
  coord n = pcen.n;

  for (int nj = 0; nj < nm; nj++){
    pcen.centroids_x[nj] = 0.0;
    pcen.centroids_y[nj] = 0.0;
    pcen.centroids_z[nj] = 0.0;

    // Compute the circulation $\Gamma = \iint \omega dA$ and position of the
    // vortex centroids $\vec{\mu} = \iint \vec{x}\omega dA$
    double var0 = 0., var1 = 0., var2 = 0., var3 = 0., var4 = 0.;
    foreach(reduction(+:var0), reduction(+:var1), reduction(+:var2), reduction(+:var3), reduction(+:var4)){
      double alpha = (_alpha - n.x*x - n.y*y - n.z*z)/Delta;
      if (fabs(alpha) > 0.87)
        continue;
      if (is_active(cell) && is_local(cell)) {
        coord p = {x,y,z};

        double sval;
        if (n.x == 1)
          sval = 0.5*(val(omega) + val(omega,1,0,0));
        else if (n.y == 1)
          sval = 0.5*(val(omega) + val(omega,0,1,0));
        else
          sval = 0.5*(val(omega) + val(omega,0,0,1));

        if (m[] == nj){
          var0 += dA ;
          var1 += dA * (sval - 2*omega_zero);
          var2 += dA * (sval - 2*omega_zero) * p.x;
          var3 += dA * (sval - 2*omega_zero) * p.y;
          var4 += dA * (sval - 2*omega_zero) * p.z;
        }
      }
    }

    pcen.Atot[nj] = var0;
    pcen.circulation[nj] = var1;
    if (pcen.circulation[nj] != 0) {
      pcen.centroids_x[nj] = var2/var1;
      pcen.centroids_y[nj] = var3/var1;
      pcen.centroids_z[nj] = var4/var1;
    }
  }

  for (int i = 0; i < nm; i++){
    for (int j = i + 1; j < nm; ++j){
      if (pcen.circulation[i] < pcen.circulation[j]){
        double a;
        sortme(pcen.circulation);
        sortme(pcen.centroids_x);
        sortme(pcen.centroids_y);
        sortme(pcen.centroids_z);
      }
    }
  }
}

v v v v v v v
=============
trace
void vorticity_moments_plane(scalar omega, scalar m, struct VorticityCentroids pcen, struct VorticityMoments pmom){

  int nm = pcen.nm;
  double omega_zero = pcen.omega_zero, _alpha=pcen._alpha;
  coord n = pcen.n;

  for (int nj = 0; nj < nm; nj++){
    pcen.centroids_x[nj] = 0.0;
    pcen.centroids_y[nj] = 0.0;
    pcen.centroids_z[nj] = 0.0;

    // Compute the circulation $\Gamma = \iint \omega dA$ and position of the
    // vortex centroids $\vec{\mu} = \iint \vec{x}\omega dA$
    double var0 = 0., var1 = 0., var2 = 0., var3 = 0., var4 = 0.;
    foreach(reduction(+:var0), reduction(+:var1), reduction(+:var2), reduction(+:var3), reduction(+:var4)){
      double alpha = (_alpha - n.x*x - n.y*y - n.z*z)/Delta;
      if (fabs(alpha) > 0.87)
        continue;
      if (is_active(cell) && is_local(cell)) {
        coord p = {x,y,z};
        double sval;
        if (n.x == 1)
          sval = 0.5*(val(omega) + val(omega,1,0,0));
        else if (n.y == 1)
          sval = 0.5*(val(omega) + val(omega,0,1,0));
        else
          sval = 0.5*(val(omega) + val(omega,0,0,1));

        if (m[] == nj){
          var0 += dA;
          var1 += dA * (sval - 2*omega_zero);
          var2 += dA * (sval - 2*omega_zero) * p.x;
          var3 += dA * (sval - 2*omega_zero) * p.y;
          var4 += dA * (sval - 2*omega_zero) * p.z;
        }
      }
    }

    pcen.Atot[nj] = var0;
    pcen.circulation[nj] = var1;
    if (pcen.circulation[nj] != 0) {
      pcen.centroids_x[nj] = var2/var1;
      pcen.centroids_y[nj] = var3/var1;
      pcen.centroids_z[nj] = var4/var1;
    }
  }

  for (int nj = 0; nj < nm; nj++){
    // Compute the circulation $\Gamma = \iint \omega dA$ and position of the
    // vortex centroids $\vec{\mu} = \iint \vec{x}\omega dA$
    double var1 = 0., var2 = 0., var3 = 0., var4 = 0., var5=0., var6=0.;
    foreach(reduction(+:var1), reduction(+:var2), reduction(+:var3), reduction(+:var4), reduction(+:var5), reduction(+:var6)){
      double alpha = (_alpha - n.x*x - n.y*y - n.z*z)/Delta;
      if (fabs(alpha) > 0.87)
        continue;
      if (is_active(cell) && is_local(cell)) {
        coord p = {x,y,z};
        double sval;
        if (n.x == 1)
          sval = 0.5*(val(omega) + val(omega,1,0,0));
        else if (n.y == 1)
          sval = 0.5*(val(omega) + val(omega,0,1,0));
        else
          sval = 0.5*(val(omega) + val(omega,0,0,1));

        if (m[] == nj){
          var1 += dA * (sval - 2*omega_zero) * sq(p.x - pcen.centroids_x[nj]);
          var2 += dA * (sval - 2*omega_zero) * sq(p.y - pcen.centroids_y[nj]);
          var3 += dA * (sval - 2*omega_zero) * sq(p.z - pcen.centroids_z[nj]);
          var4 += dA * (sval - 2*omega_zero) * (p.x - pcen.centroids_x[nj]) * (p.y - pcen.centroids_y[nj]);
          var5 += dA * (sval - 2*omega_zero) * (p.y - pcen.centroids_y[nj]) * (p.z - pcen.centroids_z[nj]);
          var6 += dA * (sval - 2*omega_zero) * (p.z - pcen.centroids_z[nj]) * (p.x - pcen.centroids_x[nj]);
        }
      }
    }

    if (n.x == 1){
      pmom.M20[nj] = var2/pcen.circulation[nj];
      pmom.M02[nj] = var3/pcen.circulation[nj];
      pmom.M11[nj] = var5/pcen.circulation[nj];
    }
    else if (n.y == 1){
      pmom.M20[nj] = var1/pcen.circulation[nj];
      pmom.M02[nj] = var3/pcen.circulation[nj];
      pmom.M11[nj] = var6/pcen.circulation[nj];
    }
    else{
      pmom.M20[nj] = var1/pcen.circulation[nj];
      pmom.M02[nj] = var2/pcen.circulation[nj];
      pmom.M11[nj] = var4/pcen.circulation[nj];
    }

    pmom.omax[nj] = interpolate (omega, pcen.centroids_x[nj], pcen.centroids_y[nj], pcen.centroids_z[nj]);
  }

  if (pid() == 0){
    @if _MPI
      MPI_Reduce (MPI_IN_PLACE, &pmom.omax, nm, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    @endif
  }
  @if _MPI
  else
    MPI_Reduce (&pmom.omax, NULL, nm, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  @endif

  for (int i = 0; i < nm; i++){
    for (int j = i + 1; j < nm; ++j){
      if (pcen.circulation[i] < pcen.circulation[j]){
        double a;
        sortme(pcen.Atot);
        sortme(pcen.circulation);
        sortme(pcen.centroids_x);
        sortme(pcen.centroids_y);
        sortme(pcen.centroids_z);
        sortme(pmom.M20);
        sortme(pmom.M02);
        sortme(pmom.M11);
        sortme(pmom.omax);
      }
    }
  }
}


*************
trace
void vorticity_moments_plane(scalar omega, scalar m, struct VorticityCentroids pcen, struct VorticityMoments pmom){

  int nm = pcen.nm;
  double omega_zero = pcen.omega_zero, _alpha=pcen._alpha;
  coord n = pcen.n;

  for (int nj = 0; nj < nm; nj++){
    pcen.centroids_x[nj] = 0.0;
    pcen.centroids_y[nj] = 0.0;
    pcen.centroids_z[nj] = 0.0;

    // Compute the circulation $\Gamma = \iint \omega dA$ and position of the
    // vortex centroids $\vec{\mu} = \iint \vec{x}\omega dA$
    double var0 = 0., var1 = 0., var2 = 0., var3 = 0., var4 = 0.;
    foreach(reduction(+:var0), reduction(+:var1), reduction(+:var2), reduction(+:var3), reduction(+:var4)){
      double alpha = (_alpha - n.x*x - n.y*y - n.z*z)/Delta;
      if (fabs(alpha) > 0.87)
        continue;
      if (is_active(cell) && is_local(cell)) {
        coord p = {x,y,z};
        double sval;
        if (n.x == 1)
          sval = 0.5*(val(omega) + val(omega,1,0,0));
        else if (n.y == 1)
          sval = 0.5*(val(omega) + val(omega,0,1,0));
        else
          sval = 0.5*(val(omega) + val(omega,0,0,1));

        if (m[] == nj){
          var0 += dA;
          var1 += dA * (sval - 2*omega_zero);
          var2 += dA * (sval - 2*omega_zero) * p.x;
          var3 += dA * (sval - 2*omega_zero) * p.y;
          var4 += dA * (sval - 2*omega_zero) * p.z;
        }
      }
    }

    pcen.Atot[nj] = var0;
    pcen.circulation[nj] = var1;
    if (pcen.circulation[nj] != 0) {
      pcen.centroids_x[nj] = var2/var1;
      pcen.centroids_y[nj] = var3/var1;
      pcen.centroids_z[nj] = var4/var1;
    }
  }

  for (int nj = 0; nj < nm; nj++){
    // Compute the circulation $\Gamma = \iint \omega dA$ and position of the
    // vortex centroids $\vec{\mu} = \iint \vec{x}\omega dA$
    double var1 = 0., var2 = 0., var3 = 0., var4 = 0., var5=0., var6=0.;
    foreach(reduction(+:var1), reduction(+:var2), reduction(+:var3), reduction(+:var4), reduction(+:var5), reduction(+:var6)){
      double alpha = (_alpha - n.x*x - n.y*y - n.z*z)/Delta;
      if (fabs(alpha) > 0.87)
        continue;
      if (is_active(cell) && is_local(cell)) {
        coord p = {x,y,z};
        double sval;
        if (n.x == 1)
          sval = 0.5*(val(omega) + val(omega,1,0,0));
        else if (n.y == 1)
          sval = 0.5*(val(omega) + val(omega,0,1,0));
        else
          sval = 0.5*(val(omega) + val(omega,0,0,1));

        if (m[] == nj){
          var1 += dA * (sval - 2*omega_zero) * sq(p.x - pcen.centroids_x[nj]);
          var2 += dA * (sval - 2*omega_zero) * sq(p.y - pcen.centroids_y[nj]);
          var3 += dA * (sval - 2*omega_zero) * sq(p.z - pcen.centroids_z[nj]);
          var4 += dA * (sval - 2*omega_zero) * (p.x - pcen.centroids_x[nj]) * (p.y - pcen.centroids_y[nj]);
          var5 += dA * (sval - 2*omega_zero) * (p.y - pcen.centroids_y[nj]) * (p.z - pcen.centroids_z[nj]);
          var6 += dA * (sval - 2*omega_zero) * (p.z - pcen.centroids_z[nj]) * (p.x - pcen.centroids_x[nj]);
        }
      }
    }

    if (n.x == 1){
      pmom.M20[nj] = var2/pcen.circulation[nj];
      pmom.M02[nj] = var3/pcen.circulation[nj];
      pmom.M11[nj] = var5/pcen.circulation[nj];
    }
    else if (n.y == 1){
      pmom.M20[nj] = var1/pcen.circulation[nj];
      pmom.M02[nj] = var3/pcen.circulation[nj];
      pmom.M11[nj] = var6/pcen.circulation[nj];
    }
    else{
      pmom.M20[nj] = var1/pcen.circulation[nj];
      pmom.M02[nj] = var2/pcen.circulation[nj];
      pmom.M11[nj] = var4/pcen.circulation[nj];
    }

    pmom.omax[nj] = interpolate (omega, pcen.centroids_x[nj], pcen.centroids_y[nj], pcen.centroids_z[nj]);
  }

  if (pid() == 0){
    @if _MPI
      MPI_Reduce (MPI_IN_PLACE, &pmom.omax, nm, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    @endif
  }
  @if _MPI
  else
    MPI_Reduce (&pmom.omax, NULL, nm, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  @endif

  for (int i = 0; i < nm; i++){
    for (int j = i + 1; j < nm; ++j){
      if (pcen.circulation[i] < pcen.circulation[j]){
        double a;
        sortme(pcen.Atot);
        sortme(pcen.circulation);
        sortme(pcen.centroids_x);
        sortme(pcen.centroids_y);
        sortme(pcen.centroids_z);
        sortme(pmom.M20);
        sortme(pmom.M02);
        sortme(pmom.M11);
        sortme(pmom.omax);
      }
    }
  }
}



















^ ^ ^ ^ ^ ^ ^

/**
## Ellipticity of the vortex dipole
The goal is to evaluate several quantities that are pertinent to describe the
flow field such as circulation, sizes, and ellipticity. When a vorticity field
possesses an elliptical shape, these quantities are sufficient to define its
geometry. In this case, we measure the ellipticity from the intersection of
each vortex with the plane $z=0$. In this case, we recover a complete rotation
of the vortex pair, the initial relaxation and transition to a viscous regime.
Results are very clean,
*/

event logfile (t += 0.05) {

  vector omega[];
  vorticity3d(u, omega);

  double ekin = 0.0, circ = 0.0, enst = 0.0, volu_tot=0.0;

  foreach(reduction(+:ekin), reduction(+:circ), reduction(+:enst), reduction(+:volu_tot)){
    ekin += dv() * (sq(u.x[]) + sq(u.y[]));
    circ += dv() * omega.z[];
    enst += dv() * sq(omega.x[]);
    volu_tot += dv();
  }

  double var0=0.0, var1=0.0, var2=0.0, var3=0.0;
  foreach_boundary (front, reduction(+:var0)){
    var0 +=  0.5*(u.y[]+u.y[ghost]) * pow(Delta, dimension - 1);
  }
  foreach_boundary (back, reduction(+:var1)){
    var1 += -0.5*(u.y[]+u.y[ghost]) * pow(Delta, dimension - 1);
  }
  foreach_boundary (top, reduction(+:var2)){
    var2 += -0.5*(u.z[]+u.z[ghost]) * pow(Delta, dimension - 1);
  }
  foreach_boundary (bottom, reduction(+:var3)){
    var3 +=  0.5*(u.z[]+u.z[ghost]) * pow(Delta, dimension - 1);
  }
  double circ_x = var0 + var1 + var2 + var3;


  var0=0.0, var1=0.0, var2=0.0, var3=0.0;
  foreach_boundary (front, reduction(+:var0)){
    var0 += -0.5*(u.x[]+u.x[ghost]) * pow(Delta, dimension - 1);
  }
  foreach_boundary (back, reduction(+:var1)){
    var1 +=  0.5*(u.x[]+u.x[ghost]) * pow(Delta, dimension - 1);
  }
  foreach_boundary (left, reduction(+:var2)){
    var2 += -0.5*(u.z[]+u.z[ghost]) * pow(Delta, dimension - 1);
  }
  foreach_boundary (right, reduction(+:var3)){
    var3 +=  0.5*(u.z[]+u.z[ghost]) * pow(Delta, dimension - 1);
  }
  double circ_y = var0 + var1 + var2 + var3;


  var0=0.0, var1=0.0, var2=0.0, var3=0.0;
  foreach_boundary (left, reduction(+:var0)){
    var0 += -0.5*(u.y[]+u.y[ghost]) * pow(Delta, dimension - 1);
  }
  foreach_boundary (right, reduction(+:var1)){
    var1 +=  0.5*(u.y[]+u.y[ghost]) * pow(Delta, dimension - 1);
  }
  foreach_boundary (top, reduction(+:var2)){
    var2 += -0.5*(u.x[]+u.x[ghost]) * pow(Delta, dimension - 1);
  }
  foreach_boundary (bottom, reduction(+:var3)){
    var3 +=  0.5*(u.x[]+u.x[ghost]) * pow(Delta, dimension - 1);
  }
  double circ_z = var0 + var1 + var2 + var3;

  var0=0.0, var1=0.0, var2=0.0;
  foreach(reduction(+:var0), reduction(+:var1), reduction(+:var2)){
    coord p = {x,y,z};
    var0 += dv() * (omega.z[] - 2.0*omega_zero) ;
    var1 += dv() * (omega.z[] - 2.0*omega_zero)  * p.x;
    var2 += dv() * (omega.z[] - 2.0*omega_zero)  * p.y;
  }
  FILE * fp = fopen("globals.asc", "a");
  fprintf (fp, "%.5g %.15g %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g \n", t, omega_zero, ekin, circ, enst, volu_tot, circ_x, circ_y, circ_z, var0, var1/var0, var2/var0);
  fclose(fp);

  scalar l2[];
  lambda2 (u, l2);

  scalar m[];
  foreach()
    m[] = l2[]<=0.0;
  int nm = tag (m) + 1;

  {
    double Atot[nm], circulation[nm], centroids_x[nm], centroids_y[nm], centroids_z[nm], M20[nm], M02[nm], M11[nm], omax[nm];
    vorticity_moments_plane(  omega.z, m, (struct VorticityCentroids){nm, (coord){0,0,1}, 0.0, omega_zero, circulation, centroids_x, centroids_y, centroids_z, Atot}, (struct VorticityMoments) {M20, M02, M11, omax});

    for (int nj = 0; nj < 2; nj++){
      double va = sqrt(M20[nj] + M02[nj]);
      double vb = sqrt((M20[nj] + M02[nj] + sqrt(4*sq(M11[nj]) + sq(M20[nj] - M02[nj]))));
      double vc = sqrt((M20[nj] + M02[nj] - sqrt(4*sq(M11[nj]) + sq(M20[nj] - M02[nj]))));
      double ve = sqrt(1 - sq(vc)/sq(vb));
      if (pid()==0) {
        fp = fopen("vortex_z0_l2.asc", "a");
        fprintf (fp, "%.5g %d %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g \n", t, nj, Atot[nj], circulation[nj], centroids_x[nj], centroids_y[nj], centroids_z[nj], M20[nj], M02[nj], M11[nj], va, vb, vc, ve, omax[nj]);
        fclose(fp);
      }
    }

    for (int i = 0; i < 2; i++){
      for (int j = i + 1; j < 2; ++j){
        if (centroids_x[i] < centroids_x[j]){
          double a;
          sortme(circulation);
          sortme(centroids_x);
          sortme(centroids_y);
          sortme(centroids_z);
        }
      }
    }

    view(camera="iso");
    squares ("m", spread = -1);
    scatter (centroids_x[0], centroids_y[0], centroids_z[0]);
    scatter (centroids_x[1], centroids_y[1], centroids_z[1]);
    save ("m.png");

    foreach(){
      coord p = {x,y,z};
      coord c1 = {centroids_x[0],centroids_y[0],z}, c2 = {centroids_x[1],centroids_y[1],z};
      if (vecdist2(p, c1) < 0.25)
        m[] = 1;
      else if (vecdist2(p, c2) < 0.25)
        m[] = 2;
      else
        m[] = 0;
    }

    squares ("m", spread = -1);
    scatter (centroids_x[0], centroids_y[0], centroids_z[0]);
    scatter (centroids_x[1], centroids_y[1], centroids_z[1]);
    save ("m2.png");

  }

  nm = 3;
  {
    double Atot[nm], circulation[nm], centroids_x[nm], centroids_y[nm], centroids_z[nm], M20[nm], M02[nm], M11[nm], omax[nm];
    vorticity_moments_plane(  omega.z, m, (struct VorticityCentroids){nm, (coord){0,0,1}, 0.0, omega_zero, circulation, centroids_x, centroids_y, centroids_z, Atot}, (struct VorticityMoments) {M20, M02, M11, omax});

    for (int nj = 0; nj < nm; nj++){
      double va = sqrt(M20[nj] + M02[nj]);
      double vb = sqrt((M20[nj] + M02[nj] + sqrt(4*sq(M11[nj]) + sq(M20[nj] - M02[nj]))));
      double vc = sqrt((M20[nj] + M02[nj] - sqrt(4*sq(M11[nj]) + sq(M20[nj] - M02[nj]))));
      double ve = sqrt(1 - sq(vc)/sq(vb));
      if (pid()==0) {
        fp = fopen("vortex_z0.asc", "a");
        fprintf (fp, "%.5g %d %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g \n", t, nj, Atot[nj], circulation[nj], centroids_x[nj], centroids_y[nj], centroids_z[nj], M20[nj], M02[nj], M11[nj], va, vb, vc, ve, omax[nj]);
        fclose(fp);
      }
    }
  }
}

event movie (t += 0.05) {
  scalar l2[];
  lambda2 (u, l2);

  vector omega[];
  vorticity3d(u, omega);

  view (camera="iso");
  squares ("u.x", linear = false);
  save ("u.mp4");

  squares ("u.y", linear = false);
  save ("v.mp4");

  squares ("l2", linear = false);
  isosurface ("l2", -1e-6);
  box();
  save ("lambda2.mp4");

  squares ("omega.z", linear = false);
  box();
  save ("omega.mp4");
}


event output (t = 20.0){
  printf( "Done!... \n" );
v v v v v v v
  //dump (list = (scalar *){u});

//
//
// /**
// # Additional outputs
//
// The vorticity and $\lambda_2$ field look something like this
//
// ![Vorticity field](corotating/omega.mp4) ![Vorticity field](corotating/lambda2.mp4)
//
// */
=============
  dump (list = (scalar *){u});
}
*************
  dump (list = (scalar *){u});
}
//
//
// /**
// # Additional outputs
//
// The vorticity and $\lambda_2$ field look something like this
//
// ![Vorticity field](corotating/omega.mp4) ![Vorticity field](corotating/lambda2.mp4)
//
// */
^ ^ ^ ^ ^ ^ ^
