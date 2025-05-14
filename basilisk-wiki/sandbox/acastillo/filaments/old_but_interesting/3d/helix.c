#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "fractions.h"
#include "lambda2.h"
#include "view.h"
#include "PointTriangle.h"

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

/**
# Simulated Navier-Stokes helical vortex

The evolution of a helical vortex is simulated as in
[Antoon's sandbox](http://basilisk.fr/sandbox/Antoonvh/helical.c). The goal
of the physical space initialisation is to map an analytically defined trefoil
vortex onto an Eulerian (static) numerical mesh. The helical trajectory
discussed in here is defined by:
$$
x = R \cos{(t)}
$$
$$
y = R \sin{(t)}
$$
$$
z = \frac{H}{2\pi} t
$$
where $t$ ranges between 0 and $2\pi$, $R$ is the radius and $H$ the pitch.

![Space-curve of a helical vortex relative to the domain size](helix/helix0.png)
*/

int n_seg;
double Gamma= 1.0, as = 0.1, Hs = 1.0, Rs = 0.5, n_turns;

#include "../filaments.h"
#include "../filaments.c"
#include "../../output_fields/output_matrix_normal.h"

#define sortme(x) a = x[i]; x[i] = x[j]; x[j] = a;
#define vecdist2(a,b) (sq((a).x - (b).x) + sq((a).y - (b).y) + sq((a).z - (b).z))
#define MINLEVEL 5

void helical_filament(vector omega, double Gamma){
  // Define the helix as the spacecurve c
  double delta_t0 = (2 + n_turns)*2*pi/((double)n_seg-1);
  double t0[n_seg];
  coord c[n_seg], dc[n_seg], d2c[n_seg];

  for (int i = 0; i < n_seg; i++){
    t0[i] = delta_t0 * (double)i - 2*pi;
    c[i].x = Rs * cos(t0[i]);
    c[i].y = Rs * sin(t0[i]);
    c[i].z = (Hs/(2*pi)) * t0[i] - L0/2 ;
  }

  view (camera="iso", fov=40);
  draw_space_curve(n_seg, c);
  box();
  save ("helix0.png");

  // Find the first and second derivatives to compute a Frenet-Serret Frame
  coord xshift = {0, 0, (2 + n_turns)*Hs}, dxshift = {0, 0, 0};
  fd_derivative(n_seg, delta_t0,  xshift,  c,  dc);
  fd_derivative(n_seg, delta_t0, dxshift, dc, d2c);

  coord tvec[n_seg], nvec[n_seg], bvec[n_seg];
  for (int i = 0; i < n_seg; i++){
    foreach_dimension(){
      tvec[i].x =  dc[i].x/sqrt(vecdot( dc[i],  dc[i]));
      nvec[i].x = d2c[i].x/sqrt(vecdot(d2c[i], d2c[i]));
    }
    bvec[i] = vecdotproduct(tvec[i], nvec[i]);
  }

  // Initialize the vorticity field
  foreach(){
    coord pcar = {x,y,z};
    coord val_omega = get_vorticity_filament(pcar, n_seg, as, t0, c, tvec, nvec, bvec, 0);
    foreach_dimension()
      omega.x[] = Gamma*val_omega.x;
  }
}

/**
For this example, we use the (compressible) Navier-Stokes equations inside a
triply periodic box. Periodicity implies that the net circulation contained
within the domain must be precisely zero, which can trigger an artificial
centrifugal instability of an isolated vortex (see,
[Pradeep and Hussain](DOI: 10.1017/S002211200400076X)).
We use the approach by [Otheguy, Chomaz, and Billant](doi:10.1017/S0022112006008901)
and work in a rotating-frame with angular velocity $\Omega^F$, such that
the circulation in the rotating frame is zero.
$$
\Omega^F = \Gamma_0/2L_xL_y
$$
*/

double omega_zero ;
int main(){
  L0 = 8;
  X0 = Y0 = Z0 = -L0/2;
  omega_zero = (-0.5*Gamma/sq(L0));

  init_grid (1 << MINLEVEL);

  //TOLERANCE = 1e-6;

  n_turns = L0/Hs;
  n_seg = 64 * n_turns;

  periodic(top);
  periodic(left);
  periodic(front);
  run();
}

event adapt (i++){
  adapt_wavelet ((scalar*){u}, (double[]){1e-4,1e-4,1e-4}, MAXLEVEL, MINLEVEL);
}


#define K0() (0.)
#define F0() (2.0 * omega_zero)
#define alpha_H 1.0
event end_timestep (i++){
  foreach(){
   coord b0 = { - K0(), - K0() }, b1 = { F0(), -F0() };
	  coord m0 = { 1. - alpha_H*dt*b0.x, 1. - alpha_H*dt*b0.y };
	  coord m1 = { - alpha_H*dt*b1.x, - alpha_H*dt*b1.y };
	  double det = m0.x*m0.y - m1.x*m1.y;
    coord r;

	  r.x = u.x[] + (1. - alpha_H)*dt*(b0.x*u.x[] + b1.x*u.y[]) + dt*g.x[];
    r.y = u.y[] + (1. - alpha_H)*dt*(b0.y*u.y[] + b1.y*u.x[]) + dt*g.y[];
    r.z = u.z[] + dt*g.z[];

    u.x[] = (m0.y*r.x - m1.x*r.y)/det - dt*g.x[];
    u.y[] = (m0.x*r.y - m1.y*r.x)/det - dt*g.x[];
    u.z[] = r.z - dt*g.x[];
  }
  boundary ((scalar *){u});
}



event output (i = 20){
  printf( "Done!... \n" );
}

/**
## Initial conditions
We consider the flow is expressed as a vector potential,
$$
\vec{u} = \nabla \times \vec{\psi}
$$
such that
$$
\nabla^2 \vec{\psi} = -\vec{\omega}
$$
which requires solving a Poisson problem for each vorticity component.

![Isocontour of $\lambda_2=0$ for the initial conditions](helix/lambda2.png)

*/

event init (t = 0){

  FILE * fp;
  char header[180];
  sprintf(header, "[1]t\t [2]tag\t [3]Atot\t [4]Gamma\t [5]mu_x\t [6]mu_y\t [7]mu_z\t [8]M20\t [9]M02\t [10]M11\t [11]a\t [12]b\t [13]c\t [14]e\t [15]maxvor \n");
  fp = fopen("vortex_x0_l2.asc", "w"); fputs(header, fp); fclose (fp);
  fp = fopen("vortex_y0_l2.asc", "w"); fputs(header, fp); fclose (fp);
  fp = fopen("vortex_x0.asc", "w"); fputs(header, fp); fclose (fp);
  fp = fopen("vortex_y0.asc", "w"); fputs(header, fp); fclose (fp);
  fp = fopen("globals.asc", "w"); fclose (fp);
  fp = fopen("closed_paths.asc", "w"); fclose (fp);

#ifdef _INIT
  printf( "Creating initial conditions... \n" );

  refine  ((sqrt(sq(x) + sq(y))) < 16*Rs && level < MAXLEVEL-3);
  refine  ((sqrt(sq(x) + sq(y))) < 8*Rs && level < MAXLEVEL-2);
  refine  ((sqrt(sq(x) + sq(y))) < 4*Rs && level < MAXLEVEL-1);
  refine  ((sqrt(sq(x) + sq(y))) < 2*Rs && level < MAXLEVEL);

  vector c[], b[];
  helical_filament(b, Gamma);

  double bsum = 0.;
  foreach(reduction(+:bsum)) {
    bsum += b.z[]*dv();
  }
  fprintf (stderr, "bsum (before): %g\n", bsum);
  fprintf (stderr, "omega_zero (before): %.15g\n", omega_zero);
  omega_zero = -0.5*bsum/pow(L0,3);

  foreach()
    foreach_dimension()
      c.x[] = 0.;

  double bsum_x = 0., bsum_y = 0., bsum_z = 0.;
  foreach(reduction(+:bsum_x), reduction(+:bsum_y), reduction(+:bsum_z)) {
    double r = sqrt(sq(x) + sq(y));
    b.z[] += 2.0*omega_zero;
    c.z[] -= sq(r)*omega_zero/2.0;
    foreach_dimension()
      bsum_x += b.x[]*dv();
  }

  fprintf (stderr, "omega_zero: %.15g\n", omega_zero);
  foreach_dimension()
    fprintf (stderr, "bsum: %g\n", bsum_x);

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

  scalar l2[];
  lambda2 (u, l2);
  squares ("l2", linear = false); isosurface ("l2", -1e-3); cells();
  save ("lambda2.png");


  printf( "Done! \n" );
#else
  printf( "Loading initial conditions... \n" );
  restore("dump0");
  printf( "Done! \n" );

  event ("logfile");
  event ("slices");
#endif
}

/**
![Evolution of the isocontour of $\lambda_2=0$](helix/lambda2.mp4)
*/


//event movie (t += 0.05) {
event movie (i += 1) {
  scalar l2[];
  lambda2 (u, l2);

  vector omega[];
  vorticity3d(u, omega);

  view (camera="iso");
  squares ("u.x", linear = false, alpha=-L0/2);
  save ("u.png");

  squares ("u.y", linear = false, alpha=-L0/2);
  save ("v.png");

  squares ("u.z", linear = false, alpha=-L0/2);
  save ("w.png");

  squares ("l2", linear = false, alpha=-L0/2);
  isosurface ("l2", -1e-2);
  box();
  save ("lambda2.png");

  squares ("omega.z", linear = false, alpha=-L0/2);
  box();
  save ("omega.png");
}




/**
## Ellipticity of the vortex cores
The goal is to evaluate several quantities that are pertinent to describe the
flow field such as circulation, sizes, and ellipticity. When a vorticity field
possesses an elliptical shape, these quantities are sufficient to define its
geometry. In this case, we measure the ellipticity from the intersection of
each vortex with the planes $x=0$ and $y=0$. We also take a slice perpendicular
to the vortex tubes, which can be used later. Here, we use the vorticity at
the vortex centroids to the define the plane intersecting the vortices.
A more refined approach like the one proposed by [Selçuk](https://www.theses.fr/2016PA066138.pdf)
would be nice.

![Intersection between the vortex tubes and the plane $y=0$](helix/m3.png)


*/




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
  double * va;
  double * vb;
  double * vc;
  double * ve;
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

  double var0[nm], var1[nm], var2[nm], var3[nm], var4[nm];
  for (int nj = 0; nj < nm; nj++){
    pcen.centroids_x[nj] = 0.0;
    pcen.centroids_y[nj] = 0.0;
    pcen.centroids_z[nj] = 0.0;
    var0[nj] = 0.0;
    var1[nj] = 0.0;
    var2[nj] = 0.0;
    var3[nj] = 0.0;
    var4[nj] = 0.0;
  }

  // Compute the circulation $\Gamma = \iint \omega dA$ and position of the
  // vortex centroids $\vec{\mu} = \iint \vec{x}\omega dA$
  foreach(reduction(+:var0[:nm]), reduction(+:var1[:nm]), reduction(+:var2[:nm]), reduction(+:var3[:nm]), reduction(+:var4[:nm])){
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

      for (int nj = 0; nj < nm; nj++){
        if (m[] == nj){
          var0[nj] += dA ;
          var1[nj] += dA * (sval - 2*omega_zero);
          var2[nj] += dA * (sval - 2*omega_zero) * p.x;
          var3[nj] += dA * (sval - 2*omega_zero) * p.y;
          var4[nj] += dA * (sval - 2*omega_zero) * p.z;
        }
      }
    }
  }

  for (int nj = 0; nj < nm; nj++){
    pcen.Atot[nj] = var0[nj];
    pcen.circulation[nj] = var1[nj];
    if (pcen.circulation[nj] != 0) {
      pcen.centroids_x[nj] = var2[nj]/var1[nj];
      pcen.centroids_y[nj] = var3[nj]/var1[nj];
      pcen.centroids_z[nj] = var4[nj]/var1[nj];
    }
  }

  for (int i = 0; i < nm; i++){
    for (int j = i + 1; j < nm; ++j){
      if (abs(pcen.circulation[i]) < abs(pcen.circulation[j])){
        double a;
        sortme(pcen.circulation);
        sortme(pcen.centroids_x);
        sortme(pcen.centroids_y);
        sortme(pcen.centroids_z);
      }
    }
  }
}

trace
void vorticity_moments_plane(scalar omega, scalar m, struct VorticityCentroids pcen, struct VorticityMoments pmom){

  int nm = pcen.nm;
  double omega_zero = pcen.omega_zero, _alpha=pcen._alpha;
  coord n = pcen.n;

  double var0[nm], var1[nm], var2[nm], var3[nm], var4[nm];
  for (int nj = 0; nj < nm; nj++){
    pcen.centroids_x[nj] = 0.0;
    pcen.centroids_y[nj] = 0.0;
    pcen.centroids_z[nj] = 0.0;
    var0[nj] = 0.0;
    var1[nj] = 0.0;
    var2[nj] = 0.0;
    var3[nj] = 0.0;
    var4[nj] = 0.0;
  }

  // Compute the circulation $\Gamma = \iint \omega dA$ and position of the
  // vortex centroids $\vec{\mu} = \iint \vec{x}\omega dA$
  foreach(reduction(+:var0[:nm]), reduction(+:var1[:nm]), reduction(+:var2[:nm]), reduction(+:var3[:nm]), reduction(+:var4[:nm])){
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

      for (int nj = 0; nj < nm; nj++){
        if (m[] == nj){
          var0[nj] += dA ;
          var1[nj] += dA * (sval - 2*omega_zero);
          var2[nj] += dA * (sval - 2*omega_zero) * p.x;
          var3[nj] += dA * (sval - 2*omega_zero) * p.y;
          var4[nj] += dA * (sval - 2*omega_zero) * p.z;
        }
      }
    }
  }

  for (int nj = 0; nj < nm; nj++){
    pcen.Atot[nj] = var0[nj];
    pcen.circulation[nj] = var1[nj];
    if (pcen.circulation[nj] != 0) {
      pcen.centroids_x[nj] = var2[nj]/var1[nj];
      pcen.centroids_y[nj] = var3[nj]/var1[nj];
      pcen.centroids_z[nj] = var4[nj]/var1[nj];
    }
  }

  double var5[nm], var6[nm];
  for (int nj = 0; nj < nm; nj++){
    var0[nj] = 0.0;
    var1[nj] = 0.0;
    var2[nj] = 0.0;
    var3[nj] = 0.0;
    var4[nj] = 0.0;
    var5[nj] = 0.0;
    var6[nj] = 0.0;
  }

  // Compute the circulation $\Gamma = \iint \omega dA$ and position of the
  // vortex centroids $\vec{\mu} = \iint \vec{x}\omega dA$
  foreach(reduction(+:var1[:nm]), reduction(+:var2[:nm]), reduction(+:var3[:nm]), reduction(+:var4[:nm]), reduction(+:var5[:nm]), reduction(+:var6[:nm])){
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

      for (int nj = 0; nj < nm; nj++){
        if (m[] == nj){
          var1[nj] += dA * (sval - 2*omega_zero) * sq(p.x - pcen.centroids_x[nj]);
          var2[nj] += dA * (sval - 2*omega_zero) * sq(p.y - pcen.centroids_y[nj]);
          var3[nj] += dA * (sval - 2*omega_zero) * sq(p.z - pcen.centroids_z[nj]);
          var4[nj] += dA * (sval - 2*omega_zero) * (p.x - pcen.centroids_x[nj]) * (p.y - pcen.centroids_y[nj]);
          var5[nj] += dA * (sval - 2*omega_zero) * (p.y - pcen.centroids_y[nj]) * (p.z - pcen.centroids_z[nj]);
          var6[nj] += dA * (sval - 2*omega_zero) * (p.z - pcen.centroids_z[nj]) * (p.x - pcen.centroids_x[nj]);
        }
      }
    }
  }

  for (int nj = 0; nj < nm; nj++){
    if (n.x == 1){
      pmom.M20[nj] = var2[nj]/pcen.circulation[nj];
      pmom.M02[nj] = var3[nj]/pcen.circulation[nj];
      pmom.M11[nj] = var5[nj]/pcen.circulation[nj];
    }
    else if (n.y == 1){
      pmom.M20[nj] = var1[nj]/pcen.circulation[nj];
      pmom.M02[nj] = var3[nj]/pcen.circulation[nj];
      pmom.M11[nj] = var6[nj]/pcen.circulation[nj];
    }
    else{
      pmom.M20[nj] = var1[nj]/pcen.circulation[nj];
      pmom.M02[nj] = var2[nj]/pcen.circulation[nj];
      pmom.M11[nj] = var4[nj]/pcen.circulation[nj];
    }
    pmom.va[nj] = sqrt(pmom.M20[nj] + pmom.M02[nj]);
    pmom.vb[nj] = sqrt((pmom.M20[nj] + pmom.M02[nj] + sqrt(4*sq(pmom.M11[nj]) + sq(pmom.M20[nj] - pmom.M02[nj]))));
    pmom.vc[nj] = sqrt((pmom.M20[nj] + pmom.M02[nj] - sqrt(4*sq(pmom.M11[nj]) + sq(pmom.M20[nj] - pmom.M02[nj]))));
    pmom.ve[nj] = sqrt(1 - sq(pmom.vc[nj])/sq(pmom.vb[nj]));
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
      if (abs(pcen.circulation[i]) < abs(pcen.circulation[j])){
        double a;
        sortme(pcen.Atot);
        sortme(pcen.circulation);
        sortme(pcen.centroids_x);
        sortme(pcen.centroids_y);
        sortme(pcen.centroids_z);
        sortme(pmom.M20);
        sortme(pmom.M02);
        sortme(pmom.M11);
        sortme(pmom.va);
        sortme(pmom.vb);
        sortme(pmom.vc);
        sortme(pmom.ve);
        sortme(pmom.omax);
      }
    }
  }
}


















void global_quantities(FILE * fp, vector omega, vector u){

  double ekin = 0.0, circ[3] = {0.0,0.0,0.0}, circ2=0.0, enst = 0.0, volu_tot=0.0, cents[3] = {0.0,0.0,0.0};
  foreach(reduction(+:ekin), reduction(+:circ[:3]), reduction(+:enst), reduction(+:volu_tot), reduction(+:cents[:3])){
    ekin     += dv() * (sq(u.x[]) + sq(u.y[]) + sq(u.z[]));
    circ[0]  += dv() * omega.x[];
    circ[1]  += dv() * omega.y[];
    circ[2]  += dv() * omega.z[];
    enst     += dv() * (sq(omega.x[])+sq(omega.y[])+sq(omega.z[]));
    volu_tot += dv();

    coord p = {x,y,z};
    cents[0] += dv() * omega.z[] * p.x;
    cents[1] += dv() * omega.z[] * p.y;
    cents[2] += dv() * omega.z[] * p.y;
  }

  fprintf (fp, "%.5f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f \n", t, omega_zero, ekin, circ[0], circ[1], circ[2], enst, volu_tot, cents[0], cents[1], cents[2]);
}

void global_quantities2(FILE * fp, vector u, double xc, double yc, double radius, int nj){

  vertex scalar phi[];
  foreach_vertex()
    phi[] = sqrt(sq(x-xc)+sq(y-yc)) - radius;
  boundary ({phi});

  scalar c[];
  face vector s[];
  fractions (phi, c, s);
  squares ("c", spread = -1);
  save ("c2.png");

  double area = 0., utang=0., unorm=0.;
  foreach (reduction(+:area), reduction(+:utang), reduction(+:unorm)){
    if (c[] > 0 && c[] < 1.) {
      coord n = interface_normal (point, c), p;
      double alpha = plane_alpha (c[], n);
      double dl = pow(Delta, dimension - 1)*plane_area_center (n, alpha, &p);
      normalize (&n);

      double uvel = interpolate (u.x, x+p.x*Delta, y+p.y*Delta, z + p.z*Delta );
      double vvel = interpolate (u.y, x+p.x*Delta, y+p.y*Delta, z + p.z*Delta );
      double tvel = ( uvel * n.y - vvel * n.x);
      double nvel = ( uvel * n.x + vvel * n.y);
      utang += tvel * dl;
      unorm += nvel * dl;
      area += dl;
    }
  }
  fprintf (fp, "%.5f %d %.15f %.15f %.15f \n", t, nj, area, utang, unorm);
}


//event logfile (t += 0.05) {
event logfile (i += 1) {

  vector omega[], u_inert[];
  vorticity3d(u, omega);

  foreach(){
    u_inert.x[] = u.x[] + y*omega_zero;
    u_inert.y[] = u.y[] - x*omega_zero;
    u_inert.z[] = u.z[];
    omega.z[] -= 2.0*omega_zero;
  }

  // Write the global quantities for verification
  FILE * fp = fopen("globals.asc", "a");
  global_quantities(fp, omega, u_inert);
  fclose(fp);

  fp = fopen("closed_paths.asc", "a");
  global_quantities2(fp, u_inert, 0.0, 0.0, 3.0, 0);
  fclose(fp);

  scalar l2[];
  lambda2 (u, l2);

  // Compute the vorticity moments on the plane x0, once using the l2 criterion
  scalar m[];
  foreach()
    m[] = ( ( ( abs(x) <= 4*Delta ) & ( y >= 0 ) ) & ( l2[]<=0.0 ) );

  int nm = tag (m) + 1;
  {
    double Atot[nm], circulation[nm], centroids_x[nm], centroids_y[nm], centroids_z[nm], omax[nm];
    double M20[nm], M02[nm], M11[nm], va[nm], vb[nm], vc[nm], ve[nm];
    vorticity_moments_plane(  omega.x, m, (struct VorticityCentroids){nm, (coord){1,0,0}, 0.0, 0.0, circulation, centroids_x, centroids_y, centroids_z, Atot}, (struct VorticityMoments) {M20, M02, M11, va, vb, vc, ve, omax});
    for (int nj = 0; nj < nm; nj++){
      if (pid()==0) {
        fp = fopen("vortex_x0_l2.asc", "a");
        fprintf (fp, "%.5g %d %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g \n", t, nj, Atot[nj], circulation[nj], centroids_x[nj], centroids_y[nj], centroids_z[nj], M20[nj], M02[nj], M11[nj], va[nj], vb[nj], vc[nj], ve[nj], omax[nj]);
        fclose(fp);
      }
    }

    view(camera="iso");
    squares ("m", spread = -1, n = {1,0,0} );
    for (int nj = 1; nj < nm; nj++){
      scatter (centroids_x[nj], centroids_y[nj], centroids_z[nj]);
    }
    save ("m.png");

    foreach(){
      m[] = 0;
      coord p = {x,y,z};
      for (int nj = 1; nj < nm; nj++){
        coord c1 = {centroids_x[nj],centroids_y[nj],centroids_z[nj]};
        if (vecdist2(p, c1) < sq(Hs/2))
          m[] = nj;
      }
    }
    squares ("m", spread = -1,  n = {1,0,0} );
    for (int nj = 1; nj < nm; nj++){
      scatter (centroids_x[nj], centroids_y[nj], centroids_z[nj]);
    }
    save ("m2.png");
  }
  // then, using a region of interest around the centroids obtained from l2
  // this avoids cutting the low intensity vorticity
  {
    double Atot[nm], circulation[nm], centroids_x[nm], centroids_y[nm], centroids_z[nm], omax[nm];
    double M20[nm], M02[nm], M11[nm], va[nm], vb[nm], vc[nm], ve[nm];

    vorticity_moments_plane(  omega.x, m, (struct VorticityCentroids){nm, (coord){1,0,0}, 0.0, 0.0, circulation, centroids_x, centroids_y, centroids_z, Atot}, (struct VorticityMoments) {M20, M02, M11, va, vb, vc, ve, omax});

    for (int nj = 0; nj < nm; nj++){
      if (pid()==0) {
        fp = fopen("vortex_x0.asc", "a");
        fprintf (fp, "%.5g %d %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g \n", t, nj, Atot[nj], circulation[nj], centroids_x[nj], centroids_y[nj], centroids_z[nj], M20[nj], M02[nj], M11[nj], va[nj], vb[nj], vc[nj], ve[nj], omax[nj]);
        fclose(fp);
      }
    }
  }
}


event slices (t += 0.5) {
  static int nslices = 0;

  vector omega[], u_inert[];
  vorticity3d(u, omega);

  foreach(){
    u_inert.x[] = u.x[] + y*omega_zero;
    u_inert.y[] = u.y[] - x*omega_zero;
    u_inert.z[] = u.z[];
    omega.z[] -= 2.0*omega_zero;
  }
  boundary ((scalar *){u_inert, omega});

  scalar l2[];
  lambda2 (u, l2);

  // Compute the vorticity centroids on the plane x0
  scalar m[];
  foreach()
    m[] = ( ( ( abs(x) <= 4*Delta ) & ( y >= 0 ) ) & ( l2[]<=0.0 ) );

  int nm = tag (m) + 1;
  double Atot[nm], circulation[nm], centroids_x[nm], centroids_y[nm], centroids_z[nm];
  vorticity_centroids_plane(  omega.x, m, (struct VorticityCentroids){nm, (coord){1,0,0}, 0.0, 0.0, circulation, centroids_x, centroids_y, centroids_z, Atot});

  for (int nj = 1; nj < nm; nj++){
    coord P = {centroids_x[nj], centroids_y[nj], centroids_z[nj]};
    coord tvec;
    tvec.x = interpolate (omega.x, centroids_x[nj], centroids_y[nj], centroids_z[nj]);
    tvec.y = interpolate (omega.y, centroids_x[nj], centroids_y[nj], centroids_z[nj]);
    tvec.z = interpolate (omega.z, centroids_x[nj], centroids_y[nj], centroids_z[nj]);
    normalize(&tvec);

    double _alpha = vecdot(tvec, P);
    coord nvec = {(_alpha/tvec.x)-P.x, 0-P.y, 0-P.z};
    normalize(&nvec);
    coord bvec = vecdotproduct(tvec, nvec);

    char name[80];
    int nres = 256;
    double len_window = Hs;
    if ((sqrt(sq(P.x)+sq(P.y)+sq(P.z)) + 1.1*len_window) < (L0/2)){

      // Slice in cartesian coordinates
      sprintf(name, "local_x0a_%3.3d_%3.3d.vts", nslices, nj);
      output_vts_normal((vector *) {u_inert, omega}, name, nres, len_window, P, tvec, nvec, bvec);

      // Slice in local coordinates
      sprintf(name, "local_x0b_%3.3d_%3.3d.vts", nslices, nj);
      output_vts_normal2((vector *) {u_inert, omega}, name, nres, len_window, P, tvec, nvec, bvec);

    }
    if (nj >= 10) break;
  }
  nslices++;
}
