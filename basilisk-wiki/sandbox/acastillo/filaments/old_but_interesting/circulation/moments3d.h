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
#define sortme(x) a = x[i]; x[i] = x[j]; x[j] = a;
#include "PointTriangle.h"
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

  @if _MPI
  MPI_Allreduce (MPI_IN_PLACE, &pmom.omax, nm, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
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

void global_quantities(char * name, vector omega, vector u){

  double ekin = 0.0, circ[3] = {0.0,0.0,0.0}, enst = 0.0, volu_tot=0.0, cents[3] = {0.0,0.0,0.0};
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

  FILE * fp;
  if (pid()==0) {
    fp = fopen(name, "a");
    fprintf (fp, "%.5f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f \n", t, omega_mean.x, omega_mean.y, omega_mean.z, ekin, circ[0], circ[1], circ[2], enst, volu_tot, cents[0], cents[1], cents[2]);
    fclose(fp);
  }
}

void global_quantities2(char * name, vector u, double xc, double yc, double radius, int nj){

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
  FILE * fp;
  if (pid()==0) {
    fp = fopen(name, "a");
    fprintf (fp, "%.5f %d %.15f %.15f %.15f \n", t, nj, area, utang, unorm);
    fclose(fp);
  }
}
