/**
   This is the extension of the [traditional Coriolis
   force](http://basilisk.fr/src/layered/coriolis.h) written by [A. Castillo
   ](http://basilisk.fr/sandbox/acastillo/filaments/circulation/coriolis3d.h)

   g is the sum of pressure gradient and advection.

 */

double f0_x = 0.;
double f0_y = 0.;
double f0_z = 0.;

double k0 = 0.;

double alpha_H = 0.5; // =1 for fully implicit, =0 for explicit

event acceleration (i++){
//  coord b0 = {-K0(),-K0(),-K0()}, b1 = { 2*F0().x, 2*F0().y, 2*F0().z};
  coord b0 = {-k0,-k0,-k0}, b1 = { f0_x, f0_y, f0_z};
  double B[3][3] = {{b0.x, b1.z, -b1.y},{-b1.z, b0.y, b1.x},{b1.y, -b1.x, b0.z}};
  double m[3][3], minv[3][3];
  for (int i = 0; i < 3; i++){
    for (int j = 0; j < 3; j++){
      m[i][j] = -alpha_H * dt * B[i][j];
      if (i == j){
        m[i][j] += 1.0;
      }
    }
  }

  double det = m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2]) - \\
               m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) + \\
               m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);

  double invdet = 1 / det;
  minv[0][0] = (m[1][1] * m[2][2] - m[2][1] * m[1][2]) * invdet;
  minv[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * invdet;
  minv[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * invdet;
  minv[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * invdet;
  minv[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * invdet;
  minv[1][2] = (m[1][0] * m[0][2] - m[0][0] * m[1][2]) * invdet;
  minv[2][0] = (m[1][0] * m[2][1] - m[2][0] * m[1][1]) * invdet;
  minv[2][1] = (m[2][0] * m[0][1] - m[0][0] * m[2][1]) * invdet;
  minv[2][2] = (m[0][0] * m[1][1] - m[1][0] * m[0][1]) * invdet;

  foreach(){
    coord r;
    r.x = u.x[] + (1. - alpha_H)*dt*(B[0][0]*u.x[] + B[0][1]*u.y[] + B[0][2]*u.z[]) + dt*g.x[];
    r.y = u.y[] + (1. - alpha_H)*dt*(B[1][0]*u.x[] + B[1][1]*u.y[] + B[1][2]*u.z[]) + dt*g.y[];
    r.z = u.z[] + (1. - alpha_H)*dt*(B[2][0]*u.x[] + B[2][1]*u.y[] + B[2][2]*u.z[]) + dt*g.z[];

    u.x[] = (minv[0][0]*r.x + minv[0][1]*r.y + minv[0][2]*r.z) - dt*g.x[];
    u.y[] = (minv[1][0]*r.x + minv[1][1]*r.y + minv[1][2]*r.z) - dt*g.y[];
    u.z[] = (minv[2][0]*r.x + minv[2][1]*r.y + minv[2][2]*r.z) - dt*g.z[];
  }
  boundary ((scalar *){u});
}
