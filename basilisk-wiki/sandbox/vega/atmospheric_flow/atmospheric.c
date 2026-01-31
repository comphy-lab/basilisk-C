/**

# Introduction

![A movie of the flow](../kolmogorov_flow/out_movie_Re16000_Fr0.25.mp4)

The first instability (near t = 30) of the forced flow along the x axis creates a tangential flow along the y axis.

And a second instability (between t = 200 and t = 300) creates thin horizontal shearing layers.

From this state, we can expect an other instability to occur looking like a Kelwin-Helmholtz instability.

*/


/**

## Run

*/

#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"

#include "layered/perfs.h"

#define drho(T) (-1.*(T))
#include "layered/dr.h"

/**

The dimensionless numbers describing the flow.

`ratio` = Lz/Lx

Lx = Ly = 2`N_osc` $\pi$

*/

double Re = 1.;
double ReH = 1.; // ReH = ReFr^2
double Fr = 100.;
double Pr = 1.;

double ratio = 0.5;

int N_osc = 2;

/**

The numerical parameters

*/

double epsilon = 1.e-3;
int N_res = 32;

char extension_phy[32];
double H;

double T_max = 100.;

/**

*/

int main () {
  cell_lim = mono_limit;
  
/**
Free slip boundary conditions at top and bottom
*/
  
  lambda_b[] = {HUGE, HUGE};
  theta_H = 0.51;
  rigid = true;
  
/**
Periodic BC in the other directions
*/
  
  L0 = 2.*pi*N_osc;
  H = ratio*L0;
    
  foreach_dimension() {
    periodic(right);
  }
  
/**

*/
  
  DT = 0.01; //should be inferior to approximatively Delta_x^2 Re / 3
  TOLERANCE = 1e-3;
  
  N = N_res;
  nl = N_res;

/**
We set two loops to iterate over the values of Re and the product of Re Fr^2 (= `ReH`)
*/
  
  for (ReH = 1000.; ReH <=1000.; ReH *= 10.) {
    for (Re = 1000.; Re <= 8000.; Re *= 4.) {
      Fr = sqrt(ReH/Re);// ReFr^2 = ReH
      snprintf(extension_phy, sizeof(extension_phy), "Re%g_Fr%g", Re, Fr);
 
      nu = 1./Re;
      G = 1/sq(Fr);
    
      run();
    } //end of Re loop
  }
  fprintf(stderr, "Simulation finished");
}

/**
We initialize the flow with the laminar profile for the velocity and a linear profile for the stratification (= flottability).
*/

event init (i = 0) {
  fprintf(stderr, "Running %s\n", extension_phy);
  foreach() {
    zb[] = 0.;
    double z = 0.;
    foreach_layer() {
      h[] = H/nl;
      z += h[]/2.;
      //initialization
      u.x[] += 1.*cos(y);  //the (theoretical) laminar profile
      T[] += +1.*(z - H/2.); //+1 : stable stratification ; -1 : unstable stratification 
      //
      z += h[]/2.;
    } //end of foreach_layer
  } //end of foreach
}

/**

We want to add a perturbation of the form $\mathbf{u'} = \boldsymbol{\nabla} \land \begin{pmatrix}
\psi_x \\
\psi_y \\
\psi_z \\
\end{pmatrix}$

*/

event perturb (i = 1) {
  
/**
We first define the 3 components of the "stream function"
*/
  
  double regul;
  scalar psi_x = new scalar[nl];
  scalar psi_y = new scalar[nl];
  scalar psi_z = new scalar[nl];
  foreach() {
    double z = zb[];
    for (int l = 0; l<nl; l++) {
      z += h[0,0,l]/2.;
      regul = exp(1/(sq((z-pi)/(1.*pi)) -1))*exp(1.); //replace pi with H/2. ?
      psi_x[0,0,l] = +0.*epsilon*cos(0.5*x)*cos(y)*regul;
      psi_y[0,0,l] = +1.*epsilon*cos(0.5*x)*cos(y)*regul;
      psi_z[0,0,l] = -1.*epsilon*sin(0.5*x)*cos(y);
//#if DR
//      T[0,0,l] += -1.*epsilon*cos(ratio*x)/sq(Fr);
//#endif
      z+= h[0,0,l]/2.;
    }
  }

/**
Thus, we take the rotational.
*/  
  foreach() {
    for (int l = 1; l<nl-1; l++) {
      u.y[0,0,l] += +(psi_x[0,0,l+1] - psi_x[0,0,l-1])/(2.*h[0,0,l]) ;
      w[0,0,l]   += -(psi_x[0,1,l]   - psi_x[0,-1,l])/(2.*Delta);

      u.x[0,0,l] += -(psi_y[0,0,l+1] - psi_y[0,0,l-1])/(2.*h[0,0,l]) ;
      w[0,0,l]   += +(psi_y[1,0,l]   - psi_y[-1,0,l])/(2.*Delta);

      u.x[0,0,l] += +(psi_z[0,1,l]   - psi_z[0,-1,l])/(2.*Delta);
      u.y[0,0,l] += -(psi_z[1,0,l]   - psi_z[-1,0,l])/(2.*Delta);
    }
  }
}

/**

*/

event acceleration (i++) {
  foreach_face(x) {
    foreach_layer() {
      ha.x[] += 1./Re*cos(y);
//      ha.x[] +=  1./Re*cos(x)*sin(y);
//      ha.y[] += -1./Re*sin(x)*cos(y);
    }
  }
}

/**

We add the horizontal diffusion for all fields and the vertical diffusion for the flottability, taking the flux-fixed boundary conditions into account.

A way to rewrite the Robin boundary condition to transform it into a Neumann BC is the following:

$$\partial_z T_b = \dfrac{T_b - T_0}{\lambda_b}$$

So we can make $\lambda_b \longrightarrow + \infty$ keeping the following ratio : $T_0/\lambda_b = -1$

*/

event viscous_term (i++) {
  horizontal_diffusion({u.x, u.y, w}, 1./Re, dt);
  horizontal_diffusion({T}, 1./(Re*Pr), dt);
  foreach()
    vertical_diffusion (point, h, T, dt, 1./(Re*Pr), 1, -1*HUGE, HUGE);
}

/**
## Outputs

We want to animate all the fields (`u.x`, `u.y`, `w` and `T`) on the vertical planes $(x=0)$, $(y=0)$ and on the horizontal plane $(z=0)$.
*/

event moviemaker (t = end) {
  if (true) {
    sleep(10.);
    char cmd[512];
    snprintf(cmd, sizeof(cmd),
             "for f in plot-*.png; do convert $f ppm:- ; rm -f $f; done | "
             "ppm2mp4 out_movie_%s.mp4", extension_phy);//, extension_num);
    system(cmd);
  }
}

/**

Because we use the multilayer solver, plotting on a horizontal plane is (approximately) the same as plotting the field at a fixed layer.

Doing so in gnuplot without formatting the data is very hard, that's why I opted for a solution where I will just plot the field point by point with squares, the color reflecting the instensity of the field.

*/

void plot_horz_field (FILE * fp, scalar f, char * field_name, double f_amp) {
  fprintf(fp, "set title \"%s\"\n", field_name);

  fprintf (fp,
//           "set pm3d map\n"// interpolate 2,2\n"
           "set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25 0 0.5647 1, 0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392, 0.625 1 0.9333 0, 0.75 1 0.4392 0, 0.875 0.9333 0 0, 1 0.498 0 0 )\n"
           "set pointsize 2\n"
           );

  if (f_amp > 0.) {
    fprintf(fp, "set cbrange [-%g:%g]\n",f_amp, f_amp);
  }

  fprintf(fp,
          "set xlabel 'x'\n"
          "set ylabel 'y'\n"
          "set xrange [0:%g]\n"
          "set yrange [0:%g]\n"
          ,L0, L0);
  fprintf(fp, "sp '-' u 1:2:3 w p linetype 5 lc palette\n");

//  int index = 0;
  foreach (serial) {
//    index += 1;
    fprintf(fp, "%g %g %g\n", x, y, f[0,0,0]);
//    fprintf(stderr, "debug - %d %g %g\n", index, x, y);
  }
  fprintf (fp, "e\n\n");
}

/**
On the other hand, plotting a field on a vertical plane can be tough in the multilayer setting. The main idea here was to test for each point if the distance between the point and the plane was inferior to the numerical resolution. That is to say:

$$- \dfrac{\Delta}{2} \leq x - x_0 < \dfrac{\Delta}{2}$$

Although not being a perfect method, it works quite well, especially when the planes are the first ones (or the last).

We use the integer `axis` to chose if we want either plot the plane $(x = 0)$ (`axis = 1`) or the plane $(y = 0)$ (`axis = 2` or any integer different from `1`).

*/

bool test_coor (double coor, double Delta) {
  return (coor < 1.1*Delta);
}

void plot_vert_field (FILE * fp, int axis, scalar f, char * field_name, double f_amp) {
  fprintf(fp, "set title \"%s\"\n", field_name);
  fprintf (fp,
           "set pm3d map\n"// interpolate 2,2\n"
           "set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25 0 0.5647 1, 0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392, 0.625 1 0.9333 0, 0.75 1 0.4392 0, 0.875 0.9333 0 0, 1 0.498 0 0 )\n"
           //"set cbrange [-%g:%g]\n", f_amp, f_amp //comment in order to have an adaptative colorbar
           );

  if (f_amp > 1e-7) {
    fprintf(fp, "set cbrange [-%g:%g]\n",f_amp, f_amp);
//    fprintf(stderr, "deb - %s \n", field_name);
  }

  else {
    fprintf(fp, "unset cbrange\n");
  }

  fprintf(fp,
          "set ylabel 'z'\n"
          "set xrange [0:%g]\n"
          "set yrange [0:%g]\n"
          ,L0, H);
  if (axis == 1) {
    fprintf(fp, "set xlabel 'y'\n");
  }
  else {
    fprintf(fp, "set xlabel 'x'\n");
  }

  fprintf(fp, "sp '-' u 1:2:3\n");
  
  bool is_good;
  double coor;

  foreach (serial) {
    if (axis == 1) {
      is_good = test_coor(x, Delta);
      coor = y;
    }
    else {
      is_good = test_coor(y, Delta);
      coor = x;
    }

    if (is_good) {
      double z = zb[];
      foreach_layer() {
        z += h[]/2.;
        fprintf (fp, "%g %g %g\n", coor, z, f[]);
        z += h[]/2.;
      }
      fprintf (fp, "\n");
    }
  }
  fprintf (fp, "e\n\n");
}

/**

Finally, we call the previous function at a regular timestep to produce the figures.

*/

event gnuplot (t+=1.) {
  static FILE * fp = popen("gnuplot 2> /dev/null", "w");

  bool xy = true;
  bool xz = true;
  bool yz = true;

  int total_index = 3;
  
  // this setup works with gnuplot v4.6 but may bug with v5.2
  fprintf(fp,
            "set term pngcairo size 1600,1400\n"
            "set output 'plot-%06d.png'\n"
            "set multiplot layout 4,%d \n", // 3 because 3 planes
            i, total_index);
  
  if (yz) {
    char field_name[32];
    snprintf(field_name, sizeof(field_name), "u_x(x = 0, t = %g)", t);
    plot_vert_field(fp, 1, u.x, field_name, 1.);
  }
  if (xz) {
    char field_name[32];
    snprintf(field_name, sizeof(field_name), "u_x(y = 0, t = %g)", t);
    plot_vert_field(fp, 2, u.x, field_name, 1.);
  }
  if (xy) {
    char field_name[32];
    snprintf(field_name, sizeof(field_name), "u_x(z = 0, t = %g)", t);
    plot_horz_field(fp, u.x, field_name, 1.);
  }

  if (yz) {
    char field_name[32];
    snprintf(field_name, sizeof(field_name), "u_y(x = 0, t = %g)", t);
    plot_vert_field(fp, 1, u.y, field_name, 1.);
  }
  if (xz) {
    char field_name[32];
    snprintf(field_name, sizeof(field_name), "u_y(y = 0, t = %g)", t);
    plot_vert_field(fp, 2, u.y, field_name, 1.);
  }
  if (xy) {
    char field_name[32];
    snprintf(field_name, sizeof(field_name), "u_y(z = 0, t = %g)", t);
    plot_horz_field(fp, u.y, field_name, 1.);
  }

  if (yz) {
    char field_name[32];
    snprintf(field_name, sizeof(field_name), "u_z(x = 0, t = %g)", t);
    plot_vert_field(fp, 1, w, field_name, 1.);
  }
  if (xz) {
    char field_name[32];
    snprintf(field_name, sizeof(field_name), "u_z(y = 0, t = %g)", t);
    plot_vert_field(fp, 2, w, field_name, 1.);
  }
  if (xy) {
    char field_name[32];
    snprintf(field_name, sizeof(field_name), "u_z(z = 0, t = %g)", t);
    plot_horz_field(fp, w, field_name, 1.);
  }

  if (yz) {
    char field_name[32];
    snprintf(field_name, sizeof(field_name), "T(x = 0, t = %g)", t);
    plot_vert_field(fp, 1, T, field_name, pi);
  }
  if (xz) {
    char field_name[32];
    snprintf(field_name, sizeof(field_name), "T(y = 0, t = %g)", t);
    plot_vert_field(fp, 2, T, field_name, pi);
  }
  if (xy) {
    char field_name[32];
    snprintf(field_name, sizeof(field_name), "T(z = 0, t = %g)", t);
    plot_horz_field(fp, T, field_name, pi);
  }


  fprintf (fp,
             "unset multiplot\n"
             "set output\n"
             ); 
}

/**
### Utils

Finally, we define some global outputs.
*/

double norm_L1 (scalar field, double slope) {
  double abssum = 0.;
  foreach(reduction(+:abssum)) {
    double z = zb[];
    foreach_layer(){
      z += h[]/2.;
      abssum += fabs(field[] - slope*(z - H/2.))*h[]*sq(Delta); //maybe *ratio / N_osc
      z += h[]/2.;
    }
  }
  return abssum;
}

double norm_L2 (scalar field, double slope) {
  double var = 0.;
  foreach(reduction(+:var)) {
    double z = zb[];
    foreach_layer(){
      z += h[]/2.;
      var += sq(field[] - slope*(z - H/2.))*h[]*sq(Delta); //maybe *ratio / N_osc
      z += h[]/2.;
    }
  }
  return var;
}

double delta (scalar field) {
  double max = -1.*HUGE;
  double min = HUGE;
  foreach(reduction(max:max) reduction(min:min)) {
    foreach_layer() {
//      fprintf(stderr, "deb - %g %g\n", l_max, field[0,0,0]);
      if (field[0,0,0] > max) {
        max = field[0,0,0];
      }
      if (field[0,0,0] < min) {
        min = field[0,0,0];
      }
    }
  }
  return max - min;
}

event snapshot (t = end) {
  static FILE * fp = fopen("snapshot", "w");
  foreach (serial) {
    double z = zb[];
    foreach_layer() {
      z += h[]/2.;
      fprintf(fp, "%g %g %g %g %g %g\n",x, y, z, u.x[], u.y[], w[]);
      z += h[]/2.;
    }
  }
}

/**
In a datafile, we output the variance of the vertical velocity and the difference between the min and the max of the flotability.
*/

event datafile (i++) {
  static FILE * fp = fopen("data", "w");
//  double Strat = norm_L1(T, 0.);
  double uz2 = norm_L2(w, 0.);
  double delta_T = delta(T);
  fprintf(fp, "%d %g %g %g\n", i, t, uz2, delta_T);
  fprintf(stdout, "deb - %g %d\n", t, i);
}
/**

*/

event stop (t = T_max) {
  return 1.;
}
