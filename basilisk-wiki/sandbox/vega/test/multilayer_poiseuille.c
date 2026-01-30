/**

We will simulate only the half of a linear poseuille flow because of the limitations of the boundary condition.

More precisely, we will use the following boundary conditions :

$$u_x(z = 0) = 0 ; \partial_z u_x (z = H) = 0$$

$$u_z(z = 0) = u_z (z = H) = 0$$

The main goal of this test case is to check the convergence in both $\Delta x$ and $dt$ of the multilayer solver with the "rigid lid" option.

*/

#include "grid/multigrid1D.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"

double Re = 1.;
double Tmax = 10.;

double H0 = 2.;

/**

We set an option to activate to do a .mp4 animation

*/

bool do_profile = false;

/**

*/

char char_extension[32];
bool correct_DT; //used for timestep debugging

/**

### Convergence of the residual (difference between theoritical and numerical solution)

*/

void plot_convergence() {
  FILE * pp = popen("gnuplot 2> /dev/null", "w");
  fprintf(pp,
          "set grid\n"
          "set logscale x\n"
          "set logscale y\n"
          "files = system(\"ls residual_*\")\n"
          "N = words(files)\n"
          //"p for [i=1:N] word(files,i) u 1:( ($4 == 1) ? $3 : NaN) t substr(word(files,i), 9, 32)\n"
          "p for [i=1:N] word(files,i) u 1:3, x**(-2)\n"
          "set term pngcairo size 1024,720\n"
          "set output 'residual.png'\n"
          "replot\n"
  );
  fflush(pp);
  //sleep(1.);
}

/**

## Run

*/

int main () {
  
/**
Definition of the boundary conditions
*/
  
  cell_lim = mono_limit;
  rigid = true;
  theta_H = 0.51; // to match with the rigid lid option
  
/**

*/
  
  //DT = 0.1;
  TOLERANCE = 1e-5;

  periodic(right);
  
  L0 = 10.;
  nu = 1./Re;
  nl = 16; // number of layer
  N = 16;  // number of points
  for (DT = 0.4; DT >= 0.1; DT /=2.) {
    for (N = 16; N <= 128; N*=2) {
      snprintf(char_extension, sizeof(char_extension), "N%d_DT%g", N, DT);
      run();
    }
  char cmd[128];
  snprintf(cmd, sizeof(cmd), "mv residual residual_DT%g", DT);
  system(cmd);
  }
  plot_convergence();
}

event stop(t=Tmax) {
  return 1.;
}

/**
The constant pressure gradient is implemented as an homogeneous force, which takes the flowing dimensionless form!

$$\mathbf{f_p} = \dfrac{1}{Re} \mathbf{e_x}$$

The prefactor is set in such way that the max of the dimensionless velocity is 1.

*/

event acceleration (i++) {
  foreach_face(x) {
    double z = zb[];
    foreach_layer() {
      z += h[]/2.;
      ha.x[] += h[]*1./Re;
      z += h[]/2.;
    }
  }
}

/**

*/

event init (i=0) {
  correct_DT = true;
  foreach() {
    zb[] = 0.;
    foreach_layer() {
      h[] = H0/nl;
      z += h[]/2.;
      u.x[] = 0.; //the flow is at rest at the begining
      z += h[]/2.;
    }
  }
}

/**

We thus define a function to compute the difference between the numerical solution and the theoretical (laminar) one.

*/

event residual (t=end) {
  double r = 0.;
  foreach() {
  double z = zb[];
  foreach_layer() {
    z += h[]/2.;
    r += fabs(u.x[] - (0.5*z*(2*H0 - z)));
    z += h[]/2.;
    }
  }
  //fprintf(stdout, "%d %g %g\n", N, DT, r/(N*nl));
  FILE * fp = fopen("residual", "a");
  fprintf(fp, "%d %g %g %d\n", N, DT, r/(N*nl), correct_DT);
  fflush(fp);
  sleep(1.);
}

/**

If set so at the begining of the code with `do_profile = true`, we generate a .mp4 movie of the horizontal component of the velocity.

*/

event animate_profile (t = end) {
  sleep(10); // wait for the buffers of profile to be written
  if (do_profile) {
    char cmd[512];
    snprintf(cmd, sizeof(cmd),
             "for f in profile_%s_*.png; do convert $f ppm:- ; rm -f $f; done | "
             "ppm2mp4 profile_movie_%s.mp4", char_extension, char_extension);
    system(cmd);
  }
}

event profile (i++) {
  if (do_profile) {
    FILE * fp = fopen("profile_log", "w");
    foreach(serial) {
      double z = zb[];
      foreach_layer() {
        z += h[]/2.;
        fprintf(fp, "%g %g %g\n", x, z, u.x[]);
        z += h[]/2.;
      }
    }
    fflush(fp);

    static FILE * pp = popen("gnuplot 2> /dev/null", "w");
    fprintf(pp,
            "reset\n"
            "set xrange [0.:%g]\n" //x lim
            "set grid\n"
            "set palette defined(0 \"red\", 1 \"blue\")\n"
            "unset colorbox\n"
            ,H0);

    char field_name[32];
    snprintf(field_name, sizeof(field_name), "ux(z, t = %g)", t);
    fprintf(pp,
            "set title \" %s \\n N = %d\"\n"
            , field_name, N);
    fprintf(pp,
    "set yrange [0:10]\n"
    "plot 'profile_log' u 2:3:1 w points lc palette, 0.5*x*(2*%g - x)\n"
    "set term pngcairo size 1024,720\n"
    "set output 'profile_%d_%06d.png'\n"
    "replot\n",
    H0, N, i);
    fflush(pp);

    system("rm profile_log");
  }
}

event errfile (i++) {
  fprintf(stderr, "%d %g\n", i, t);
  if ((fabs(dt - DT) > 0.1*DT) && (i > 2)) {
    correct_DT = false;
  }
}
