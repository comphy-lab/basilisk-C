/**

# Introduction

# The Code

## How to execute?

## Set Up

*/

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "ppm2mp4"
#include "utils.h"
#include "embed.h"
#include "view.h"

#define dimension 2

double Re = 1.;

double Tmax = 40;

int NX = 10;
int NY = 1;
int N_osc = 1;
double ratio; // aspect ratio

double epsilon = 0.001;

bool do_profile = false;
bool do_movie = false;

int N_cells = 128;

char simulation_type[16];
char extension_1[64];
char extension_2[64];
char extension_3[64];
face vector av[];
  /**
  


*/
void gnu_title (FILE * fp, char * title) {
  char gnu_command[512];
  snprintf(gnu_command, sizeof(gnu_command),
           "set title \" %s %s \\nRe = %g, {/Symbol a} = %g\\n {/Symbol e} = %g, n_res = %d \"\n",
           title, simulation_type, Re, ratio, epsilon, N_cells);
  fprintf(fp, gnu_command);
}

  /**
  


*/

void plot_stability() {
  FILE * fp = popen("gnuplot 2> /dev/null", "w");
  fprintf(fp,
          "set ylabel '{/Symbol a}'\n"
          "set xlabel 'Re'\n"
          "set grid\n"
          "set title \" Stability Diagram\\nn_res = %d, DT = %g\"\n",
          N_cells, DT
  );

  fprintf(fp,
          "files = system(\"ls sigma_*_%s\")\n"
          "N = words(files)\n"
          "set palette defined (0 \"blue\", 1 \"red\")\n"
          "set cbrange [0:1]\n"
          "set pointsize 2\n"
          "unset colorbox\n"
          "p for [i=1:N] word(files,i) u 2:( ($5 < 0.01) ? $1 : NaN):(($3 > 0) ? 1 : 0) w p linetype 7 lc palette t '', '../Rec' u 1:2 w l t 'Theory Periodic - G.M.', '../Rec_freeslip' u 1:2 w l t 'Theory FreeSlip - G.M.'\n"
          ,extension_3
  );
  fprintf(fp,
          "set term pngcairo size 1024,720\n"
          "set output 'ReCRatio_%s.png'\n"
          "replot\n"
          , extension_3);
  fflush(fp);

}

  /**
  


*/

void plot_sigma_all (FILE * fp) {
  gnu_title(fp, "");
  fprintf(fp,
          "set xlabel 'Re'\n"
          "set ylabel '{/Symbol s}'\n"
          "set grid\n"
         );
  fprintf(fp,
          "files = system(\"ls sigma_*_%s\")\n"
          "N = words(files)\n"
          "p for [i=1:N] word(files,i) u 2:3:4 w errorbars\n"
         , extension_3);
  fprintf(fp,
          "set term pngcairo size 1024,720\n"
          "set output 'sigma_out_%s.png'\n"
          "replot\n"
          , extension_3);
  fflush(fp);
}

  /**
  


*/

int main () {
  X0 = 0.;//-L0/2.;
  TOLERANCE = 1e-5; 

  DT = 0.1; // need to be non-dimensional in order to work with non-dimensional system
  CFL = 0.5;


  a = av;
  
  for (N_cells = 32; N_cells <= 256; N_cells*=2) {
    N = N_cells;
    for (DT = 0.01; DT >= 0.0001; DT/=10.) {

      snprintf(extension_3, sizeof(extension_3), "DT%g_N%d", DT, N_cells);

      for (NY = 7; NY<= 7; NY+=1) {
        ratio = 1.*NY/NX;
        L0 = 2.*pi/ratio;//length of the box

        dimensions(nx=NX, ny=NY);
        if (false) { // full periodic
          foreach_dimension() {
            snprintf(simulation_type, sizeof(simulation_type), "periodic");
            periodic(right);
          } //end of foreach
        } //endif
        else {
          snprintf(simulation_type, sizeof(simulation_type), "freeslip");
          periodic(right);
          u.t[top] = neumann(0.);
          u.t[bottom] = neumann(0.);
          u.n[top] = dirichlet(0.);
          u.n[bottom] = dirichlet(0.);
        } //endelse

        snprintf(extension_2, sizeof(extension_2), "alpha%g", ratio);

        for (Re = 5.; Re <= 5.; Re += 1.) {
          mu[] = {1./Re,1./Re};
  
          snprintf(extension_1, sizeof(extension_1), "Re%g_alpha%g", Re, ratio);
      
          run();
        } //end of Re loop
      } //end of ratio loop
    } //end of DT loop
  } //end of N loop
  FILE * fp = popen("gnuplot 2> /dev/null", "w");
  plot_sigma_all(fp);
  plot_stability();
}

event init (i = 0) {
  foreach() {
    u.x[] = 1.*cos(N_osc*y);
  }
  fprintf(stderr, "running %s %s\n", extension_1, extension_3);
}

event noise (t=10) {
  scalar psi[];
  double regul;
  foreach() {
    regul = exp(1/(sq((y-pi)/(1.*pi)) -1))*exp(1.);
    psi[] = epsilon*cos(ratio*x)*regul;
  }
  
  foreach() {
    u.x[] += +(psi[0,1] - psi[0,-1])/(2.*Delta);
    u.y[] += -(psi[1,0] - psi[-1,0])/(2.*Delta);
  }
}

event acceleration (i++) {
  foreach_face(x) {
    av.x[] += 1./Re*cos(N_osc*y);
  }
}

  /**
  


*/

//// Utilities ////



//// OUTPUT ////

event animate_profile (t = end) {
  if (do_profile) {
    sleep(10); // wait for the buffers of profile to be written
    char cmd[512];
    snprintf(cmd, sizeof(cmd),
             "for f in profile_%s_%s*.png; do convert $f ppm:- ; rm -f $f; done | "
             "ppm2mp4 profile_movie_%s_%s.mp4", extension_1, extension_3, extension_1, extension_3);
    system(cmd);
  }
}


event profile (i++) {
  if (do_profile) {
    FILE * fp = fopen("profile_log", "w");
    foreach(serial) {
        fprintf(fp, "%g %g %g\n", x, y, u.x[]);
    }
    fflush(fp);

    char cmd[256];
    snprintf(cmd, sizeof(cmd), "mv profile_log profile_log_%s_%s", extension_1, extension_3);
    system(cmd);
    static FILE * pp = popen("gnuplot 2> /dev/null", "w");
    fprintf(pp,
            "reset\n"
            "set xrange [0.:%g]\n" //x lim
            "set grid\n"
            "set palette defined(0 \"red\", 1 \"blue\")\n"
            "unset colorbox\n"
            ,2*pi
            );
    char field_name[32];
    snprintf(field_name, sizeof(field_name), "ux(z, t = %g)", t);
    gnu_title(pp, field_name);
    fprintf(pp,
//    "set yrange [-1.2:1.2]\n"
    "plot 'profile_log_%s_%s' u 2:3:1 w points lc palette, cos(x)\n"
    "set term pngcairo size 1024,720\n"
    "set output 'profile_%s_%s_%06d.png'\n"
    "replot\n",
    extension_1, extension_3, extension_1, extension_3, i);
    fflush(pp);
    fprintf(stderr, "Written %d\n", i);
  }
}

void setup_view () {
  view(tx=-1., sx = 2., ty=-0.5, sy = 2.
      //, width = 200.*NX, height = 200.*NY
  );
}

event plot_field (t+=1.) {
  if (do_movie) {
  setup_view();
  squares ("u.x", spread = -1, min = -1, max = 1, linear = false);
  char name1[256];
  snprintf(name1, sizeof(name1), "ux_%s_%s.mp4", extension_1, extension_3);
  save(name1);

  setup_view();
  squares ("u.y", spread = -1,
           // min = -1, max = 1,
            linear = false);
  char name2[256];
  snprintf(name2, sizeof(name2), "uy_%s_%s.mp4", extension_1, extension_3);
  save(name2);

  }
}

double FT_1D_mode (scalar field, double k) {
  scalar Ek[];
  scalar Ek_cos[];
  scalar Ek_sin[];
  foreach() {
    Ek[] = sq(field[]*cos(k*x)) + sq(field[]*sin(k*x));
    Ek_cos[] = field[]*cos(k*x);
    Ek_sin[] = field[]*sin(k*x);
  }
  stats s_Ek = statsf(Ek);
  stats s_field = statsf(field);
  return sqrt(s_Ek.sum);
}

event datafile (i++) {
  static FILE * fp = fopen("data", "w");
  double Ek = FT_1D_mode(u.y, ratio);
  fprintf(fp, "%d %g %g\n", i, t, Ek);
  fflush(fp); //write the buffer immediately, so that growth will always have the complete data set
  //if (Ek > 10.) { // stop condition
    //fprintf(stderr, "Non linear state reached\n");
    //return 1.;
  //}
}

event growth (t = Tmax) {
  
  char cmd[256];
  snprintf(cmd, sizeof(cmd), "mv data data_%s_%s", extension_1, extension_3);
  system(cmd);

  static FILE * fp = popen("gnuplot 2> /dev/null", "w");
  gnu_title(fp, "");

  fprintf(fp,
          "set fit errorvariables\n"
          "FIT_CONVERGED = 0\n"
          //"a = 0\n"
          //"b = 0.0001\n"
          "f(x) = a + b*x\n"
          "fit [15:][-7*log(10):-2*log(10)] f(x) 'data_%s_%s' u 2:(log($3)) via a,b\n"
          "if (!FIT_CONVERGED) {a = NaN; b = NaN; b_err = NaN; FIT_STDFIT = NaN}\n"
          "set print 'sigma_%s_%s' append\n"
          "print %g, %g, b, b_err, FIT_STDFIT, %d, %g\n"
          "unset print\n"
          , extension_1, extension_3, extension_2, extension_3, ratio, Re, N, DT
         );
  fprintf(fp,
          "set logscale y\n"
          "set xlabel 't [L_0/U_lam]'\n"
          "set ylabel 'abs(ux_k)'\n"
          "set grid\n"
          "plot 'data_%s_%s' u 2:3 w p, exp(a+b*x) t 'linear fit'\n"
          "set term pngcairo size 1024,720\n"
          "set output 'amp_growth_%s_%s.png'\n"
          "replot\n"
          , extension_1, extension_3, extension_1, extension_3);
  fflush(fp);
}

event snapshot(t = end) {
  static FILE * fp = fopen("snapshot", "w");
  foreach(serial) {
    fprintf(fp, "%g %g %g %g\n", x, y, u.x[], u.y[]);
  }
}

event logfile (i+=10) {
  fprintf(stderr, "%d %g\n", i, t);
}
