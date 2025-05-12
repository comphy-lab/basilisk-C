/** Inspired by **[WENO scheme in rajarshi's sandbox](http://basilisk.fr/sandbox/rajarshi/WENO_CODES/Saint_Venant_Test/)**
*/

/** Note that wrong included files would lead to weird results*/
#include "grid/multigrid1D.h"
// #include "green-naghdi.h"
// #include "grid/bitree.h"
#include "saint-venant_modified.h"

#define MAXLEVEL 8
#define MINLEVEL 8
#define MAXMAXLEVEL 13

// problem-sepcific parameters
double So = 0.05011;
double normalDepth = 0.00798;
double normalVelocity = 1.0377;
double gravityCoeff = 9.81;
double disMag = 0.05;
double disPeriod = 0.933;
double simTime = 120.0;
double cf; // = gravityCoeff * So * 2. * normalDepth / (sq(normalVelocity));
double Qin;
double froudeNum=3.71;

/**
Reproduction of Brock (1967)'s experiment. */

int main()
{
     N = 2048;
     L0 = 40.;
     G = gravityCoeff;
     CFL = 0.30; // CFL number should be sufficiently small
//      theta = 2.0;
//      gradient = NULL;
     //nu = 5;
     //nl = 8;
     run();
}

/**
We use Dirichlet (supercritical) inlet and zero-order extrapolation outlet.*/

// h[left] = dirichlet(normalDepth + disMag * normalDepth * (t < disPeriod ? 1. : 0.));


event init(i = 0)
{
    Qin = normalDepth*normalVelocity;
     // for 2-D case
     // mask (y > 0.05 ? top : none);
//     h[left] = dirichlet(normalDepth + disMag * normalDepth * (t < disPeriod ? 1. : 0.));
    h[left] = dirichlet(normalDepth + disMag * normalDepth * sin(2. * pi * t / disPeriod));
//     h[left] = dirichlet((normalDepth + 1.0 * normalDepth));
     u.n[left] = dirichlet(froudeNum*(sqrt(gravityCoeff*(normalDepth + disMag * normalDepth * sin(2. * pi * t / disPeriod)))));

     // u.n[right] = neumann(0.);
     // // h[right] = neumann(0.);

     // u.n[right] = radiation(normalVelocity);
     // h[right] = radiation(normalDepth);
     u.n[right] = neumann(0.);
     h[right] = neumann(0.);
     
     cf = gravityCoeff * So * 2. * normalDepth / (sq(normalVelocity));

     foreach ()
     {
          // zb[] = -So * x;
          zb[] = 0.;
          h[] = normalDepth;
          u.x[] = normalVelocity;
     }
}

static double chezyBedFriction(double u, double h, double cf, double g, double So)
{
     double rhs;
     rhs = -(cf / 2.) * u * fabs(u) / h + g*So;
     return rhs;
}

// Quadratic bottom friction
event friction(i++)
{
//     double uMed;
     
     foreach ()
     {
          double a = h[] < dry ? HUGE : 1. + (cf / (2.)) * dt * norm(u) / h[];
          // double a = 1. + (cf / (2.)) * dt * u.x[] / h[];
          foreach_dimension(){
//                   u.x[] /= a;
              u.x[] = (u.x[] + G * So * dt) / (1. + (cf / (2.)) * dt * norm(u) / (h[]));
          
//           uMed = u.x[] + dt * chezyBedFriction(u.x[], h[], cf, G, So);
//           uMed = (3. / 4.) * u.x[] + (1. / 4.) * uMed + (1. / 4.) * dt * chezyBedFriction(uMed, h[], cf, G, So);
//           u.x[] = (1. / 3.) * u.x[] + (2. / 3.) * uMed + (2. / 3.) * dt * chezyBedFriction(uMed, h[], cf, G, So);
              
        }
     }
     boundary((scalar *){u.x}); // note that the input should be a list (at least for 1d)
     
     fprintf (stderr, "%g\n", t);
}

/**
We use gnuplot to visualise the wave profile as the simulation
runs and to generate a snapshot at $t=40$. 

![Snapshot of waves. The top of the bar is seen in white.](bar/snapshot.png)
*/

void plot_profile(double t, FILE *fp)
{
     fprintf(fp,
             "set term pngcairo enhanced size 800,600 font \",10\"\n"
             "set output 't%.0f.png'\n"
             "set title 't = %.2f'\n"
             "set xlabel 'x(m)'\n"
             "set ylabel 'eta(m)'\n"
             "plot [0:][0:]'-' u 1:2 w l lw 2\n",
             t, t);
     foreach ()
     {
          fprintf(fp, "%g %g\n", x, h[]);
     }
     fprintf(fp, "e\n\n");
     fflush(fp);
}

event gnuplot(t = 0; t <= simTime; t += 2)
{
     static FILE *fp = popen("gnuplot 2> /dev/null", "w");
     plot_profile(t, fp);
     // fprintf(fp,
     //         "set term pngcairo enhanced size 800,600 font \",10\"\n"
     //         "set output 't%.0f.png'\n"
     //         "set title 't = %.2f'\n"
     //         "set xrange [0:40]\n"
     //         "plot u 1:2 w l t\n",
     //         t, t);
     // fprintf(fp, "\n");
     // foreach ()
     //      fprintf(fp, "%g %g\n", x, h[]);
     // fprintf(fp, "e\n\n");
     // fflush(fp);
     // fprintf(stderr, "%.3f %.3f\n", t, statsf(h).max); // uncomment if needed
}

event output(t = 0; t <= simTime; t += 0.50)
{
     char name[80];
     sprintf(name, "out-%.1f", t);
     FILE *fp = fopen(name, "w");
     foreach ()
          fprintf(fp, "%g %g %g \n", x, h[], u.x[]);
     fprintf(fp, "\n");
     fclose(fp);
     fprintf(stderr, "%g\n", dt);
}

