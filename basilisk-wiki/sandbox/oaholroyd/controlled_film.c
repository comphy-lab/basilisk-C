/**
# Falling liquid film with proportional control

A liquid film, governed by the 2D Navier-Stokes equations, falls down an
inclined plane and is controlled towards the uniform film (or Nusselt) solution.
The control consists of a finite number of proportional actuators ($M$), which
inject or remove fluid at the base.
*/
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "embed.h" // embedded boundaries (supersedes mask)
#include "tension.h" // surface tension
#include "heights.h" // interfacial height

#define LEVEL 7
#define N (1 << (LEVEL))
#define DX (LX/N)

/**
The dynamics of the film are governed by a number of dimensionless parameters
*/
// domain parameters
const double LX = 30.0; // domain length
const double LY = 8; // domain height
const double THETA = 1.0; // inclination angle

// fluid parameters
const double RE = 15.0; // Reynolds number
const double CA = 0.02; // capillary number
const double RHO_RATIO = 1000.0; // density ratio
const double MU_RATIO = 1000.0; // viscosity ratio

// control parameters
const double W = 0.1; // control width
#define M 5 // number of controls
#define ALPHA 1.0 // control strength
#define ETA -1.0 // control shift

#define TSTART 100.0 // time controls turn on
#define TMAX 250.0 // final time




/**
## Interfacial measurement

Define a function to measure the interfacial height $h$ at a given x-coordinate.
 */
#define NP 10
vector hei[]; // heights
double interfacial_height(double xp) {
  if (xp < 0) {
    xp += LX;
  } else if (xp > LX) {
    xp -= LX;
  }

  double dh[NP];
  double yh[NP];
  double yp;
  double y0 = 0.0, y1 = 2.0;

  /* try and find range of possible heights */
  for (int i = 0; i < NP; i++) {
    yp = y0 + i*(y1-y0)/(NP-1);
    Point point = locate(xp,yp);

    if (hei.y[] != nodata) {
      yh[i] = y + height(hei.y[])*Delta;
      dh[i] = fabs(y-yh[i]);
      return yh[i];
    } else {
      yh[i] = -1000;
      dh[i] = 1000;
    }
  } // i end

  /* find the closest one */
  int j = 0;
  for (int i = 1; i < NP; i++) {
    if (dh[i] < dh[j]) { j = i; }
  } // i end

  return yh[j];
}

/**
Utilities to allow conversion between grid index and x-coordinate, and actuator
index and x-coordinate.
 */
#define ITOX(i) (DX*(i+0.5))
#define ILOC(x) ((0.5+i) * (LX/M) + (ETA))


/**
## Boundary conditions and control mechanism

At the lower boundary we use no-slip alongside the applied controls. The
controls consist of the sum of M actuators:
$$
f(x) = \sum_{i=1}^M A_i B \exp\left[\frac{\cos(2\pi x / L_x)-1}{\omega^2}\right],
$$
where $A_i$ is the actuator strength and $B$ a normalising constant. The actuator
strength is governed by the proportional relationship
$$
A_i = -\alpha (h(x_i + \eta) - 1),
$$
where $x_i$ is the location of the $i$th actuator and $\eta$ is an up- or
downstream shift.
 */
double NORM; // control normaliser
double actuator(double x) {
  if (x < -LX/2) {
    x += LX;
  } else if (x > LX/2) {
    x -= LX;
  }

  return NORM*exp((cos(2*M_PI*x/LX)-1.0)/(W*W));
}

double A[M]; // actuator strengths
double control(double xp) {
  double f = 0.0;

  for (int i = 0; i < M; i++) {
    f += A[i]*actuator(xp - ILOC(i));
  }

  return f;
}

u.n[bottom] = dirichlet(control(x));
u.t[bottom] = dirichlet(0.0);

/**
We use [embed](http://basilisk.fr/src/embed.h) for the top boundary, applying a
free-outflow condition.
 */
u.n[embed] = neumann(0.0);
p[embed] = neumann(0.0);
pf[embed] = neumann(0.0);


face vector av[]; // acceleration vector field
double G[2]; // gravity
int main(int argc, char const *argv[]) {
  periodic(right);

  init_grid(1 << (LEVEL));
  L0 = LX;
  X0 = 0.0;
  Y0 = 0.0;

  TOLERANCE = 1e-4;
  DT = 5e-2;

  rho1 = RE;
  rho2 = rho1 / RHO_RATIO;
  mu1 = 1.0;
  mu2 = mu1 / MU_RATIO;
  f.sigma = 1.0/CA;

  a = av;
  G[0] = 2.0/RE;
  G[1] = -2.0/tan(THETA)/RE;

  /* control normaliser */
  NORM = 1.0;
  double integral = 0.0;
  for (int i = 0; i < N; i++) {
    integral += actuator(DX*i - LX/2);
  } // i end
  NORM = 1.0/(DX*integral);

  run();
}


/**
## Initial condition

The film height is initialised with a small perturbation to trigger an
instability, and the pressure and velocity fields start at values given by the
Nusselt solution.
 */
event init(i=0) {
  solid(cs, fs, LY + y);

  fraction(f, 1.0-y+0.05*sin(1.0*(2.0/(LX))*M_PI*(x+10)));

  foreach () {
    u.x[] = f[]*y*(2.0-y) + (1.0-f[]);
    u.y[] = 0.0;
    p[] = f[]*2*cos(THETA)/sin(THETA)*(1-y);
  }

  boundary({u,f,p});
  heights(f,hei);
}


/**
## Static grid refinement

Static grid refinement makes it simpler to sample the interface while also avoiding
wasted computational effort in the largely inconsequential gas layer.
 */
event refinement(i=0) {
  refine(level < LEVEL);
  unrefine(y > 3 && level > LEVEL-2);
  unrefine(y > 5 && level > LEVEL-3);
  unrefine(y > 7 && level > LEVEL-4);
}


/**
We impose acceleration due to gravity in the fluid.
 */
event acceleration(i++) {
  foreach_face (x) {
    av.x[] += f[]*G[0];
  }
  foreach_face (y) {
    av.y[] += f[]*G[1];
  }
}


/**
## Controls

At each time step we measure the interface at the site of each actuator and
scale the actuator output proportionally. The controls are not switched on
immediately to allow a travelling wave to develop.
 */
event controls(i++) {
  heights(f,hei);

  for (int i = 0; i < M; i++) {
    if (t > TSTART) {
      A[i] = - ALPHA * (interfacial_height(ILOC(i))-1.0);
    } else {
      A[i] = 0.0;
    }

  } // i end

  boundary({u.x, u.y});
}


/**
## Output and Results

At each time step we plot the interface and the controls with gnuplot, before
converting this into a movie at the end with ffmpeg.

The film is successfully controlled:

![](controlled_film/mov.mp4)

with perturbations damped exponentially.

![](controlled_film/lines.png)
 */
FILE *gp_pipe;
event init (t = 0) {
  gp_pipe = popen("gnuplot", "w");
  fprintf(gp_pipe,
        "set term pngcairo\n"
        "set xr [0: %.1lf]\n"
        "set yr [%.1lf: %.1lf]\n"
        "set key bottom left\n"
        "set grid\n"
        "set xlabel 'x'\n"
        "set ylabel 'h, f'\n", LX, -2.0, 2.0);

  FILE *fp = fopen("out.dat", "a");
  fprintf(fp, "t h ");
  for (int i = 0; i < M; i++) {
    fprintf(fp, "A_%d ", i);
  } // i end
  fprintf(fp, "\n");
  fclose(fp);
}

int frame = 0;
event output (t += 1) {
  fprintf(gp_pipe, "set output 'plot%d.png'\n", frame);
  fprintf(gp_pipe, "set title 'Interfacial height (t = %.1lf)'\n", t-TSTART);
  fprintf(gp_pipe, "plot \
          '-' w l lw 3 lt rgb \"blue\" t 'h', \
          '-' w l lw 3 lt rgb \"red\" t 'f'\n");
  for (int i = 0; i < N; i++) {
    fprintf(gp_pipe, "%lf %lf \n", ITOX(i), interfacial_height(ITOX(i)));
  } // i end
  fprintf(gp_pipe, "e\n");
  for (int i = 0; i < N; i++) {
    fprintf(gp_pipe, "%lf %lf \n", ITOX(i), control(ITOX(i)));
  } // i end
  fprintf(gp_pipe, "e\n");
  frame++;

  FILE *fp = fopen("out.dat", "a");
  double h2 = 0.0;
  for (int i = 0; i < N; i++) {
    double h = interfacial_height(ITOX(i))-1.0;
    h2 += h*h;
  } // i end

  fprintf(fp, "%lf %lf ", t-TSTART, DX*sqrt(h2));
  for (int i = 0; i < M; i++) {
    fprintf(fp, "%lf ", fabs(A[i]));
  } // i end
  fprintf(fp, "\n");
  fclose(fp);
}

event stop (t = TMAX) {
  system("rm mov.mp4");
  system("ffmpeg -r 25 -f image2 -i plot%d.png -c:v libx264 -vf format=yuv420p -y mov.mp4");
  system("rm plot*");

  fprintf(gp_pipe,
        "set term pngcairo\n"
        "set xr [%lf: %lf]\n"
        "set yr [%lf: %lf]\n"
        "set logscale y 10\n"
        "set key bottom left\n"
        "set title 'Interfacial height/Actuator strength'\n"
        "set grid\n"
        "set xlabel 't'\n"
        "set ylabel 'h, a'\n", -TSTART, TMAX-TSTART, 0.001, 1.0);
  fprintf(gp_pipe, "set output 'lines.png'\n");
  fprintf(gp_pipe, "plot for [col=2:%d] 'out.dat' u 1:col w l t columnheader", M+2);
  // system("rm out.dat");
  return 1;
}
