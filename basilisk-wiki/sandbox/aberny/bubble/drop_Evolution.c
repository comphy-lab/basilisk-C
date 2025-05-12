/** With this small code, we want to obtain the evolution of a drop
position when knowing its initial position and velocity.

The idea behing this code is to compare a theoretical evolution of a drop with
the one computed by basilisk*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/**
We will solve the following differential equation:

$$ m_{drop}\frac{\partial v}{\partial t} + \frac{1}{2}\rho_{air}Sv^2C_D = 0 $$

We choose the following model for the viscous drag coefficient:

$$ c_D = \frac{24}{Re}\left(1+0.15\text{Re}^{0.687}\right)$$

The droplet masse is obtain by it's radius with the following formula:

$$ m_{drop} = \frac{4}{3}\pi R_{drop}^3 \rho_{water} $$

Then, by replacing the term, and after simplification, the differential equation
we need to solve is:

$$ \frac{\partial v}{\partial t} = -\frac{9 \mu_{air}v}{R_{drop}^2}\left(1+0.15\left(\frac{\rho R_{drop} v}{\mu}\right)^{0.687}\right)$$

For that we will use a simple euler algorithm.

The function f is computing the value of the derivation of $v$, with respect to $t$
*/

double f (double v, double mu, double R) {
  double term1 = -0.15*v*pow( ( (1./998.)*R*v/mu), (0.687) )*(9.*mu/(R*R));
  // double term1 = 0;
  double term2 = -(9.*mu/(R*R))*v;
  return (term1 + term2);
}

/** The function vn1 is computing the velocity $v$ at time-step $t + dt$, when
knowing the velocity at $t$.

The computation is the following:

$$ v(t+dt) = v(t) + \text{dt}f\left(v\left(t\right),t\right)$$

*/

double vn1 (double vn, double h, double mu, double R) {
  return (vn + h*f(vn, mu, R));
}

/**
To obtain the position of the drop, we just perform a numerical integration of the velocity.*/

double xn1 (double xn, double h, double vn) {
  return (xn + h*vn);
}
/**
Initialisation of the external parameters.

 * La is the Laplace number which give us the viscosity of the drop.
 * x0 is the initial position of the drop
 * v0 is the initial velocity of the drop
 * R is the radius of the drop
 * t0 is the initial time
*/

double La = 5000;
double x0 = 0;
double v0 = 0;
double R = 1;
double t0 = 0;

double h = 0.0001;

int main(int argc, char const *argv[]) {

  /**
  All the external parameters can be change*/
  
  if (argc >= 2){
    La = atof (argv[1]);
  }

  if (argc >= 3){
    x0 = atof (argv[2]);
  }

  if (argc >=4) {
    v0 = atof (argv[3]);
  }

  if (argc >=5) {
    R = atof (argv[4]);
  }

  if (argc >=6) {
    t0 = atof (argv[5]);
  }

  /**
  The air viscosity is define by the following equation:

$$ \mu_{air} = \frac{1}{\mu_{ratio}\sqrt{\text{La}}}$$
  
  where $\mu_{ratio}$ is the ratio between the liquid viscosity and the air
  viscosity. Here this value is equal to 55*/

  double mu = 1./(55*sqrt(La));

  /**
  Declaration of the table for x and t*/

  int size = 1/h;

  double x[size];
  double t[size];
  double v[size];
  
  /**
  Initialisation of x and t*/

  int i = 0;

  x[i] = x0;
  t[i] = t0;
  v[i] = v0;

  fprintf(stdout, "%f %f %f\n", t[i], x[i], v[i]);

  /**
  Resolution*/

  while (t[i]<=1. && x[i]<=5.) {
    i++;
    v[i] = vn1(v[i-1], h, mu, R);
    x[i] = xn1(x[i-1], h, v[i-1]);
    t[i] = t[i-1] + h;

    /**     
    Once the equation is solved, we output the position and the velocity
    of the droplet*/     

    fprintf(stdout, "%f %f %f\n", t[i], x[i], v[i]);   }

  return 0;
}
