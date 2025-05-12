/** 
With this small code, we want to obtain the evolution of a drop when knowing
its initial size and velocity. This code is taking into account the
evaporation of the droplets.

Concerning the evolution of the velocity, we have to take into account the drag
of the air arround the droplet. We also take into account the gravity. Then,
the general form of the equation driving the velocity of the droplet is:

$$ 
m_{droplet}\frac{\partial v}{\partial t} = -\frac{1}{2}\rho_{air}S_{droplet}C_xv-m_{droplet}g
$$

We choose the following model for the drag coefficient: 

$$
C_x = \frac{24}{Re}\left(1+0.15Re^{0.687}\right)
$$

Concerning the radius of the droplet, it will change du to the evaporation of
the droplet during its evolution. The evaporation should compensate the amount
of water in the air that diffuse away from the droplet to the rest of the air.

We have: $$ \text{div}\left(\rho_{vap}v\right) = -\text{div}\left(j_{vap}\right)$$

We define the dimensionless number $Y_{vap} = \rho_{vap}/\rho_{air}$

We applie the Fick's law and the divergence theorem. We assume that the density
of the gaz is a constant. We then obtain:

$$
4\pi r^2_{droplet}\rho_{air}vY_{vap}-4\pi r^2D\rho\frac{dY_{vap}}{dr} = cst
$$

where $D$ is the duffision coefficient of water in the air.

The constant in this equation correspond to the mass variation $\dot{m}$. This
is explain by the fact that the divergence theorem is apply on a flux. This
flux is the variation of mass coming from the droplet. Assuming that we are
working with incompressible gazes, we can express the mass variation as:
$\dot{m} = 4\pi r^2\rho_{gaz}v$

Then the diffusion equation looks like:

$$
\dot{m} = \dot{m}Y_{vap}-4\pi r^2D\rho_{gaz}\frac{dY_{vap}}{dr}
$$

This droplet is moving. We assume that the diffusive boundary layer is a
perfect circle arround the droplet. Let $\delta M$ being the thickness of this
layer. We integrate the equation between $r_{droplet}$ and $r_{droplet}  +\delta
M$, and we get:

$$
\dot{m} = 4\pi D \rho_{gaz}r_{droplet}\left(\frac{r_{droplet}}{\delta M} +1\right)\text{ln}\left(B_M+1\right)
$$

where $B_M$ is the Spalding number. $B_m = \frac{Y_{vap}^{surf}-Y_{vap}^{\infty}}{1-Y_{vap}^{surf}}$

We know need to replace $\delta M$ by the actual size of the boundary layer. We choose the following model for the ratio $r_{droplet}/\delta M = 0.3\text{Sc}^{1/3}\text{Re}^{1/2}/F(B_M)$

By replacing everithing, we finally get the equation driving the evolution of
the radius of the droplet:

$$ \frac{dr_{droplet}}{dt} = -\frac{1}{r_{droplet}}\frac{\rho_{air}}{\rho_l}D\left(1+\frac{0.3\text{Sc}^{1/3}}{F(\text{Bm})} \left(\frac{\rho_{air} v r_{droplet}}{\mu_{air}}\right)^{1/2}\right)\text{ln}\left(1+\text{Bm}\right)
$$

where Sc is the Schmidt number, $\text{Sc} = \frac{\mu_{air}}{\rho_{air}D}$

We will compute the evolution of the velocity of the droplet but also the
evcolution of the droplet size. For that, we will solve a system coople system
of 2 differential equations. The system is:

$$ \frac{dr_{droplet}}{dt} = -\frac{1}{r_{droplet}}\frac{\rho_{air}}{\rho_l}D\left(1+\frac{0.3\text{Sc}^{1/3}}{F(\text{Bm})} \left(\frac{\rho_{air} v r_{droplet}}{\mu_{air}}\right)^{1/2}\right)\text{ln}\left(1+\text{Bm}\right)
$$

$$ \frac{\partial v}{\partial t} = -\frac{9 \mu_{air}v}{r_{droplet}^2}\left(1+0.15\left(\frac{\rho_{air} r_{droplet} v}{\mu_{air}}\right)^{0.687}\right) - g$$

We want to have a dimensionless system. Indeed, the data we want to use, such as
the droplets parameters or the liquid parameters are dimensionless.

As a characteristic length scale, we use the radius of the bubble that generated
the liquid jet and the droplet we are working on.

For the fluid property, we will use the Laplace and the Bond numnber. They
express  as $\rho_{liq}\gamma R/\mu_{liq}^2$ and $\rho_{liq}g
R^2/\gamma$ respectivly.

We assume that we work with water and air. For that, we then fix the viscosity
and the density ratio. We have:

$$
\frac{\rho_{air}}{\rho_{liq}} = \frac{1}{998} \quad \frac{\mu_{air}}{\mu_{liq}}=\frac{1}{55}
$$

Finally, concerning the timescale, we will use the same time-scale as the one
used for the bursting bubble problem. This problem is driven by surface tension
forces. This imply that we are working arround the capillary length scale. At
the capillary scale, there is propagation of capilary wave. We will then use 
the inverse of the frequency for our characteristiq time-scale. This time is 
also called the capillary time-scale. It is define as $t_c = \sqrt{\frac{\rho 
R_{Bubble}^3}{\gamma}}$.

With all that said, the velocity must be dimensionless, as the derivation of 
the position with respect to the time (both of them are dimensionless).

Since we are playing with the liquid parameters, we choose to have the liquid 
density as a characteristic density, which fix it to 1, and fix $\gamma$ to 1, 
so the only fluid parameters that is changing is the viscosity. We set the 
liquid viscosity a the "main" viscosity, since the Laplace number is using the 
liquid property.

With all that set, the liquid viscosity can be express as: $1/\sqrt{La}$.

The gaz properties are then define with the liquid parameters:

$$
\rho_{air}= \rho_{liq}\frac{\rho_{air}}{\rho_{liq}} \quad
\mu_{air} = \frac{1}{\sqrt{La}}\frac{\mu_{air}}{\mu_{liq}}
$$

To transform the capillary number used for the velocity, with $Ca =
\mu_{liq}/\gamma v$ to a velocity that is in the code "units" system, we just
multiply Ca by $\gamma\sqrt{La}$.

The droplet radius $r_{droplet}$ is in fact the ratio between the original
radius of the bubble generating the droplet and the droplet. We will keep
noting it $r_{droplet}$. The same rule applies to $v$, since the velocity is
express in the "bursting bubble unit system" (see the simulation for more
information).

The bond number $\rho R^2 g/\gamma$ is just defining the gravity, since all the
other parameters are charectistic parameters of the problem.

We note $\rho_r = \rho_{air}/\rho_{liq}$ and $\mu_r = \mu_{air}/\mu_{liq}$

The equation system can know be transform into a dimensionless system, assuming
that $\rho_{liq}$ is a dimensionless parameter, equal to 1.

$$
\frac{dr_{droplet}}{dt} = -\frac{1}{r_{droplet}}\rho_r
\frac{1}{Sc} \sqrt{\frac{\frac{1}{\sqrt{La}}\mu_r}{\rho_l\rho_r}}
\left(\sqrt{\frac{\frac{1}{\sqrt{La}}\mu_r}{\rho_l\rho_r}}+
\frac{0.3\text{Sc}^{1/3}}{F(\text{Bm})} \sqrt{v r_{droplet}}\right)\text{ln}\left(1+\text{Bm}\right)
$$

$$ 
\frac{\partial v}{\partial t} = 
-\frac{9 \frac{1}{\sqrt{La}}\mu_{r}v}{r_{droplet}^2}
\left(1+0.15\left(\frac{\rho_{liq}\rho_r }{\frac{\mu_r}{\sqrt{La}}}
r_{droplet} v\right)^{0.687}\right) - \text{Bo}
$$

The liquid condition will be set as input parameters (with the use of the
Laplace number).

For our equation, we will also need an initial condition in velocity and in
size. This will be set with an input parameters for this code.

Finaly, there is also the external condition, impacting the evaporation of the
liquid. We will work with water, in the air. The external pressure will be 1013
hPa. We will set the air and the water temperature at 20°C. The humidity rate
will be set to 80%, since we are evaporating water droplet above a huge surface
of water. 
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "runge-kutta.h"

#define Gravity 1

/**
Initialisation of the external parameters.

 * La is the Laplace number which give us the viscosity of the drop.
 * Bo is the Bond number
 * v0 is the initial velocity of the drop
 * R is the radius of the drop (dimensionless)
*/

double La = 50000;
double Bo = 0.01;

double v0 = 0.1;
double R = 0.13;

/** 
We define a correction factor. The air is not always at rest. The correction
factor take that into account. At 0.9 it means that the air velocity is about
10% of the drop velocity. This correction is apply during the ascending part of
the droplet.
*/

double Correction = 1;
/**
We also setup some physical parameters:

 * Bm is the Spalding number (transfer rate of mass), $Bm = \frac{Y_{vap}^{surf}-Y_{vap}^{\infty}}{1-Y_{vap}^{surf}}$
 * Sc is the Schmidt number (momentum diffusivity vs mass diffusivity), $\text{Sc} = \frac{\mu_{air}}{\rho_{air}D}$
*/

// double Bm = 0.0004629324985; //Water in air at 20°C, for 80% of humidity
double Bm = 0.0019;
// double Bm = 0.0288; // Ethanol in air at 20°C, for no ethanol outside
// double Sc = 1.85/(2.4*1.20); // in the air, at 20°C

/** 
To match with the model of Veron, we can compute the diffusion coefficient
based on the droplet size.

Dstar is a corrected diffusion coefficient.*/

double Dstar (double Bo, double R) {
  double dv = 104e-9;
  double gamma = 72e-3;
  double rho = 998;
  double g = 9.81;
  double D = 2.4e-5;
  double Rb = sqrt(gamma*Bo/(rho*g));
  
  double alpha = 0.036;
  double beta = 145.5;

  double term1 = 1./(1+(dv/(R*Rb)));
  double term2 = D/(R*Rb*alpha*beta);

  double sol = D/(term1+term2);
  return sol;
}

/**
Since we take into account the possible change in D, we compute the Schmidt 
number

We use the constant macro to switch between the maxwell model and the one from
Prupacher*/

#define CONSTANT 1


double Sc(double Bo, double R) {
  double mu = 1.85e-5;
  double rho = 1.2;
  double D;
  #if CONSTANT
  D = 2.4e-5;
  #else
  D = Dstar(Bo, R);
  #endif
  // fprintf(stderr, "%g\n", D);
  return (mu/(rho*D));
}

/** 
The function F compute a specific function needed for the evaporation of the
droplet*/

double F (double Bm) {
  return (log (1+Bm) / (Bm*pow(1+Bm,0.7)) );
  // return 1;
}

double ratio(void) {
  double rho = 998;
  double gamma = 72e-3;
  double g = 9.81;
  double muG = 1.e-3/55.;

  double Rb = sqrt(gamma*Bo/(rho*g) );
  double muL = sqrt(rho*gamma*Rb/La);
  return (muL/muG);

  // return 55;
}

/**
The function Cd will compute the drag coefficient of the droplet*/

double Cd(double v, double r) {


  /**
  When the drop is rising, we take the correction on the velocity into account.*/

  double Re;

  if (v>0.){
    Re = Correction*1.204/998.*fabs(v)*2*r/(1./(ratio()*sqrt(La)));
  }
  else{
    Re = 1.204/998.*fabs(v)*2*r/(1./(ratio()*sqrt(La)));
  }
  double sol;

  if (Re <0) return 1;

  sol =(1+0.15*pow(Re,0.687));
  // fprintf(stderr, "%g %g\n",Re, sol);
  return (sol);
}

/**
The function deriv compute both of the derivative
*/

static void deriv (scalar * rl, double t, scalar * kl) {

  /**
  First, we define the variables of the 2 equations*/

  scalar v = rl[0], r = rl[1];
  scalar kv = kl[0], kr = kl[1];

  /**
  We then define the integer variables (to perform an order 4 integration)*/

  scalar x = rl[2];
  scalar kx = kl[2];

  double mu = 1./(sqrt(La));
  double muG = mu/ratio();
  double ratioD = 1.204/998.; // rho_g/rho_l
  double rhoMu = 1./ratioD*muG; // mu_g/rho_g
  foreach() {

    /**
    Velocity equation*/


    if (v[]>0.)
      kv[] = -9./2.*muG*v[]/sq(r[])*Cd(v[], r[])*Correction;
    else
      kv[] = -9./2.*muG*v[]/sq(r[])*Cd(v[], r[]);
    // kv[] = -0.15*v[]*pow((2*r[]*fabs(v[])*1./rhoMu ),(0.687))*(9./2.*(muG)/((r[]*r[])));
    // kv[] += -(9./2.*muG/(r[]*r[]))*(v[]);
    
    #if Gravity
    kv[] += -Bo;
    #endif

    /**
    Evaporation equation*/

    if (r[]>=0){
      kr[] = -1./(r[])*(ratioD)*1./Sc(Bo, r[])*rhoMu;
      // kr[] *=(1+0.3*pow(Sc,1./3.)/F(Bm)*sqrt(fabs(2*v[])*r[]*1./rhoMu));


      if (v[]>0.)
        kr[] *=(1+0.3*pow(Sc(Bo, r[]),1./3.)*sqrt(fabs(2*v[]*Correction)*r[]*1./rhoMu));
      else
        kr[] *=(1+0.3*pow(Sc(Bo, r[]),1./3.)*sqrt(fabs(2*v[])*r[]*1./rhoMu));
      kr[] *= log(1+Bm);
    }

    /**
    Integration of velocity*/

    kx[] = v[];

    /**     
    In some case, the radius of the drop can be negative. It means that
    there is no more droplet. We filter that in this function to avoid the
    error. The derivative part of everything goes to 0.*/

    if (r[] <0) {
      kr[] = 0;
      kv[] = 0;
      kx[] = 0;
    }
  }
}

double stokesVelo(double r) {
  double mu = 1./(sqrt(La));
  double muG = mu/ratio();
  double ratioD = 1.204/998.; // rho_g/rho_l
  double rhoMu = 1./ratioD*muG;// mu_g/rho_g

  double sol = 2./9.1/muG*Bo*sq(r);
  return sol;
}

/** We need to compute the position of the droplet. We will use it as stop
condition.*/

double xn1 (double xn, double h, double vn) {
  return (xn + h*vn);
}

/** This code is use to compute the evaporation of the droplet (and there
evolution). To know the amount of water that goes in the air, we just need to
evluate the evolution of the volume of the drop.

r1 is the initial radius, r2 is the second radius*/

double evaporation (double r1, double r2) {
  double vol1 = 4/3.*pi*pow(r1,3);
  double vol2 = 4/3.*pi*pow(r2,3);
  return vol1-vol2;
}

/**
We want to compare the evaporation of the droplet with the one coming from a surface of water. We fix the following parameters:

 * The surface of evaporation is a disc of radii 5 times the bubble radius
 * The high of integration (between the water surface and the infinity) is 
 take equal to 1 meter
 * The time of integration is equal to the final time of the droplet movements

*/

double waterSurface (double tFinal) {
  double muR = ratio();
  double L = 1;
  double BoUnit = 998*9.81*sq(L)/72e-3;
  double r = sqrt(BoUnit/Bo);
  double S = sq(5)*pi;
  double volume = -tFinal*1./(sqrt(La)*muR)*S/r*1./Sc(1,1)*(0.0115-0.013488786);
  return volume;
}

/**
We set the time steps at $1.10^{-1}$*/

double dt = 0.01;
int order = 4;
double x0 = 0.000001;

int main(int argc, char const *argv[]) {

  /**
  All the external parameters can be change*/
  
  if (argc >= 2){
    La = atof (argv[1]);
  }

  if (argc >=3) {
    Bo = atof (argv[2]);
  }

  if (argc >=4) {
    v0 = atof (argv[3]);
  }

  if (argc >=5) {
    R = atof (argv[4]);
  }

  if (argc >=6) {
    Correction = atof (argv[5]);
  }

  /**
  The air viscosity is define by the following equation:

  $$ \mu_{air} = \frac{1}{\mu_{ratio}\sqrt{\text{La}}}$$
  
  where $\mu_{ratio}$ is the ratio between the liquid viscosity and the air
  viscosity. Here this value is equal to 55*/

  // double mu = 1./(55*sqrt(La));

  /**
  We will use 4 tables for our resolution. One is for the radius, the second 
  one for the time, the third one for the velocity, and the last one is for the
  position.*/

  // fprintf(stderr, "%g\n", ratio());

  init_grid(1);

  scalar velocity[];
  scalar rayon[];
  scalar position[];

  
  // double r[size];
  // double t[size];
  // double v[size];
  // double x[size];

  /**
  We define the lost of volume and the original volume of the drop*/

  double volLostRK = 0;
  double volDrop = 4./3.*pi*pow(R,3);

  /**   
  If the drop lost 99.5% of its volume, then we assume that all the drop
  has evaporated*/

  double epsilonDrop = 0.1e-2;

  
  /**

  We initialise the variables. The initial radius correspond to the input
  radius. The initial velocity is the input velocity with the basilisk
  dimension. Indeed, we input the velocity as a dimensionless parameter, the
  Capillary number. We change it back to the "code" velocity. */

  // long int i = 0;

  // double ri = R;
  // double vi = v0*sqrt(La);
  // double xi = x0;
  double ti = 0.;


  double riRK = R;
  double viRK = v0*sqrt(La);
  double xiRK = x0;

  foreach() {
    velocity[] = viRK;
    rayon[] = riRK;
    position[] = xiRK;
  }

  // fprintf(stdout, "%f %f %f %f %f %f %f 0 0\n", ti, xi, xiRK, vi, viRK, ri, riRK);
  fprintf(stdout, "%f %f %g %g 0\n", ti, xiRK, viRK, riRK);

  /**
  Resolution*/

  while (xiRK>=0. && ((volDrop - volLostRK)/(volDrop)>epsilonDrop ) && riRK>0.) {
    // t[i] = t[i-1] + dt;
    ti += dt;

    double ri_1RK = riRK;

    // double vi_1 = vi;
    // double ri_1 = ri;
    // double xi_1 = xi;

    runge_kutta ({velocity, rayon, position}, (ti-dt), dt, deriv, order);

    /**
    We update our list v. We also "reset" the velocity scalar before the next 
    resolution for the radius*/

    // foreach(){
    //   // v[i] = velo[];
    //   viRK = velocity[];
    //   // velo[] = v[i-1];
    //   velocity[] = vi_1RK;
    //   rayon[] = ri_1RK;
    // }

    /**
    For a couple system, the parameter we are working one must be the first.*/

    // runge_kutta({rayon, velocity}, (ti-dt), dt, dray, order);

    // foreach() {
    //   // r[i] = ray[];
    //   riRK = rayon[];
    //   // velo[] = v[i];
    //   velocity[] = viRK;
    // }

    /**
    We update the parameters for the output and the previous step*/

    foreach() {
      viRK = velocity[];
      riRK = rayon[];
      xiRK = position[];
    }

    // vi += dt * velo(vi_1, ri_1);
    // ri += dt * ray(vi_1, ri_1);

    // x[i] = xn1(x[i-1], dt, v[i-1]);
    // xiRK = xn1(xi_1RK, dt, vi_1RK);
    // xi = xn1(xi_1, dt, vi_1);


    /**
    We also compute the lost in volume*/

    // volLost += evaporation(r[i-1], r[i]);
    volLostRK += evaporation(ri_1RK, riRK);
    // volLost += evaporation(ri_1, ri);


    /**
    Once the equations are solved, we output the time, the
    variation in position, the velocity, the radius and the total lost since 
    the beginning*/

    //fprintf(stdout, "%f %f %f %f %f %f %f %g %g\n", ti, xi, xiRK, vi, viRK, ri, riRK, volLost, volLostRK);
    fprintf(stdout, "%f %f %g %g %g %g\n", ti, xiRK, viRK, riRK, volLostRK, evaporation(ri_1RK, riRK));

    // fprintf(stderr, "%f %g %g %g\n", ti, viRK, -stokesVelo(riRK), viRK + stokesVelo(riRK) );

  }

  if (fabs((volDrop - volLostRK)/(volDrop))<epsilonDrop)
    volLostRK = volDrop;

  if (volDrop - volLostRK <0.)
    volLostRK = volDrop;

  /**
  Since we are interested in the lost of volume, we output it as a final results*/


  fprintf(stderr, "%g\n", volLostRK);

  /**
  We also output the flying time and the volume loss with a surface of
  water during the flying time*/

  double waterLoss = waterSurface(ti);

  fprintf(stderr, "%g %g\n", ti, waterLoss);

  return 0;
}

/**
We can observe the evolution of the position of the drop:
~~~gnuplot
plot 'out' u 1:2 w l t 'Evolution of the position'
~~~

We can also plot the radius evolution

~~~gnuplot
plot 'out' u 1:4 w l t 'Evolution of the droplet radius'
~~~

Lastly, we can observe the evolution of the "lost" volume (or, the total amount of water lost since the beginning of the travel)

~~~gnuplot
plot 'out' u 1:5 w l t 'Evolution of the lost volume'
~~~
*/