/**
# Elasto-viscoplastic 2D droplet in a Couette Newtonian shear flow

We model the shear deformation of an elasto-viscoplastic (Saramito's models (2007) and (2009)) droplet in a Newtonian matrix. After a time, the imposed deformation stops, and the droplet is left to retract.

The polymeric viscosity of the droplet ($\mu_{p1}$) is calculated as a function of the polymeric stress tensor ($\mathbf{\tau_p}$). In turn, $\mu_{p1}$ is used to calculate $\tau_p$ at the next time step and the relaxation time ($\lambda$) at the current time step. If we multiply $\mu_{p1}$ by the volume fraction, ($f$), we change $\mu_{p1}$ on the interface and, consequently, we change $\mathbf{\tau_p}$. To work around this problem, we calculate $\mu_{p1}$ in the cells on the interface as the average of  $\mu_{p1}$ in the neighboring cells (not using neighboring cells also on the interface).

## Problem nondimensionalization

The Reynolds number ($Re$) is based on the outer fluid (2) properties:

$Re = \rho_2 D U / \mu_2$

$D$ is the droplet diameter and $U$ is the velocity difference at the droplet poles at $t = 0.0 s$. The characteristic stran rate is $\dot{\gamma}_c = U/D = 1 s^{-1}$.

The Weber number ($We$) is defined as:

$We = \rho_2 D U^2 / \sigma$


Where, $\sigma$ is the surface tension coefficient. The Deborah number ($De$) is defined as:

$De = t_r / t_f = \lambda_c \dot{\gamma}_c$

Where, $t_r$ and $t_f$ are the relaxation and flow characteristic times, respectively. The relaxation time ($\lambda$), is given as:

$\lambda = \mu_{p1} / G$

Where, $G$ is the elastic modulus, which is a constant. Therefore, $\lambda$ varies with $\mu_{p1}$, which is given as (saramito's models):

Bingham fluid:
$\mu_{p1} = |\tau_{pd}| / (|\tau_{pd}| - \tau_y)$

Herschel-Bulkley fluid:
$\mu_{p1} = (K|\tau_{pd}|^n / (|\tau_{pd}| - \tau_y))^{1/n}$

Here, the subscript $d$ stands for deviatoric, $\tau_y$ is the yield stress, $K$ is the consistency index, and $n$ is the flow index.

The characteristic viscosity of the droplet is $\mu_{c1}$, and it is composed of a solvent contribution ($\mu_{s1}$) and a polymeric contribution ($\mu_{pc1}$).

Bingham fluid:
$\mu_{c1} = \mu_{s1} + \mu_{pc1} = \mu_{s1} + \tau_y / \dot{\gamma}_c + \mu_b$

Herschel-Bulkley fluid:
$\mu_{c1} = \mu_{s1} + \mu_{pc1} = \mu_{s1} + \tau_y / \dot{\gamma}_c + K\dot{\gamma}_c^{n-1}$

Also, $\mu_{s1} = \beta\mu_{c1}$ and $\mu_{pc1} = (1-\beta)\mu_{c1}$. The viscosity and density ratios are defined as:

$\mu_r = \mu_2 / \mu_{c1}$ and $\rho_r = \rho_2 / \rho_1$

We define the characteristic relaxation time:

$\lambda_c = \lambda \frac{\mu_{pc1}}{\mu_{p1}} = \frac{\mu_{p1}}{G} \frac{\mu_{pc1}}{\mu_{p1}} = \frac{\mu_{pc1}}{G}$

Finally, the Plastic number ($\Pl$) is defined as:

$Pl = \tau_y / (\mu_{pc1}\dot{\gamma}_c)$

## Code
*/

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "log-conform.h"
#include "tension.h"
#include "view.h"

#define Pl 0.5      // Plastic number
#define Re 0.5      // Reynold number
#define We 0.5     // Weber number
#define Deb 5.0     // Deborah number
#define MUr 1.      // ratio of outer(matrix) to inner(droplet) viscosity
#define RHOr 1.     // ratio of outer to inner density
#define Beta 0.5    // ratio of the solvent visc. to the total viscoelastic visc.

int MAXLEVEL = 9;
int MINLEVEL = 6;

double tau_y, K, epsilon, gamma_c, muc1, mupc1;
double D = 2.;       // Droplet diameter
double U = 2;        // Characteristic velocity U = V * D / H
double nnn = 0.8;    // Flow index
double Nreg = 1.e1;  // Dimensionless regularization parameter.
double tEnd = 5.;

scalar mupd[], lam[], yielded[], tau_norm[];

/**
The top and bottom boundary conditions are those of a Couette flow. After $t = 3.0 s$ the shear deformation stops and the droplet retracts */

u.t[top]    = dirichlet (t <= 3. ? y : 0);
u.t[bottom] = dirichlet (t <= 3. ? y : 0);

/** 
The domain is a 16x16 box which will be masked later to a become a 16x8 size box. */

int main() {
  L0 = 16.;
  origin (-L0/2, -L0/4.);
  periodic (right);
//  DT = 1e-3;

  /**
  We set the viscosities, densities, surface tension and visco-elastic
  parameters. We also calculate $\tau_y$ and $K$. $\mu_{p1}$ and $\lambda$ are defined in the properties event. */
  
  gamma_c = U/D;
  muc1 = 1/(MUr*Re);        // Droplet characteristic viscosity
  mu1 = Beta*muc1;          // Solvent viscosity of the droplet
  mupc1 = (1-Beta)*muc1;     // Polymeric characteristic viscosity of the droplet
  tau_y = Pl*mupc1*gamma_c;  // Yield stress
  K = (mupc1-tau_y/gamma_c)/pow(gamma_c, nnn-1);  // Powe law index
  epsilon = tau_y/(Nreg * muc1);  // Regularization parameter
  mu2 = 1./Re;      // Outer phase Newtonian phase viscosity
  rho1 = 1/RHOr;    // Droplet density
  rho2 = 1.;        // Outer phase density
  f.sigma = 1./We;  // Surface tension

  lambda = lam;   
  mup = mupd;

  init_grid (1 << 8);
  run();
}


/**
The droplet radius is equal to 1.
*/
event init (i = 0) 
{
  mask (y > 4 ? top : none);

  fraction (f, 1. - (sq(x) + sq(y)));
  foreach()
    u.x[] = y;
  boundary ((scalar *){u});
}



event properties (i++)
{
  double ff = 0.99;

  foreach()
  {
    tau_norm[] = sqrt(sq((tau_p.x.x[]-tau_p.y.y[])/2) + sq(tau_p.x.y[]));  // Norm of the deviatoric part of the polymeric stress tensor

     if (tau_norm[] <= tau_y + epsilon && f[] > ff)
    {
      mupd[] = max((tau_y + epsilon)/(epsilon), 0.);  // Bingham fluid
  //    mupd[] = pow( (K * pow(tau_norm[] + epsilon, nnn)/(epsilon)), 1/nnn);   // Herschel-Bulkley fluid
    }
    else if (f[] < 1. - ff)
      mupd[] = 0.;
    else if (f[] > 1. - ff && f[] < ff) // We calculate the polymeric viscosity on the interface as the average of the neighboring cells
    {
      int c1, c2, c3, c4;
      c1 = c2 = c3 = c4 = 1;
      if (f[-1,0] > 1. - ff && f[-1,0] < ff)
       c1 = 0;
      if (f[1,0] > 1. - ff && f[1,0] < ff)
       c2 = 0;
      if (f[0,-1] > 1. - ff && f[0,-1] < ff)
       c3 = 0;
      if (f[0,1] > 1. - ff && f[0,1] < ff)
       c4 = 0; 

      mupd[] = (mupd[-1,0]*c1 + mupd[1,0]*c2 + mupd[0,-1]*c3 + mupd[0,1]*c4) * f[] / (c1 + c2 + c3 +c4);
    }
    else{
      mupd[] = max(tau_norm[]/(tau_norm[] - tau_y), 0.);     // Bingham fluid
  //    mupd[] = pow( (K * pow(tau_norm[], nnn)/(tau_norm[] - tau_y)), 1/nnn);    // Herschel-Bulkley fluid
    }

    lam[] = max(Deb*mupd[]/mupc1,0); 
  }
  boundary ({tau_norm, mupd, lam});

 foreach()  // We identify the yield and unyield regions
 {
    if(tau_norm[] > tau_y*f[] || f[] < ff)
    {
      yielded[] = 1.;
    }
    else
    {
      yielded[] = 0.; 
    }
  }
  boundary ({yielded});
}

/**
The mesh is adapted according to the errors on volume fraction,
velocity, yield surface and polymeric viscosity. */

event adapt (i++) {
  adapt_wavelet ({f, u, yielded, mupd}, (double[]){1e-3, 1e-2, 1e-2, 1e-2, 1e-2}, MAXLEVEL, MINLEVEL);
}


/**
Here we generate some videos
*/
event movies (i += 10)
{
  clear(); 
  view(fov = 6., tx = 0.0, ty = 0.0);
  
  scalar vel[];
  foreach()
    vel[] = pow(sq(u.x[])+sq(u.y[]), 0.5);
  boundary({vel});
  draw_vof ("f");
  squares("vel", linear = true);
  save("movie_vel.mp4");

  draw_vof ("f");
  squares("mupd", linear = true);
  save("movie_mu.mp4");

  scalar tau_pp[];
  foreach()
  {
      tau_pp[] = pow(0.5*(sq(tau_p.x.x[])+sq(tau_p.y.y[])+sq(tau_p.x.y[])+sq(tau_p.y.x[])),0.5);
  }
  boundary({tau_pp});
  draw_vof ("f");
  squares("tau_pp", linear = true);
  save("movie_tau_p.mp4");

  draw_vof ("f");
  squares("yielded");
  save("movie_yield.mp4");
}

event end (t = tEnd) {}


/**
## Results
Results for the Bingham fluid case.


### Velocity field
![Velocity field](shear_droplet2/movie_vel.mp4)

### Polymeric stress field
![Polymeric stress](shear_droplet2/movie_tau_p.mp4)

### Polymeric viscosity ($\mu_p$) field
![Pressure field](shear_droplet2/movie_mu.mp4)

### Yielded region
![Yielded region](shear_droplet2/movie_yield.mp4)

## References
*/