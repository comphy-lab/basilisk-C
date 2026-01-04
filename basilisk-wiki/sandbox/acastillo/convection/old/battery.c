 /**
 
# THE BATTERY PROBLEM

A layer of two (or more) liquids is located between two horizontal parallel
planes of extension *W* in *x*, where coordinate *y* is directed upwards.

- Liquid A is located at the top between *y=$H_E$/2* and *y=H/2*
- Liquid B is located at the bottom between *y=-H/2* and *y=-$H_E$/2*,
- Liquid E is located at the top of Liquid 1 between *y=-$H_E$/2* and *y=$H_E$/2*.

Where the densities of each fluid is such that $\rho_B < \rho_E < \rho_A$.

Each fluid is characterized by a density $\rho$, a thermal expansion coefficient $\beta$, a dynamic viscosity $\mu$,
a thermal conductivity $\lambda$, and an specific heat capacity $C_p$.

The top and bottom walls have an imposed temperature $T_0$, the sidewalls are adiabatic and no-slip conditions
are imposed on all boundaries.
An electrical current $J_0$ is imposed on the system, which in turn creates a *Joule heating* effect due to the
electrical ressitivity of the fluids. This heating effect takes the form of an internal heat source per unit of volume, $Q_0=J_0^2/\sigma$,
with $\sigma$ being the electrical conductivity of the fluid.

![Schematic representation of three layers of fluid contained between parallel planes](/battery.png)


## Characteristic Scales

We define $\rho_c$, $\beta_c$, $\mu_c$, $\lambda_c$, $C_{pc}$, $\alpha_c = \frac{\lambda_c}{\rho_c C_{pc}}$ and
$\nu_c = \frac{\mu_c}{\rho_c}$ as characteristic values for our system, and we take these values from the middle layer
which correspond to a molten salt electrolyte (LiCl-KCl). */

    #define RHOC (1.63e3)
    #define LAMBDAC (0.42)
    #define MUC (1.15e-3)
    #define BETAC (2.9e-4)
    #define CPC (1.21e3)
    #define ALPHAC (LAMBDAC/(RHOC*CPC))
    #define NUC (MUC/RHOC)
    #define Q0C (588235.)
    #define G (9.81)

/**
We can write the properties for each fluid as a function of these characteristic values in dimensionless form as :
- Density : $\tilde{\rho}_0 = \frac{\rho_0}{\rho_c}$
*/

    #define RHO0A (0.491e3/RHOC)
    #define RHO0B (9.23e3/RHOC)
    #define RHO0E (1.63e3/RHOC)

/**
- Coefficient of thermal expansion of the fluid : $\tilde{\beta} = \frac{\beta}{\beta_c}$
*/

    #define BETAA (1.8e-4/BETAC)
    #define BETAB (1.82e-4/BETAC)
    #define BETAE (2.9e-4/BETAC)

/**
- Dynamic viscosity : $\tilde{\mu} = \frac{\mu}{\mu_c}$
*/

    #define MUA (3.5e-4/MUC)
    #define MUB (1.29e-3/MUC)
    #define MUE (1.15e-3/MUC)


/**
- Thermal Conductivity : $\tilde{\lambda}=\frac{\lambda}{\lambda_c}$
*/

    #define LAMBDAA (52./LAMBDAC)
    #define LAMBDAB (16./LAMBDAC)
    #define LAMBDAE (0.42/LAMBDAC)


/**
- Heat Capacity : $\tilde{C_p} = \frac{C_p}{C_{p,c}}$
*/

    #define CPA (4.24e3/CPC)
    #define CPB (0.188e3/CPC)
    #define CPE (1.21e3/CPC)


/**
- Thermal diffusivity : $\tilde{\alpha} = \tilde{\lambda}/\tilde{\rho(T_0)}\tilde{C_p}$
*/

    #define ALPHAA (LAMBDAA/(RHO0A*CPA))
    #define ALPHAB (LAMBDAB/(RHO0B*CPB))
    #define ALPHAE (LAMBDAE/(RHOE*CPE))


/**
- Internal Heating : $\tilde{Q_0} = \frac{Q_0}{Q_c}$
*/

    #define Q0A (0./Q0C)
    #define Q0B (0./Q0C)
    #define Q0E (588235./Q0C)


/**
where values for fluid *A* correspond to a light metal (Molten Lithium) and values for fluid B to a heavy metalloid (Pb-Li17) at a reference
temperature of 723K.



## Dimensionless Parameters

We define a characteristic length, temperature and velocity scales,

- $[L] = H$
- $[\Theta]= Q_0[L]^2/\lambda_c$
- $[U] = \sqrt{\beta_c g [\Theta][L]}$


From these values we define the following dimensionless numbers,

- A *Reynolds* number, $Re_c = [U][L]/\nu_c$
- A *Prandtl* number,  $Pr_c = \nu_c/\alpha_c = 3.31$
- A *Rayleigh* number, $Ra_c = g \beta_c [\Theta][L]^3/{\nu_c \alpha_c}$ = $Pr_c Re_c^2$
- A dimensionless volume change, $B = \beta_c [\Theta]$

There exist a relation between $B$ and $Re_c$,
$$
B = \pi_1 Re_c^{4/5}
$$
where $\pi_1$ is a dimensionless group,
$$
\pi_1 (Q_0) \equiv  \left[{\frac{\beta^3_c Q^3_0 \nu^4_c}{g^2 \lambda_c^3}}\right]^{1/5}
$$
such that we can rewrite the parameter $B$ as a function of $\pi_1$ and $Re_c$,
*/

    #define REYNOLDS (5000.)
    #define PRANDTL (3.31)
    #define PI1 pow(pow(BETAC,3.)*pow(Q0C,3.)*pow(MUC,4.)/pow(RHOC,4.)/pow(G,2.)/pow(LAMBDAC,3.),1./5.)

    double Re = REYNOLDS, Pr = PRANDTL, Ra, B;

/**
# System of Equations

The resulting equation system is of the form
$$
\nabla \cdot \vec{u} = 0
$$
$$
\partial_t \vec{u} + \nabla \cdot \left(\vec{u} \otimes \vec{u} \right) =
	\frac{1}{\tilde{\rho}} \left[-\nabla p  + \frac{1}{Re_c} \cdot \left( \tilde{\mu} \nabla \vec{u} \right) \right] +
   \frac{1}{B}\left[\frac{1}{\tilde{\rho}}-1\right] \vec{e}_y
$$
$$
  \tilde{\rho} \tilde{C_p} \left[ \partial_t T + \nabla \cdot \left( \vec{u} T \right) \right] = \frac{1}{Re_c Pr_c} \left[  \nabla \cdot
   \left[ \tilde{\lambda} \nabla T \right] + \tilde{Q}_0  \right]
$$

where $\tilde{\rho} = \tilde{\rho}_0[1-\tilde{\beta} B T]$

For this problem we use a quad-tree grid and the centered Navier-Stokes solver.
*/

    #include "navier-stokes/centered.h"


/**
## Interface between fluids
We use Volume-of-Fluid advection and define the interfaces 12 and 23. We define $Z_E$ as the vertical position of the center of the middle layer, and $H_E$ as its depth.
*/

    #include "vof.h"
    scalar f12[], f23[];
    scalar * interfaces = {f12,f23};
    #define ZE 0.0
    #define HE 0.3

/**
The properties of the fluid ($\tilde{\rho},\tilde{\beta},...$) are evaluated as functions of the volume fraction of each fluid.

  $$ \tilde{\rho} = (1-f_{12})(1-f_{23})\tilde{\rho}_B + f_{12}~(1-f_{23})~\tilde{\rho}_E + f_{12}~f_{23}~\tilde{\rho}_A $$
  $$ \tilde{\beta} = (1-f_{12})(1-f_{23})\tilde{\beta}_B + f_{12}~(1-f_{23})~\tilde{\rho}_E + f_{12}~f_{23}~\tilde{\beta}_A $$
  $$ \vdots $$

*/

    #define FLUID(f12,f23,x1,x2,x3) (((1-(f12))*(1-(f23))*x1) + ((f12)*(1-f23)*x2) + ((f12)*(f23)*x3))
    #define RHO(RHO0,BETA,T) (RHO0*(1-BETA*B*(T)))
    #define RHOTILDE(f12,f23,T) 	FLUID(f12,f23,RHO(RHO0B,BETAB,T),RHO(RHO0E,BETAE,T),RHO(RHO0A,BETAA,T))
    #define BETATILDE(f12,f23) 	FLUID(f12,f23,BETAB,BETAE,BETAA)
    #define MUTILDE(f12,f23) 	FLUID(f12,f23,MUB,MUE,MUA)
    #define LAMBDATILDE(f12,f23) 	FLUID(f12,f23,LAMBDAB,LAMBDAE,LAMBDAA)
    #define CPTILDE(f12,f23) 	FLUID(f12,f23,CPB,CPE,CPA)
    #define Q0TILDE(f12,f23) 	FLUID(f12,f23,Q0B,Q0E,Q0A)

/**
# Solution Approach
We need to solve an advection--diffusion equation for the temperature.
We use the corresponding predefined solvers. */

    #include "tracer.h"
    #include "diffusion.h"


/**  We define a scalar field for the temperature. The temperature is the only
scalar field transported by the flow i.e. it is the only field in the list
of tracers required by [tracer.h](/src/tracer.h). */

    scalar T[];
    scalar * tracers = {T};


/** We set the total length of the simulation. */

    double EndTime = 3.;

/** and store the statistics on the Poisson solver for the diffusion
   step in `mgT`. */

    mgstats mgT;

/**
## Boundary Conditions
We impose no-slip boundary conditions on all the walls, and impose
a temperature zero top and bottom plates.
Side-walls are considered to be periodic. */

    T[top] 	    = dirichlet(0.);
    u.t[top]    = dirichlet(0.);
    u.n[top]    = dirichlet(0.);

    T[bottom]   = dirichlet(0.);
    u.t[bottom] = dirichlet(0.);
    u.n[bottom] = dirichlet(0.);

    u.t[left]   = periodic();
    u.n[left]   = periodic();
    T[left] = periodic();


/**
## Base state and initial conditions
The base state of the system is in a configuration with zero velocity and a temperature profile
that satisfies the following relation
$$
\frac{d}{dz}\left[\lambda(z) \frac{dT}{dz} \right] = -S_q(z), \quad \quad T(\pm H/2) = T_0 \quad \quad \mbox{  for } t=t_0
$$
which takes the form of a piecewise function.

And the initial condition correspond to the base state plus an initial perturbation.
*/

    #define R1  (1./LAMBDAB)
    #define R2  (1./LAMBDAE)
    #define R3  (1./LAMBDAA)
    #define C0  (Q0E*HE*(HE*(R2-R3) + R3)/(HE*(2*R2 - R1 - R3) + (R1+R3)))
    #define ZMAX = C0/Q0E - HE/2.;

    #define DELTA1 (R1*C0*(1.-HE)/2.)
    #define DELTA2 (R2*(C0*HE - Q0E*pow(HE,2.)/2.))
    #define DELTA3 (R3*(C0 - HE*Q0E)*(1.-HE)/2.)
    #define FDELTA(x) (pow(x,2.) + (x)*HE)

    #define SOLUTION1(y) (C0*R1*(y+1./2.))
    #define SOLUTION2(y) (DELTA1 + C0*(y+HE/2.)*R2 - Q0E/2.*R2*(FDELTA(y) FDELTA(-HE/2.)))
    #define SOLUTION3(y) (-DELTA3+(C0-HE*Q0E)*(y -HE/2.)*R3)

    event init (t = 0)
    {
      scalar phi12[];
      foreach_vertex()
        phi12[] = 0.00*cos(2.*pi*x) + y - (ZE-HE/2.);
      fractions (phi12, f12);

      scalar phi23[];
      foreach_vertex()
        phi23[] = 0.00*cos(2.*pi*x) + y - (ZE+HE/2.);
      fractions (phi23, f23);

      foreach(){
        T[] = FLUID(f12[],f23[],SOLUTION1(y),SOLUTION2(y),SOLUTION3(y)) + (6e-5)*noise();
        u.x[] = (6e-5)*noise();
        u.y[] = (6e-5)*noise();
      }
      boundary ({T,u});

      mu = new face vector;
      alpha = new face vector;
      a = new face vector;

    }


/**
## Update fluid properties
Fluid properties are going to depend on the volume fraction for each fluid and
on the temperature field.
*/

    scalar alphac[];
    event properties (i++) {

      face vector muv = mu;
      face vector alphav = alpha;

      foreach_face() {
        double f12m = (f12[] + f12[-1,0])/2.;
        double f23m = (f23[] + f23[-1,0])/2.;
        double Tm = (T[] + T[-1,0])/2.;
        alphav.x[] = 1./RHOTILDE(f12m,f23m,Tm);
        muv.x[] = (1./sqrt(Re))*MUTILDE(f12m,f23m);
      }

      foreach()
        alphac[] = 1./RHOTILDE(f12[],f23[],T[]);
    }



/**
## Buoyancy Force
The buoyancy force is added trough the acceleration term *a*, which depends on the temperature field $T$.
Hence, we solve the equation for the temperature fields, and then we add the corresponding forcing term.
*/
    event tracer_diffusion (i++) {
      face vector D[];
      foreach_face(){
        double f12m = (f12[] + f12[-1,0])/2.;
        double f23m = (f23[] + f23[-1,0])/2.;
        D.x[] = LAMBDATILDE(f12m,f23m);
      }

      scalar r[];
      foreach()
        r[] = Q0TILDE(f12[],f23[]);

      scalar beta[];
      foreach()
        beta[] = 0.;

      scalar g[];
      foreach()
        g[] = Re*Pr*CPTILDE(f12[],f23[])/alphac[];

      mgT = diffusion (T, dt, D, r, beta, g);
    }

    event acceleration (i++) {
    	face vector av = a;
    	foreach_face(y)
        av.y[]  = 1./B*(alpha.y[]-1.);
    }

    #define MINLEVEL 8
    #define MAXLEVEL 10
    #if TREE
    event adapt (i++) {
      adapt_wavelet ((scalar *){T}, (double[]){1e-3}, MAXLEVEL);
    	refine (level < MINLEVEL);
    }
    #endif

/**
# Simulation Parameters
Finally we do consecutive runs by changing the $Re_c$, $Pr_c$ and *B*, and then
calling the run() method of the Navier–Stokes solver.
We set the size of the domain `L0` and select an initial and maximum
time-steps, the CFL criteria, and the grid size */

    int main() {
      L0 = 1.;
      X0 = Y0 = -0.5;
      DT = 0.01;
      TOLERANCE = 1e-6;

      Re = REYNOLDS ; Ra = Pr*pow(Re,2.); B = PI1*pow(Re,4./5.); N=1<<MINLEVEL;
      run();
    }

/**
## Outputs
*/

    event logfile (i++)
    {
      fprintf (stderr, "%d %g %g %d %d %d \n", i, t, dt, mgT.i,mgp.i,mgu.i);
    }

    #include "output_fields/output_vtu_foreach.h"
    void backup_fields (scalar T, vector u, int nf)
    {
    	char name[80];
    	FILE * fp ;

      scalar lvl[];
      foreach()
        lvl[] = level;

    	nf > 0 ? sprintf(name, "battery_%6.6d_n%3.3d.vtu", nf,pid()) : sprintf(name, "battery_n%3.3d.vtu",pid());
    	fp = fopen(name, "w");
    	output_vtu_bin_foreach ((scalar *) {T, lvl, f12,f23}, (vector *) {u}, N, fp, false);
    	fclose (fp);
    }

    event logfile (t = EndTime) {
    	backup_fields(T,u,0);
    }

    event export (t += 1.0 ; t <=EndTime)
    {
      static int nf = 0;
      backup_fields(T,u,nf);
      nf++
    }


/**
We record the Nusselt number and kinetic energy and enstropy every
X time units.
*/

    event time_series (t += 0.05; t <= EndTime)
    {
    	char cmd[1024];
    	FILE * fp;
    	sprintf(cmd, "battery_%g_%g_%g_%g_%i.asc", Re,Pr,Ra,B,N);
    	fp = i > 0 ? fopen(cmd, "a") : fopen(cmd, "w");

      double ekin_a=0., ekin_e=0., ekin_b=0. ;
      foreach(reduction(+:ekin_a),reduction(+:ekin_e),reduction(+:ekin_b))
        foreach_dimension()
        {
          ekin_a += f23[]*sq(u.x[]);
          ekin_e += f12[]*(1-f23[])*sq(u.x[]);
          ekin_b += (1-f12[])*(1-f23[])*sq(u.x[]);
        }

    	fprintf (fp, "%g %g %g %g \n", t, ekin_a, ekin_b, ekin_e);
      fprintf (stderr, "%g %g %g %g \n", t, ekin_a, ekin_b, ekin_e);

    	fclose (fp);

    }
