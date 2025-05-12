/** 
# Introduction to turbulence
There is no precise definition of turbulence, but it can be described as an unsteady, aperiodic motion in which all three velocity components and the other transported quantities fluctuate in space and time. Most of natural and engineering flows are turbulent, this is why modeling turbulence is this important to predict and study fluid motion. Many turbulence models have been created to describe the different physical problems.
According to the simulation, the appropriate turbulence model is chosen: RANS, DNS, LES,… 


Despite the advanced models known to date as LES and DNS, RANS still the most turbulence model used in industry for being the less expensive the most efficient one to describe the mean flow so we can understand what happening without calculating the fluctuating terms.


### 1. RANS Modelling
 Reynolds averaged Navier Stokes equations are based on the decomposition of the flow variables (like velocity $u$) into the mean (time-averaged) component $\overline{u}$ and the fluctuating component $u'$ :

$$ \mathbf{u}(\mathbf{x}, t) = \overline{\mathbf{u}}(\mathbf{x}) + \mathbf{u'}(\mathbf{x}, t)$$

where $$ \mathbf{u} = (u, v, w)$$ and $$\mathbf{x} = (x, y, z) $$

Starting from the Navier–Stokes equations of motion for an incompressible Newtonian fluid :
$$
 \frac{\partial u_i}{\partial x_i} = 0 
$$
$$
\frac{\partial u_i}{\partial t} + u_j \frac{\partial u_i}{\partial x_j}
= f_i 
- \frac{1}{\rho} \frac{\partial p}{\partial x_i}
+ \nu \frac{\partial^2 u_i}{\partial x_j \partial x_j}.
$$


Substituting,

$$ u_i = \bar{u_i} + u_i^\prime, p = \bar{p} + p^\prime$$
and taking a time-average of these equations yields,

$$ \frac{\partial \bar{u_i}}{\partial x_i} = 0$$

$$ 
 \bar{u_j}\frac{\partial \bar{u_i} }{\partial x_j}
+ \overline{u_j^\prime \frac{\partial u_i^\prime }{\partial x_j}}
= \bar{f_i}
- \frac{1}{\rho}\frac{\partial \bar{p}}{\partial x_i}
+ \nu \frac{\partial^2 \bar{u_i}}{\partial x_j \partial x_j}.  $$

On further manipulations this yields,
$$
\rho \frac{\partial \bar{u_j} \bar{u_i} }{\partial x_j}
= \rho \bar{f_i}
+ \frac{\partial}{\partial x_j} 
\left[ - \bar{p}\delta_{ij} 
+ 2\mu \bar{S_{ij}}
- \rho \overline{u_i^\prime u_j^\prime} \right ]$$
where, 
$$\bar{S_{ij}} = \frac{1}{2}\left( \frac{\partial \bar{u_i}}{\partial x_j} + \frac{\partial \bar{u_j}}{\partial x_i} \right)$$
is the mean rate of strain tensor, and $- \rho \overline{u_i^\prime u_j^\prime}$ constitute the Reynolds stress.
 
Considering the Boussinesq approximation, we assume that Reynolds stress definded above, is proportionnal to $\frac{\partial \bar{u}}{\partial y}$ as :
$$- \rho \overline{u_i^\prime u_j^\prime} = \mu_t \frac{\partial \bar{u}}{\partial y}  $$

where $\mu_t$ is the Eddy viscosity, an artificial viscosity that controls the strength of the diffusion.

Eddy viscosity models are a class of turbulence models used to calculate the Reynolds stress by modellinzing the turbulent viscosity using algebraic models (zero equation models), If we solve one or more
turbulence transport equations, then we refer to the models as one, two, etc.,
equation models, based on the number of equations solved.

In our case, we will focus on implementing the Prandtl mixing length model; an algebraical model using an empirical definition of the turbulent viscosity. 

### 2. Prandtl Mixing Length 
Algebraic models are known to be efficient with simple physical problems like channel flow without very fine meshes.

### 3. Channel flow

As a first approach of turbulence in Basilisk, the channel flow seemed to be the appropriate case in which we can test the viability of the turbulent model implemented, as it’s a fundamental physical problem with a basic geometry and the advantages that this case involve in term of simplifications associated with the symmetries of the problem.

We study the unsteady 2D flow between two parallel plane plates. As the flow is fully developed : 
$$\frac{\partial u}{\partial x} = 0  ; v = 0 ; w = 0$$

which implies
$$\frac{\partial p}{\partial y} = 0 ; \frac{\partial p}{\partial x} = cst$$

The no-slip boundary condition at the wall imposes $$ u(y = +/- \frac{H}{2}) = 0 $$

Before considering the turbulent case, let's test the laminar one.

#### 3.1 Laminar case


At laminar state, we expect at the end of the simulation a parabolic velocity profile (Poiseuille profile). As it was demonstrated long before :
According to the NS equations,
$$ \frac{\partial u}{\partial x} = 0 $$
$$ \frac{\partial p}{\partial x} = \mu  \frac{\partial ²u}{\partial x²} = cst$$
After integration, we obtain :
$$ u(y) = U_{max} (1 – (\frac{y}{H/2})²)$$ with $$U_{max} = \frac{H²}{8\mu}(- \frac{\partial p}{\partial x})$$ 


#### 3.2 Turbulent case
*/


/** 
Turbulence model on a channel flow configuration
*/
#include "embed.h"
#include "navier-stokes/centered.h"
#include "view.h"

//Defining the first and second derivative operators

#define dd_ddy_2(a) (fs.x[] && fs.x[1] ? (a[1] - 2*a[] + a[-1])/sq(Delta) : \
			    fs.x[1] ? (a[2] - 2*a[1] + a[])/sq(Delta) :		    \
			    fs.x[]  ? (a[-2] - 2*a[-1] + a[])/sq(Delta) : 0.)

#define d_dy_4(a)   (fs.x[] && fs.x[1] ? (a[-2] - 8*a[-1] + 8*a[1] - a[2] )/(12*Delta) : \
                            fs.x[1] ? (-a[2] + 4*a[1] - 3*a[])/(2*Delta) :             \
                            fs.x[]  ? (a[-2] - 4*a[-1] + 3*a[])/(2*Delta) : 0.)


face vector av[];
face vector muf[];
face vector dtau [];
scalar lm[];
scalar un[];

double k      = 0.41; // Karman constant
double H      = 2.;    // Total channel height
double Re     = 2*0.32500000E+04*27/20.;
double Tau;
double Re_target = 180.;


#define N_max                         7
#define N_min                         3
#define Xplot                         0.0
#define Umax                       27/20.
#define mu1                      (Umax*H/Re)
#define dpdx      8.*sq(mu1)*Re/(H*H*H)
#define domain_size     4.* pow(2. , N_max) / (pow(2. , N_max) - 2.)

const face vector unit[] = {dpdx,0.};


int main()
{
  L0 = domain_size;
  origin (-L0/2., -L0/2.);
  periodic (right);
  DT = 0.1;
//  stokes = true;
  TOLERANCE = 1e-7;

  N = 128;
  run();

}

event init (t = 0) {

  a = av;
  mu = muf;

  /**
  The channel geometry is defined using Constructive Solid Geometry. */

  solid (cs, fs,intersection(
			  -y+H/2. ,// y positive halve
			  +y+H/2. //  y negative halve
			  ));

  /**
  The boundary condition is zero velocity on the embedded boundaries. */

  u.n[embed] = dirichlet(0);
  u.t[embed] = dirichlet(0);


  /**
  We use "third-order" [face flux interpolation](/src/embed.h). */

  for (scalar s in {u})
    s.third = true;


  // Initialize velocity field
  foreach()
    u.x[] = cs[] > 0. ? Umax*(1. - sq(y)) : 0. ;

  foreach()
    un[] = u.x[];
}


event properties (i++)
{
  foreach_face()
  {
    muf.x[] = mu1*fs.x[];
  }
  foreach()
  if (cs[] > 0. && cs[] < 1.)
     u.x[] = 0.;

}

scalar f[];
scalar mut[];
vector gradf[];
vector gradff[];
vector dmut[];
face vector muft[];

event compute_gradient (i++)
{
   foreach()
     f[] = u.x[];
   foreach()
   {
     foreach_dimension()
        gradf.x[] = d_dy_4(f);
     foreach_dimension()
        gradff.x[] = dd_ddy_2(f);
     foreach_dimension()
        dmut.x[] = d_dy_4(mut);
    }
}


// iterative loop to maintain a constant flow rate
coord flag = {1.,0.,0.};
double Q0 = 1.8;
double Q;

event flowrate (t += 0.001) {

    FILE * fp1 = fopen("Qout", "a");
    Q=0;

    double dy=L0/N;
    for (double y = -1.; y <= 1.; y += dy)
        Q+=interpolate (u.x, Xplot, y);
     Q=Q*dy;

    fprintf (fp1, "%6.4g %g \n", t, Q);
    fclose (fp1);
    if (Q < Q0)
    {
       flag.x += 0.01;

    }
    if (Q > Q0)
    {
       flag.x -= 0.01;
    }
}


// Computing the acceleration vector to combine the pressure gradient and the turbulent term
event acceleration (i++)
{
   foreach()
      lm[] = y <= 0. ? k*(1.+ y) : k*(1.- y);

   foreach()
      mut[] = sq(lm[])*fabs(gradf.y[]);


   foreach_face()
   {
     muft.x[] = mut[]*fs.x[];
     dtau.y[] = (muft.x[]*face_value(gradff.y, 0) + face_value(dmut.y, 0)*face_value(gradf.y, 0));
     av.x[]   = flag.x * unit.x[] + dtau.y[];
   }
}

// Computing the turbulent parameters

event Tau (i++)
{
  double bin     = 0.;
  double count   = 0.;
  double Re_tau;

  FILE *fp2 = fopen("Tau", "a");


  foreach()
    if (cs[] > 0. && cs[] < 1.)
    {
      bin     += fabs(gradf.y[]);
      count   += 1.;
    }

  bin   /= count;
  Tau    = mu1*bin;
  Re_tau = H*sqrt(Tau)/mu1;
  if (Re_tau < Re_target)
      Re += 10.;
  else
      Re -= 10.;

  fprintf (fp2 ,"%+6.5e \t %+6.5e \t %+6.5e \t %+6.5e \t %+6.5e \n",t, bin ,Tau , Re_tau, Re);
  fclose (fp2);
}


scalar u_star[];
scalar y_star[];

event star (i++)
{

   double u_tau = sqrt(Tau);
   double delta_mu = mu1/u_tau;
   foreach()
   {
     u_star[] = y <= fabs(H) ? u.x[]/u_tau : 0.;
     y_star[] = y <= fabs(H) ? (y + 1.)/delta_mu : 0.;

   }
}




event profile  (t = end)
{
  printf ("\n");
  for (double y = -L0/2. ; y <= 0.0 ; y += 0.001)
      {
      fprintf (stdout,"%+6.5e \t %+6.5e \t %+6.5e \t %+6.5e \t %+6.5e \t %+6.5e \n", t, Xplot, y, interpolate(u.x, Xplot, y),
                           interpolate(y_star, Xplot, y),
                           interpolate(u_star, Xplot, y));
      }


}



/**
The mesh is adapted following the wavelet-based estimation for the
discretization error.
 */
event adapt (i++) {
  adapt_wavelet ({cs, u}, (double[]){0.001, 1e-3, 1e-3, 1e-3}, N_max,N_min);
}



event logfile (t += 0.01; i <= 1e9)
{

  double du = change (u.x, un);
  FILE *fp3 = fopen("Residus", "a");
  fprintf (fp3 , " %d %e \n", i , du);
  if (i > 0 && du < 1e-5)
    return 1; /* stop */
  fclose(fp3);

}


