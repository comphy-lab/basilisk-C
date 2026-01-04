/**
# Gouy-Chapman Debye layer

The [Debye
layer](http://en.wikipedia.org/wiki/Double_layer_%28interfacial%29) is
the ionic concentration and potential distribution structure that
appears on the surface of a charged electrode in contact with solvents
in which are dissolved ionic species.  Louis Georges Gouy and David
Leonard Chapman at the beginning of the XX century proposed a model of
the Debye layer resulting from the combined effect of its thermal
diffusion and its electrostatic attraction or repulsion.  In effect,
in a stationary situation and assuming fluid at rest, the
Poisson-Nernst-Planck equations are,

$$
0 = \nabla \cdot (e \omega_i Z_i c_i \nabla \phi) + \nabla \cdot 
(\omega_i k_B T \nabla c_i) \quad \mathrm{with} \quad 
\nabla \cdot (\epsilon \nabla \phi) = \sum_i e c_i
$$

where $\phi$ is the electric potential and $c_i$ is the number of
$i$-ions per volume. $\omega_i$ and $Z_i$ are the $i$-ion mobility and
valence.  $k_B$ is the Boltzmann constant, $e$ is the electron charge,
$\epsilon$ the electrical permittivity and $T$ the temperature.

The above equations, written in dimensionless form, reduces in the
case of a fully dissolved binary system in a planar geometry to,

$$
\hat{c}_+ = exp (-\hat{\phi}), \, \hat{c}_- = exp (\hat{\phi}) 
\quad \mathrm{with} \quad (\hat{\phi})_{xx} = 2 \sinh (\hat{\phi}).   
$$
*/

#include "grid/multigrid1D.h"
#include "diffusion.h"
#include "run.h"
#include "ehd/pnp.h"

#define Volt 1.0
#define DT 0.01

/**
We assume a fully dissolved binary system labelling the positive ion as $Cp$ 
and the counterion as $Cm$. The valence is one, ($|Z|=1$). */

scalar phi[];
scalar Cp[], Cm[];
int Z[2] = {1,-1};
scalar * sp = {Cp, Cm};

/**
Ions are repelled by the electrode due to its positive volume
conductivity while counterions are attracted (negative
conductivity). */

#if 1
const face vector kp[] = {1., 1.};
const face vector km[] = {-1., -1.};
vector * K = {kp, km};
#endif

/**
On the left the charged planar electrode is set to a constant
potential $\phi =1$. The concentrations of the positive and negative
ions depend exponentially on the voltage electrode. */

phi[left] = dirichlet(Volt);
Cp[left]  = dirichlet (exp(-Volt));
Cm[left]  = dirichlet (exp(Volt));

/**
In the bulk of the liquid, on the right boundary, the electrical
potential is zero and the ion concentrations match the bulk
concentration i.e */

phi[right] = dirichlet (0.);
Cp[right]  = dirichlet (1.);
Cm[right]  = dirichlet (1.);

/**
Initially, we set the ion concentration to their bulk values together
with a linear decay of the electric potential $\phi$. */
 
event init (i = 0)
{
  foreach() {
    phi[] = Volt*(1.-x/5.);
    Cp[] = 1.0; 
    Cm[] = 1.0;
  }
}

event integration (i++) {
  dt = dtnext (DT);

  /**
  At each instant, the concentration of each species is updated taking into
  account the ohmic transport. */

#if 1
  ohmic_flux (sp, Z, dt, K);
#else
  ohmic_flux (sp, Z, dt); // fixme: this does not work yet
#endif

  /**
  Then, the thermal diffusion is taken into account. */

  for (scalar s in sp)
    diffusion (s, dt);

  /**
  The electric potential $\phi$ has to be re-calculated since the net
  bulk charge has changed. */

  scalar rhs[];
  foreach() {
    int i = 0;
    rhs[] = 0.;
    for (scalar s in sp)
      rhs[] -= Z[i++]*s[];
  }
  poisson (phi, rhs);
}

event result (t = 3.5) {
  foreach()
    fprintf (stderr, "%g %g %g %g \n", x, phi[], Cp[], Cm[]);
}

/**
## Results

We compare the numerical results (symbols) with the analytical
solution (lines).

~~~gnuplot Profiles of electric potential and concentrations
set xlabel 'x'
gamma = tanh(0.25)
fi(x) = 2*log((1+gamma*exp(-sqrt(2)*x))/(1-gamma*exp(-sqrt(2)*x)))
nplus(x) = exp(-fi(x))
nminus(x) = exp(fi(x))
plot 'log' u 1:2 notitle, fi(x) t '{/Symbol f}',\
     'log' u 1:3 notitle, nplus(x) t 'n+',\
     'log' u 1:4 notitle, nminus(x) t 'n-' lt 7
~~~
*/

int main() {
  N = 32;
  L0 = 5;
  TOLERANCE = 1e-4;
  run();
}
