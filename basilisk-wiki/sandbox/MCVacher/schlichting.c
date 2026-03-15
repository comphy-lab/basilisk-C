/**
# Schlichting 2D Jet

We try to recover the famous problem of the laminar 2D jet in an immersed domain, also called the Schlichting jet or the Bickley jet. This analytical solution can be found in H. Schlichting's book "Boundary layer theory" (McGraw-Hill) in section IX.f. It was derived by H. Schlichting (1933) himself and W. Bickley (1939).

It is a monophasic problem depending only on the Reynolds number $Re$. Here we take $Re = 10$.

*/

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "navier-stokes/perfs.h"
#include "view.h"
#define BVIEW 1


scalar s[];
scalar * tracers = {s};

double U0;
double R_d;

double Re;

FILE * fpmax; 

face vector muv[];

/**
We choose to impose left-right symmetry by simulating only one half of the domain (using slip boundary conditions at the right). Left and top boundaries are outflow conditions, and imposed flow at the bottom. 
*/

int main() {

  Re=10;

  R_d=0.01; //Initial radius for the jet
  L0=1;  
  U0=1;

  TOLERANCE = 1e-3 [*];

  u.n[bottom] = dirichlet (U0*(x > L0/2-R_d));
  u.t[bottom] = dirichlet(0.);

  u.n[top] = u.n[] > 0. ? neumann(0) : dirichlet(0);
  p[top] = dirichlet(0.);

  s[bottom] = dirichlet (U0*(x > L0/2-R_d));

  u.n[left] = neumann(0);
  p[left] = dirichlet(0);

  u.n[right] = dirichlet(0.);
  p[right]=neumann(0.);

  N=256;
  origin (-L0/2, 0);
  init_grid(N);
  
  char param_dim[80];
  sprintf (param_dim, "param_dim.txt");
  FILE * fparam = fopen(param_dim, "w");
  fprintf (fparam, "%g %g %g %d\n",L0,R_d,U0,N);
  fclose (fparam);

  mu=muv;

  fpmax =  fopen("log.dat", "w"); 

  run();
}

/**
We choose the viscosity so that it matches the Re.
*/

event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]*U0*R_d/Re;
}


/**
At each iteration, we calculate residuals to see if the solution has converged. They are saved in "log.dat".
*/

double sum_u = 0.;
double sum_u_t = 0.;
double res_u = 0.;
double sum_p = 0.;
double sum_p_t = 0.;
double res_p = 0.;

event logfile (i++) { 
  sum_u_t=0.;
  sum_p_t=0.;

  foreach_face(){
    sum_u_t += sqrt(sq(u.x[]) + sq(u.y[]));
  }
  foreach(){
    sum_p_t += p[];
  }
  res_u = fabs(sum_u_t - sum_u);
  res_p = fabs(sum_p_t - sum_p);

  sum_u=sum_u_t;
  sum_p=sum_p_t;

  fprintf (stderr, "%d %g %g %g\n", i, t, res_u,res_p); 
  fprintf (fpmax, "%d %g %g %g\n", i, t, res_u,res_p);
}

/**
~~~pythonplot Residuals
import matplotlib.pyplot as plt
import numpy as np

list_t=[]
list_res_u=[]
list_res_p=[]

with open('log.dat', encoding='utf8') as File:
    for line in File:
        if line.isspace()==0:
            temp=line.split()
            t=float(temp[1])
            list_t.append(t)
            norm_u=float(temp[2])
            list_res_u.append(norm_u)
            p=float(temp[3])
            list_res_p.append(p)
    File.close()

plt.figure()
plt.semilogy(list_t,list_res_p,'g-',label= r'Res. $p$')
plt.semilogy(list_t,list_res_u,'r-',label= r'Res. $|u|$')
plt.legend()
plt.savefig('residuals.png')
~~~

Residuals are converged around t=800.
*/

/**
We plot the vertical velocity as well as a passive tracer to observe the establishment of the flow.
*/

event ppm_output (t = 0; t += 0.05; t <= 200) {
  
    char name1[80];
    sprintf (name1, "uY.mp4");
    output_ppm (u.y, file = name1, n = 512, min = -U0, max = +U0, linear = true);

    // optionally tracer
    char name2[80];
    sprintf (name2, "s.mp4");
    output_ppm (s, file = name2, n = 512, min = 0., max = U0, linear = true);
}

event movie_streamlines (t += 0.05; t <= 200)
{
  scalar omega[];
  vertex scalar stream[];
  scalar psi[];

  vorticity(u, omega);

  psi[bottom] = dirichlet (0.);
  psi[top]    = neumann (0.);
  psi[left]   = neumann (0.);
  psi[right]  = dirichlet (0.);

  poisson (psi, omega);
  boundary ({psi});

  foreach_vertex()
    stream[] = (psi[0,-1] + psi[-1,-1] + psi[] + psi[-1])/4.;

  clear();
  view(
    width = L0/2,
    height = L0/2,
    tx     = 0,       // center horizontally
    ty     = -1.5*L0/4     // shift up so bottom aligns with y=0
);
  isoline("stream", n = 30);
  box();
  save ("streamlines.mp4");
}

/**
![Passive tracer](schlichting/s.mp4)
![Vertical velocity](schlichting/uY.mp4)
![Streamlines](schlichting/streamlines.mp4)

We can clearly see that the Neumann condition at the top for the velocity creates a slow down of the flow, because the condition does not respect the flow solution : *equations to be written soon*.

*/

/**
We save the different fields at the last iteration.
*/

event profile (t = end) {
  char name[80];
  sprintf (name, "res_end.txt");
  FILE * fpres = fopen(name, "w");
  foreach()
    fprintf (fpres, "%g %g %g %g %g %g \n", x, y, u.x[], u.y[], p[],s[]);
  fclose(fpres);
  
  printf ("-----END-----\n");
}

/**
## Comparison with theory

We consider the boundary layer equations (also called Prandtl model) :

$$\frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} = 0\\
        u\frac{\partial v}{\partial x} + v\frac{\partial v}{\partial y} = \nu\frac{\partial^{2}v}{\partial x^{2}} $$

The flow has to recover zero velocity at infinity. We write this set of equation with the stream function as varaible:

$$u(x,y) = \frac{\partial \psi}{\partial y} \\ v(x, y) = -\frac{\partial \psi}{\partial x}$$

Mass conservation:

$$\frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} = \frac{\partial^{2} \psi}{\partial x \partial y} - \frac{\partial^{2} \psi}{\partial y \partial x} = 0$$

Momentum:

$$-\frac{\partial\psi}{\partial y}\frac{\partial^{2} \psi}{\partial x^{2}} + \frac{\partial\psi}{\partial x}\frac{\partial^{2} \psi}{\partial y \partial x} = -\nu \frac{\partial^{3}\psi}{\partial x^{3}}$$

We look for a self-similar solution for the stream function $\psi$ that would be satisified far from the nozzle: 

$$\psi(x,y) = y^{p}F(\frac{x}{y^{p}})$$

To find $p$ and $q$, we refer to the two following conditions. First, the total y-momentum fluw has to be constant (see Schlichting's book for details):

$$ J = \rho \int_{-\infty}^{\infty} v^{2}\mathrm{d}x = cst~~~~(=\mathrm{2R_dU_0^{2}~in~this~case})$$
$$\Rightarrow y^{2(p-q)+q } \sim 1\\ \Rightarrow 2p - q = 0$$

Second, inertial is balanced by viscous forces:

$$v\frac{\partial v}{\partial y} \sim \nu\frac{\partial^{2}v}{\partial x^{2}} $$

$$\Rightarrow p+q = 1$$

Hence:

$$\psi(x,y) = y^{1/3}F(\frac{x}{y^{2/3}})$$

This is equivalent to: 

$$\psi(x,y) = \nu^{1/2}y^{1/3}F(\frac{1}{3\nu^{1/2}}\eta)$$

with $\eta = \frac{x}{y^{2/3}}$. Using this into the previous momentum equation, it gives that F is solution of the following ODE:

$$ (F')^{2} + FF'' + F''' = 0$$

It can be shown that the solution of this ODE is: 

$$ F'(\frac{1}{3\nu^{1/2}}\eta) = 2\alpha \times \mathrm{sech^{2}}(\alpha\frac{1}{3\nu^{1/2}} \eta)$$

where $\alpha$ is deduced from $J$ (see Schlichting's book). In the end, we can write :

$$ y^{1/3} \times v(x,y) = A~\mathrm{sech^{2}}(B \eta) $$

with $A$ and $B$ constants.

This is verified just below for most of the positions, far enough from the boundaries.                                    
*/



/**
~~~pythonplot Self-similarity check
import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# Load numerical data
# ============================================================
data = np.loadtxt("res_end.txt")

x = data[:, 0]
y = data[:, 1]
u_y = data[:, 3]

# ============================================================
# Simulation parameters
# ============================================================

Re = 10
V_0 = 1
rho = 1
R_d = 0.01
nu = V_0 * R_d / Re
D = 2 * R_d

# ============================================================
# Define usable y/D range (avoid inlet & boundary)
# ============================================================

yD = y / D

yD_min = 10    # avoid inlet
yD_max = 40     # avoid external boundary

mask = (yD >= yD_min) & (yD <= yD_max)

x = x[mask]
y = y[mask]
u_y = u_y[mask]
yD = yD[mask]

# ============================================================
# Build structured 2D matrix U_y[y,x]
# ============================================================
x_unique = np.unique(x)
y_unique = np.unique(y)

nx = len(x_unique)
ny = len(y_unique)

U_y = np.zeros((ny, nx))

for i, yi in enumerate(y_unique):
    for j, xj in enumerate(x_unique):
        idx = np.where((x == xj) & (y == yi))[0]
        if len(idx) > 0:
            U_y[i, j] = u_y[idx[0]]
        else:
            U_y[i, j] = np.nan

# ============================================================
# Self-similar collapse from numerical data
# ============================================================

eta_all = []
f_all = []
yD_all = []

for i, yi in enumerate(y_unique):

    if yi <= 0:
        continue

    v_slice = U_y[i, :]

    if np.max(np.abs(v_slice)) < 1e-10:
        continue

    # shifted x
    x_shift = x_unique - 1/2

    # similarity variable
    eta = x_shift / (yi**(2/3))

    # rescaled velocity (multiply by y^(2/3))
    f = v_slice * (yi**(1/3))

    eta_all.extend(eta)
    f_all.extend(f)
    yD_all.extend([yi/D]*len(eta))

eta_all = np.array(eta_all)
f_all = np.array(f_all)
yD_all = np.array(yD_all)

# ============================================================
# Theoretical profile: A y^{-1/3} sech²(B x / y^{2/3})
# In similarity variables:
# y^{1/3} v = A sech²(B eta)
# ============================================================

J=rho*V_0**2*0.01
alpha=0.8255*(J/(rho*nu**(1/2)))**(1/3)

A=2/3*alpha**2 #Theory
B=alpha*1/(3*nu**(1/2)) #Theory

A = 0.3 #best fit according to me ^^
B = 7.1 #best fit according to me ^^

eta_theory = np.linspace(-1, 1, 500)
f_theory = A * (1 / np.cosh(B * eta_theory))**2

# ============================================================
# Gaussian comparison
# ============================================================

sigma = 0.12
f_gauss = A * np.exp(-(eta_theory**2) / (2 * sigma**2))

# ============================================================
# Plot
# ============================================================

plt.figure(figsize=(7,6))

sc = plt.scatter(
    eta_all,
    f_all,
    c=yD_all,
    cmap='viridis',
    s=5,
    alpha=0.6,
)

# symmetry
plt.scatter(
    -eta_all,
    f_all,
    c=yD_all,
    cmap='viridis',
    vmin=yD_min,
    vmax=yD_max,
    s=6,
    alpha=0.7
)

plt.plot(
    eta_theory,
    f_theory,
    'r',
    linewidth=2,
    label=r'Theory: $A\,\mathrm{sech}^2(B\eta)$'
)

plt.plot(
    eta_theory,
    f_gauss,
    'k--',
    linewidth=2,
    label='Gaussian'
)

plt.xlabel(r'$\eta = \frac{x}{y^{2/3}}$')
plt.ylabel(r'$y^{1/3} v(x,y)$')
plt.xlim([-0.5,0.5])
plt.ylim([-0.025,0.4])
plt.legend()
plt.grid(True)

cbar = plt.colorbar(sc)
cbar.set_label(r'$\mathrm{Vertical~position~~~y/D}$')

plt.tight_layout()
plt.savefig('self_sim.png')
~~~
*/

/**
## See also

See also the axisymmetric case [here](/sandbox/MCVacher/schlichting_axi.c).

*/