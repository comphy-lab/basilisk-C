/**
# "U" Shaped flow around a cylinder

##Equation on the stream function

In this section, the goal is to create a particular flow around a cylinder. The input flow is located jsut beneath the cylinder, and the output flow is located above the cylinder. 

Here we chose to resolve the Poisson equation for the stream function $\psi$. 

Just like [Potential flow around a half cylinder](http://basilisk.fr/sandbox/mgarnier/PotentialCylinder/HalfCylinder.c), we just shift the origin to the right in order to obtain a half cylinder. Globaly, the code remains the same as my code for [Potential flow around a half cylinder](http://basilisk.fr/sandbox/mgarnier/PotentialCylinder/HalfCylinder.c). We only modify the boundary conditions. 
*/

#include "embed.h"
#include "poisson.h"
#include "view.h"

double R = 1.; //Rayon du bout du biseau
double Vitesse = 1.;

static
double dirichlet_homogeneous_bc (Point point, Point neighbor,
				 scalar s, void * data) {
  return 0.;
}

int main()
{
  for (N = 512; N <= 512; N *= 2) {
    size (6. [0]); // dimensionless
    origin (-L0, -L0/2);
    init_grid (N);

/**
The domain of computation remains the same.
*/
    
    vertex scalar phi[];
    foreach_vertex() {
     double r = sqrt (sq(x) + sq(y));

    phi[] = r - R ;

    }
    fractions (phi, cs, fs);  
#if TREE
    cs.refine = cs.prolongation = fraction_refine;
#endif
    restriction ({cs,fs});

    cm = cs;
    fm = fs;
    
/**
##Boudary conditions

Let's talk about the boundary conditions. We choose the following ones : 

#####Top, left, bottom boundaries

 On the left, top, bottom (i.e wall) the normal velocity is null. If note the velocity $\vec{U} = \left( u , v \right)$, we have for the top and bottom boundary : 
 
$$
v = 0 = -\frac{\partial \psi}{\partial x} = \frac{\partial \phi}{\partial y}
$$

So $\psi$ is steady on the top and bottom boundary. Note that this condition is equivalent to $\frac{\partial \phi}{\partial \vec{n}}=0$

On the left boundary we have : 

$$
u = 0 = \frac{\partial \psi}{\partial y} = \frac{\partial \phi}{\partial x}
$$

Like before, $\psi$ is steady on the top and bottom boundary. Note that his condition is equivalent to $\frac{\partial \phi}{\partial \vec{n}}=0$

Because the stream function is continous, the stream function is equal to the same value on the left, top and bottom boundaries. In our simulation, we impose $\psi = 0$ ont theses boundaries. 

#####Right boundary

On this boundary we need to impose the input and output flow. For the input flow (beneath the cylinder), we have $u=-C$ on the right boundary (where $C$ is a positive constant). So, on the boundary : 

$$
u = -C = \frac{\partial \psi_{input}}{\partial y}
$$

then, 

$$
\psi_{input} \left( y \right) = -C y + B 
$$

We can simplifiy this expression because the stream function is continuous  $\psi \left( y = -L/2 \right) = 0$ with L the size of the domain

$$
\psi_{input} \left( y \right) = -C( y + L/2)
$$
\
\

We can do the same for the output flow (above the cylinder). This time $u=C$. Following the same steps as for the input flow we obtain : 

$$
\psi_{output} \left( y \right) = C y + D 
$$

With the condition of continuity $\psi_{output} \left( y = L/2 \right) = 0$. We arrive to the following expression. 

$$
\psi_{output} \left( y \right) = C \left( y - L/2 \right) 
$$

Note that for the other cases where the solid is not a simple cylinder, or where the origin is not in the middle theses expressions are not valid. 

#####Embedded boundary

For the embedded boundary, we suppose that the normal velocity is null. It means that on the solid (the cylinder) $\psi$ is steady. With the conditions of continuity, we know that : 

$$
\psi_{embedded} = \psi_{input} \left( y = -R \right) = \psi_{output} \left( y = R \right) = -C \left( R - L/2 \right)
$$

#####Note on the code

Just below this section you can find the code to implement the boundary conditions. Note that if you use the dirichlet function, you must not use the 
"dirichlet_homogeneous_bc" function. Otherwise the computation will not be stable. 
*/
   
    scalar psip[], b[];

psip[left] = dirichlet(0.);
//psip.boundary_homogeneous[left] = dirichlet_homogeneous_bc;
psip[right] = (y <= -R) ? dirichlet(-Vitesse*(y+3)) : (y >= R) ? dirichlet(Vitesse*(y-3)) : dirichlet(0.); 
//psip.boundary_homogeneous[right] = dirichlet_homogeneous_bc;
psip[top] = dirichlet(0.);
//psip.boundary_homogeneous[top] = dirichlet_homogeneous_bc;
psip[bottom] = dirichlet(0.);
//psip.boundary_homogeneous[bottom] = dirichlet_homogeneous_bc;
psip[embed] = dirichlet(Vitesse*(R-3));
//psip.boundary_homogeneous[embed] = dirichlet_homogeneous_bc;

foreach(){
  if(cs[]==0)
  psip[] = nodata;
}
    
    foreach() {
      double xc = x, yc = y;
      if (cs[] > 0. && cs[] < 1.) {
	coord n = facet_normal (point, cs, fs), p;
	double alpha = plane_alpha (cs[], n);
	line_center (n, alpha, cs[], &p);
	xc += p.x*Delta, yc += p.y*Delta;
      }
      b[] = 0;
    }

    /**
    The Poisson equation is solved. */
    
    struct Poisson p;
    p.alpha = fs;
    p.lambda = zeroc;
    p.embed_flux = embed_flux;
    scalar resdir[];
    double maxpdir = residual ({psip}, {b}, {resdir}, &p), maxfdir = 0.;
    foreach()
      if (cs[] == 1. && fabs(resdir[]) > maxfdir)
	maxfdir = fabs(resdir[]);
    fprintf (stderr, "maxres %d %.3g %.3g\n", N, maxfdir, maxpdir);

    // FIXME: need to set minlevel to 4
    timer t = timer_start();
    mgstats s = poisson (psip, b, alpha = fs, tolerance = 1e-6, minlevel = 4);
    double dt = timer_elapsed (t);
    printf ("%d %g %d %d\n", N, dt, s.i, s.nrelax);

    FILE * fp1 = fopen ("U_cyl.dat","w");
    foreach() 
      fprintf(fp1,"%g %g %g \n", x, y, psip[]);
    fflush (fp1);

    
/**
##Equation on the potential function

Another possibility is to resolve the Poisson equation on the potential function $\phi$. For this case, we need to modify the boundary conditions. I took the code of [Quasi 2D incompressible flow in a recorder window](http://basilisk.fr/sandbox/hajczak/recorderWindow.c) by  A.Hajczak. I replaced his geometry by mine. The goal is to compare the stream function obtained by resolving the Poisson on the stream function and the stream function obtained by resolving the Poisson equation on the potential. 

*/
    
scalar phi1[]; 


phi1[bottom] = neumann(0.);
phi1[left] = neumann(0.);
phi1[top] = neumann(0.);
phi1[right] = (y < -R) ? -Vitesse : (y > R) ? Vitesse : 0.; 
phi1.boundary_homogeneous[right] = dirichlet_homogeneous_bc; 
phi1[embed] = neumann(0.);


phi1.third = true;


foreach(){
  if(cs[]==0)
  phi1[] = nodata;
}


  struct Poisson p1;
  p1.alpha = fs;
  p1.lambda = zeroc;
  p1.embed_flux = embed_flux;
  scalar res1[]; 
  double maxp1 = residual ({phi1}, {b}, {res1}, &p1), maxf1 = 0.;  
  foreach()  
    if (cs[] == 1. && fabs(res1[]) > maxf1)  
      maxf1 = fabs(res1[]);  
  fprintf (stderr, "maxres %d %.3g %.3g\n", N, maxf1, maxp1);  
  timer t1 = timer_start();
  mgstats s1 = poisson (phi1, b, alpha = fs, tolerance = 1e-4, minlevel = 4);
  double dt1 = timer_elapsed (t1);
  printf ("%d %g %d %d\n", N, dt1, s1.i, s1.nrelax);

  vector u[];
  gradients ({phi1}, {u});
  
  scalar psi[];
  psi[bottom] = 0.;
  psi[top] = 0.;
  
  foreach()
    {
      psi[] = -u.y[]*Delta + psi[-1,0];
    }
  
  foreach(){
  if(cs[]==0)
  psi[] = nodata;
}

  FILE * fp = fopen ("streamlines.dat","w");
  foreach() 
    fprintf(fp,"%g %g %g %g %g \n", x, y, u.x[], u.y[], psi[]);
  fflush (fp);

  }
}





/**

With the figure below we can compare the two solutions. 

~~~pythonplot

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

L = 6
l = L / 2
R = 1

data_U = np.loadtxt('U_cyl.dat')
data_streamline = np.loadtxt('streamlines.dat')

x_U = data_U[:, 0]
y_U = data_U[:, 1]
psi_U = data_U[:, 2]

x_streamline = data_streamline[:, 0]
y_streamline = data_streamline[:, 1]
psi_streamline = data_streamline[:, 4]

xi_U = np.linspace(x_U.min(), x_U.max(), 5000)
yi_U = np.linspace(y_U.min(), y_U.max(), 5000)
xi_U, yi_U = np.meshgrid(xi_U, yi_U)
zi_U = griddata((x_U / l, y_U / l), psi_U, (xi_U, yi_U), method='linear')

xi_streamline = np.linspace(x_streamline.min(), x_streamline.max(), 5000)
yi_streamline = np.linspace(y_streamline.min(), y_streamline.max(), 5000)
xi_streamline, yi_streamline = np.meshgrid(xi_streamline, yi_streamline)
zi_streamline = griddata((x_streamline / l, y_streamline / l), psi_streamline, (xi_streamline, yi_streamline), method='linear')

fig, ax = plt.subplots(figsize=(10, 5))

contour_U = ax.contour(xi_U, yi_U, zi_U, levels=np.arange(-2.02, 0, 0.25), linewidths=2., colors='r', linestyles='solid')
contour_streamline = ax.contour(xi_streamline, yi_streamline, zi_streamline, levels=np.arange(-2.02, 0, 0.25), linewidths=2., colors='b', linestyles='dashed')

theta = np.linspace(0, 2 * np.pi, 100)
x_circ = R * np.cos(theta)
y_circ = R * np.sin(theta)
ax.fill(x_circ / l, y_circ / l, color='black')

ax.axis('square')
ax.set_xlim(left=-1.75, right=0)
ax.set_ylim(bottom=-0.875, top=0.875)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('Comparative Isolines of Stream Functions from different computations')

stream_line, _ = contour_U.legend_elements()
potential_line, _ = contour_streamline.legend_elements()
ax.legend([stream_line[0], potential_line[0]], ['Stream equation', 'Potential equation'], loc='upper left')

plt.tight_layout()
plt.savefig('comparison_isolines.svg')
~~~

*/ 
