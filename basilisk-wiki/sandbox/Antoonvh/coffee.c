/**
![Some people prefer to mix some milk in their coffee. Image via [Antony Hall](http://www.antonyhall.net/).](https://i1.wp.com/antonyhall.net/images/coffee-cup-2015/4.jpg)

# Mixing milk with Coffee

A liquid layer with depth $h_0$ and density $\rho_1$ is hosted in a less dense atmosphere 
($\rho_2$) and kept in a mug-like geometry with radius $R$ by
the acceleration due to gravity ($g$). It is brought into a solid-body
rotation with an angular velocity $v_{\theta} = \omega$. At some point
the mug becomes stationary and due the the viscous boundary layers at
the no slip walls a secondary (Ekman) circulation will be
generated. This warrants to introduce both fluid's viscosities
($\mu_1$ and $\mu_2$), giving rise to the following dimensional groups,
with typical values for coffee mixing:

$$\frac{h_0}{R}= 0.5$$

$$\frac{\rho_1}{\rho_2} = 1000$$

$$\frac{\mu_1}{\mu_2} = 50$$

$$\frac{\omega^2 R}{g} = 0.2$$

$$\frac{\rho_1 \omega R^2}{\nu_1} = 15000$$

Using the approximate properties of water and air, $R=5\mathrm{cm}$, $\omega
= 2\pi\mathrm{Rad/s}$, $g = 9.81 \mathrm{m}^2/\mathrm{s}$ and a near-empty cup of coffee $h_0 = 2.5\mathrm{cm}$. For the dimensionless set up we
choose $g = 1$, $R = 1$ and set $\rho_1 = 1$, then: $h_0 = 0.5R$, $\omega = 0.45\sqrt{g/R}$,
$\nu_1 = 3\times 10^{-5}\rho_1 \omega R^2$, $\nu_2 = 6\times
10^{-7} \rho_1 \omega R^2$ and $\rho_2 = 1/1000 \rho_1$. Noting that we do not consider the
surface tension associated with the liquid-atmosphere interface
($\sigma$) an important parameter for the dynamics ($Bo
\rightarrow \infty$).

## Numerical Set up

We assume the swirling flow to be axisymmetric, this reduces the
computational costs considerably and makes it *easy* to implement the
circular mug boundaries. 

*/
#include "grid/multigrid.h" // I do not understand axi + adaptivity
#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "reduced.h"
#include "view.h"
#include "particles.h"
/**
## Swirl

By default the axisymmetric flow solver only considers $u_z$ and
$u_r$, we need to extend it with the angular velocity (in Rad/s)
$u_{\theta}$. We advect and diffuse the angular velocity field (`ut`),
naively assuming this works well with the axi-symmetric
metric. Furtheremore, we take into account the centrifugal force due
to the non-zero swirl.
*/
#include "tracer.h"
#include "diffusion.h"
scalar ut[];
scalar * tracers = {ut}; 
face vector av[];

event acceleration(i++){
  coord f = {0, 1};
  foreach_face(y){
    if (y > Delta/2.)
      av.y[] += f.y*y*sq((ut[0, -1] + ut[])/2.);
    else
      av.y[] += 0.;
  }
  DT = (min(DT*1.01, 0.1));
}

event tracer_diffusion(i++)
  diffusion(ut, dt, mu);

#define BC (0.45*(t < tf))
double tf = 1;
uf.n[top] = 0;// This solved a prominent mass-conservation issue
ut[top] = dirichlet(BC);
ut[left] = dirichlet(BC);

u.t[left] = dirichlet(0.);
u.t[top] = dirichlet(0.);
f[bottom] = neumann(0.);

/**
We initialize and set some more properties before we start the
`run()`.
 */
int main(){
  Y0 = 0.;
  init_grid(128);
  a = av;
  mu1 = 3E-5;
  mu2 = 6E-7;
  G.x = -1;
  G.y = 0.;
  rho1 = 1.;
  rho2 = 1./1000.;
  run();
}
/**
## Mixing of milk

In order to check if the aforementioned secondary circularion is an
efficient mixing mechanism we add tracer particles to the flow.
 */
event init( t = 0){
  init_particles_2D_square_grid(50, 0.1, 0.10, 0.1); 
  int j = 0;
  while (j < n_part)
    loc[j++].z = 0; // This is not axisymmetric
  fraction(f, (0.5 + 0.5*sq(BC)*sq(y)) - x);
  TOLERANCE = 10E-5;
  DT = 0.001;
  foreach()
    ut[] = BC;
  boundary(all);
}
/**
Apart from the particle advection in the *r* and *z* direction, the $\theta$
direction needs to be threated as well.
 */
event advance_particles(i++){
  int j = 0;
  while (j < n_part){
    loc[j].z += dt*interpolate(ut, loc[j].x, loc[j].y);
    j++;
  }
}
/**
## Generating the visialization

If we wish to view the paths of the tracers, some
extra care is required not to let the 2D Basilisk view functionalities limit our
artistic freedom.
 */

struct _scatter {
  coord * loc;    // Coordinates
  float s, pc[3], coefs[3]; // point size, colour and distance attenuation coefs.
};


trace
void scatter3Din2D(struct _scatter p){// Modified from scatter.h
  bview * view = draw();
  if (!p.coefs[0]){ // A guess:
    p.coefs[0] = 0.01;
    p.coefs[1] = 0.2;
    p.coefs[2] = 0.5;
  }
  if (p.pc)
    glColor3f(p.pc[0], p.pc[1], p.pc[2]);
  if (!p.s)
    p.s = 10;
  glPointSize(p.s);
  glBegin (GL_POINTS);
  for (int j = 0; j < n_part; j++)
    glvertex3d (view, p.loc[j].x, p.loc[j].y*cos(loc[j].z), p.loc[j].y*sin(loc[j].z));
  glEnd();
  view->ni++; 
}
/**
This is the movie-making event.
*/
event movie( t += 0.05; t <= 100){
  view(phi = 0.5);
  scalar omega[];
  vorticity(u, omega);
  view (fov = 30, quat = {-0.170586,0.0694095,-0.696274,0.693747},
	tx = -0.0209548, ty = -0.204598, width = 650, height = 400);
  draw_vof ("f", lw=6);
  box();
  scatter3Din2D(loc = loc);
  mirror ({0,-1}){
    draw_vof ("f", lw=6);
    box(notics = true);
  }
  if (fabs(t-tf) < 1)
    draw_string("Switch BCs", size = 25);
  save ("mov.mp4");
}

/**
## Movie 

Remember that the flowfield $\overrightarrow{u}$ is 2D: ($\overrightarrow{u}(r,z)$). 

![The bounding domain (box), the free surface (black line) and flow tracers (dots)](coffee/mov.mp4)

It appears that the tracers are well mixed by the secondary circulation after a few 
rotations, a nice example of laminar mixing. Consider using this mechanism next time you wish to mix milk in
your coffee.
 */
