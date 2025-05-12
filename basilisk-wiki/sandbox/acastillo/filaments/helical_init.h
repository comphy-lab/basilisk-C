/**
## Initialisation
For the initialisation process, we proceed as follows.

1. The position of the vortex segments is given by a space-curve $c$
parametrized as function of t0.

2. We require the first and second derivatives
of $c$, as to define a local Serret-Frenet frame $(\vec{t},\vec{n},\vec{b})$.
$$
\vec{t} = \frac{d\vec{c}}{dt_0}/\|\frac{d\vec{c}}{dt_0}\|
$$
$$
\vec{n} = \frac{d^2\vec{c}}{dt_0^2}/\|\frac{d^2\vec{c}}{dt_0^2}\|
$$
$$
\vec{b} = \vec{t} \times \vec{n}
$$

3. Each position $\vec{p}$ is projected into the local Frenet-Serret frame
to obtain a set of local coordinates, such that:
$$
(\vec{p} - \vec{c}(t_0)) \cdot \vec{t} = 0
$$
$$
(\vec{p} - \vec{c}(t_0)) \cdot \vec{n} = x_n
$$
$$
(\vec{p} - \vec{c}(t_0)) \cdot \vec{b} = x_b
$$
which is done through a minization process.

4. We use the local coordinates $(x_n, x_b)$ to define a radial coordinate
$\rho$ required to compute the vorticity of a Lamb-Oseen vortex as
$$
\vec{\omega} = \Gamma/(\pi a^2) \exp(-\rho^2/a^2) \cdot \vec{t}
$$
where $\Gamma$ is the circulation and $a$ the core size.

5. We consider the flow is expressed as a vector potential,
$$
\vec{u} = \nabla \times \vec{\psi}
$$
such that
$$
\nabla^2 \vec{\psi} = -\vec{\omega}
$$
which requires solving a Poisson problem for each vorticity component.

6. We use the vector potential to evaluate the velocity field.
*/

#include "filaments.h"
#include "filaments.c"

void helical_filament(vector u){
  // Step 1. Define the helix as the spacecurve c
  double delta_t0 = (2 + n_turns)*2*pi/((double)n_seg-1);
  double t0[n_seg];
  coord c[n_seg], dc[n_seg], d2c[n_seg];

  for (int i = 0; i < n_seg; i++){
    t0[i] = delta_t0 * (double)i - 2*pi;
    c[i].x = (Rs + delta_init * sin(t0[i]/(2*pi*k_init))) * cos(t0[i]);
    c[i].y = (Rs + delta_init * sin(t0[i]/(2*pi*k_init))) * sin(t0[i]);
    c[i].z = (Hs/(2*pi)) * t0[i] - L0/2 + delta_init * sin(t0[i]/(2*pi*k_init));
  }

  view (camera="iso", fov=40);
  draw_space_curve(n_seg, c);
  box();
  save ("helical0.png");

  // Step 2. Find the first and second derivatives to compute a Frenet-Serret Frame
  coord xshift = {0, 0, (2 + n_turns)*Hs}, dxshift = {0, 0, 0};
  fd_derivative(n_seg, delta_t0,  xshift,  c,  dc);
  fd_derivative(n_seg, delta_t0, dxshift, dc, d2c);

  coord tvec[n_seg], nvec[n_seg], bvec[n_seg];
  for (int i = 0; i < n_seg; i++){
    foreach_dimension(){
      tvec[i].x =  dc[i].x/sqrt(vecdot( dc[i],  dc[i]));
      nvec[i].x = d2c[i].x/sqrt(vecdot(d2c[i], d2c[i]));
    }
    bvec[i] = vecdotproduct(tvec[i], nvec[i]);
  }

  // Steps 3 and 4. Initialize the vorticity field
  foreach(){
    coord pcar = {x,y,z};
    coord val_omega = get_vorticity_filament(pcar, n_seg, as, t0, c, tvec, nvec, bvec, 0);
    foreach_dimension()
      omega.x[] = -val_omega.x;
  }
  boundary ((scalar*){omega});

  // Step 5. Solve a Poisson problem for psi
  vector psi[];
  foreach()
    foreach_dimension()
      psi.x[] = 0.;
  boundary ((scalar*){psi});
  poisson (psi.x, omega.x);
  poisson (psi.y, omega.y);
  poisson (psi.z, omega.z);

  // Step 6. Compute the velocity field from the potential field
  foreach(){
    u.x[] = ((psi.z[0,1,0] - psi.z[0,-1,0]) - (psi.y[0,0,1] - psi.y[0,0,-1]))/(2.*Delta);
    u.y[] = ((psi.x[0,0,1] - psi.x[0,0,-1]) - (psi.z[1,0,0] - psi.z[-1,0,0]))/(2.*Delta);
    u.z[] = ((psi.y[1,0,0] - psi.y[-1,0,0]) - (psi.x[0,1,0] - psi.x[0,-1,0]))/(2.*Delta);
  }
  boundary ((scalar *){u});
  if (pid()==0)
    printf( " Velocity fields from filament solutions ... Done! \n" );
}

event init (t = 0) {

  printf( "Initial conditions from filamentary solutions \n" );
  refine  (RAD < Rs + 3*as   && RAD > Rs - 3*as   && level < MAXLEVEL - 1);
  refine  (RAD < Rs + 1.5*as && RAD > Rs - 1.5*as && level < MAXLEVEL);

  helical_filament(u);

  lambda2 (u, l2);

  view (camera="iso", fov=40);
  box();
  isosurface ("l2", 0);
  cells (n = {0, 0, 1});
  save ("helical1.png");

  adapt_wavelet ((scalar*){u}, (double[]){1e-4, 1e-4, 1e-4}, MAXLEVEL, MINLEVEL);
  printf( "Adaptative refinement ... Done! \n" );

  view (camera="iso", fov=40);
  box();
  isosurface ("l2", 0);
  cells (n = {0, 0, 1});
  save ("helical2.png");

  view (camera="front");
  squares ("omega.z", n = {0, 0, 1});
  save ("helical3.png");
}
