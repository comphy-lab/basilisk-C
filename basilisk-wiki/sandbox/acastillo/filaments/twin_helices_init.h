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

void helical_pair_filaments(vector omega, char * name){

  double Rs, Hs, as, circulation, axialvel;

  // Step 1. Define the helix as the spacecurve c
  FILE * fp = fopen(name, "r");
  fscanf(fp, "%d",  &n_seg);
  fscanf(fp, "%lf", &Rs);
  fscanf(fp, "%lf", &as);
  fscanf(fp, "%lf", &Hs);
  fscanf(fp, "%lf", &circulation);
  fscanf(fp, "%lf", &axialvel);

  double t0[n_seg];
  coord c[n_seg], dc[n_seg], d2c[n_seg];
  for (int i = 0; i < n_seg; i++){
    fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &t0[i], &c[i].x, &c[i].y, &c[i].z, &dc[i].x, &dc[i].y, &dc[i].z, &d2c[i].x, &d2c[i].y, &d2c[i].z);
    c[i].z -= L0/2;
  }
  fclose (fp);

  view (camera="iso", fov=40);
  box();
  cells( alpha=-L0/2 );
  cells( n = {1,0,0} );
  draw_space_curve(n_seg, c);
  char subname[80];
  sprintf(subname, "%s.png", name);
  save (subname);

  // Step 2. Compute a Frenet-Serret Frame
  coord tvec[n_seg], nvec[n_seg], bvec[n_seg];
  for (int i = 0; i < n_seg; i++){
    foreach_dimension(){
      tvec[i].x =  dc[i].x/sqrt(vecdot( dc[i],  dc[i]));
      nvec[i].x = d2c[i].x/sqrt(vecdot(d2c[i], d2c[i]));
    }
    bvec[i] = vecdotproduct(tvec[i], nvec[i]);
  }

  // Steps 3. Initialize the vorticity field
  foreach(){
    coord pcar = {x,y,z};
    coord val_omega = get_vorticity_filament2(circulation, axialvel, pcar, n_seg, as/1.36, t0, c, tvec, nvec, bvec, 0);
    foreach_dimension()
      omega.x[] -= val_omega.x;
  }
}

event init (t = 0) {

  printf( "Initial conditions from filamentary solutions \n" );

  // Refine around the region of interest
  char name[80];
  sprintf (name, "./base_hext.asc");
  FILE * fp = fopen(name, "r");
  fscanf(fp, "%d",  &n_seg);
  fscanf(fp, "%lf", &Rs1);
  fscanf(fp, "%lf", &as);
  fscanf(fp, "%lf", &Hs1);
  fclose (fp);

  sprintf (name, "./base_hint.asc");
  fp = fopen(name, "r");
  fscanf(fp, "%d",  &n_seg);
  fscanf(fp, "%lf", &Rs2);
  fscanf(fp, "%lf", &as);
  fscanf(fp, "%lf", &Hs2);
  fclose (fp);

  if (MAXLEVEL <= 8){
    refine  (RAD < Rs1 + 16*as && RAD > Rs2 - 16*as && level < MAXLEVEL - 2);
    refine  (RAD < Rs1 +  8*as && RAD > Rs2 -  8*as && level < MAXLEVEL - 1);
    refine  (RAD < Rs1 +  4*as && RAD > Rs2 -  4*as && level < MAXLEVEL);
  }
  else {
#ifdef _BETA
    sprintf (name, "./case%d/level%d/0_init/dump", _BETA, MAXLEVEL-1);
    printf(name);
    restore(name);
#else
    restore("dump");
#endif
    adapt_wavelet ((scalar*){u}, (double[]){1e-3, 1e-3, 1e-3}, MAXLEVEL, MINLEVEL);
  }

  // Initialize the vorticity fields
  vector omega[];
  foreach()
    foreach_dimension()
      omega.x[] = 0.;

  helical_pair_filaments(omega, "base_hext.asc");
  helical_pair_filaments(omega, "base_hint.asc");

  view (camera="iso");
  box();
  cells();
  squares ("omega.x", linear = false, n = {1,0,0} );
  save ("twin_helices0.png");

  // Solve a Poisson problem for each vorticity component
  vector psi[];
  foreach()
    foreach_dimension()
      psi.x[] = 0.;
  boundary ((scalar*){psi,omega});
  poisson (psi.x, omega.x);
  poisson (psi.y, omega.y);
  poisson (psi.z, omega.z);

  // Compute the velocity from the vector potential field
  foreach(){
    u.x[] = ((psi.z[0,1,0] - psi.z[0,-1,0]) - (psi.y[0,0,1] - psi.y[0,0,-1]))/(2.*Delta);
    u.y[] = ((psi.x[0,0,1] - psi.x[0,0,-1]) - (psi.z[1,0,0] - psi.z[-1,0,0]))/(2.*Delta);
    u.z[] = ((psi.y[1,0,0] - psi.y[-1,0,0]) - (psi.x[0,1,0] - psi.x[0,-1,0]))/(2.*Delta) - 1.0;
  }
  boundary ((scalar *){u});

  // Refine the grid
  adapt_wavelet ((scalar*){u}, (double[]){1e-3, 1e-3, 1e-3}, MAXLEVEL, MINLEVEL);

  scalar l2[];
  lambda2 (u, l2);

  view (camera="iso");
  box();
  cells();
  squares ("u.x", linear = false, n = {1,0,0} );
  isosurface ("l2", 0);
  save ("twin_helices1.png");

  dump (list = (scalar *){u});
}
