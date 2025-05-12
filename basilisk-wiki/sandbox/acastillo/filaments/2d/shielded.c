
#include "navier-stokes/centered.h"
#define MAXLEVEL 8
#define MINLEVEL 4

int main(){
  L0 = 8;
  X0 = Y0 = -L0/2;
  init_grid (1 << MAXLEVEL);
  periodic(top);
  periodic(left);
  //double reynolds= 2500;
  //const face vector muc[] = {1./reynolds,1./reynolds};
  //mu = muc;
  run();
}

event adapt (i++) {
  adapt_wavelet ((scalar *){u}, (double[]){1e-6,1e-6}, MAXLEVEL, MINLEVEL);
}

#define RAD (sqrt(sq(x) + sq(y)))
event init (t = 0){
  scalar psi[], omega[];

  double alpha=2.0;
  foreach() {
    omega[]  = (1 - 0.2*0.5*alpha*pow(RAD, alpha)) * exp(-pow(RAD,alpha));
    psi[] = 0.;
  }
  boundary ({psi,omega});

  poisson (psi, omega);
  coord f = {-1.,1.};
  foreach()
    foreach_dimension()
      u.x[] = f.x*(psi[0,1] - psi[0,-1])/(2.*Delta);
  boundary ((scalar *){u});

  FILE * fp = fopen("vortex.asc", "w");
  fclose(fp);
}


#include "view.h"
event movie (t += 0.2) {
  scalar omega[];
  vorticity (u, omega);

  view (fov=25);
  squares ("omega", linear = false);
  box();
  save ("omega.mp4");
}

event output (t = 100) {
  scalar omega[];
  vorticity (u, omega);

  scalar psi[];
  boundary ({psi, omega});
	poisson (psi, omega);

  view (fov=25);
  squares ("omega", linear = false);
  isoline("psi", n=11);
  box();
  save ("omega.png");
}
