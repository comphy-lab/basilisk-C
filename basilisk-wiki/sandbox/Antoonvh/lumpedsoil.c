/**
# A Lumped-parameter soil underneath an atmosphere

![](lumpedsoil/mp2.mp4)

See more:

![Note the surface cooling tendency towards the asomptotic scale $b_{surf}= b_{\Lambda} = -0.4$](lumpedsoil/mov.mp4)

We prescribe a radiative loss term.
*/
#define QFLX (-0.005)
/**
In the lowest part of the domain we want to evaluate $b_{surf}$ from our solution. For that we adopt a linear-profile subgrid-scale model. The macro is to be evaluated in the cells that are at the `bottom` surface. 
*/
#define BSURF (1.5*b[] - 0.5*b[0, 1])
/**
Using this, we van evaluate the feedback flux with respect to a lumped parameter ($\Lambda$) and a reference buoyancy scale. 
*/
#define GFLX (-Lambda*(BSURF - bd))
double bd = 0, Lambda = 0.0125;
/**
We also assume an buoyancy profile that, at the surface, corresponds to the asymptotic buoyancy scale $b_{\Lambda}=Q/\Lambda + b_d$  
*/
#define STRAT (log(y + a1+ 1.) - log(a1 + 1.) + (QFLX/Lambda + bd))
double a1 = 1;
/**
We setup the the case wit the usual suspects.
*/
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"
#include "view.h"
#include "profile5c.h"

#define RAD (pow(pow((x - xo), 2)+(pow((y - yo), 2)), 0.5))
#define ST (-(x - xo)/RAD)

scalar b[];
scalar * tracers = {b};
b[top] = dirichlet(STRAT);
/** 
Since the ghost-cell values of the buoyancy at the `bottom` boundary are used for interpolations when making profiles and define refined-cell values, we set the boundary condition consistent with our subgrid-scale model.  
*/
b[bottom] = dirichlet(BSURF);

u.t[bottom] = dirichlet(0.);
double xo = 7.7, yo = 9.6;
double temp = 30;
double Re = 1000;
const face vector muc[] = {1/Re, 1/Re};

int main(){
  L0 = 15.;
  init_grid (1 << (8));
  foreach_dimension()
    u.x.refine = refine_linear;
  mu = muc;
  run();
}
/**
We open a pipeline for the varying-profile movie and initialize the buoyancy field. Furthermore we place a downward-targeted, self-advecting, vortex structure in the warm part of the domain.
*/
FILE * gnuplotPipe;
event init (t = 0){
  foreach()
    b[] = STRAT;
  refine (RAD < 2.0 && level <= 8);
  refine (RAD < 1.0 && level <= 9);
  scalar psi[];
  double k = 3.83170597;
  foreach() 
    psi[] = ((RAD > 1)*((1/RAD))*ST) + ((RAD < 1)*((-2*j1(k*RAD)*ST/(k*j0(k))) + (RAD*ST)));
  boundary({psi});
  foreach() {
    u.x[] = -((psi[0, 1] - psi[0, -1])/(2*Delta));
    u.y[] = (psi[1, 0] - psi[-1, 0])/(2*Delta);
  }
  boundary(all);
  gnuplotPipe = popen ("gnuplot", "w");
  fprintf(gnuplotPipe,
	  "set term pngcairo\n"
	  "set xr [-0.5 : 2.3]\n"
	  "set yr [0 : 15]\n"
	  "set key off\n"
	  "set grid\n"
	  "set title 'Buoyancy profile'\n"
	  "set xlabel 'buoyancy'\n"
	  "set ylabel 'height'\n"
	  "set size square\n");
}
/**
Our buoyancy field is diffusive with the same diffusivity als the momentum ($Pr = 1$). However, since we calculate the flux according to our own lumped-parameter model we need to ensure there is no flux over the `bottom` boundary. This is done by setting the buoyancy diffusiviy to zero at the surface.   
*/
face vector muz[];
event tracer_diffusion(i++){
  scalar r[];
  foreach_face(){
    if (y < Y0 + 1E-10)
      muz.x[] = 0.;
    else
      muz.x[] = muc.x[];
  }
  /**
  Furthermore, we set a tendency term in the lowest cells due to the adopted *soil-buoyancy balance*. We also diagnose the resulting surface-averaged flux and the mean surface buoyancy.   
  */
  foreach(){
    r[] = 0;
    if (y < Delta)
      r[] = (QFLX + GFLX)/Delta;
  }
  double flx = 0, bt = 0;
  foreach_boundary(bottom, reduction(+:flx) reduction(+:bt)){
    flx += (QFLX + GFLX) * Delta;
    bt += BSURF * Delta;
  }
  bt /= L0;
  flx /= L0;
  fprintf(stderr, "%g %g %g %d\n", t, flx, bt, i);  
  diffusion(b, dt, muz, r = r);
}

event adapt(i++)
  adapt_wavelet((scalar *){u, b}, (double []){0.05, 0.05, 0.02}, 9); 

event bviewer(t += 0.1; t <= temp){
  scalar omega[];
  view(fov = 25, tx = -0.5, ty = -0.4, width = 1200, height = 500);
  vorticity(u, omega);
  squares("omega", map = cool_warm);
  translate(x = L0){
    squares("b", min = -0.2, max = 1);
  }
  translate(x = -L0){
    cells();
  }
  draw_string(" Cells, Vorticity and Temperature", size = 35);
  save("mp2.mp4");
}
/**
We output many snapshots of the profiles in a `.png` format.
*/
int frame = 0;
event profs(t += 0.1){
  fprintf(gnuplotPipe, "set output 'plot%d.png'\n", frame);
  fprintf(gnuplotPipe, "plot '-' w l lw 5\n");
  double yp = Y0 + 1.E-8;
  boundary ({b});
  double by[1];
  while (yp < L0 + Y0){
    double dy = average_over_yp({b}, by, yp);
    fprintf(gnuplotPipe, "%g %g\n", by[0], yp);
    yp += dy;
  }
  yp = Y0 + L0 - 1.E-8;
  average_over_yp({b}, by, yp);
  fprintf(gnuplotPipe, "%g %g\n", by[0], yp);
  fprintf(gnuplotPipe, "e\n");
  frame++;
}
/**
These `.png`s are concatenated into a `.mp4` movie and then removed from the disk.
*/
event stop (t = end){
  system("rm mov.mp4");
  system("ffmpeg -r 25 -f image2 -i plot%d.png -c:v libx264 -vf format=yuv420p -y mov.mp4");
  system("rm plot*");
  return 1;
}