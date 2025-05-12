/**
## 2D planar vorticity near a bubble
We use the velocity/fraction field form [code](http://basilisk.fr/sandbox/yonghui/vortex/2dvortwrite.c), and run the simulation in a 2D case.
It's equivalent to a planar simulation with $k=0$
*/
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "tag.h"  
//#include "vtknew.h"
#include "view.h"
/**
Parameters
*/
#define LVM 10
#define R0 0.125
#define Re 5000. // Re_sim * R_0
#define We 100. // We_sim * R_0
#define mur 0.02 //ratio vis
#define tend 5.0

FILE * fp1, * fp2;
scalar omega[];
/**
## main 
We set periodic bc at both sides to avoid the border effect
*/
int main() {
  size(2.);
  init_grid (512);
  origin(-L0/2.,-L0/2.);
  rho1 = 1.e-3, mu1 = mur*R0/Re;
  rho2 = 1., mu2 = R0/Re, f.sigma = R0/We;
  TOLERANCE = 1e-5;
  periodic(left);
  periodic(top);
  //DT = 0.001;
  fp1 = fopen("bub1","w");
  fp2 = fopen("bub2","w");
  run();
}


event init (t = 0) {
  int numr = 140000;
  double dons[numr][6];
  /**
We read the file contains the output fields 
*/
  FILE *file2 ;
  file2 = fopen("data_t220", "r");
  for(int i=0; i < numr; i++){
    for(int j=0; j < 6; j++){
      fscanf(file2, "%lf",&dons[i][j]);
    }
  }
  fclose(file2);

  double rlim=0.2;
  scalar omega0[], psi[];
  psi[left]   = dirichlet(0);
  psi[right]  = dirichlet(0);
  psi[top]    = dirichlet(0);
  psi[bottom] = dirichlet(0);

  /**
We read the fraction field (bubble shape) and the vorticity field from the file.
Copied from [vortex](http://basilisk.fr/src/test/vortex.c)
{
We need to convert the initial vorticity field into the velocity field. 
To do so we first declare the streamfunction $\psi$ and vorticity $\omega$ fields. }
*/

  foreach(){
    double rcell = sqrt(sq(y)+sq(x));
    if (rcell < rlim ){
      for (int ii=0; ii < numr ; ii++){
        if ( sq(x - dons[ii][0]) + sq(-y - dons[ii][1]) < sq(0.5*Delta) ){
          omega0[] = rcell < rlim*0.9 ? dons[ii][2]:0.;									
          f[] = dons[ii][3];
          break;			
        }
      }//end find ii
    } //end if r < 0.1
  } //end foreach
  boundary ({psi,omega0});

  // ============The readed vorticity field ============
  view (fov = 24.,bg = {1,1,1}, width = 1024, height = 1024, samples = 1);
  clear();
  box();
  squares("omega0", max=10., min=-10.);
  save ("omega_read.png");

  /**
We then solve the Poisson equation
  $$
  \nabla^2 \psi = \omega
  $$
  and compute the centered velocity components by differentation of the
  streamfunction i.e.
  $$
  u_x = - \partial_y \psi ,
	\quad 
  u_y = \partial_x \psi
  $$
*/

  poisson (psi, omega0);
  coord ff = {-1.,1.};
  foreach()
    foreach_dimension()
      u.x[] = ff.x*(psi[0,1] - psi[0,-1])/(2.*Delta);
  boundary ((scalar *){u});
  vorticity (u, omega);

  // ============The recoverd vorticity field============
  view (fov = 3.,bg = {1,1,1}, width = 1024, height = 1024, samples = 1);
  box();
  draw_vof("f",lw=5);
  squares ("omega", linear = true,max=30.,min=-30.);
  save ("omega_recoverd.png");
}

/***/
event mydata (i++){
  fprintf (stderr, "%d %d %d %g %g %g\n", \
           i, mgp.i, mgu.i, t, dt, statsf(f).sum);
}

event myvid (i += 20; t <= tend){

  char legensd[2000]; //time
  sprintf(legensd, "t=%0.2g", t);
  vorticity (u, omega);
  // ============video 1============
  view (fov = 15.,bg = {1,1,1}, width = 1024, height = 1024, samples = 1);
  box();
  draw_vof("f",lw=2);
  draw_string(legensd, 0, size = 60.,lw = 5.);
  squares ("omega", linear = true,max=30.,min=-30.);
  save ("omega.mp4");
}

#if 0
event mypic (t += tend/10.; t <= tend){
  char picname[800];
  sprintf (picname, "2dbubt%03g.png",t*100);
  vorticity (u, omega);
  // ============pic 1============
  view (fov = 3.,bg = {1,1,1}, width = 1024, height = 1024, samples = 1);
  box();
  draw_vof("f",lw=5);
  squares ("omega", linear = true,max=30.,min=-30.);
  save (picname);
}
#endif



/**
We use the adapted mesh.
*/

event adapt (i++) {
  double uemax = 0.01;
  adapt_wavelet ({f,u},(double[]){0.001,uemax,uemax}, LVM, 5);
}

/**
# Results

![Read](planarvort/omega_read.png)
![recovered](planarvort/omega_recovered.png)

![Bubbble's evolution](planarvort/omega.mp4)(loop)
*/
