/**
#Regular bubble entrapment

Bubbles are entrapped in the crater produced by the impacting drops only in a sharply delimited region of parameter space, which we shall refer to as the bubble region. If, at the moment of impact, the drop is assumed to be spherical, then two parameters are sufficient to completely characterize an impact event. They can be conveniently chosen to be the impact velocity $U$ and the drop diameter $D$ or, alternatively, the Froude number $Fr=\frac{U^2}{gD}$ and the Weber number $We=\frac{\rho D U^2}{\sigma}$

The inertia, viscous force and capillary force together decide the regular bubble entrapment. which has been extensively studied by experiment. To validate previous experimental result and study some extreme case such as immiscible flow, this numerical method is needed.
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "view.h"
#include "reduced.h"

#define diameter 1. // Numerical diameter equals to 50% of the box size
#define length 10. // Numerical size of the domain
// Define all physical parameters

#define D 3.0e-3 //       m    diameter of the droplet
#define hdrop 7.7*D //  m    position of the droplet's center
#define wheight 7.*D //  m    water layer's thickness
#define V 2.04    //    m/s  Initial imapct velocity

#define Mu_w 0.001 //     Pa.s
#define Mu_g 1.6e-5 //    Pa.s
#define rho_w 1000. //    kg/m**3
#define rho_g 1.6 //       kg/m**3
#define sigma12 72.e-3 //  Surface tension water/air
#define G 9.8//  	   m**2/s

//Dimensionless numbers

#define Re (rho_w*V*D/Mu_w)
#define We (rho_w*sq(V)*D/sigma12)

#define vrho_w 1.
#define vrho_g (rho_g/rho_w)
#define vmu_w (1./Re)
#define vmu_g ((1./Re)*(Mu_g/Mu_w))
#define SIGMA 1./We
#define FR V/sqrt(G*D)
//convert ratio between physical and numerical values

#define ratio diameter/D // convert all geometrical lengths in numerical value
// the physical domain of calcul is define by length/ratio

// Convert all real values in numerical values

#define wheight_num wheight*ratio
#define hdrop_num hdrop*ratio

// Geometry set up

#define droplet(x,y) (-sq(x-hdrop_num)-sq(y)+sq(diameter/2))
#define pool(x,y) wheight_num-x
//#define epsilon 0.95e-4*ratio
#define epsilon 0.1*D
#define initwater(x,y) (-sq(x-hdrop_num)-sq(y)+sq(diameter/2+epsilon))

// Boundaries conditions

//u.t[bottom] = neumann(0.);
//u.n[bottom] = neumann(0);

// Numerical simulation

int maxlevel = 13;
int minlevel = 6;

double dtuser = 0.1;
double tmax = 20.;
double uem = 0.01;
double fem = 0.01;
double tc = 0.;// restart point


int frame = 0; // Frame counter

double tframe = 0.;// time counter 

//define a global varaible for vorticity
scalar omega[];
//define a global varaible for drop  //it was in the event init before

int main() {
  size (length);
  init_grid (1 << minlevel);
  origin(0.,0.);
  rho2 = vrho_g;
  rho1 = vrho_w;
  mu2 = vmu_g;
  mu1 = vmu_w;
  f.sigma = SIGMA;
  //size(length);
  run();     
}


event init (t = tc) {
  
  if (!restore (file = "restart")) {
  printf("we are going through the init");
  scalar f5[];
  refine(sq(x-hdrop_num)+sq(y)<sq(diameter/2+0.2) && sq(x-hdrop_num)+sq(y) > sq(diameter/2-0.2) && level < maxlevel);
  refine(x<wheight_num+0.1 && x>wheight_num+0.1 && level < maxlevel);
  fraction(f, max(pool(x, y), droplet(x, y)));
  fraction(f5, droplet(x,y));
  foreach() {
    u.x[] = -1. * f5[];
  }
}
else{
 printf("we restore the dumpfile");
 }
}


#if !REDUCED
event acceleration (i++) {
  face vector av = a;
  foreach_face(x)
    av.x[] -= 1./sq(FR);
}
#endif

event adapt (i++) {
  adapt_wavelet({u, f}, (double[]){uem, uem, fem}, maxlevel, minlevel);
}





  //log the voulume velocity and position of bubble
  event logfile (i += 10) {
  double xb = 0., yb = 0., zb = 0., sb = 0.;
  double vbx = 0., vby = 0., vbz = 0.;
 
  foreach(reduction(+:xb) reduction(+:yb) reduction(+:zb)
	  reduction(+:vbx) reduction(+:vby) reduction(+:vbz)
	  reduction(+:sb)) {
    double dv = (1. - f[])*dv();
   if(x<5.5){
   xb += x*dv;
   yb += y*dv;
   zb += z*dv;
   vbx += u.x[]*dv;
   vby += u.y[]*dv;
   vbz += u.z[]*dv;
   sb += dv;
   }
    
          }
 fprintf (stderr,
	  "%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n", 
	   t, sb,
	   xb/sb, yb/sb, zb/sb,
	   vbx/sb, vby/sb, vbz/sb);
    fflush (stderr);
  }



event viewing (t = tc; t += dtuser; t <= tmax) {
  view (fov = 10, tx = 0, ty = -0.7, quat = {0, 0, -0.707, 0.707}, width=2000,height=2000);

  clear();
  
    draw_vof (c = "f", edges = true, larger = 1, filled = 1, max = 1, fc = {0.20784313725490197,0.5176470588235295,0.8941176470588236}, lw = 2);
    draw_vof ("f", lw = 5, lc = {0, 0, 0});

  box(notics = true);
  mirror({0,1}) {
    draw_vof (c = "f", edges = true, larger = 1, filled = 1, max = 1, fc = {0.20784313725490197,0.5176470588235295,0.8941176470588236}, lw = 2);
    draw_vof ("f", lw = 5, lc = {0, 0, 0});
    //squares("u.y", linear=true);
    box(notics=true);
  }
  //char ppm_filename[100];
  char png_filename[100];
  // tframe = 0.1 * frame;
  //sprintf(ppm_filename, "time-%04.1f.ppm", tframe);
  sprintf(png_filename, "time-%04.1f.png", t);

  save(png_filename);

  //output_ppm (f, file = "tframe.png", n = 1024, min = -1, max = 1);
 
  frame++; 
}

event movie (t =tc; t +=dtuser; t<= tmax) {
  
  vorticity (u, omega);
  view(fov = 10, tx = 0, ty = -0.7,
       quat = {0, 0, -0.707, 0.707},width=2000,height=2000);
  clear();
  //draw_vof ("f");
  squares ("omega", linear = true, spread = 10);

  //draw_vof (c = "f", edges = true, larger = 1, filled = 1, max = 1, fc = {0.20784313725490197,0.5176470588235295,0.8941176470588236}, lw = 3);
  //squares ("u.x", linear = true);

  box(notics = true);
  mirror({0,1}) {
  //draw_vof (c = "f", edges = true, larger = 1, filled = 1, max = 1, fc = {0.20784313725490197,0.5176470588235295,0.8941176470588236}, lw = 3);
  //squares ("u.x", linear = true);
 
  //draw_vof ("f");
  squares("omega", linear=true,spread=10);

  box(notics=true);
  }
  save("movie.mp4");
}


event snapshots (t =tc; t +=dtuser; t<= tmax) {
  char name[100];
  sprintf (name, "dump-%0.1f", t);
  dump (name);
}
