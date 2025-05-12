/**
   We set up a 2D plane Poiseuille flow. This is in part an investigation of the mask function. This problem file is courtesy of J. Eggers. */

#include "navier-stokes/centered.h"

/**
   Define Reynolds number, channel length, half-width, total time, and resolution. */
#define Re (1.) // Reynolds number 
#define L (2.) // length of channel
#define hw (L / 11.0) // half-width of channel
#define tfinal (5.) // length of run
#define LEVEL (9)

/**
The domain is eight units long, centered vertically. */

int main() {
  L0 = L;
  origin (-hw, -L0/2.);
  init_grid(1 << LEVEL);
  /**
     Set viscosity to unity. */
  const face vector muc[] = {1.,1.};
  mu = muc;
  run(); 
}

/**
   Set the pressure drop over channel length, and the velocity amplitude. */

double dp = 2.*L*Re/hw; 
double U = Re / hw;     
/**
   Set the inlet velocity at the left according to known solution. An
outflow condition is used on the right boundary, with the exit
pressure specified. No-slip is enforced. */
u.n[left]  = dirichlet(U*(sq(hw)-sq(y)));
p[left]    = neumann(0.);   
pf[left]   = neumann(0.);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

u.t[top]   = dirichlet(0.);
u.t[bottom]  = dirichlet(0.);


event init (t = 0) {
  /**
     Define the masks at the channel halfwidths. */
  mask (y >  hw ? top :
	y < -hw ? bottom :
  	none);

  /**
     Initialize the interior flow with a Poiseuille profile, although
     a quiescent initial flow should also work. */
  foreach()
    {
      u.x[] = U*(sq(hw) - sq(y));
    }
  boundary((scalar *){u});
}

/**
   Output statistics. */

event logfile (i++)
  fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);

/**
   Output movies for pressure, horizontal and vertical velocity, and
   local refinement level. */
event movies (i ++; t <= tfinal) {
  scalar omega[];
  vorticity (u, omega);
  static FILE * fp = fopen ("p.ppm", "w");
  static FILE * fp1 = fopen ("ux.ppm", "w");
  static FILE * fp2 = fopen ("uy.ppm", "w");
  static FILE * fp3 = fopen ("l.ppm", "w");
  scalar l[];
  foreach(){
    l[] = level;
  }
  output_ppm (p, fp, box= {{-hw,-hw},{L-hw,hw}},
	      spread = 2, linear = true);
  output_ppm (u.x, fp1, box = {{-hw,-hw},{L-hw,hw}},
  	      linear = true, spread = 2);  
  output_ppm (u.y, fp2, n=1024, box= {{-hw,-hw},{L-hw,hw}},
	      linear = true, spread = 2);
  output_ppm (l, fp3, n=1024, box= {{-hw,-hw},{L-hw,hw}} );
	      
}
/**
   Output slices through the velocity field; time steps of 1 until 
   final time tfinal. */


event output_velocity_profiles (t += 1.; t<=tfinal) {
   static FILE * fux = fopen ("poiseuille_profiles.dat", "w");
   static FILE * fp = fopen ("poiseuille_pressure.dat", "w");
   fprintf (fux, "\n");
   fprintf (fux, "# time = %.4f\n",t);
   fprintf (fp, "\n");
   fprintf (fp, "# time = %.4f\n",t);
   for (double yy = -hw; yy <= hw; yy += hw/100.0) {
     double xslice = 1.;
     fprintf (fux, "%.4f %.4f\n", yy, interpolate (u.x, xslice, yy));
   }
   for (double xx = -hw; xx <= L-hw; xx += 0.1) {
     fprintf (fp, "%.4f %.4f\n", xx, interpolate (p, xx, 0.));
   }
}  

/**
   Adapt on the velocity field. */
event adapt (i++) {
  adapt_wavelet ({u}, (double[]){1e-2,1e-2,1e-2}, LEVEL, LEVEL-4);
}
