/**
# Forward step

Here we compare the standard adherence condition u.t = 0. on the embedded boundary.

![Animation of the vorticity field.](adherentWall/vort_0.mp4)(loop)

![Animation of the vorticity field, adherent wall.](adherentWall/vort_1.mp4)(loop)

There is a huge difference in the first instants.

*/

#include "embed.h"
#include "navier-stokes/centered.h"

double Reynolds = 160.;
int maxlevel = 9;
face vector muv[];

double eps = 1.e-6;
double mygeom(double x, double y){
  return max(-x+eps,y+L0/2-L0/20.-eps); // marche montante
}


int boundary_type;

int main() {
  L0 = 8.;
  origin (-L0/2., -L0/2.);
  N = 512;
  mu = muv;
  for ( boundary_type = 0; boundary_type < 2; boundary_type++)
  {
    fprintf(stderr,"OK %d boundary type", boundary_type);
    run(); 
  }

}


event properties (i++)
{
  if(boundary_type == 1){
    foreach(){
      if(cs[]*(1-cs[]) != 0.)
        u.x[] = 0.;
    }
  }

  foreach_face()
    muv.x[] = fm.x[]*0.125/Reynolds;
}

u.n[left]  = dirichlet(1.);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);

u.n[right] = neumann(0.);

u.t[bottom] = dirichlet(0.);

p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);

event init (t = 0)
{

  vertex scalar distn[];
 
  foreach_vertex() {
    distn[] = mygeom(x,y);
    }
  boundary({distn});
  fractions(distn,cs,fs);
  fractions_cleanup(cs,fs);

  
  foreach()
    u.x[] = cs[] ? 1. : 0.;
}

event logfile (i++)
  fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);


event movies (i += 4; t <= 4.)
{
  scalar omega[], m[];
  vorticity (u, omega);
  foreach()
    m[] = cs[] - 0.5;
  char name[80];
  sprintf (name, "vort_%d.mp4", boundary_type);
  output_ppm (omega, file = name,
	      min = -10, max = 10, linear = true, mask = m);
}

event adapt (i++) {
  adapt_wavelet ({cs,u}, (double[]){1e-2,3e-2,3e-2}, maxlevel, 4);
}

