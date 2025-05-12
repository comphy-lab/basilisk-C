#include "grid/octree.h"
#include "navier-stokes/centered.h"
#define TWO_WAY 1
#include "stokes-particles.h"
#include "view.h"
#include "scatter3.h"  

//air flow
u.n[front] = dirichlet (0.7);
u.r[front]  = neumann(0);
u.r[front]  = neumann(0);
p[front]   = neumann(0);
pf[front]  = neumann(0);


u.n[back] = neumann(0);
u.t[back] = neumann(0);
u.r[back] = neumann(0);
p[back]   = dirichlet(0);

Particles drop;

#define MUZ 1.8e-5      // air vicosity
#define U0 27.6 //m/s   // particle velocity
#define rout 1e-4 //m.  //nozzle outer diameter

double rp[1000000], p_new, mu_new, sigma_new;
int  pdn = 30, n[1000000], n_particles, j, k, N0 = 1e3;   


int main() 
{
  const face vector muc[] = {MUZ, MUZ, MUZ}; //mu Water
  mu = muc;
  G.y = -9.81;
  L0 = 1.1;
  X0 = Y0 = Z0 = -L0/2;
  mu = muc;
  N = 64;
  DT = 0.00006;
  run();
}


event init (t = 0) 
{
  const scalar rhof[] = 1.2;
  rho = rhof;
  drop = new_inertial_particles (0);
}


//adding particles after every 20 steps
event add_drops( i += 20)
{
  double rin = rout/3;
  double delt = rout -rin, rav = 0.5*(rin + rout);
  n_particles = 0;
  for ( int j = 0; j < pdn; j++ )
  {
    rp[j] = 0.5*(55e-6 + 55e-6*noise());

    //intermediate
    p_new = 0.3208137453;
    mu_new = 3.0615978666;
    sigma_new = 0.0704145049;          
    n[j] = p_new * N0 * exp( -pow((log(2*rp[j]*1e6) - mu_new),2) /(2*sigma_new))/(sqrt(2*pi) * 2*rp[j]*1e6 * sqrt(sigma_new));
    n_particles = n_particles + n[j];
  }
  j = 0;
  k = 0;
  //defining particle properties
  for ( int l = 0; l < n_particles; l++ )
  {
    double r = rav + 0.5*delt*noise();
    double th = pi*noise();
    if ( j < n[k])
      {
      particle p = { .x = r*cos(th), .y = r*sin(th), .z = 0.15,  
        .u.x = U0*(p.x/sqrt(sq(p.x)+sq(p.y)))*sin(pi/6 + pi*12.*noise()/180), 
        .u.y = U0*(p.y/sqrt(sq(p.x)+sq(p.y)))*sin(pi/6 + pi*12.*noise()/180), 
        .u.z = -U0*cos(pi/6), .u2.x = 1000, .u2.y = rp[k], .tag = n_particles*t};
      add_particle (p, drop);
      }  
    else
      {
      k++;
      j = 0;
      particle p = { .x = r*cos(th), .y = r*sin(th), .z = 0.15, 
        .u.x = U0*(p.x/sqrt(sq(p.x)+sq(p.y)))*sin(pi/6 + pi*12.*noise()/180), 
        .u.y = U0*(p.y/sqrt(sq(p.x)+sq(p.y)))*sin(pi/6 + pi*12.*noise()/180), 
        .u.z = -U0*cos(pi/6), .u2.x = 1000, .u2.y = rp[k], .tag = n_particles*t};
      add_particle (p, drop);
      }
    j++;
    
  }
  
}

//removing particle if it is going out of boundary
event remover (i++) 
{
  remove_particles (drop, z > L0/2 || y > 0.15 || x > 0.175 || z < -0.40 || y < -0.2 || x < -0.175);
}
//z<-0.292


event adapt (i++) 
  adapt_wavelet ((scalar*){u}, (double[]){0.01, 0.01, 0.01}, 9, 5);

//making movie
event mov (t += .00006) 
{
  view (theta = 0*10. - 1.4, bg = {1,1,1});
  foreach_particle_in(drop)
  {
    scatter (drop, s = p().u2.y*2e6, pc = {0, 0, 1}, x = p().x, y = p().y, z = p().z);
  }
  box();
  translate (z = -L0/2)
  //cells();
  save ("locs.mp4");
}

event end (t = 10);

//storing particle data
event particle_data (i=500; i++)
{
  char name[80];
  double vol=0;
  sprintf (name, "data.dat");
  FILE * fp = fopen (name, "a");
  foreach_particle_in(drop)
    {
    fprintf (fp, "%g %g %g %g %g %g %g %g\n", t, p().x, p().y, p().z, p().u.x, p().u.y, p().u.z, p().u2.y);
    }
  fclose (fp); 
}