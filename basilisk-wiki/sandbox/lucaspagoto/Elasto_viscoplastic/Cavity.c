/**
# Flow a elasto-viscoplastic material over a cavity

We perform the simulation of the an elasto-viscoplastic fluid (Saramito's model) over a cavity. The yielded region assimetry shear region asymmetry is observed.
*/

#include "embed.h"
#include "navier-stokes/centered.h"
#include "log-conform.h"
#include "view.h"

#define Re 1.0  // Reynolds number  
#define Pl 0.5  // Plastic number
#define WI 1.0  // Weissenberg number
#define U 0.1   // Inlet velocity
#define D 0.1001  // Channel height 
#define H 0.2001  // Cavity depth
#define h 0.201   // Cavity length
#define HEIGHT 1. // Domain size
#define DT_MAX 1e-4 // time step size
#define nnn 0.7   // Flow index
#define BETA 0.5  // Solvent to total viscosity ratio
#define MU0 D*U/Re  // Total viscosity
#define MUP ((1. - BETA)*MU0) // Polimeric characteristi viscosity
#define MUS (BETA*MU0)  // Solvent viscosity

double tEnd = 5.;

face vector muv[], muv2[];
scalar rhoc[], shear[], eta[], yielded[], un[], le[];
scalar mupc[], fnorm[], srp[], lamb[], deviatoric[];

int LEVEL = 6;

double tau_y, K, etac, etamx, epsilon, gamma_c;
double Nreg = 1e2; // regularization parameter

u.n[bottom] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);

u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);

u.n[left] = dirichlet(y < H && y > H - D ? U : 0);
p[left] = neumann(0);
pf[left] = neumann(0);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);



int main()
{

  DT = DT_MAX;
  init_grid (1 << LEVEL);
  size(HEIGHT);

  gamma_c = U/(D);
  etac = D*U*(1-BETA)/Re;
  tau_y = Pl*U*U*(1-BETA)/Re;
  K = (MUP*gamma_c - tau_y)/pow(gamma_c,nnn);
  etamx = Nreg * etac;
  epsilon = tau_y/etamx;

  mu = muv;
  rho = rhoc;
  lambda = lamb;
  mup = mupc;

  TOLERANCE = 1E-4;

 run();
}


//int DEVIATORIC = 1;
event properties (i++)
{
  foreach()
  {
    //Frobernoius norm of the polymeric stress tensor
  //    fnorm[] = sqrt(0.5)*sqrt(sq(tau_p.x.x[]/2) + sq(tau_p.x.y[]) + sq(tau_p.y.x[]) + sq(tau_p.x.x[]/2)); 
 fnorm[] = sqrt(sq((tau_p.x.x[]-tau_p.y.y[])/2) + sq(tau_p.x.y[]));

    // strain rate on the plastic component of the mechanical analog
    if (fnorm[] <= tau_y)
      srp[] = 0.;
    else
      srp[] = pow((fnorm[] - tau_y)/K, 1/nnn);

    if (fnorm[] <= tau_y + epsilon)
     deviatoric[] = tau_y + epsilon;
    else
      deviatoric[] = fnorm[];
  }
  boundary ({fnorm, srp, deviatoric});

  foreach()
  {
    mupc[] = cm[]*pow( (K * pow(deviatoric[], nnn)/(deviatoric[] - tau_y)), 1/nnn);  //Elasto-viscoplastic
  //  mupc[] =cm[]*MUP;   // Viscoelastic
    rhoc[]= cm[];
  }
  boundary ({mupc, rhoc});

  foreach()
  {
     lamb[] = cm[] * WI * mupc[] / MUP;  //Elasto-viscoplastic
   //  lamb[] = cm[];   // Viscoelastic
  }
  boundary ({lamb});

  foreach_face()
  {
    muv.x[] = fm.x[]*MUS;
  }
  boundary ((scalar *){muv});

 foreach()
  {
    if(fnorm[] > tau_y)
    {
      yielded[] = 1.;
    }
    else
    {
      yielded[] = 0.; 
    }
  }
  boundary ({yielded});
}




event init (i = 0) 
{
//  refine (y < H + 0.01 && level < MAXLEVEL);
  
  vertex scalar phi[];
  foreach_vertex()
  {
  	phi[] = intersection (x - 0.401, 0.601 - x);
    phi[] = union (y - H + D, phi[]);
    phi[] = intersection (H - y, phi[]);
  }
  boundary ({phi});
  fractions (phi, cs, fs);

  foreach()
  {
    u.x[] = 0.1;
    u.y[] = 0.;
  }
}



event adapt (i++)
{
   adapt_wavelet ({yielded, u}, (double[]){0.01, 0.01, 0.01}, 7, 6);
}


/*
event uprofile2 (i++) {
  char cen[50];
  sprintf (cen, "center.txt");
  static  FILE * cenn = fopen (cen, "w");
  fprintf (cenn, "%.5f %.5f %.5f %.5f %.5f\n", t, interpolate(u.x, 0.9, 0.1505), interpolate(tau_p.x.y, 0.9, 0.2), interpolate(p, 0.1, 0.1505), interpolate(p, 0.9, 0.01505));
  fflush (cenn);
}*/


event movie_vel (t += 0.05)
{
  scalar vel[];
  foreach()
  {
      vel[] = pow(sq(u.x[])+sq(u.y[]),0.5);
  }
  boundary({vel});
  char movie_vel[80];
  sprintf (movie_vel, "movie_vel.mp4");
  clear(); 
  view(tx = -0.50, ty = 0.0);
  draw_vof ("cs");
  squares("vel", linear = true);
  save(movie_vel);
}

event movie_tau_p (t += 0.05)
{
  scalar tau_pp[];
  foreach()
  {
      tau_pp[] = pow(0.5*(sq(tau_p.x.x[])+sq(tau_p.y.y[])+sq(tau_p.x.y[])+sq(tau_p.y.x[])),0.5);
  }
  boundary({tau_pp});
  char movie_vel[80];
  sprintf (movie_vel, "movie_tau_pp.mp4");
  clear(); 
  view(tx = -0.50, ty = 0.0);
  draw_vof ("cs");
  squares("tau_pp", linear = true);
  save(movie_vel);
}

event movie_p (t += 0.05)
{
  char movie_vel[80];
  sprintf (movie_vel, "movie_p.mp4");
  clear(); 
  view(tx = -0.50, ty = 0.0);
  draw_vof ("cs");
  squares("p", linear = true);
  save(movie_vel);
}

event movie_yield (t += 0.05)
{
  char movie_vel[80];
  sprintf (movie_vel, "movie_yield.mp4");
  clear(); 
  view(tx = -0.50, ty = 0.00);
  draw_vof ("cs");
  squares("yielded");
  save(movie_vel);
}


event end (t = tEnd) {}


/**
## Results
### Velocity field
![Velocity field](Cavity/movie_vel.mp4)

### Polymeric stress
![Polymeric stress](Cavity/movie_tau_pp.mp4)

### Pressure field
![Pressure field](Cavity/movie_p.mp4)

### Yielded region
![Yielded region](Cavity/movie_yield.mp4)
*/