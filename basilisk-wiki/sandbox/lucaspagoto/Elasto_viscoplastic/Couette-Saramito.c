/**
# Saramito's Elasto-viscoplastic Model in a Couette Flow
This code is based on the Elasto-viscoplastic model of Saramito (2009) and it modifies the log-conformation model to include a variable polymeric viscosity, $\eta_p$, in order to include plasticity.

## Non-dimensionalization

$Re = \frac{8\rho U^2}{\tau_c} = \frac{8\rho U^2}{\eta \dot{\gamma_c}}$ 

$\eta_c = \eta_s + \eta_{pc}$ 

$Pl = \frac{\tau_y}{\eta_{pc} \dot{\gamma_c}} = \frac{\tau_y}{\left(\frac{\tau_y}{\dot{\gamma_c}} + K \dot{\gamma_c}^{n-1} \right)\dot{\gamma_c}} = \frac{\tau_y}{\tau_y + K\dot{\gamma_c}^n}$

$\lambda = \dot{\gamma_c} Wi$

## Code
*/

// Maximum time step
#define DT_MAX 1e-3
// Reynolds number
#define Re 1.
// Plastic number
#define Pl 0.2
// Characteristic velocity
#define U 1.
// The domain size
#define HEIGHT 1.
// Characteristic strain rate
#define gamma_c U/HEIGHT
// Density
#define density 1.
// Power index
#define nnn 0.7
// Ratio of the solvent viscosity to the total viscosity
#define BETA 0.5
// Total viscoelastic viscosity 
#define MU0 8*density*U*U/(Re*gamma_c)
// Polymer viscosity
#define MUP ((1. - BETA)*MU0)
// Solvent viscosity 
#define MUS (BETA*MU0)
// Weissberger number
#define WI 1.0
// The wall velocity
#define uwall(t) (1.)

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "log-conform.h"
#include "view.h"

double tEnd = 10.;

scalar yielded[];
scalar mupc[];  
scalar fnorm[];
scalar srp[];
scalar lamb[];
scalar deviatoric[];
//const scalar lamb[] = WI;

double tau_y, K, epsilon;
double Nreg = 1e2;
int DEVIATORIC = 1;

int main()
{
  periodic (right);
  DT = DT_MAX;
  init_grid (1 << 8);
  size(HEIGHT);

  
  lambda = lamb;
//  const scalar mupc[] = MUP;
  mup = mupc;

  const face vector mus[] = {MUS,MUS};
  mu = mus;
  const scalar rhoc[] = density;
  rho = rhoc;

  tau_y = Pl*8*density*U*U*(1-BETA)/Re;
  K = (MUP*gamma_c - tau_y)/pow(gamma_c,nnn);
  epsilon = tau_y/(Nreg * MUP * gamma_c);

  stokes = true ;
//  TOLERANCE = 1E-4;
  run();
  
//  run();
}

u.t[top]    = dirichlet(uwall(t));
u.t[bottom] = dirichlet(0.);
//u.n[left] = neumann(0.);
//u.n[right] = neumann(0.);

/* Initialization */
event init (i = 0) 
{
  foreach()
    u.x[] = y;
  boundary ((scalar *){u});
}


event logfile(i+=1)
{
  int j = 0;
  scalar le[];
  foreach()
  {
    le[] = level;
    j++;
  }
  stats s = statsf(le);
  printf ("i = %d t = %g, #cells = %d, le.min = %g, le.max = %g DEV = %d WI = %.3f tau_y = %.3f K = %.3f n = %.3f BETA = %.3f\n", i,t,j, s.min, s.max, DEVIATORIC, WI, tau_y, K, nnn, BETA);
  fflush(stdout);
}

/* Calculating $\eta_p$ in the event properties */
event properties (i=0; i++)
{
  foreach()
  {
   //mupc[] = MUP;   // this line reduces the code to the Oldroyd-B model
    if (DEVIATORIC == 1)   
      fnorm[] = sqrt(0.5)*sqrt(sq((tau_p.x.x[]-tau_p.y.y[])/2) + sq(tau_p.x.y[]) + sq(tau_p.y.x[]) + sq((tau_p.y.y[]-tau_p.x.x[])/2));   //Frobernoius norm of the deviatoric polymeric stress tensor
    if (DEVIATORIC == 0)
      fnorm[] = tau_p.x.y[]; //Frobernoius norm of the polymeric stress tensor shear components

    // strain rate on the plastic component of the mechanical analog
    if (fnorm[] - tau_y <= 0)
      srp[] = 0.;
    else
      srp[] = pow((fnorm[] - tau_y)/K, 1/nnn);

    if (fnorm[] - tau_y <= 0)
     deviatoric[] = tau_y + epsilon;
    else
      deviatoric[] = fnorm[];
  }
  boundary ({fnorm});
  boundary ({srp});
  boundary ({deviatoric});
  
  foreach()
  {
    mupc[] = pow( (K * pow(deviatoric[], nnn)/(deviatoric[] - tau_y)), 1/nnn);
  }
  boundary ({mupc});

  foreach()
  {
     lamb[] = 1. * mupc[] / MUP;
  }
  boundary ({lamb});

/*
  foreach()
  {
      if (srp[]*mupc[]/(gamma_c*MUP) > Pl)
      {
        yielded[] = 1.;
      }
      else
      {
        yielded[] = 0.; 
      }
  }
  boundary ({yielded});*/
  
  // Variables in the middle of the domain (0.5, 0.5)
  char nameu2[50];
  if (DEVIATORIC == 1) 
    sprintf (nameu2, "pos_dev.txt");
  if (DEVIATORIC == 0)
    sprintf (nameu2, "pos_she.txt");
  static FILE * fname2 = fopen (nameu2, "w");
  fprintf (fname2, "%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n", t, interpolate(u.y, 0.5, 0.5), interpolate(u.x, 0.5, 0.5), interpolate(tau_p.x.y, 0.5, 0.5), interpolate(tau_p.y.x, 0.5, 0.5), interpolate(tau_p.x.x, 0.5, 0.5), interpolate(tau_p.y.y, 0.5, 0.5), interpolate(fnorm, 0.5, 0.5), interpolate(srp, 0.5, 0.5), interpolate(mupc, 0.5, 0.5), interpolate(lamb, 0.5, 0.5));
  fflush (fname2);
}

event end (t = tEnd) {DEVIATORIC = 0;}

/** 
## Outputs

~~~gnuplot
set yrange [0:10]
set xlabel 'Time'
set ylabel 'Stress Components, Plastic Viscosity and Plastic Strain Rate'
set ytics 1,1,10
set xtics 1,1,10
set size square
set grid
set key outside right center
plot [0:10][0:8] 'pos_dev.txt' u 1:4 w l lw 3 lc 1 dt 1 t "tau-p.x.y", '' u 1:6 w l lw 3 lc 2 dt 1 t "tau-p.x.x", '' u 1:9 w l lw 3 lc 3 dt 1 t "Plastic Strain Rate", '' u 1:10 w l lw 3 lc 4 dt 1 t "Plastic Viscosity eta_p", '' u 1:11 w l lw 3 lc 5 t "lambda", '' u 1:(($10)/($11)) w l lw 3 lc 6t "G"
~~~
![Deviatoric Stress Tensor]

~~~gnuplot
set yrange [0:10]
set xlabel 'Time'
set ylabel 'Stress Components, Plastic Viscosity and Plastic Strain Rate'
set ytics 1,1,10
set xtics 1,1,10
set size square
set grid
plot [0:10][0:8] 'pos_she.txt' u 1:4 w l lw 3 lc 1 dt 1 t "tau-p.x.y", '' u 1:6 w l lw 3 lc 2 dt 1 t "tau-p.x.x", '' u 1:9 w l lw 3 lc 3 dt 1 t "Plastic Strain Rate", '' u 1:10 w l lw 3 lc 4 dt 1 t "Plastic Viscosity eta_p", '' u 1:11 w l lw 3 lc 5 t "lambda", '' u 1:(($10)/($11)) w l lw 3 lc 6t "G"
~~~
![Only the Shear Components]

## References
Saramito, P. (2009). A new elastoviscoplastic model based on the Herschel–Bulkley viscoplastic model. Journal of Non-Newtonian Fluid Mechanics, 158(1-3), 154-161.
*/