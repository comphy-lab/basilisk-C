/**
#  Retraction of a viscoplastic liquid sheet
The base code for this problem is the same as in [column_SCC.c](/sandbox/M1EMN/Exemples/column_SCC.c). We have included the [tension.h](/src/tension.h) to incorporate the effect of surface tension.

We consider an elongated liquid sheet of initial width $h_0$ and length $L_0$. We have used the Bingham model to mimic the viscoplastic material. We consider that the liquid sheet has yield stress $\tau_y$ and plastic viscosity $\mu_p$. The liquid sheet has surface tension $\sigma$ and density $\rho_1$ and the surrounding Newtonian fluid has density $\rho_2$ and viscosity $\mu_2$. The constitutive equations are given in 
[column_SCC.c](/sandbox/M1EMN/Exemples/column_SCC.c).

##  Non-dimensionalization 
Using the scaling
$\bar{x}_i = {x_i}/{h_0}, \hspace{1mm} 
\bar{u}_i = {u_i}/{U_{TC}},
\hspace{1mm} (\bar{\mu}_1, \bar{\mu}_2) = ( \mu_1, \mu_2)/{\mu_c} ,
\hspace{1mm} (\bar{\rho}_1, \bar{\rho}_2) = ( \rho_1, \rho_2)/{\rho_1},\hspace{1mm} \bar{P} = {h_0 P}/{\sigma},  \hspace{1mm} \bar{t} = {t}/\sqrt{\sigma/\rho_1 h_0^3},$

we can write the Navier-Stokes equations in the form 


$\bar{\rho} \left[ \frac{\partial \bar{u_i}}{\partial \bar{t}} + \bar{u_i} \frac{\partial \bar{u_j}}{\partial \bar{x_j}} \right] = \frac{\partial \bar{P}}{\partial \bar{x_i}} +  Oh \frac{\partial}{\partial \bar{x_j}} \left[ \bar{\mu} \left( \frac{\partial \bar{u_i}}{\partial \bar{x_j}} + \frac{\partial \bar{u_j}}{\partial \bar{x_i}}\right)\right] + \bar{\kappa} \hat{\bf n} \delta_s.$

Here, $\mu_c$ the characteristic viscosity of the liquid sheet given as $\mu_c \left(\equiv \mu_p + \tau_y/\dot{\gamma}_c \right)$  and $\dot{\gamma}_c$ is the characteristic strain rate and in the present study, it is the inverse of the capillary time ($\sqrt{\sigma/\rho_1 h_0^3}$). The viscosity of the viscoplastic material is calculated as
$\bar{\mu}_1 = \mu_1/\mu_c= (1-Pl)\left[ 1+ {Pl}/{((1-Pl)\bar{\dot{\gamma}})}  \right]$
where $Pl\left( \equiv {\tau_y}/(\tau_y+\mu_p \dot{\gamma}_c) \right)$  is Plastic number and is a representation of the dimensionless yield stress of the material.

## Code
*/
#include "navier-stokes/centered.h"
# include "vof.h"
#include "tension.h"

scalar f[], * interfaces = {f};
scalar eta[], shear[], s_reg[];
double rho1 = 1., mu1 = 1., rho2 = 1e-3, mu2 = 1e-3;

double tmax;
double etamx;

/**The density inside the droplet is $\rho_1$ and outside $\rho_2$ */
# define rho(f) (clamp(f,0.,1.)*(rho1 - rho2) + rho2)
#define LEVEL 10
#define INLEVEL 6
#define MINLEVEL 6

/**We are simulating only one forth of the liquid sheet.*/
#define L0 15
#define L 9

/**Parameters */
#define R 1.0
#define Oh 0.1
#define Pl 0.1

double Bi = Pl/(1.-Pl);
double etac =  1/(1-Pl);  //mu1/(1-Pl);
double tC = 1; //R*R*Oh*(1.-Pl);
double utc = 1./((1-Pl)*R*Oh);//velocity scale for dimensional form
double tEnd = 100;
/**
Some "smearing" of the density/viscosity jump. */
scalar cf[];
face vector alphav[];
scalar rhov[];
face vector muv[];

int main() {
  printf ("Pl = %g Bi = %g etac = %g\n", Pl, Bi, etac);
  /**
  The domain */
  size (L0);
  init_grid (1 << INLEVEL);
  TOLERANCE = 1e-5;
  NITERMAX = 200;

  /**The density and viscosity are variable. */
  alpha = alphav;
  rho = rhov;
  mu = muv;

    etamx=100000;    
    tmax = 199.99;
    DT = 0.025;

  f.sigma = 1.0;

  run();
}
/**
The viscosity is defined at each timestep via the *properties()* event
declared by the Navier--Stokes solver. */

event properties (i++) {

    trash ({alpha});
    foreach() {
        eta[] = etamx;
        double D2 = 0.;
        foreach_dimension() {
            double dxx = u.x[1,0] - u.x[-1,0];
            double dxy = (u.x[0,1] - u.x[0,-1] + u.y[1,0] - u.y[-1,0])/2.;
            D2 += sq(dxx) + sq(dxy);
        }
        if (D2 > 0.) {
            D2 = sqrt(D2)/(2.*Delta);
            eta[] = (Bi/(sqrt(2)*D2 + Pl/etamx) + 1)*mu1 ;
            shear[] = D2;

        }
    }
    boundary ({eta});
    boundary ({shear});
    scalar fa[];
    // filtering density twice
    foreach()
    fa[] =  (4.*f[] +
             2.*(f[-1,0] + f[1,0] + f[0,-1] + f[0,1]) +
             f[1,1] + f[-1,1] + f[1,-1] + f[-1,-1])/16.;
    boundary ({fa});
   if (i>1){
 foreach() {
   if (fa[] > 0){
   if ((fa[]*shear[]*(sqrt(2))) > ((Pl*Oh)/(muv.x[]))) {
    s_reg[] = 0;
     }
   else {
    s_reg[] = 1; 
     }
    }
   else
    s_reg[] = 0;
    }
    }
    foreach_face() {
        double fm = (fa[] + fa[-1,0])/2.;
        muv.x[] =  Oh*(1-Pl)*(fm*(eta[] + eta[-1,0])/2. + (1. - fm)*mu2/(1-Pl));
//        muv.x[] =  Oh*(fm*(eta[] + eta[-1,0])/2. + (1. - fm)*mu2);
        //        muv.x[] = 1./(2.*fm/(eta[] + eta[-1,0]) + (1. - fm)/mu2);
        alphav.x[] = 1./rho(fm);
    }
    foreach()
    rhov[] = rho(fa[]); 
    boundary ({muv,alphav,rhov});
}

double geometry (double x, double y) {
  double Line = y-R;
  double Line_vert = x-L;
  double Circle = sq(x-L)+sq(y)-sq(R);
  double right_part = max(-Circle,-Line_vert);
  return min(-Line, right_part);
}

       /**
       We initialise the geometry of the interface.*/

event init (t = 0) {
    if (!restore (file = "dump-xxxx")){
  refine (x < (L + 2*R) && y < (6*R) && level < LEVEL);
  refine (x < (L + 3*R) && y > (6*R) && y < (7*R) && level < LEVEL-1.);
  refine (x > (L + 2*R) && x < (L + 3*R) && y > 0. && y < (7*R) && level < LEVEL-1.);
  refine (x < (L + 4*R) && y > (7*R) && y < (8*R) && level < LEVEL-2.);
  refine (x > (L + 3*R) && x < (L + 4*R) && y > 0. && y < (8*R) && level < LEVEL-2.);
  refine (x < (L + 5*R) && y > (8*R) && y < (9*R) && level < LEVEL-3.);
  refine (x > (L + 4*R) && x < (L + 5*R) && y > 0. && y < (9*R) && level < LEVEL-3.);
            fraction(f, geometry(x,y));
       #ifndef cf
         foreach()
           cf[] = f[];
       #endif
    }

  scalar le[];
  foreach(){
    le[] = level;
  }
  char name_grid[100];
  sprintf (name_grid, "grid.png");
  FILE* fgrid = fopen (name_grid, "w");
  output_ppm (le, file = name_grid,min=MINLEVEL,max=LEVEL,n = 1000,box = {{0.,0.},{L0,L0}});
  fclose (fgrid);

}


/**
We output : the interface position and scalar / vector fields . */

event saveDatas (t += tC; t <= tEnd) {

  /**
  We only output the interface in this function. */

  char name[100];
  sprintf (name, "interface_%g_%d.dat", t,pid());

  FILE * fp = fopen (name, "w");
  scalar pid[], ff[];
  foreach() {
    pid[] = fmod(pid()*(npe() + 37), npe());
    ff[] = f[] < 1.0e-6 ? 0 : f[] > 1. - 1.0e-6 ? 1. : f[];
  }
  boundary ({pid,ff});

  output_facets (ff, fp);

  fclose (fp);  


  char s[80];
  sprintf (s, "field_%g_%d.dat", t,pid());
  FILE * fpmu = fopen (s, "w");
//  scalar pid[], ff[];
  foreach() {
    pid[] = fmod(pid()*(npe() + 37), npe());
    ff[] = f[] < 1.0e-6 ? 0 : f[] > 1. - 1.0e-6 ? 1. : f[];
  }
  boundary ({pid,ff});

  //output_facets (ff, fp);
  foreach(){
  fprintf (fpmu," %f %f %g %g %g\n",x,y,shear[]*sqrt(2),muv.x[]/Oh, s_reg[]);
 // output_field ({f, muv}, fpf, linear = true);
  }
  fclose (fpmu); 

}

/**
We output the maximum x position of the liquid finger. Indeed this
position should linearly evolve in time, with a few variations coming
from the capillary wave. */

double xPrev = -1, tPrev = -1, veloPrev = 1;

event extractPosition (i ++) {

  /**
  We define a new vector field, h. */

  vector h[];

  /**
  We reconstruct the height function field and take the corresponding 
  maximum along the x axis. */

  heights (f, h);
  double xMax = -HUGE;;
  foreach(reduction(max:xMax))
    if (h.x[] != nodata) {
      double xi = x + height(h.x[])*Delta;
      if (xi > xMax)
  xMax = xi;
  }
  /**
  We also output the velocity of the end of the bulge.*/
  //if (fabs(xPrev - xMax) >= 0.01){
  double veloTip = tPrev >= 0 ? (xPrev - xMax)/(t - tPrev) : 0.;
  char name[80];
  sprintf (name, "velocity_pos.dat");
  static FILE * fp = fopen (name, "w");

  fprintf (fp," %g %g %g\n", t, xMax, veloTip);

  fflush (fp);
  //fprintf(stderr, "%g %g %g\n", t, xMax, veloTip);
  //fflush (stderr);
  tPrev = t, xPrev = xMax, veloPrev=abs(veloTip);
 //}
/**
We output, in the standard output file, the step with the corresponding
time. */

}

event logfile(i++) {
  printf ("i = %d t = %g\n", i,t);
  //printf ("utc = %g", utc);
  fflush(stdout);
}

event end (t = tEnd) {}