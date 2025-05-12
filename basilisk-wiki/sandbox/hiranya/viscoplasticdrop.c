/**
#  Retraction of a viscoplastic drop
This is a modification of the [two dimensional code](/sandbox/hiranya/viscoplasticsheet.c) to study the retraction of slender axisymmetric drops. 

We consider a slender drop of radius $R_0$ and length $L_0$. The Bingham model is used to mimic the viscoplastic properties. The slender drop has yield stress $\tau_y$, plastic viscosity $\mu_p$, surface tension $\sigma$ and density $\rho_1$. The constitutive equations are given in 
[column_SCC.c](/sandbox/M1EMN/Exemples/column_SCC.c).

##  Non-dimensionalization 
We use the same scaling as mentioned in the [viscoplastic sheet](/sandbox/hiranya/viscoplasticsheet.c) case.

## Code
*/ 

#include "axi.h"
#include "navier-stokes/centered.h"
# include "vof.h"
#include "tension.h"

/**
The interface is represented by the volume fraction field *f*. */

scalar f[], * interfaces = {f};
scalar eta[], shear[], s_reg[];
        scalar D3[] ; 
double rho1 = 1., mu1 = 1., rho2 = 1e-3, mu2 = 1e-3;

double tmax;
double etamx;

/**
The density inside the droplet is rho1 and outside rho2 */
#ifndef rho
# define rho(f) (clamp(f,0.,1.)*(rho1 - rho2) + rho2)
#endif

#define LEVEL 10 
#define INLEVEL 5
#define MINLEVEL 5

/**
We are simulating only one half of the liquid sheet.*/
#define L0 15
#define L 9

/**Parameters */
#define R 1.0
#define Oh 0.8
#define Pl 0.8

double Bi = Pl/(1.-Pl);
double etac =  1/(1-Pl);  
double tC = 1; 
double utc = 1./((1-Pl)*R*Oh);
double tEnd = 10.1;
/**
Some "smearing" of the density/viscosity jump. */
scalar cf[];

/**
The density and viscosity are variables. */

face vector alphav[];
scalar rhov[];
face vector muv[];

int main() {
  printf ("Pl = %g Bi = %g etac = %g rho2=%g\n", Pl, Bi, etac, rho2);
  /**
  The domain */
  size (L0);
  init_grid (1 << INLEVEL);
  //mass conservation tolerance : very important parameter
  TOLERANCE = 1e-5;
  NITERMAX = 200;


  /**
  The density and viscosity are variable. */
  alpha = alphav;
  rho = rhov;
  mu = muv;
 
    etamx=100000*etac;
    tmax = 199.99;
    DT = 0.025;
    f.sigma = 1.0;
  run();
}

/**The density is defined at each timestep via the *properties()* event
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
            D2 = sqrt((D2/(4.*Delta*Delta)) + sq(u.y[]/y));
            eta[] = (Bi/(sqrt(2)*D2 + mu1*Bi/etamx) + 1)*mu1 ;
            shear[] = D2;
        }
    }
    boundary ({eta});
    boundary ({shear});
    scalar fa[];
    foreach()
    fa[] =  (4.*f[] +
             2.*(f[-1,0] + f[1,0] + f[0,-1] + f[0,1]) +
             f[1,1] + f[-1,1] + f[1,-1] + f[-1,-1])/16.;

    boundary ({fa});

  /**We seperate yielded/unyielded regions within the drop */
 foreach() {
   if (fa[] > 0.5){
   if ((fa[]*shear[]*(sqrt(2))) > ((Pl*Oh)/(muv.x[]/fm.x[]))) {
    s_reg[] = 0;
     }
   else {
    s_reg[] = 1; 
     }
    }
   else
    s_reg[] = 0;
    }
/** and update the viscosity */
  foreach_face() {
    double ff = (fa[] + fa[-1])/2.;
    alphav.x[] = fm.x[]/rho(ff);
      muv.x[] = fm.x[]*Oh*(1-Pl)*(ff*(eta[] + eta[-1,0])/2. + (1. - ff)*mu2/(1-Pl));
  }
  foreach()
    rhov[] = cm[]*rho(fa[]);

    boundary ({muv,alphav,rhov});
}

/** We define the viscoplastic drop by the variavle geometry*/
double geometry (double x, double y) {
  double Line = y-R;
  double Line_vert = x-L;
  double Circle = sq(x-L)+sq(y)-sq(R);
  double right_part = max(-Circle,-Line_vert);
  return min(-Line, right_part);
}

      /**
       We allocate the volume fraction to the drop defined above.*/

event init (t = 0) {
    if (!restore (file = "dump-xxxx")){

  refine (x < (L + 2*R) && y < (5*R) && level < LEVEL);
  refine (x < (L + 3*R) && y > (5*R) && y < (6*R) && level < LEVEL-1.);
  refine (x > (L + 2*R) && x < (L + 3*R) && y > 0. && y < (6*R) && level < LEVEL-1.);
  refine (x < (L + 4*R) && y > (6*R) && y < (7*R) && level < LEVEL-2.);
  refine (x > (L + 3*R) && x < (L + 4*R) && y > 0. && y < (7*R) && level < LEVEL-2.);
  refine (x < (L + 5*R) && y > (7*R) && y < (8*R) && level < LEVEL-3.);
  refine (x > (L + 4*R) && x < (L + 5*R) && y > 0. && y < (8*R) && level < LEVEL-3.); 
  refine (x < (L + 6*R) && y > (8*R) && y < (9*R) && level < LEVEL-4.);
  refine (x > (L + 5*R) && x < (L + 6*R) && y > 0. && y < (9*R) && level < LEVEL-4.); 
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
  /** We print the initial grid */
  char name_grid[100];
  sprintf (name_grid, "grid.png");
  FILE* fgrid = fopen (name_grid, "w");
  output_ppm (le, file = name_grid,min=MINLEVEL,max=LEVEL,n = 1000,box = {{0.,0.},{L0,L0}});
  fclose (fgrid);

}


/**
We output the interface position and scalar / vector fields . */

event saveDatas (t += tC; t <= tEnd) {

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

  char frac2[80];
  sprintf (frac2, "vof-%f_%d.dat", t, pid());
  FILE * fpf = fopen (frac2, "w");
  foreach() {
    pid[] = fmod(pid()*(npe() + 37), npe());
    ff[] = f[] < 1.0e-6 ? 0 : f[] > 1. - 1.0e-6 ? 1. : f[];
  }
  boundary ({pid,ff});
  foreach(){
  fprintf(fpf, "%f %f %g\n",x,y,f[]);
  }
 // output_field ({p}, fpf, linear = true);
  fclose (fpf);

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
  fprintf (fpmu," %f %f %g %g %g\n",x,y,shear[]*sqrt(2),muv.x[]/(fm.x[]*Oh),s_reg[]);
 // output_field ({f, muv}, fpf, linear = true);
  }
  fclose (fpmu); 

}

/**
We output the tip position of the slender drop. */

double xPrev = -1, tPrev = -1, veloPrev = 1;

event extractPosition (i ++) {
  vector h[];
  
  /**
  We reconstruct the height function field and take the corresponding 
  maximum along the $x$-axis. */
  
  heights (f, h);
  double xMax = -HUGE;;
  foreach(reduction(max:xMax))
    if (h.x[] != nodata) {
      double xi = x + height(h.x[])*Delta;
      if (xi > xMax)
  xMax = xi;
  }
  /**We also output the velocity at the tip of the drop.*/
  
  double veloTip = tPrev >= 0 ? (xPrev - xMax)/(t - tPrev) : 0.;
  char name[80];
  sprintf (name, "velocity_pos.dat");
  static FILE * fp = fopen (name, "w");

  fprintf (fp," %g %g %g\n", t, xMax, veloTip);
 
  fflush (fp);
  tPrev = t, xPrev = xMax, veloPrev=abs(veloTip);
 }
/**
We output, in the standard output file, the step with the corresponding
time. */
event logfile(i++) {
  printf ("i = %d t = %g\n", i,t);
  //printf ("utc = %g", utc);
  fflush(stdout);
}


/* Movies post processing */
/*
event movies (t += 0.01; t <= tEnd){
  scalar vort[];
  vorticity (u, vort);
  scalar m[];
  foreach() {
    if (f[]<0.95 && f[] >0.05)
      m[] = -10;
    else
      m[] = 0.;
  }*/
//  static FILE * fo = popen ("ppm2gif > vort.gif", "w");
//  output_ppm (vort, fo, mask =m,min=-1.,max=1., n = 1000,box = {{0.,0.},{10.,2.}});
/*
  char name_vort[100];
  sprintf (name_vort, "vort-%g.png", t);
  FILE* fvort = fopen (name_vort, "w");
  double om = (mu1 /Oh)*1./rho1/R/R;
  output_ppm (vort, fvort, mask =m,min=-om,max=om, map=cool_warm, n = 1000,box = {{0.,0.},{L0,3*R}});
  fclose (fvort);


}*/

event end (t = tEnd) {}

